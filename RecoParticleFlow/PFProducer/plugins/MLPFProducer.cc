#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoParticleFlow/PFProducer/interface/MLPFModel.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"


using namespace cms::Ort;

//use this to switch on detailed print statements in MLPF
//#define MLPF_DEBUG

class MLPFProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
  explicit MLPFProducer(const edm::ParameterSet&, const ONNXRuntime*);

  void produce(edm::Event& event, const edm::EventSetup& setup) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // static methods for handling the global cache
  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ONNXRuntime*);

private:
  const edm::EDPutTokenT<reco::PFCandidateCollection> pfCandidatesPutToken_;
  const edm::EDGetTokenT<edm::View<reco::GsfElectron>> gsfElectrons_;
  const edm::EDGetTokenT<reco::PFBlockCollection> inputTagBlocks_;
};

MLPFProducer::MLPFProducer(const edm::ParameterSet& cfg, const ONNXRuntime* cache)
    : pfCandidatesPutToken_{produces<reco::PFCandidateCollection>()},
      gsfElectrons_{consumes<edm::View<reco::GsfElectron>>(edm::InputTag("gedGsfElectronsTmp"))},
      inputTagBlocks_{consumes<reco::PFBlockCollection>(cfg.getParameter<edm::InputTag>("src"))} {}

void MLPFProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  using namespace reco::mlpf;

  const auto& blocks = event.get(inputTagBlocks_);
  const auto& all_elements = getPFElements(blocks);

  const auto& gsfElectrons = event.get(gsfElectrons_);

  std::vector<const reco::PFBlockElement*> selected_elements;
  unsigned int num_elements_total = 0;
  for (const auto* pelem : all_elements) {
    if (pelem->type() == reco::PFBlockElement::PS1 || pelem->type() == reco::PFBlockElement::PS2 || pelem->type() == reco::PFBlockElement::BREM) {
      continue;
    }
    num_elements_total += 1;
    selected_elements.push_back(pelem);
  }
  const auto tensor_size = num_elements_total;

#ifdef MLPF_DEBUG
  std::cout << "tensor_size=" << tensor_size << std::endl;
#endif

  //Fill the input tensor (batch, elems, features) = (1, tensor_size, NUM_ELEMENT_FEATURES)
  //std::vector<TritonInputContainer<std::vector<float>>> inputs;
  std::vector<std::vector<float>> inputs;
  inputs.push_back(std::vector<float>(NUM_ELEMENT_FEATURES * tensor_size, 0.0));
  //inputs.push_back(TritonInputContainer<std::vector<float>>(NUM_ELEMENT_FEATURES * tensor_size, 0.0));
  inputs.push_back(std::vector<float>(tensor_size, 0.0));
  unsigned int ielem = 0;
  for (const auto* pelem : selected_elements) {
    if (ielem > tensor_size) {
      continue;
    }

    const auto& elem = *pelem;

    //prepare the input array from the PFElement
    const auto& props = getElementProperties(elem, gsfElectrons).as_array();

    //copy features to the input array
    for (unsigned int iprop = 0; iprop < NUM_ELEMENT_FEATURES; iprop++) {
      inputs[0][ielem * NUM_ELEMENT_FEATURES + iprop] = normalize(props[iprop]);
    }
    //mask
    inputs[1][ielem] = 1.0;
    ielem += 1;
  }

  //run the GNN inference, given the inputs and the output.
  const auto& outputs = globalCache()->run({"Xfeat_normed", "mask"}, inputs, {{1, tensor_size, NUM_ELEMENT_FEATURES}, {1, tensor_size}});
  const auto& output_binary = outputs[0];
  const auto& output_pid = outputs[1];
  const auto& output_p4 = outputs[2];
  std::cout << "check-point 94_"<< tensor_size << std::endl;
  std::cout << "check-point 95_"<< output_binary[0] <<"_____________________________________" << output_binary[1]<< std::endl;
 // std::cout << "check-point 96_"<< output_pid[0] <<"_____________________________________" << output_pid[1]<< std::endl;
 // std::cout << "check-point 97_"<< output_p4[0] <<"_____________________________________" << output_p4[1]<< std::endl;
    
#ifdef MLPF_DEBUG
  std::cout << "output_binary=" << output_binary.size() << std::endl;  
  assert(output_binary.size() == tensor_size * 2);

  std::cout << "output_pid=" << output_pid.size() << std::endl;  
  assert(output_pid.size() == tensor_size * NUM_OUTPUT_FEATURES_CLS);

  std::cout << "output_p4=" << output_p4.size() << std::endl;  
  assert(output_p4.size() == tensor_size * NUM_OUTPUT_FEATURES_P4);
#endif

  std::vector<reco::PFCandidate> pOutputCandidateCollection;
  for (size_t ielem = 0; ielem < num_elements_total; ielem++) {
    std::vector<float> pred_id_probas(pdgid_encoding.size(), 0.0);
    const reco::PFBlockElement* elem = selected_elements[ielem];

#ifdef MLPF_DEBUG
    std::cout << "ielem=" << ielem << " inputs:";
    for (unsigned int iprop = 0; iprop < NUM_ELEMENT_FEATURES; iprop++) {
      std::cout << iprop << "=" << inputs[0][ielem * NUM_ELEMENT_FEATURES + iprop] << " ";
    }
    std::cout << std::endl;
#endif

    const auto logit_no_ptcl = output_binary[ielem * 2 + 0];
    const auto logit_ptcl = output_binary[ielem * 2 + 1];
#ifdef MLPF_DEBUG
    std::cout << "binary: " << logit_no_ptcl << " " << logit_ptcl << std::endl;
#endif
   
    // Check if the binary classifier of the model predicted a particle 
    int pred_pid = 0;
    if (logit_ptcl > logit_no_ptcl) {
      for (unsigned int idx_id = 0; idx_id < pred_id_probas.size(); idx_id++) {
        auto pred_proba = output_pid[ielem * NUM_OUTPUT_FEATURES_CLS + idx_id];
#ifdef MLPF_DEBUG
        std::cout << "pid proba: " << pred_proba << std::endl;
        assert(!std::isnan(pred_proba));
#endif
        pred_id_probas[idx_id] = pred_proba;
      }
  
      auto imax = argMax(pred_id_probas);
  
      //get the most probable class PDGID
      pred_pid = pdgid_encoding.at(imax);
    }

#ifdef MLPF_DEBUG
    std::cout << "p4: "
      << output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + 0] << " "
      << output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + 1] << " " 
      << output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + 2] << " " 
      << output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + 3] << " " 
      << output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + 4] << std::endl;
#endif

    //a particle was predicted for this PFElement, otherwise it was a spectator
    if (pred_pid != 0) {
      //muons and charged hadrons should only come from tracks, otherwise we won't have track references to pass downstream
      if (((pred_pid == 13) || (pred_pid == 211)) && elem->type() != reco::PFBlockElement::TRACK) {
        pred_pid = 130;
      }

      float pred_charge = 0.0;
      if (elem->type() == reco::PFBlockElement::TRACK) {
        const auto* eltTrack = dynamic_cast<const reco::PFBlockElementTrack*>(elem);
	//for now, just take the charge from the track
	if (eltTrack->trackRef().isNonnull()) {
          pred_charge = eltTrack->trackRef()->charge();
        }

        //a track with no muon ref should not produce a muon candidate, instead we interpret it as a charged hadron here
        if (pred_pid == 13 && eltTrack->muonRef().isNull()) {
          pred_pid = 211;
        }
       
        //taus are reconstructed downstream based on other criteria, instead we interpret it as a charged hadron here
	if (pred_pid == 15) {
          pred_pid = 211;
        }

        //tracks from displaced vertices need reference debugging downstream as well, so we just treat them as neutrals for the moment
        if ((pred_pid == 211) && (eltTrack->isLinkedToDisplacedVertex())) {
          pred_pid = 130;
        }
      }

      //get the predicted momentum components from the model
      float pred_pt = output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_PT];
      pred_pt = exp(pred_pt) * inputs[0][ielem * NUM_ELEMENT_FEATURES + 1]; 
      float pred_eta = output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_ETA];
      float pred_sin_phi = output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_SIN_PHI];
      float pred_cos_phi = output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_COS_PHI];
      float pred_e = output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_ENERGY];
      pred_e = exp(pred_e) * inputs[0][ielem * NUM_ELEMENT_FEATURES + 5];

      auto cand = makeCandidate(pred_pid, pred_charge, pred_pt, pred_eta, pred_sin_phi, pred_cos_phi, pred_e);
      setCandidateRefs(cand, selected_elements, ielem);
      pOutputCandidateCollection.push_back(cand);

#ifdef MLPF_DEBUG
      std::cout << "ielem=" << ielem << " pred: pid=" << cand.pdgId() << " E=" << cand.energy() << " pt=" << cand.pt()
                << " eta=" << cand.eta() << " phi=" << cand.phi() << " charge=" << cand.charge() << std::endl;
#endif
    }
  }  //loop over PFElements

  event.emplace(pfCandidatesPutToken_, pOutputCandidateCollection);
}

std::unique_ptr<ONNXRuntime> MLPFProducer::initializeGlobalCache(const edm::ParameterSet& params) {
  auto session_options = ONNXRuntime::defaultSessionOptions(params.getParameter<bool>("use_cuda") ? Backend::cuda : Backend::cpu);
  return std::make_unique<ONNXRuntime>(params.getParameter<edm::FileInPath>("model_path").fullPath(), &session_options);
}

void MLPFProducer::globalEndJob(const ONNXRuntime* cache) {}

void MLPFProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("particleFlowBlock"));
  desc.add<edm::FileInPath>("model_path",
                            edm::FileInPath("RecoParticleFlow/PFProducer/data/mlpf/"
                                            "mlpf_5M_attn2x3x256_bm5_relu_checkpoint5_1xa100_fp32_fused.onnx"));
  desc.add<bool>("use_cuda", false);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(MLPFProducer);
