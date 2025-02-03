#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoParticleFlow/PFProducer/interface/MLPFModel.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"

#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonData.h"

// using namespace cms::Ort;
class MLPFSONICProducer : public TritonEDProducer<> {
public:
  explicit MLPFSONICProducer(const edm::ParameterSet &);
  ~MLPFSONICProducer() override;
  
  void acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput) override;
  
  void produce(edm::Event &iEvent, edm::EventSetup const &iSetup, Output const &iOutput) override;
  static void fillDescriptions(edm::ConfigurationDescriptions & );

  // static methods for handling the global cache
   // static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  //      static void globalEndJob(const ONNXRuntime*);

private:
  const edm::EDPutTokenT<reco::PFCandidateCollection> pfCandidatesPutToken_;
  const edm::EDGetTokenT<edm::View<reco::GsfElectron>> gsfElectrons_;
  const edm::EDGetTokenT<reco::PFBlockCollection> inputTagBlocks_;
};

class SelectedElementsManager {
public:
    static SelectedElementsManager& getInstance() {
        static SelectedElementsManager instance; // Single instance for the program
        return instance;
    }

    void fill(const std::vector<const reco::PFBlockElement*>& all_elements) {
        selected_elements_.clear();

        for (const auto* pelem : all_elements) {
            if (pelem->type() == reco::PFBlockElement::PS1 || 
                pelem->type() == reco::PFBlockElement::PS2 || 
                pelem->type() == reco::PFBlockElement::BREM) {
                continue;
            }
            selected_elements_.push_back(pelem);
        }
    }

    const std::vector<const reco::PFBlockElement*>& get() const {
        return selected_elements_;
    }

private:
    SelectedElementsManager() = default;
    ~SelectedElementsManager() = default;

    SelectedElementsManager(const SelectedElementsManager&) = delete;
    SelectedElementsManager& operator=(const SelectedElementsManager&) = delete;

    std::vector<const reco::PFBlockElement*> selected_elements_;
};



MLPFSONICProducer::MLPFSONICProducer(const edm::ParameterSet &iConfig)    
    : TritonEDProducer<>(iConfig),
      pfCandidatesPutToken_{produces<reco::PFCandidateCollection>()},
      gsfElectrons_{consumes<edm::View<reco::GsfElectron>>(edm::InputTag("gedGsfElectronsTmp"))},
      inputTagBlocks_{consumes<reco::PFBlockCollection>(iConfig.getParameter<edm::InputTag>("src"))} {}


// // Define the function
// std::vector<const reco::PFBlockElement*> filterSelectedElements(
//     const std::vector<const reco::PFBlockElement*>& all_elements) {
    
//     std::vector<const reco::PFBlockElement*> selected_elements;

//     for (const auto* pelem : all_elements) {
//         // Skip elements with type PS1, PS2, or BREM
//         if (pelem->type() == reco::PFBlockElement::PS1 || 
//             pelem->type() == reco::PFBlockElement::PS2 || 
//             pelem->type() == reco::PFBlockElement::BREM) {
//             continue;
//         }
//         selected_elements.push_back(pelem);
//     }

//     return selected_elements;
// }

MLPFSONICProducer::~MLPFSONICProducer() {}
void MLPFSONICProducer::acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput){

  using namespace reco::mlpf;
  const auto& blocks = iEvent.get(inputTagBlocks_);
  const auto& all_elements = getPFElements(blocks);

  const auto& gsfElectrons = iEvent.get(gsfElectrons_);

  
  SelectedElementsManager::getInstance().fill(all_elements); // Fill data once
  std::cout << "filled selected_elements." << std::endl;
  const auto& selected_elements = SelectedElementsManager::getInstance().get();
  // Total Number of selected_elements 
  unsigned int num_elements_total = selected_elements.size();


  const auto tensor_size = num_elements_total;

// #ifdef MLPF_DEBUG
//   assert(num_elements_total < NUM_MAX_ELEMENTS_BATCH);
//   //tensor size must be a multiple of the bin size and larger than the number of elements
//   assert(tensor_size <= NUM_MAX_ELEMENTS_BATCH);
//   assert(tensor_size % LSH_BIN_SIZE == 0);
// #endif

#ifdef MLPF_DEBUG
  std::cout << "tensor_size=" << tensor_size << std::endl;
#endif

  //Fill the input tensor (batch, elems, features) = (1, tensor_size, NUM_ELEMENT_FEATURES)
  std::vector<std::vector<float>> inputs;
  //std::vector<TritonInputContainer<std::vector<float>>> inputs;
  //inputs.push_back(TritonInputContainer<std::vector<float>>(NUM_ELEMENT_FEATURES * tensor_size, 0.0));
  //inputs.push_back(TritonInputContainer<std::vector<float>>(tensor_size, 0.0));
    
  inputs.push_back(std::vector<float>(NUM_ELEMENT_FEATURES * tensor_size, 0.0));
  inputs.push_back(std::vector<float>(tensor_size, 0.0));  // mask
  unsigned int ielem = 0;
  const char* mask = "mask";
  const char* Xfeat = "Xfeat_normed";
    
  auto &data1 = iInput.at(mask);
  auto &data2 = iInput.at(Xfeat);
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
    inputs[1][ielem] = 1.0;
    ielem += 1;
    }
    auto wrapped_input1 = std::vector<std::vector<float>>{inputs[1]};
    auto input_ptr1 = std::make_shared<TritonInput<float>>(wrapped_input1);
    auto wrapped_input0 = std::vector<std::vector<float>>{inputs[0]};
    auto input_ptr0 = std::make_shared<TritonInput<float>>(wrapped_input0);
    TritonInputContainer<float> input_container1 = input_ptr1;
    TritonInputContainer<float> input_container0 = input_ptr0;
    data1.toServer(input_container1);
    data2.toServer(input_container0);
//     for (const auto& input : inputs) {
//         auto input_ptr = std::make_shared<TritonInput<float>>(input);
//         TritonInputContainer<float> input_container = input_ptr;
//         data1.toServer(input_container);
// }

}

void MLPFSONICProducer::produce(edm::Event &iEvent,
                                const edm::EventSetup &iSetup,
                                Output const &iOutput){
    using namespace reco::mlpf;
    const auto& output1 = iOutput.at("bid");
    const auto& output2 = iOutput.at("id");
    const auto& output3 = iOutput.at("momentum");
    const auto & output_binary = output1.fromServer<float>();
    const auto & output_pid = output2.fromServer<float>();
    const auto & output_p4 = output3.fromServer<float>();
    const auto& selected_elements = SelectedElementsManager::getInstance().get();
    // Total Number of selected_elements 
    unsigned int num_elements_total = selected_elements.size();
    unsigned int tensor_size = num_elements_total;
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

      const auto logit_no_ptcl = output_binary[ielem * 2 + 0].begin();
      const auto logit_ptcl = output_binary[ielem * 2 + 1].begin();
      // Check if the binary classifier of the model predicted a particle 
      int pred_pid = 0;
      if (logit_ptcl > logit_no_ptcl) {
        for (unsigned int idx_id = 0; idx_id < pred_id_probas.size(); idx_id++) {
            auto pred_proba = output_pid[ielem * NUM_OUTPUT_FEATURES_CLS + idx_id].begin();
  #ifdef MLPF_DEBUG
        std::cout << "pid proba: " << pred_proba << std::endl;
        assert(!std::isnan(pred_proba));
  #endif
            pred_id_probas[idx_id] = *pred_proba;
        }
  
        auto imax = argMax(pred_id_probas);
  
        //get the most probable class PDGID
        pred_pid = pdgid_encoding.at(imax);
      }
    



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
        float pred_pt = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_PT].begin();
       // pred_pt = exp(pred_pt) * inputs[0][ielem * NUM_ELEMENT_FEATURES + 1]; 
        float pred_eta = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_ETA].begin();
        float pred_sin_phi = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_SIN_PHI].begin();
        float pred_cos_phi = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_COS_PHI].begin();
        float pred_e = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_ENERGY].begin();
      //  pred_e = exp(pred_e) * inputs[0][ielem * NUM_ELEMENT_FEATURES + 5];
      
        auto cand = makeCandidate(pred_pid, pred_charge, pred_pt, pred_eta, pred_sin_phi, pred_cos_phi, pred_e);
        setCandidateRefs(cand, selected_elements, ielem);
        pOutputCandidateCollection.push_back(cand);
      }
    }//end loop of elements

    iEvent.emplace(pfCandidatesPutToken_, pOutputCandidateCollection);
}
    

// std::unique_ptr<ONNXRuntime> MLPFSONICProducer::initializeGlobalCache(const edm::ParameterSet& params) {
//   return std::make_unique<ONNXRuntime>(params.getParameter<edm::FileInPath>("model_path").fullPath());
// }

void MLPFSONICProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("particleFlowBlock"));
  // desc.add<edm::FileInPath>("model_path",
  //                           edm::FileInPath("RecoParticleFlow/PFProducer/data/mlpf/"
  //                                           "dev.onnx"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(MLPFSONICProducer);

