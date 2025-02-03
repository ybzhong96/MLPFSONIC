#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoParticleFlow/PFProducer/interface/MLPFModel.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"

#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonData.h"

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



MLPFProducer::MLPFSONICProducer(const edm::ParameterSet& cfg, const ONNXRuntime* cache)
    : pfCandidatesPutToken_{produces<reco::PFCandidateCollection>()},
      gsfElectrons_{consumes<edm::View<reco::GsfElectron>>(edm::InputTag("gedGsfElectronsTmp"))},
      inputTagBlocks_{consumes<reco::PFBlockCollection>(cfg.getParameter<edm::InputTag>("src"))} {}


// Define the function
std::vector<const reco::PFBlockElement*> filterSelectedElements(
    const std::vector<const reco::PFBlockElement*>& all_elements) {
    
    std::vector<const reco::PFBlockElement*> selected_elements;

    for (const auto* pelem : all_elements) {
        // Skip elements with type PS1, PS2, or BREM
        if (pelem->type() == reco::PFBlockElement::PS1 || 
            pelem->type() == reco::PFBlockElement::PS2 || 
            pelem->type() == reco::PFBlockElement::BREM) {
            continue;
        }
        selected_elements.push_back(pelem);
    }

    return selected_elements;
}


void MLPFSONICProducer::acquire(edm::Event& event, const edm::EventSetup& setup,
                                Input &iInput){

  using namespace reco::mlpf;
  const auto& blocks = event.get(inputTagBlocks_);
  const auto& all_elements = getPFElements(blocks);

  const auto& gsfElectrons = event.get(gsfElectrons_);

  
  SelectedElementsManager::getInstance().fill(all_elements); // Fill data once
  std::cout << "filled selected_elements." << std::endl;
  const auto& selected_elements = SelectedElementsManager::getInstance().get();
  // Total Number of selected_elements 
  unsigned int num_elements_total = selected_elements.size();


  const auto tensor_size = LSH_BIN_SIZE * std::max(2u, (num_elements_total / LSH_BIN_SIZE + 1));

#ifdef MLPF_DEBUG
  assert(num_elements_total < NUM_MAX_ELEMENTS_BATCH);
  //tensor size must be a multiple of the bin size and larger than the number of elements
  assert(tensor_size <= NUM_MAX_ELEMENTS_BATCH);
  assert(tensor_size % LSH_BIN_SIZE == 0);
#endif

#ifdef MLPF_DEBUG
  std::cout << "tensor_size=" << tensor_size << std::endl;
#endif

  //Fill the input tensor (batch, elems, features) = (1, tensor_size, NUM_ELEMENT_FEATURES)
  std::vector<std::vector<float>> inputs(1, std::vector<float>(NUM_ELEMENT_FEATURES * tensor_size, 0.0));
  unsigned int ielem = 0;
  
  auto &input = iInput #   
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
    ielem += 1;
  }
  input.toServer(inputs)

}

void MLPFSONICProducer::produce(edm::Event &iEvent,
                                const edm::EventSetup &iSetup,
                                Output const &iOutput){
    const auto &output1 = iOutput[0].begin()->second;
    const auto &output = output1.fromServer<float>();
   
    const auto& selected_elements = SelectedElementsManager::getInstance().get();
    // Total Number of selected_elements 
    unsigned int num_elements_total = selected_elements.size();
    
    std::vector<reco::PFCandidate> pOutputCandidateCollection;
    for (size_t ielem = 0; ielem < num_elements_total; ielem++) {
      std::vector<float> pred_id_probas(pdgid_encoding.size(), 0.0);
      // when to use that info
      const reco::PFBlockElement* elem = selected_elements[ielem];

      for (unsigned int idx_id = 0; idx_id < pred_id_probas.size(); idx_id++) {
        auto pred_proba = output[ielem * NUM_OUTPUT_FEATURES + idx_id];
#ifdef MLPF_DEBUG
        assert(!std::isnan(pred_proba));
#endif
        pred_id_probas[idx_id] = pred_proba;
      }
    }
    auto imax = argMax(pred_id_probas);

    //get the most probable class PDGID
    int pred_pid = pdgid_encoding.at(imax);

    //a particle was predicted for this PFElement, otherwise it was a spectator
    if (pred_pid != 0) {
      //muons and charged hadrons should only come from tracks, otherwise we won't have track references to pass downstream
      if (((pred_pid == 13) || (pred_pid == 211)) && elem->type() != reco::PFBlockElement::TRACK) {
        pred_pid = 130;
      }

      if (elem->type() == reco::PFBlockElement::TRACK) {
        const auto* eltTrack = dynamic_cast<const reco::PFBlockElementTrack*>(elem);

        //a track with no muon ref should not produce a muon candidate, instead we interpret it as a charged hadron
        if (pred_pid == 13 && eltTrack->muonRef().isNull()) {
          pred_pid = 211;
        }

        //tracks from displaced vertices need reference debugging downstream as well, so we just treat them as neutrals for the moment
        if ((pred_pid == 211) && (eltTrack->isLinkedToDisplacedVertex())) {
          pred_pid = 130;
        }
      }

      //get the predicted momentum components
      float pred_pt = output[ielem * NUM_OUTPUT_FEATURES + IDX_PT];
      float pred_eta = output[ielem * NUM_OUTPUT_FEATURES + IDX_ETA];
      float pred_sin_phi = output[ielem * NUM_OUTPUT_FEATURES + IDX_SIN_PHI];
      float pred_cos_phi = output[ielem * NUM_OUTPUT_FEATURES + IDX_COS_PHI];
      float pred_e = output[ielem * NUM_OUTPUT_FEATURES + IDX_ENERGY];
      float pred_charge = output[ielem * NUM_OUTPUT_FEATURES + IDX_CHARGE];
    
                                }
    }//end loop of elements

std::unique_ptr<ONNXRuntime> MLPFProducer::initializeGlobalCache(const edm::ParameterSet& params) {
  return std::make_unique<ONNXRuntime>(params.getParameter<edm::FileInPath>("model_path").fullPath());
}

void MLPFProducer::globalEndJob(const ONNXRuntime* cache) {}

void MLPFProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("particleFlowBlock"));
  desc.add<edm::FileInPath>("model_path",
                            edm::FileInPath("RecoParticleFlow/PFProducer/data/mlpf/"
                                            "dev.onnx"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(MLPFProducer);

