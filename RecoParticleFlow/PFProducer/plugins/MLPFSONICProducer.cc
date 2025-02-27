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
  std::vector<std::string> input_names_;
  std::vector<std::string> output_names_;
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
      inputTagBlocks_{consumes<reco::PFBlockCollection>(iConfig.getParameter<edm::InputTag>("src"))},
      input_names_(iConfig.getParameter<std::vector<std::string>>("input_names")),
      output_names_(iConfig.getParameter<std::vector<std::string>>("output_names"))
      {}

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

  //Fill the input tensor (batch, elems, features) = (1, tensor_size, NUM_ELEMENT_FEATURES)
  std::vector<std::vector<float>> inputs;
  inputs.push_back(std::vector<float>(tensor_size*NUM_ELEMENT_FEATURES, 0.0));
  inputs.push_back(std::vector<float>(tensor_size, 0.0));  // mask

  unsigned int ielem = 0;
  const auto &mask_name = input_names_[0];
  const auto &X_name = input_names_[1];
    
  auto &data1 = iInput.at(mask_name);
  auto &data2 = iInput.at(X_name);  
  data1.setShape(0, tensor_size);
  data2.setShape(0, tensor_size);
  auto tdata1 = data1.allocate<float>(true);
  auto tdata2 = data2.allocate<float>(true);
  for (const auto* pelem : selected_elements) {
    if (ielem > tensor_size) {
      continue;
    }

    const auto& elem = *pelem;

    //prepare the input array from the PFElement
    const auto& props = getElementProperties(elem, gsfElectrons).as_array();

    //copy features to the input array
    for (unsigned int iprop = 0; iprop < NUM_ELEMENT_FEATURES; iprop++) {
      inputs[0][ielem*NUM_ELEMENT_FEATURES+ iprop]=normalize(props[iprop]);
    }
    inputs[1][ielem] = 1.0;
    ielem += 1;
  }
  auto &vdata1 = (*tdata1)[0]; 
  auto &vdata2 = (*tdata2)[0];
  vdata1 = inputs[1];
  vdata2 = inputs[0];
  data1.toServer(tdata1); 
  data2.toServer(tdata2);
  std::cout << "check-point Producer-143_tensorsize_"<< tensor_size << std::endl;
  // for (unsigned int i = 0; i < tensor_size; i++) {
  //   std::cout << "check-point Producer-151_mask_ "<< (*tdata1)[0][i] << std::endl;
  //   std::cout << "check-point Producer-152_elem_"<< (*tdata2)[0][i] << std::endl; 
  // } Non-empty
}
void MLPFSONICProducer::produce(edm::Event &iEvent,
                                const edm::EventSetup &iSetup,
                                Output const &iOutput){
    using namespace reco::mlpf;
    const auto &bid_name = output_names_[0];
    const auto &id_name = output_names_[1];
    const auto &momentum_name = output_names_[2];
    const auto &output1 = iOutput.at(bid_name);
    const auto &output2 = iOutput.at(id_name);
    const auto &output3 = iOutput.at(momentum_name); 
    const auto &output_binary = output1.fromServer<float>();
    const auto &output_pid = output2.fromServer<float>();
    const auto &output_p4 = output3.fromServer<float>();
    // const auto &output_from_server = output0.fromServer<float>();
   // std::cout << "check-point Producer-163_(should be x2)"<< output_binary[0].size()<<std::endl;
   // std::cout << "check-point Producer-164_(should be x9)"<< output_pid[0].size()<<std::endl;
   // std::cout << "check-point Producer-165_(should x5)"<< output_p4[0].size()<<std::endl;
    std::cout << "check-point Producer-166_(should be x2)"<< output_binary[0][0] << "______"<< output_binary[0][1] <<"____"<< output_binary[0][2]<< "____"<< output_binary[0][3]<<std::endl;
    const auto &selected_elements = SelectedElementsManager::getInstance().get();
    // Total Number of selected_elements 
    unsigned int num_elements_total = selected_elements.size();
    unsigned int tensor_size = num_elements_total;
    //std::cout << "check-point Producer-171_"<<tensor_size<<std::endl; 
    
    std::vector<reco::PFCandidate> pOutputCandidateCollection;
    for (size_t ielem = 0; ielem < num_elements_total; ielem++) {
      std::vector<float> pred_id_probas(pdgid_encoding.size(), 0.0);
      const reco::PFBlockElement* elem = selected_elements[ielem];
      const auto logit_no_ptcl = output_binary[0][ielem * 2 + 0]; //.begin();
      const auto logit_ptcl = output_binary[0][ielem * 2 + 1]; //.begin();
     
      // Check if the binary classifier of the model predicted a particle 
      int pred_pid = 0;
      if (logit_ptcl > logit_no_ptcl) {
        for (unsigned int idx_id = 0; idx_id < pred_id_probas.size(); idx_id++) {
            auto pred_proba = output_pid[0][ielem * NUM_OUTPUT_FEATURES_CLS + idx_id]; //.begin();
            pred_id_probas[idx_id] = pred_proba; //*pred_proba;
        }
  
        auto imax = argMax(pred_id_probas);
        //get the most probable class PDGID
        pred_pid = pdgid_encoding.at(imax);
      }
      //a particle was predicted for this PFElement, otherwise it was a spectator
     // if (pred_pid != 0) {
      if(false){   
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
      // std::cout << "check-point Producer-224" << std::endl;
      //get the predicted momentum components from the model
       float pred_pt = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_PT].begin();
       //pred_pt = exp(pred_pt) * inputs[0][ielem * NUM_ELEMENT_FEATURES + 1]; 
       float pred_eta = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_ETA].begin();
       float pred_sin_phi = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_SIN_PHI].begin();
       float pred_cos_phi = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_COS_PHI].begin();
       float pred_e = *output_p4[ielem * NUM_OUTPUT_FEATURES_P4 + IDX_ENERGY].begin();
       //pred_e = exp(pred_e) * inputs[0][ielem * NUM_ELEMENT_FEATURES + 5];
      
       auto cand = makeCandidate(pred_pid, pred_charge, pred_pt, pred_eta, pred_sin_phi, pred_cos_phi, pred_e);
       setCandidateRefs(cand, selected_elements, ielem);
       pOutputCandidateCollection.push_back(cand);
     }
   }  //end loop of elements
    std::cout << "check-point Producer-239" << std::endl;
    iEvent.emplace(pfCandidatesPutToken_, pOutputCandidateCollection);
}
    

void MLPFSONICProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<edm::InputTag>("src", edm::InputTag("particleFlowBlock"));
  desc.add<std::vector<std::string>>("input_names", {"mask", "Xfeat_normed"});
  desc.add<std::vector<std::string>>("output_names", {"bid", "id", "momentum"});
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(MLPFSONICProducer);

