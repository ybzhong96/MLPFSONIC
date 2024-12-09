#include <ostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "DataFormats/Math/interface/Vector3D.h"

#include "IOMC/ParticleGuns/interface/FlatEtaRangeGunProducer.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

edm::FlatEtaRangeGunProducer::FlatEtaRangeGunProducer(const edm::ParameterSet& params)
    : FlatRandomEGunProducer(params),
      nParticles_(params.getParameter<ParameterSet>("PGunParameters").getParameter<int>("nParticles")),
      exactShoot_(params.getParameter<ParameterSet>("PGunParameters").getParameter<bool>("exactShoot")),
      randomShoot_(params.getParameter<ParameterSet>("PGunParameters").getParameter<bool>("randomShoot")),
      minDr_(params.getParameter<ParameterSet>("PGunParameters").getUntrackedParameter<double>("minDr", -1.)),
      debug_(params.getUntrackedParameter<bool>("debug")) {
}

void edm::FlatEtaRangeGunProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &(rng->getEngine(event.streamID()));

  if (debug_) {
    LogDebug("FlatEtaRangeGunProducer") << " : Begin New Event Generation" << std::endl;
  }

  // create a new event to fill
  auto* genEvent = new HepMC::GenEvent();

  // determine the number of particles to shoot
  int n = 0;
  if (exactShoot_) {
    n = (int)fPartIDs.size();
    std::cout << "exactShoot n=" << nParticles_ << std::endl;
  } else if (randomShoot_) {
    n = CLHEP::RandFlat::shoot(engine, 1, nParticles_ + 1);
    std::cout << "randomShoot n=" << nParticles_ << std::endl;
  } else {
    n = nParticles_;
    std::cout << "else n=" << nParticles_ << std::endl;
  }

  std::vector<math::XYZVector > previousp4;
  int particle_counter = 0;
  // shoot particles
  for (int i = 0; i < 2 * n; i++) {  //n for positive and n for negative eta
    // create a random deltaR

    // obtain kinematics
    int id = fPartIDs[exactShoot_ ? particle_counter : CLHEP::RandFlat::shoot(engine, 0, fPartIDs.size())];
    std::cout << "generating id=" << id << std::endl;
    particle_counter++;
    if (particle_counter >= n)
      particle_counter = 0;

    const HepPDT::ParticleData* pData = fPDGTable->particle(HepPDT::ParticleID(abs(id)));
    double eta = CLHEP::RandFlat::shoot(engine, fMinEta, fMaxEta);
    if (i < n)
      eta *= -1;
    double phi = CLHEP::RandFlat::shoot(engine, fMinPhi, fMaxPhi);



    double e = CLHEP::RandFlat::shoot(engine, fMinE, fMaxE);
    double m = pData->mass().value();
    double p = sqrt(e * e - m * m);
    math::XYZVector pVec = p * math::XYZVector(cos(phi), sin(phi), sinh(eta)).unit();


    //check
    if(minDr_>0){
        bool isgood=true;
        for(const auto ppvec: previousp4){
            double drsq = reco::deltaR2(pVec,ppvec);
            if(drsq<minDr_){
                isgood=false;
                break;
            }
        }
        if(!isgood)
            continue;
    }

    previousp4.push_back(pVec);

    HepMC::GenVertex* vtx = new HepMC::GenVertex(HepMC::FourVector(0, 0, 0, 0));

    // create the GenParticle
    HepMC::FourVector fVec(pVec.x(), pVec.y(), pVec.z(), e);
    HepMC::GenParticle* particle = new HepMC::GenParticle(fVec, id, 1);
    particle->suggest_barcode(i + 1);

    // add the particle to the vertex and the vertex to the event
    vtx->add_particle_out(particle);
    genEvent->add_vertex(vtx);

    if (debug_) {
      vtx->print();
      particle->print();
    }
  }

  // fill event attributes
  genEvent->set_event_number(event.id().event());
  genEvent->set_signal_process_id(20);

  if (debug_) {
    genEvent->print();
  }

  // store outputs
  std::unique_ptr<HepMCProduct> BProduct(new HepMCProduct());
  BProduct->addHepMCData(genEvent);
  event.put(std::move(BProduct), "unsmeared");
  auto genEventInfo = std::make_unique<GenEventInfoProduct>(genEvent);
  event.put(std::move(genEventInfo));

  if (debug_) {
    LogDebug("FlatEtaRangeGunProducer") << " : Event Generation Done " << std::endl;
  }
}
