#ifndef FlatRandomEGunProducer_H
#define FlatRandomEGunProducer_H

/** \class FlatRandomEGunProducer
 *
 * Generates single particle gun in HepMC format
 * Julia Yarba 10/2005 
 ***************************************/

#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

namespace edm {

  class FlatRandomEGunProducer : public BaseFlatGunProducer {
  public:
    FlatRandomEGunProducer(const ParameterSet& pset);
    ~FlatRandomEGunProducer() override;

    void produce(Event& e, const EventSetup& es) override;

  protected:
    // data members

    double fMinE;
    double fMaxE;
  };
}  // namespace edm

#endif
