import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8PtGun",
                         PGunParameters = cms.PSet(
        MaxPt = cms.double(1000.0),
        MinPt = cms.double(0.7),
        ParticleID = cms.vint32(-13),
        AddAntiParticle = cms.bool(False),
        MaxEta = cms.double(3.0),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-3.0),
        MinPhi = cms.double(-3.14159265359) ## in radians
        ),
                         Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
                         psethack = cms.string('single muminus pt 0.7 to 1000'),
                         firstRun = cms.untracked.uint32(1),
                         PythiaParameters = cms.PSet(parameterSets = cms.vstring())
                         )
