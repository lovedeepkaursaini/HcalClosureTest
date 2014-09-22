import FWCore.ParameterSet.Config as cms

#calcrespcorrdijets = cms.EDProducer(
calcrespcorrdijets = cms.EDAnalyzer(
    'CalcRespCorrDiJets',
    caloJetCollName     = cms.string('ak5CaloJets'),
    caloJetCorrName     = cms.string('ak5CaloL2L3'),
    pfJetCollName       = cms.string('ak5PFJets'),
    pfJetCorrName       = cms.string('ak5PFL2L3'),
    genJetCollName      = cms.string('ak5GenJets'),
    genParticleCollName = cms.string('genParticles'),
    genEventInfoName    = cms.string('generator'),
    hbheRecHitName      = cms.string('hbhereco'),
    hfRecHitName        = cms.string('hfreco'),
    hoRecHitName        = cms.string('horeco'),
    rootHistFilename    = cms.string('dijettree.root'),
    maxDeltaEta         = cms.double(1.5),
    minTagJetEta        = cms.double(0.0),
    maxTagJetEta        = cms.double(5.0),
    minSumJetEt         = cms.double(10.), #10.
    minJetEt            = cms.double(5.0), #5.0
    maxThirdJetEt       = cms.double(100.), #100.
    maxJetEMF           = cms.double(0.9),
    doCaloJets          = cms.bool(True),
    doPFJets            = cms.bool(True),
    doGenJets           = cms.bool(True),
    debug               = cms.untracked.bool(False)
    )
