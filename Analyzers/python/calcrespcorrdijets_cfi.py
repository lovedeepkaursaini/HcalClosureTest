import FWCore.ParameterSet.Config as cms

#calcrespcorrdijets = cms.EDProducer(
calcrespcorrdijets = cms.EDAnalyzer(
    'CalcRespCorrDiJets',
#    jetCollName      = cms.string('sisCone5CaloJets'),
    jetCollName      = cms.string('ak5CaloJets'),
#    jetCollName      = cms.string('ak5PFJets'),
#    jetCorrName      = cms.string('L2L3JetCorrectorSC5Calo'),
    jetCorrName      = cms.string('ak5CaloJetsL2L3'),
    genJetCollName   = cms.string('ak5GenJets'),
    rootHistFilename = cms.string('dijettree.root'),
    maxDeltaEta      = cms.double(1.5),
    minTagJetEta     = cms.double(0.0),
    maxTagJetEta     = cms.double(5.0),
    minSumJetEt      = cms.double(10.),
    minJetEt         = cms.double(5.0),
    maxThirdJetEt    = cms.double(100.),
    maxJetEMF        = cms.double(0.9),
    debug            = cms.untracked.bool(False)
    )
