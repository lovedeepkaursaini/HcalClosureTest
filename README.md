HcalClosureTest
===============

HCAL Calibration

From https://github.com/johnpaulchou/usercode/tree/master/HcalClosureTest

Running Analyzers/test/testRespCorrDiJets_cfg.py creates a skim tree
Running DataFormat/test/runCaloJetCorr.C will take that tree and produce the corrections

To shut off CaloJets or PFJets, include in the config file:
   process.calcrespcorrdijets.doCaloJets = cms.bool(False)
or
   process.calcrespcorrdijets.doPFJets = cms.bool(False)
