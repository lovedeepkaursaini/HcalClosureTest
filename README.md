HcalClosureTest
===============

Code to calibrate CMS's HCAL in eta using dijet balance. Created by J. P. Chou (Brown) https://github.com/johnpaulchou/usercode/tree/master/HcalClosureTest and updated by David G. Sheffield (Rutgers).

# Creating the Tree

Running Analyzers/test/testRespCorrDiJets_cfg.py creates a tree containing the information needed for dijet balancing. It is set by default to run over CaloJets, PFJets, and GenJets. Any of these can be shut off by adding the lines

```
process.calcrespcorrdijets.doCaloJets = cms.bool(False)
process.calcrespcorrdijets.doPFJets   = cms.bool(False)
process.calcrespcorrdijets.doGenJets  = cms.bool(False)
```

When using PFJets, the dataset being run over must be RECO as AOD does not save the necessary PF block information. When using GenJets, the dataset must contain SIM information.

# Getting response corrections

The macro DataFormat/test/runCaloJetCorr.C takes the tree and determines the response corrections.
