import FWCore.ParameterSet.Config as cms
process = cms.Process('ANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag=autoCond['startup']

#process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.kt6PFJets = process.kt6PFJets.clone(rParam = 0.6, doRhoFastjet = True)          

#load the response corrections calculator
process.load('HcalClosureTest.Analyzers.calcrespcorrphotonplusjet_cfi')

process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')

# run over files
process.calcrespcorrphotonplusjet.rootHistFilename = cms.string('PhoJet_tree.root')
process.calcrespcorrphotonplusjet.doCaloJets = cms.bool(False)
process.calcrespcorrphotonplusjet.doPFJets = cms.bool(True)
process.calcrespcorrphotonplusjet.doGenJets = cms.bool(True)

##process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/relval/CMSSW_5_3_16/RelValPyquen_GammaJet_pt20_2760GeV/GEN-SIM-RECO/PU_STARTHI53_LV1_mar03-v2/00000/20FE26F4-65A3-E311-B12C-0025904C6378.root'))

process.source = cms.Source("PoolSource", 
fileNames = cms.untracked.vstring('file:/uscms/home/lovedeep/eos/RelValPhotonJets_Pt_10_CMSSW_5_3_12_patch2_A4609359-9E2B-E311-B331-0025905964A6.root')

##fileNames = cms.untracked.vstring(
##    '/store/mc/Summer12_DR53X/G_Pt-170to300_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v1/00000/5846302F-1A18-E211-A060-00266CF2AE10.root',
##    '/store/mc/Summer12_DR53X/G_Pt-170to300_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v1/00000/586126E2-0F18-E211-9323-0030487D864B.root',
##    '/store/mc/Summer12_DR53X/G_Pt-170to300_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v1/00000/A80FB82E-1018-E211-B444-0025904B130E.root',
##    '/store/mc/Summer12_DR53X/G_Pt-170to300_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v1/00000/5809400A-F917-E211-8D3D-0030487F1C51.root',
##    '/store/mc/Summer12_DR53X/G_Pt-170to300_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v1/00000/A40C5492-F917-E211-AB13-002481E0DC82.root'
##    )

)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(1000)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# timing
#process.Timing = cms.Service('Timing')

process.p = cms.Path(process.calcrespcorrphotonplusjet)