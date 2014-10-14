import FWCore.ParameterSet.Config as cms
process = cms.Process('ANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag=autoCond['startup']

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
#fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/G_Pt-0to15_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v1/00000/002EC875-901A-E211-B893-0030487D5EB5.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(1000)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# timing
#process.Timing = cms.Service('Timing')

process.p = cms.Path(process.calcrespcorrphotonplusjet)
