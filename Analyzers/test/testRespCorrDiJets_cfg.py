import FWCore.ParameterSet.Config as cms

process = cms.Process('ANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag=autoCond['startup']

#load the response corrections calculator
process.load('HcalClosureTest.Analyzers.calcrespcorrdijets_cfi')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')

# run over files

process.calcrespcorrdijets.rootHistFilename = cms.string('Trees/QCD_Pt-15to3000_0030487D5E5F_recHitFractions')
#process.calcrespcorrdijets.doCaloJets = cms.bool(False)
#process.calcrespcorrdijets.doPFJets = cms.bool(False)
#process.calcrespcorrdijets.doGenJets = cms.bool(False)
#process.calcrespcorrdijets.debug = cms.untracked.bool(True)

#import FWCore.Python.FileUtils as FileUtils
#readFiles = cms.untracked.vstring( FileUtils.loadListFromFile ('Pion_Pt-50.list') )
#process.source = cms.Source ("PoolSource",fileNames = readFiles)
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'root://cmsxrootd-site.fnal.gov//store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/001288D3-2AE3-E111-A943-0030487D5E5F.root'
    'file:001288D3-2AE3-E111-A943-0030487D5E5F.root'
    )
                            )
#print readFiles

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(100)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# timing
#process.Timing = cms.Service('Timing')

process.p = cms.Path(process.calcrespcorrdijets)
