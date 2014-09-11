import FWCore.ParameterSet.Config as cms

process = cms.Process('ANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
##from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.GlobalTag.globaltag = 'START53_V19D::All' #'GR_R_52_V9::All' # not sure what tag to use
process.GlobalTag.globaltag = 'START53_V7A::All'
##process.GlobalTag.globaltag=autoCond['startup']
#load the response corrections calculator
process.load('HcalClosureTest.Analyzers.calcrespcorrdijets_cfi')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')
#process.load('JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff')

#from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *

# run over files

#process.calcrespcorrdijets.rootHistFilename = cms.string('test.root')
process.calcrespcorrdijets.rootHistFilename = cms.string('Trees/QCD_Pt-15to3000_0030487D5E5F_severity.root')
#process.calcrespcorrdijets.rootHistFilename = cms.string('Trees/Pion_Pt-50.root')
#process.calcrespcorrdijets.rootHistFilename = cms.string('Trees/test_QCD_Pt-120to170_TuneZ2star.root')
#process.calcrespcorrdijets.rootHistFilename = cms.string('Trees/tmpQCD.root')
#process.calcrespcorrdijets.doGenJets = cms.bool(False)
#process.calcrespcorrdijets.debug = cms.untracked.bool(True)

#import FWCore.Python.FileUtils as FileUtils
#readFiles = cms.untracked.vstring( FileUtils.loadListFromFile ('Pion_Pt-50.list') )
#process.source = cms.Source ("PoolSource",fileNames = readFiles)
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'file:/eos/uscms/store/user/maravin/HCAL/singlePiPt50_1_AODSIM.root'
    'file:001288D3-2AE3-E111-A943-0030487D5E5F.root'
    #'root://cmsxrootd-site.fnal.gov//store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/001288D3-2AE3-E111-A943-0030487D5E5F.root'
    #'/store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v2/0000/002D1EF6-33E3-E111-96E1-0030487E54B7.root'
    #'/store/data/Run2012A/MultiJet/AOD/22Jan2013-v1/20000/0036C47E-0B74-E211-B992-00266CF32684.root'
    #'file:/uscmst1b_scratch/lpc1/3DayLifetime/dgsheffi/0022FFCF-2F19-E211-AD61-00266CFFA68C.root'
    #'root://cmsxrootd-site.fnal.gov//store/mc/Summer12_DR53X/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v1/00000/0022FFCF-2F19-E211-AD61-00266CFFA68C.root'
    #'root://cmsxrootd-site.fnal.gov//store/mc/Summer12_DR53X/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6/GEN-SIM-RECO/PU_S10_START53_V7A-v3/00000/0292D61E-180C-E211-8027-003048F0E3D2.root'
    )
                            )
#print readFiles

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(100)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.load("RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff")
process.p = cms.Path(
    #process.particleFlowCluster+
	process.calcrespcorrdijets
    )
# timing
#process.Timing = cms.Service('Timing')

#process.p2 = cms.Path(process.calcrespcorrdijets)

