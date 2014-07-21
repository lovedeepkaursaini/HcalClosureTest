import FWCore.ParameterSet.Config as cms

process = cms.Process('ANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')

#load the response corrections calculator
process.load('HcalClosureTest.Analyzers.calcrespcorrdijets_cfi')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
#process.load('JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff')

# run over files

process.calcrespcorrdijets.rootHistFilename = cms.string('plots.root')

###import FWCore.Python.FileUtils as FileUtils
#readFiles = cms.untracked.vstring( FileUtils.loadListFromFile ('QCD_Pt30.list') )
#process.source = cms.Source ("PoolSource",fileNames = readFiles)
process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
    #'/store/mc/Summer12/QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v1/0000/009E226E-F499-E111-B2BD-0025904B12FC.root'
    #'file:SinglePiPt100_cfi_py_GEN_FASTSIM_HLT_VALIDATION.root'
    'file:/eos/uscms/store/user/maravin/HCAL/singlePiPt50_2_AODSIM.root'
    )
                            )
#print readFiles

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(1000)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# timing
#process.Timing = cms.Service('Timing')

process.p = cms.Path(process.calcrespcorrdijets)
