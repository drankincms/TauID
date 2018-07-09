import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (1000 ) )

process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(10000),
)

#load in source file

readFiles = cms.untracked.vstring()

process.source = cms.Source( "PoolSource",
                            fileNames = readFiles,
                            inputCommands = cms.untracked.vstring('keep *'),
                            )

readFiles.extend( [
'/store/mc/RunIIFall17MiniAOD/WToTauNu_M-200_TuneCP5_13TeV-pythia8-tauola/MINIAODSIM/94X_mc2017_realistic_v11_ext1-v1/90000/C8FFB21A-2E34-E811-975F-FA163ECF9D5A.root'
#'/store/mc/RunIIFall17MiniAOD/BulkGravToWW_narrow_M-4500_13TeV-madgraph/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/1C2EFCEA-DD60-E811-B008-FA163E7BB6F5.root'
#'/store/mc/RunIIFall17MiniAOD/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/F67C7CED-490D-E811-9676-FA163E103CB2.root'
] );

#run the producer
process.GenNtupler = cms.EDAnalyzer("GenNtupler",
			generator = cms.InputTag("prunedGenParticles"),
			genjets = cms.InputTag("slimmedGenJets"),
                        minGenParticlePt = cms.double(0.),
                        minGenJetPt = cms.double(10.),
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('GenNtuple.root')
                                   )

process.out = cms.OutputModule("PoolOutputModule",
                                       fileName = cms.untracked.string('Output.root') ,
                                       outputCommands = cms.untracked.vstring("drop *"),
                                       )

process.p = cms.Path(process.GenNtupler)

#process.e = cms.EndPath(process.out)
