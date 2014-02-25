import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
              
#66 PU files: weighted eta, phi
      'file:crab/crab_0_1211104_66PU/res/E10_19_1_Xfn.root',
      'file:crab/crab_0_1211104_66PU/res/E10_22_1_3zv.root',
      'file:crab/crab_0_1211104_66PU/res/E10_24_1_hGE.root',
      'file:crab/crab_0_1211104_66PU/res/E10_25_1_uiX.root',
          ),
   skipEvents = cms.untracked.uint32(0)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Tree.root')
)
#
#process.Timing=cms.Service("Timing")
#
#process.MessageLogger = cms.Service("MessageLogger",
#                                    qemloosenewlog = cms.untracked.PSet(
#    threshold = cms.untracked.string('INFO')
#    ),
#                                    destinations = cms.untracked.vstring('qemloosenewlog')
#                                    )
#


 
process.demo = cms.EDAnalyzer('Marte',
                                   #RecoVertices = cms.InputTag("offlinePrimaryVertices"),
                                   FilteredCircle8= cms.InputTag("calibJetProducer"),
                                   extra = cms.VInputTag("l1extraParticles:Central",
                                                         "l1extraParticles:Tau"
                                                         ),
                                   CaloJets = cms.InputTag("cAlCAloProducer"),
                                   RhoCaloJets = cms.InputTag("ak5CaloJets:rho"),
                                   RhoTowerJets = cms.InputTag("calibJetProducer:Rho"),
)           
            
            
process.dump=cms.EDAnalyzer('EventContentAnalyzer')
            
process.p = cms.Path(process.demo)
                             
                                                                                                
