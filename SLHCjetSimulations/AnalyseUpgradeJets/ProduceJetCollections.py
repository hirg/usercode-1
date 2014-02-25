import FWCore.ParameterSet.Config as cms

process = cms.Process("RCTofflineTEST")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


## Create a TFileService to allow modules to write output to a common TFile
#process.TFileService = cms.Service("TFileService",
#                                                                          fileName = cms.string("output.root")
#                                                                          )


process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring(
		#run locally at IC
    'root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/data/Run2012C/ZeroBias/RECO/PromptReco-v1/000/198/609/0E478DDC-26CD-E111-90D2-001D09F2B30B.root',
   ),

   # RAW files
   secondaryFileNames = cms.untracked.vstring(
   	#run locally at IC
    'root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/data/Run2012C/ZeroBias/RAW/v1/000/198/609/F0C83922-74CA-E111-909B-003048673374.root',
			# 'root://eoscms//eos/cms///	
   )
)



process.o1 = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('JetCollections.root'),
    #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('trigger_')),
    outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                           'keep *_SLHCL1ExtraParticles_*_*',
                                           'keep *_PUsubAK5CaloJetProducer_*_*',
                                           'keep *_L1CaloTowerProducer_*_*',
                                           'keep *_L1RingSubtractionProducer_*_*',
                                           'keep *_L1TowerJetPUEstimator_*_*',
                                           'keep *_L1TowerJetPUSubtractedProducer_*_*',
                                           'keep *_L1CalibFilterTowerJetProducer_*_*',
                                           'keep *_ak5CaloJets_*_RECO',
                                           'drop *_ak5CaloJets_rho_RECO',
                                           'keep *_ak5CaloJets_rho_RCTofflineTEST',
                                           'keep *_ak5PFJets_*_RECO',
                                           'drop *_ak5PFJets_rho_RECO',
                                           'keep *_ak5PFJets_rho_RCTofflineTEST',
 					   'keep *_l1extraParticles_*_*',
                                           'keep *_offlinePrimaryVertices_*_*',
                                           'keep *_ak5JetID_*_*',
),
)


process.maxEvents = cms.untracked.PSet(
    # restrict number of events to 10
    input = cms.untracked.int32(1000)
    # run over all events
#    input = cms.untracked.int32(-1)
)

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'POSTLS161_V12::All'

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
process.load("RecoParticleFlow.PFProducer.particleFlow_cff") 
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")
process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi")
process.load("Configuration.StandardSequences.Services_cff")



#create collection of PU corrected AK5 Calo Jets
process.PUsubAK5CaloJetProducer = cms.EDProducer('CalCaloProducer',
                                   CaloJets = cms.InputTag("ak5CaloJets"),
                                   JetCorrector = cms.string('ak5CaloL1FastL2L3Residual'), #jet corrections residual for data
)

from RecoJets.JetProducers.JetIDParams_cfi import *

process.ak5JetID = cms.EDProducer('JetIDProducer', JetIDParams,
        src = cms.InputTag('ak5CaloJets')

)

process.ak5PFJets.doRhoFastjet = cms.bool(True)

process.ak5CaloJets.doRhoFastjet = cms.bool(True)

process.p1 = cms.Path(
                      process.RawToDigi+
                      process.SLHCCaloTrigger+
                      process.ak5CaloJets+
                      process.ak5PFJets+
                      process.PUsubAK5CaloJetProducer+
                      process.ak5JetID  
)
process.outpath = cms.EndPath(process.o1)

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )
    
