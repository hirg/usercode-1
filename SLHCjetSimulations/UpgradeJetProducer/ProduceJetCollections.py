import FWCore.ParameterSet.Config as cms

process = cms.Process("RCTofflineTEST")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
 
process.source = cms.Source("PoolSource",    
fileNames = cms.untracked.vstring(
  'root://eoscms//eos/cms/store/data/Run2012D/ZeroBias25ns1/RECO/PromptReco-v1/000/209/148/6CAC87B9-5549-E211-AECA-BCAEC518FF62.root',
  #'root://eoscms//eos/cms/store/data/Run2012D/ZeroBias25ns1/RECO/PromptReco-v1/000/209/148/42B9854A-4D49-E211-BE25-003048D37580.root',          
  #'root://eoscms//eos/cms/store/data/Run2012D/ZeroBias25ns1/RECO/PromptReco-v1/000/209/148/E670C73E-4D49-E211-9293-001D09F295A1.root',  
  #'root://eoscms//eos/cms/store/data/Run2012D/ZeroBias25ns1/ALCARECO/TkAlMinBias-v1/000/209/148/FE509FCD-5749-E211-B268-003048D2BC30.root',              
  ),

#   # RAW files
  #fileNames = cms.untracked.vstring(# 
    secondaryFileNames = cms.untracked.vstring(
 'root://eoscms//eos/cms//store/data/Run2012D/ZeroBias25ns1/RAW/v1/000/209/148/B2783F9D-9F47-E211-B15B-0025901D625A.root',
 'root://eoscms//eos/cms//store/data/Run2012D/ZeroBias25ns1/RAW/v1/000/209/148/BCCC3530-9F47-E211-BBD3-003048F118E0.root' ,
 'root://eoscms//eos/cms//store/data/Run2012D/ZeroBias25ns1/RAW/v1/000/209/148/AA7A205F-A347-E211-AEBC-003048F1C832.root' ,
 'root://eoscms//eos/cms//store/data/Run2012D/ZeroBias25ns1/RAW/v1/000/209/148/AC01CFA3-9A47-E211-8F92-5404A63886B0.root' ,
 'root://eoscms//eos/cms//store/data/Run2012D/ZeroBias25ns1/RAW/v1/000/209/148/AEC2097D-9B47-E211-B5D5-001D09F29146.root' , 
   )
)
process.L1TowerJetPUcalcMedian = cms.EDProducer("L1TowerJetPUcalc",
    src = cms.InputTag("L1TowerJetProducer"),
    PUSubtractionType = cms.string("median") # "mean", "median", "mode", "truncated mean"
)
        
            
process.o1 = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('JetCollections.root'),
    #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('trigger_')),
    outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                           'keep *_SLHCL1ExtraParticles_*_*',
                                           'keep *_PUsubAK5CaloJetProducer_*_*',
                                           'keep *_calibTowerJetProducer_*_*',
                                           'keep *_ak5CaloJets_*_RECO',
                                           'drop *_ak5CaloJets_rho_RECO',
                                           'keep *_ak5CaloJets_rho_RCTofflineTEST',
                                           'keep *_ak5PFJets_*_*',
                                           'keep *_l1extraParticles_*_*',
                                           'keep *_offlinePrimaryVertices_*_*',
                                           'keep *_L1TowerJetProducer_*_*', 
                                           'keep *_L1TowerJetPUcalc_*_*',
                                           'keep *_L1TowerJetFilter2D_*_*',
                                           'keep *_L1TowerFwdJetFilter2D_*_*',
                                           'keep *_ak5JetID_*_*',
                                           'keep *_L1TowerJetPU*_*_*' 
),
)

#process.TFileService=cms.Service("TFileService",
    #fileName=cms.string('E11.root')
#)

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'GR_R_53_V2::All'
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')

# UNCOMMENT THIS LINE TO RUN ON SETTINGS FROM THE DATABASE
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource', 'GlobalTag')

process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#unpacking
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.Reconstruction_Data_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

#process.load("RecoParticleFlow.PFProducer.particleFlow_cff") ##added particle flow cff
 
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTriggerAnalysisOnData_cfi")

process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi")


process.PUsubAK5CaloJetProducer = cms.EDProducer('CalCaloProducer',
                                     CaloJets = cms.InputTag("ak5CaloJets"),
                                   JetCorrector = cms.string('ak5CaloL1FastL2L3Residual'), #jet corrections residual for data
)

#process.load("ProduceJetCollections.CalibTowerJetProducer.calibtowerjetproducer_cfi")
    
        
process.calibTowerJetProducer = cms.EDProducer('CalibTowerJetProducer',
                                   FilteredCircle8= cms.InputTag("L1TowerJetFilter2D"),
                                   FilteredFwdCircle8=cms.InputTag("L1TowerFwdJetFilter2D"),
                                   inPtdata_file = cms.FileInPath(
                                       "ProduceJetCollections/CalibTowerJetProducer/data/looselogval.log"),
                                   inRhodata_file = cms.FileInPath(
                                       "ProduceJetCollections/CalibTowerJetProducer/data/rho_lookup.txt"),
                                   inMVA_weights_file = cms.FileInPath(   
                                        "ProduceJetCollections/CalibTowerJetProducer/data/TMVARegression_BDT.weights.xml"),
                                                   
                                   )           


process.hltHighLevel = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring(#'AlCa_LumiPixels_v8',
                           'HLT_ZeroBias_v*'
                           ),           # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),
    throw = cms.bool(True)    # throw exception on unknown path names
)

process.ak5PFJets.doRhoFastjet = cms.bool(True)

process.ak5CaloJets.doRhoFastjet = cms.bool(True)

process.p1 = cms.Path(
                      process.RawToDigi+
                      process.SLHCCaloTrigger+
											process.L1TowerJetPUcalcMedian+

                      process.ak5CaloJets+
                      process.ak5PFJets+
                      process.PUsubAK5CaloJetProducer+  
                      process.calibTowerJetProducer
)




process.outpath = cms.EndPath(process.o1)


# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

#CALO TRIGGER CONFIGURATION OVERRIDE
process.load("L1TriggerConfig.RCTConfigProducers.L1RCTConfig_cff")
process.RCTConfigProducers.eMaxForHoECut = cms.double(60.0)
process.RCTConfigProducers.hOeCut = cms.double(0.05)
process.RCTConfigProducers.eGammaECalScaleFactors = cms.vdouble(1.0, 1.01, 1.02, 1.02, 1.02,
                                                      1.06, 1.04, 1.04, 1.05, 1.09,
                                                      1.1, 1.1, 1.15, 1.2, 1.27,
                                                      1.29, 1.32, 1.52, 1.52, 1.48,
                                                      1.4, 1.32, 1.26, 1.21, 1.17,
                                                      1.15, 1.15, 1.15)
process.RCTConfigProducers.eMinForHoECut = cms.double(3.0)
process.RCTConfigProducers.hActivityCut = cms.double(4.0)
process.RCTConfigProducers.eActivityCut = cms.double(4.0)
process.RCTConfigProducers.jetMETHCalScaleFactors = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0)
process.RCTConfigProducers.eicIsolationThreshold = cms.uint32(6)
process.RCTConfigProducers.etMETLSB = cms.double(0.25)
process.RCTConfigProducers.jetMETECalScaleFactors = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0)
process.RCTConfigProducers.eMinForFGCut = cms.double(100.0)
process.RCTConfigProducers.eGammaLSB = cms.double(0.25)


                                                     


 
