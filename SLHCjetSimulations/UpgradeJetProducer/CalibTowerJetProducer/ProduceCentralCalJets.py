import FWCore.ParameterSet.Config as cms

process = cms.Process("RCTofflineTEST")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://eoscms//eos/cms//store/data/Run2012D/ZeroBias25ns1/RAW/v1/000/208/931/20A5CC56-EE44-E211-87AE-5404A63886AA.root'
    ),
   #skipEvents = cms.untracked.uint32(4500)
)

    
            
process.o1 = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('testJetCollections.root'),
    #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('trigger_')),
    outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                           'keep *_SLHCL1ExtraParticles_*_*',
                                           'keep *_calibTowerJetProducer_*_*',
                                           'keep *_l1extraParticles_*_*',
                                           'keep *_L1TowerJetProducer_*_*', 
                                           'keep *_L1TowerJetPUcalc_*_*',
                                           'keep *_L1TowerJetFilter2D_*_*',
                                           'keep *_L1TowerFwdJetFilter2D_*_*',
                                           'keep *_L1TowerJetPU*_*_*' 
),
)

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'GR_R_53_V2::All'
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')

#Dont need these
#process.load("Configuration.StandardSequences.Services_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")    
#process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
    

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTriggerAnalysisOnData_cfi")

    
process.calibTowerJetProducer = cms.EDProducer("CalibTowerJetProducer",
    inRhodata_file = cms.FileInPath('ProduceJetCollections/CalibTowerJetProducer/data/rho_lookup.txt'),
    inMVA_weights_file = cms.FileInPath('ProduceJetCollections/CalibTowerJetProducer/data/TMVARegression_BDT.weights.xml'),
    inPtdata_file = cms.FileInPath('ProduceJetCollections/CalibTowerJetProducer/data/looselogval.log'),
    FilteredCircle8 = cms.InputTag("L1TowerJetFilter2D"),
    FilteredFwdCircle8 = cms.InputTag("L1TowerFwdJetFilter2D")
)
process.p1 = cms.Path(
                      process.RawToDigi+
                      process.SLHCCaloTrigger+
                      process.calibTowerJetProducer
)


process.outpath = cms.EndPath(process.o1)


# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

