import FWCore.ParameterSet.Config as cms

process = cms.Process("RCTofflineTEST")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'root://eoscms//eos/cms/store/data/Run2012C/MinimumBias1/RAW/v1/000/198/609/A2E65187-72CA-E111-87D8-003048D374F2.root', 
    'root://eoscms//eos/cms/store/data/Run2012C/MinimumBias1/RAW/v1/000/198/609/F0A0334F-75CA-E111-9833-003048D2C0F4.root'
    #'file:/data/pioppi/work/SLHC/CMSSW_5_3_1/src/66PU/AE31726D-74CA-E111-8A26-5404A6388692.root',
    #'file:/data/pioppi/work/SLHC/CMSSW_5_3_1/src/66PU/B0E72E71-74CA-E111-88F4-003048F11CF0.root',
    #'file:/data/pioppi/work/SLHC/CMSSW_5_3_1/src/66PU/B04B461E-74CA-E111-947A-0025901D6268.root',
    #'file:/data/pioppi/work/SLHC/CMSSW_5_3_1/src/66PU/B06FBDB6-72CA-E111-95AA-0025B320384C.root',
    ),
   #skipEvents = cms.untracked.uint32(4500)
)


process.o1 = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('E10.root'),
    outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                           'keep *_calibJetProducer_*_*',
                                           'keep *_SLHCL1ExtraParticles_*_*',
                                           'keep *_cAlCAloProducer_*_*',
                                           'keep *_ak5CaloJets_*_*',
                                           'keep *_l1extraParticles_*_*',
                                           'keep *_offlinePrimaryVertices_*_*',     
                                           'keep *_L1TowerJetFilter2D_*_*'                         
),
)

#process.TFileService=cms.Service("TFileService",
    #fileName=cms.string('E11.root')
#)

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_P_V32::All'
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
#process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load("Configuration.StandardSequences.Reconstruction_Data_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)


process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

#process.load("RecoParticleFlow.PFProducer.particleFlow_cff") ##added particle flow cff
 
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTriggerAnalysisOnData_cfi")

process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi")

process.calibJetProducer = cms.EDProducer('L1TowerJetFilter2D',
        src = cms.InputTag("L1TowerJetFilter1D"),
	ComparisonDirection = cms.string("phi"), # "eta" or "phi"
	NumOfOutputJets = cms.uint32(12)
    )

process.cAlCAloProducer = cms.EDProducer('CalCaloProducer',
                                     CaloJets = cms.InputTag("ak5CaloJets"),
                                   JetCorrector = cms.string('ak5CaloL1FastL2L3Residual'), #jet corrections residual for data
)


process.demo = cms.EDAnalyzer('Marte',
                                   FilteredCircle8= cms.InputTag("calibJetProducer:CalibJets"),
                                   extra = cms.VInputTag("l1extraParticles:Central",
                                                         "l1extraParticles:Tau"
                                                         ),
                              CaloJets = cms.InputTag("cAlCAloProducer"),
                              RhoCaloJets = cms.InputTag("ak5CaloJets:rho"),
                              RhoTowerJets = cms.InputTag("calibJetProducer:Rho")
)

process.hltHighLevel = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring(#'AlCa_LumiPixels_v8',
                           'HLT_ZeroBias*'
                           ),           # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),
    throw = cms.bool(True)    # throw exception on unknown path names
)




process.p1 = cms.Path(
                      process.hltHighLevel+
                      process.RawToDigi+
                      process.SLHCCaloTrigger+
                      process.l1extraParticles+
    process.trackerlocalreco+
    process.calolocalreco+
#                      process.localreco+
#                      process.globalreco+
process.offlineBeamSpot+
                          process.recopixelvertexing+
                          process.trackingGlobalReco+
                          process.hcalGlobalRecoSequence+
#                          process.particleFlowCluster+
                          process.ecalClusters+
                          process.caloTowersRec+
                          process.vertexreco+
                          process.jetGlobalReco+
                      process.recoJetAssociations+
                      process.cAlCAloProducer+
                      process.calibJetProducer
                      #process.demo
)




process.outpath = cms.EndPath(process.o1)


process.ak5CaloJets.doRhoFastjet = cms.bool(True)
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


                                                     

