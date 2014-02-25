import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        ##'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_20_1_hmg.root ',
        ##'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_22_1_Jof.root ',
        ##'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_23_1_7wR.root ',
        ##'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_25_1_aix.root ',
        ##'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_26_1_IBZ.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_27_1_Pay.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_28_1_77j.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_29_1_HAH.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_30_1_IVM.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_32_1_l73.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_33_1_klR.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_34_1_Fvc.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_35_1_UsR.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_36_1_dGh.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_37_1_Kzd.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_3_1_dWa.root  ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_41_1_3yM.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_45_1_1qv.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_47_1_0sT.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_4_1_pSd.root  ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_51_1_S8S.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_52_1_lZQ.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_53_1_keJ.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_54_1_frL.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_57_1_P1L.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_58_1_6p9.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_59_1_J7k.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_5_1_EFG.root  ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_62_1_7ba.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_63_1_VLz.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_64_1_fZt.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_65_1_ido.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_66_1_9BQ.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_67_1_Iz4.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_68_1_JZA.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_69_1_Bzi.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_6_1_Khx.root  ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_70_1_Dor.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_71_1_iUs.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_73_1_vom.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_74_1_Oft.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_75_1_VFs.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_77_1_aB6.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_79_1_0Wh.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_7_1_2wL.root  ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_80_1_G0B.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_81_1_kqN.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_83_1_kuP.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_84_1_KTj.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_86_1_YYa.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_88_1_oZB.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_89_1_u7Q.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_8_1_rgO.root  ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_91_1_7C7.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_92_1_0tA.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_93_1_VOh.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_95_1_GWS.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_96_1_sQ3.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_97_1_bjw.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_99_1_ahD.root ',
        'file:/afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/SinMuRootFiles/JetCollections_SinMu/JetCollections_9_1_8FZ.root  ',
            
                
#these are the HIGH PU FILES
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_10_1_hjn.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_11_1_UPA.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_12_1_YMZ.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_13_1_sOc.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_14_1_9Mj.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_15_1_CcO.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_16_1_c8w.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_17_1_ZZu.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_18_1_1a8.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_19_1_Jw6.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_1_1_C0u.root ',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_20_1_WTN.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_21_1_9Y6.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_22_1_0tQ.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_23_1_OK2.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_24_1_6S2.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_25_1_ydb.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_26_1_REY.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_27_0_xFv.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_28_0_vaG.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_29_0_IKI.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_2_1_XeV.root ',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_30_0_pio.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_31_0_FTA.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_32_0_KGL.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_33_0_L2o.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_34_0_AIR.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_35_0_Uz0.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_37_0_ZzE.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_38_0_Vsp.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_39_0_C2f.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_3_1_5oc.root ',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_40_0_d50.root',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_4_1_8ZI.root ',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_5_1_j1Q.root ',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_6_1_ZAv.root ',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_7_1_Y45.root ',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_8_1_1e2.root ',
#'file://afs/cern.ch/work/r/rlucas/L1TriggerUpgrades/cvsTesting/ZBRootFilesHighPU/JetCollections_9_1_a4K.root ',


   ),
   #skipEvents = cms.untracked.uint32(160000)

)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('sinMu_30.root')
)

process.demo = cms.EDAnalyzer('CalibrateJets',
                                   FilteredCircle8= cms.InputTag("L1TowerJetFilter2D"),
                                   CalibCircle8= cms.InputTag("calibTowerJetProducer:CalibCenJets"),
                                   extrajet = cms.VInputTag("l1extraParticles:Central",
                                                         "l1extraParticles:Tau",
                                            "l1extraParticles:Forward",
                                                            ),
                                   extraem =  cms.VInputTag("l1extraParticles:Isolated",
                                                         "l1extraParticles:NonIsolated",
                                                         ),
                                   CaloJets = cms.InputTag("PUsubAK5CaloJetProducer"),#cAlCAloProducer"),     #PUsubAK5CaloJetProducer
                                   RhoCaloJets = cms.InputTag("ak5CaloJets:rho"),   
                                   RhoTowerJets = cms.InputTag("calibTowerJetProducer:Rho"),
                                   jetIDHelperConfig = cms.InputTag("ak5JetID"),
                                   L1extraHT = cms.InputTag("l1extraParticles:MHT"),

                                   inRhodata_file = cms.FileInPath(
                                       "ProduceJetCollections/CalibTowerJetProducer/data/rho_lookup.txt"),
																			# on cvs this file is here: /UserCode/rlucas/SLHCjetSimulations/UpgradeJetProducer/CalibTowerJetProducer/data/rho_lookup.txt
                                   inPtdata_file = cms.FileInPath(
                                       "CalibrateJets/CalibrateJets/data/loglookup.log"),
																			# on cvs this file is here:/UserCode/rlucas/SLHCjetSimulations/UpgradeJetSimulations/Marte/data/loglookup.log
                                   inMVA_weights_file =  cms.FileInPath(
                                       "CalibrateJets/weights/TMVARegression_BDT.weights.xml"),
																			#on cvs this file is here: /UserCode/rlucas/SLHCjetSimulations/UpgradeJetSimulations/weights/CentralJets/TMVARegression_BDT.weights.xml

)                                          


process.p = cms.Path(process.demo)   
