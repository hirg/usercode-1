import FWCore.ParameterSet.Config as cms

process = cms.Process("Analyze")

process.load("FWCore.MessageService.MessageLogger_cfi")

# Stop stupid output after every event
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

    # Robyn's sample file
    #'file:JetCollections.root'



    # ****************************************
    # *   Test file                          *
    # ****************************************

    #  'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_7_1_elP.root'
    #  'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_23_1_h44.root'
    #  'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_1_1_lnP.root'
    'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_1_1_v1T.root'    

    # ****************************************
    # *   Zero bias PU45 dataset - 18 Apr    *
    # ****************************************


    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_7_1_elP.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_19_1_8MX.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_16_1_mZV.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_2_1_9gw.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_8_1_kR3.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_17_1_8Do.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_35_1_Kba.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_37_1_izD.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_15_1_6pi.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_6_1_66L.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_24_1_fbL.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_30_1_98L.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_14_1_RUv.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_39_1_eak.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_4_1_PgQ.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_20_1_ws1.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_9_1_0UK.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_18_1_gsY.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_26_1_SnJ.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_25_1_hfp.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_27_1_kv9.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_29_1_J8f.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_12_1_CQr.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_3_1_2gw.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_21_1_fFg.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_34_1_E4z.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_36_1_R5L.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_10_1_hHd.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_31_1_2mz.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_22_1_xqU.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_13_1_HPe.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_5_1_qN1.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_1_1_hZt.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_23_1_Sr6.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_11_1_TIP.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_38_1_a3Z.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_33_1_mct.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_28_1_gIk.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_40_2_nAX.root',
    #_18_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr/JetCollections_32_2_624.root',



    # **************************************************
    # *   Zero bias PU45 dataset - 18 Apr JetSatFix    *
    # **************************************************

    
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_6_1_E7m.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_21_1_PfJ.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_7_1_5LG.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_17_1_tZI.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_33_1_lIV.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_3_1_ukZ.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_24_1_5g8.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_26_1_83E.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_10_1_b1e.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_32_1_xM0.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_39_1_H1C.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_27_1_ByW.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_12_1_vuE.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_5_1_PAN.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_14_1_iS7.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_34_1_FUO.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_18_1_fuU.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_19_1_znR.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_4_1_ODK.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_30_1_3nD.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_28_1_231.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_1_1_LpD.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_13_1_gXN.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_31_1_2sl.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_9_1_3sH.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_8_1_nn4.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_38_1_Pzg.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_2_1_n6K.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_22_1_dN3.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_11_1_8mX.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_23_1_45A.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_20_1_97f.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_16_1_kuu.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_15_1_aYs.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_25_1_gKB.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_29_1_LEi.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_36_2_KD7.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_37_2_8QL.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_35_2_cWn.root',
    #_18_Sat_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_18Apr_JetSaturationFix/JetCollections_40_2_JGx.root',
 
    # ********************************************************
    # *   Zero bias PU45 dataset - 19 Apr JetSatAndEtaFix    *
    # ********************************************************

#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_35_1_aq4.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_38_1_h8b.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_29_1_5RV.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_10_1_9Gr.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_3_1_hPt.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_8_1_M5B.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_28_1_1WB.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_40_1_vJ6.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_37_1_Q10.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_33_1_IVt.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_23_1_dO9.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_2_1_1ck.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_6_1_rVr.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_39_1_W0E.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_5_1_Yqp.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_30_1_fh4.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_31_1_kMb.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_7_1_73l.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_11_1_wgE.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_36_1_Nui.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_14_1_1xW.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_34_1_ybA.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_9_1_p2p.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_15_1_rai.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_13_1_PU5.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_16_1_4SQ.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_27_1_8MR.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_4_1_Fox.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_1_1_2L5.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_12_1_9Ld.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_26_1_uQQ.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_21_1_cDV.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_22_1_aEX.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_19_1_OHO.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_18_1_VGK.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_24_1_rMt.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_17_1_5Ty.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_25_1_Y2d.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_20_1_wTX.root',
#_19_Sat_Eta#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_19Apr_JetSaturationPlusAreaFix/JetCollections_32_2_WZN.root',





    # ********************************************************
    # *   Zero bias PU45 dataset - 26 Apr MARKIBugFix        *
    # ********************************************************



#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_23_1_h44.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_31_1_5ak.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_37_1_God.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_6_1_h5r.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_18_1_xzI.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_33_2_PkD.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_5_1_RHS.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_19_1_MPT.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_9_1_EtQ.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_8_1_PD0.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_24_1_Fez.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_16_1_L0b.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_39_1_ZSJ.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_15_1_I4p.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_38_1_KMI.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_30_1_9Kk.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_34_1_pkq.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_12_1_zMd.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_32_1_6lg.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_3_1_GE2.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_28_1_Ab6.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_26_1_0QW.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_4_1_0eK.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_21_1_ywo.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_1_1_p8T.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_14_1_E4X.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_27_1_Psk.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_29_1_Us8.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_25_1_gB8.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_40_1_vgQ.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_20_1_uGQ.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_13_1_f3x.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_17_1_KLr.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_22_1_jhe.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_2_1_i93.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_7_1_MFN.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_11_1_FLn.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_10_1_0Gq.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_36_1_QG0.root',
#_26_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Apr_MARKI_BugFixes/JetCollections_35_1_bfC.root',







    # ********************************************************
    # *   Zero bias PU45 dataset - 30 Apr MARKIBugFix        *
    # ********************************************************


#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_29_2_EBE.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_31_2_2Dn.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_36_2_Nlz.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_35_2_ZK5.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_1_1_lnP.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_28_1_X18.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_27_1_ZHl.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_21_1_zem.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_33_1_Bl6.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_32_1_MqV.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_30_1_XNw.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_34_1_JNw.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_38_1_t2y.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_39_1_DK6.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_37_1_IMM.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_19_1_ZOu.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_25_1_v29.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_40_1_dUi.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_12_1_dy3.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_7_1_puk.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_14_1_fgM.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_24_1_385.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_15_1_V5l.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_13_1_7jw.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_18_1_fJE.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_22_1_qbD.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_6_1_5AF.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_4_1_8SU.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_8_1_3xS.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_20_1_8Wp.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_5_1_jAH.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_2_1_HME.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_10_1_fhG.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_17_1_5Ed.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_9_1_elk.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_16_1_4EA.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_11_1_HBF.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_23_1_wBX.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_26_1_TH3.root',
#_30_MARKI#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_30Apr_MARKI_BugFixes/JetCollections_3_1_Oqo.root',





    # **************************************************************************************
    # *                       Zero bias PU45 dataset - 06 May Mk1                          *
    # **************************************************************************************



    # *******************************************************************
    # *                              Calib                              *
    # *******************************************************************

#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_39_2_BOL.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_40_2_Osc.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_1_1_3E2.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_28_1_lQX.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_21_1_etr.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_31_1_1BW.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_37_1_Zd7.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_27_1_Cy3.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_22_1_zSD.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_20_1_zUt.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_32_1_PgI.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_11_1_WL3.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_10_1_T6P.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_26_1_FDk.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_24_1_nBJ.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_33_1_UEl.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_25_1_LKM.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_35_1_sA5.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_29_1_4MK.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_14_1_zQ0.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_19_1_KbL.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_18_1_3o4.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_17_1_4Nx.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_23_1_CUV.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_9_1_FwS.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_13_1_CYI.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_36_1_xyg.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_5_1_efr.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_16_1_vdZ.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_34_1_7hG.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_15_1_myD.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_8_1_bLw.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_12_1_LV2.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_38_1_WUK.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_30_1_D2N.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_7_1_vI4.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_2_1_Los.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_3_1_7SL.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_6_1_oLC.root',
#06Mk1Calib'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_Calib/JetCollections_4_1_NHO.root',


    # *******************************************************************
    # *                            5GeVSeed                             *
    # *******************************************************************





#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_1_1_GMm.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_28_1_Uj9.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_21_1_foG.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_31_1_fWi.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_32_1_qOz.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_27_1_6Ct.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_29_1_cRf.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_39_1_FIf.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_40_1_xpT.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_36_1_XeJ.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_11_1_rUQ.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_18_1_i8q.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_34_1_Ywt.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_37_1_IYf.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_33_1_yIP.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_35_1_Xlg.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_12_1_L1Y.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_25_1_fuf.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_17_1_TfP.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_6_1_rn3.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_14_1_F0e.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_19_1_L4c.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_4_1_r29.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_15_1_bpX.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_13_1_8VB.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_23_1_KBd.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_20_1_yTE.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_8_1_3LX.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_10_1_WjP.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_9_1_evD.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_30_1_uPm.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_2_1_gAq.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_3_1_Nny.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_7_1_lCh.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_5_1_LVF.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_16_1_mbJ.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_38_1_IIR.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_22_1_Hsa.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_24_1_gLr.root',
#06Mk15GevSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_5GeVSeed/JetCollections_26_1_mJ1.root',


    # *******************************************************************
    # *                        RingSub_5GeVSeed                         *
    # *******************************************************************





#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_1_1_XXE.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_28_1_Yrp.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_21_1_amr.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_27_1_obu.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_17_1_DGD.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_22_1_bkO.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_33_1_4Tx.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_34_1_pMj.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_29_1_czo.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_36_1_B29.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_37_1_vuc.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_32_1_lEB.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_38_1_duP.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_31_1_EBo.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_40_1_MgQ.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_35_1_zL2.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_39_1_imc.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_19_1_JBO.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_25_1_zRq.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_20_1_TVV.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_14_1_p0m.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_15_1_sSk.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_23_1_wfc.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_9_1_m0O.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_7_1_lLO.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_18_1_4k9.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_10_1_duY.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_16_1_JcI.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_8_1_71n.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_30_1_0sn.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_5_1_x3i.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_13_1_Fyx.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_6_1_AzM.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_4_1_lhx.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_2_1_Pu8.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_11_1_X0z.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_3_1_I4W.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_24_1_LgT.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_26_1_oDp.root',
#06Mk1RingSub5GeVSeed'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed/JetCollections_12_1_wvu.root',



    # *******************************************************************
    # *                            RingSub                              *
    # *******************************************************************




#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_1_1_v1T.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_27_1_xRP.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_21_1_Q62.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_28_1_GMI.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_25_1_k8Y.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_33_1_C5a.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_18_1_fOu.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_39_1_SAe.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_35_1_1fr.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_37_1_sIK.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_29_1_X3C.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_36_1_lmp.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_32_1_Iot.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_24_1_KtA.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_34_1_E2z.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_38_1_SOu.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_31_1_xdd.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_23_1_7L5.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_17_1_lwv.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_22_1_ahY.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_14_1_3Ix.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_7_1_rqs.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_19_1_IGp.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_16_1_Ngz.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_26_1_tV6.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_20_1_z7m.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_15_1_23j.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_11_1_ASV.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_13_1_ZMI.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_30_1_Fli.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_9_1_tEM.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_8_1_ykB.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_4_1_3ZN.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_5_1_i16.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_40_1_wPe.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_10_1_ubI.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_3_1_Xay.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_6_1_fQX.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_2_1_VrQ.root',
#06Mk1RingSub'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_06May_Mk1_RingSub/JetCollections_12_1_KQi.root',



    # **************************************************************************************
    # **************************************************************************************
    



    # ****************************************
    # *   Zero bias PU45 dataset - 01 Mar    *
    # ****************************************
    #_#  'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_26_1_58B.root',

    #######################################
    ## Commented out for debugging
    ####################################### 
    
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_9_1_ENp.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_19_1_7Ei.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_3_1_yoI.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_10_1_O2i.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_7_1_51m.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_6_1_O8r.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_5_1_qHT.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_8_1_lij.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_25_1_Hhb.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_17_1_MCa.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_2_1_vhs.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_13_1_5Y6.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_12_1_Q44.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_11_1_9vL.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_20_1_tgl.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_15_1_wtp.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_16_1_yp4.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_1_1_00H.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_4_1_P7W.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_22_1_Wit.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_23_1_BPm.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_18_1_SwA.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_14_1_hp4.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_24_1_DN2.root',
    #_#'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_01Mar/JetCollections_21_1_cXE.root'


    #####################################################################################################################################
    
    # Zero bias PU45 dataset - 26 Feb
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_2_1_Aaj.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_2_1_Aaj.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_11_1_Yss.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_18_1_EAf.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_15_1_3aG.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_8_1_lsV.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_3_1_6Dh.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_17_1_6Vz.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_21_1_GSC.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_25_1_W2h.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_7_1_4as.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_14_1_Paj.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_12_1_1hd.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_26_1_Gug.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_19_1_XMb.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_22_1_ujE.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_4_1_iaX.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_1_1_oiF.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_20_1_3V3.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_5_1_TqY.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_6_1_lm9.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_10_1_jam.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_24_1_7lp.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_13_1_cgM.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_9_1_Oi6.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_16_1_83t.root',
    #'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/ZeroBias45PU_26Feb/JetCollections_23_1_wgQ.root'
    ),
   skipEvents = cms.untracked.uint32(0)
)



process.TFileService = cms.Service("TFileService",

#    fileName = cms.string('ZeroBias45PU_SLHCJetAnalysis_30Apr_MARKI_TEST.root')
#    fileName = cms.string('ZeroBias45PU_SLHCJetAnalysis_18Apr.root')
#    fileName = cms.string('ZeroBias45PU_SLHCJetAnalysis_18Apr_JetSatFix.root')
#    fileName = cms.string('ZeroBias45PU_SLHCJetAnalysis_19Apr_JetSatEtaFix.root')
#    fileName = cms.string('ZeroBias45PU_SLHCJetAnalysis_26Apr_MARKI.root')
#    fileName = cms.string('ZeroBias45PU_SLHCJetAnalysis_30Apr_MARKI.root')
     fileName = cms.string('ZeroBias45PU_06May_Mk1_TEST_RS.root')
#     fileName = cms.string('ZeroBias45PU_06May_Mk1_Calib.root')
#     fileName = cms.string('ZeroBias45PU_06May_Mk1_5GeVSeed.root')
#     fileName = cms.string('ZeroBias45PU_06May_Mk1_RingSub_5GeVSeed.root')
#     fileName = cms.string('ZeroBias45PU_06May_Mk1_RingSub.root')
                                   
)


##For timing information:
#process.Timing=cms.Service("Timing")
#
#process.MessageLogger = cms.Service("MessageLogger",
#                                    qemloosenewlog = cms.untracked.PSet(
#    threshold = cms.untracked.string('INFO')
#    ),
#                                    destinations = cms.untracked.vstring('qemloosenewlog')
#                                    )
#


process.analyzer = cms.EDAnalyzer('AnalyzeJets',
                                   CalibCircle8MHT     = cms.InputTag("L1CalibFilterTowerJetProducer:TowerMHT"),
                                   CalibCircle8l1extra = cms.InputTag("L1CalibFilterTowerJetProducer:Cen8x8"),   
                                   CalibCircle8        = cms.InputTag("L1CalibFilterTowerJetProducer:CalibCenJets"),

                                   TowerHT             = cms.InputTag("L1CalibFilterTowerJetProducer:TowerHT"),
                                   UpgrCenJet          = cms.InputTag("L1CalibFilterTowerJetProducer:CenJets"),
                                  # PrePUSubUpgrCenJet  = cms.InputTag("L1CalibFilterTowerJetProducer:PrePUSubCenJets"),
                                   PrePUSubUpgrCenJet  = cms.InputTag("L1TowerJetPUSubtractedProducer:PrePUSubCenJets"),
                                  
                                  # CaloRho             = cms.InputTag("L1CalibFilterTowerJetProducer:Rho"),
                                  # CaloRho             = cms.InputTag("L1TowerJetPUEstimator:Rho"),
                                  
                                  
                                   extrajet = cms.VInputTag("l1extraParticles:Central","l1extraParticles:Tau",),
                                   L1extraMHT = cms.InputTag("l1extraParticles:MHT"),
                                   RecoVertices = cms.InputTag("offlinePrimaryVertices"),
                                   PUsubCaloJets = cms.InputTag("PUsubAK5CaloJetProducer"),
                                   jetIDHelperConfig = cms.InputTag("ak5JetID"),
                                   RhoCaloJets = cms.InputTag("ak5CaloJets:rho"),
jetCollection = cms.InputTag("ak5CaloJets"),
jetIDMap = cms.InputTag("ak5JetID"),

)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.p = cms.Path(process.analyzer)

