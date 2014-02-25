// -*- C++ -*-
//
// Package:    CalibrateJets
// Class:      CalibrateJets
// 
/**\class CalibrateJets CalibrateJets.cc CalibrateJets/CalibrateJets/src/CalibrateJets.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Robyn Elizabeth Lucas,510 1-002,+41227673823,
//         Created:  Sun Dec 16 12:49:20 CET 2012
// $Id$
//
//

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <fstream>
#include "DataFormats/JetReco/interface/JetID.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

using namespace l1slhc;
using namespace edm;
using namespace std;
using namespace reco;

bool myfunction (float i,float j) { return (i>j); }

bool sortTLorentz (TLorentzVector i,TLorentzVector j) { return ( i.Pt()>j.Pt() ); }

//
// class declaration
//

class CalibrateJets : public edm::EDAnalyzer {
   public:
      explicit CalibrateJets(const edm::ParameterSet&);
      ~CalibrateJets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      double Median(vector<double> aVec);
      float get_rho(double L1rho);
      float val_pt_cal(float l1pt, float l1eta);
      std::pair< double, int> Leading_Match(vector<TLorentzVector> offlineJets,vector<TLorentzVector> L1Jets );
      void look_up_table();
      void read_lookup_pt();
      float get_cal_pt( float pt , float eta );
      void DoRateCalc(TH1D* h1, TH1D *h2, int preScale, int nLumis);
      void TMVA_calibration();

      // ----------member data ---------------------------


  ParameterSet conf_;
   edm::FileInPath inRhoData_edm;
   edm::FileInPath inPtData_edm;
   edm::FileInPath inMVAweights_edm;
  TTree* tree_;
  TTree* rhoTree_;
  float sub_;
  float oPhi_;
  float orho_;
  float oEta_;
  float l1Pt_;
  float l1Phi_;
  float l1Eta_;
  float oPt_;
  float ratio_;
  float rhoraw_;
  float mva_rho_;
  float rhorawCal_;
  float rhoraw0_;
  float rhoratio_;
  TMVA::Reader *reader;
  TMVA::Reader *rhoreader;
  float l1Pt, l1Eta, l1Phi, RhoL1, RhoRaw_in, MVA_Rho;

  ifstream indata;
  ifstream inrhodata;
  float pt_, eta_, phi_, cal_pt_;
  vector < pair<double, double> > rho_cal_vec;

  int binsEta;
  int binsPhi;
  int binsLogPt;
  float corrFactors[100][100];

  vector<int> lumis;
  TH1F* NVTX;


      TH1F*    DeltaEta        ;
      TH1F*    CurrDeltaEta    ;
      TH1F*    DeltaPhi        ;
      TH1F*    CurrDeltaPhi    ;
      TH1F*    DeltaPt         ;
      TH1F*    CurrDeltaPt     ;

      TH2F* OffPt_L1Pt  ; 
      TH2F* OffPt_L1Pt_1;
      TH2F* OffPt_L1Pt_2;

      TH1F* DeltaPt_1 ;
      TH2F* DeltaPt_1_pt ;
      TH2F* CurrDeltaPt_pt ;
      TH2F* DeltaPt_2_pt ;
      TH1F* DeltaPt_2;

      TH2F*   DeltaPt_1_nvtx   ;
      TH2F*   CurrDeltaPt_nvtx ;
      TH2F*   DeltaPt_2_nvtx   ;

      TH2F*   DeltaPt_1_eta   ;
      TH2F*   CurrDeltaPt_eta ;
      TH2F*   DeltaPt_2_eta   ;

      TH3F*    OfflinePt3DIC40  ;
      TH3F*    OfflinePt3DIC45  ;
      TH3F*    OfflinePt3DIC50  ;
      TH3F*    OfflinePt3DIC60  ;
      TH3F*    OfflinePt3DIC80  ;
      TH3F*    OfflinePt3DIC100 ;
      TH3F*    OfflinePt3DIC120 ;
      TH3F*    OfflinePt3DIC140 ;
      TH3F*    OfflinePt3DICDen ;


      TH3F*    OfflinePt3DIC40_1  ;
      TH3F*    OfflinePt3DIC45_1  ;
      TH3F*    OfflinePt3DIC50_1  ;
      TH3F*    OfflinePt3DIC60_1  ;
      TH3F*    OfflinePt3DIC80_1  ;
      TH3F*    OfflinePt3DIC100_1 ;
      TH3F*    OfflinePt3DIC120_1 ;
      TH3F*    OfflinePt3DIC140_1 ;
      TH3F*    OfflinePt3DICDen_1 ;
                              
                              
      TH3F*    OfflinePt3DCurr40  ;
      TH3F*    OfflinePt3DCurr45  ;
      TH3F*    OfflinePt3DCurr50  ;
      TH3F*    OfflinePt3DCurr60  ;
      TH3F*    OfflinePt3DCurr80  ;
      TH3F*    OfflinePt3DCurr100 ;
      TH3F*    OfflinePt3DCurr120 ;
      TH3F*    OfflinePt3DCurr140 ;
      TH3F*    OfflinePt3DCurrDen ;


      TH2F*  pt_value_2D;
      TH3F*  pt_value;
      TH3F*  large_value;



      //Rates
      TH1D *ICjets_hist[4];
      TH1D *ICjets_rate[4];
    
      TH1D *CurrL1jets_hist[4];
      TH1D *CurrL1jets_rate[4];

      //HT turn on
      TH1F*  IC_HT_40 ; 
      TH1F*  IC_HT_50 ; 
      TH1F*  IC_HT_60 ; 
      TH1F*  IC_HT_75 ; 
      TH1F*  IC_HT_100; 
      TH1F*  IC_HT_120; 
      TH1F*  IC_HT_150; 
      TH1F*  IC_HT_175; 
      TH1F*  IC_HT_200; 
      TH1F*  IC_HT_300; 
      TH1F*  HT_IC_all; 

      TH1F*  Curr_HT_40 ;
      TH1F*  Curr_HT_50 ;
      TH1F*  Curr_HT_60 ;
      TH1F*  Curr_HT_75 ;
      TH1F*  Curr_HT_100;
      TH1F*  Curr_HT_120;
      TH1F*  Curr_HT_150;
      TH1F*  Curr_HT_175;
      TH1F*  Curr_HT_200;
      TH1F*  Curr_HT_300;
      TH1F*  Curr_HT_all;
      
      TH1F*  Offline_HT;
      TH1F*  upgrade_HT;
      TH1F*  current_HT;

      TH2F*  Offline_HT_vs_upgrade_HT; 
      TH2F*  Offline_HT_vs_current_HT; 


      int preScale;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CalibrateJets::CalibrateJets(const edm::ParameterSet& iConfig):conf_(iConfig)

{

  inRhoData_edm = iConfig.getParameter<edm::FileInPath> ("inRhodata_file");
  inPtData_edm = iConfig.getParameter<edm::FileInPath> ("inPtdata_file");
  inMVAweights_edm = iConfig.getParameter<edm::FileInPath> ("inMVA_weights_file");

  Service<TFileService> fs;


  NVTX            = fs->make<TH1F>("NVTX"           ,"NVTX"          , 100,0.0,100.0);

  DeltaEta        = fs->make<TH1F>("L1_DeltaEta"    ,  "L1_DeltaEta"    , 100,-0.5,0.5);
  CurrDeltaEta    = fs->make<TH1F>("CurrL1_DeltaEta", "CurrL1_DeltaEta" , 100,-0.5,0.5);
  DeltaPhi        = fs->make<TH1F>("L1_DeltaPhi"    ,  "L1_DeltaPhi"    , 100,-0.5,0.5);
  CurrDeltaPhi    = fs->make<TH1F>("CurrL1_DeltaPhi", "CurrL1_DeltaPhi" , 100,-0.5,0.5);
  DeltaPt         = fs->make<TH1F>("L1_DeltaPt"     , "L1_DeltaPt"      , 100,-2,2);
  CurrDeltaPt     = fs->make<TH1F>("CurrL1_DeltaPt" ,  "CurrL1_DeltaPt" , 100,-2,2);

  OffPt_L1Pt      = fs->make<TH2F>("OffPt_L1Pt"    , "OffPt_L1Pt"    , 300,0,600, 300,0,600 );
  OffPt_L1Pt_1    = fs->make<TH2F>("OffPt_L1Pt_1"  , "OffPt_L1Pt_1"  , 300,0,600, 300,0,600 );
  OffPt_L1Pt_2    = fs->make<TH2F>("OffPt_L1Pt_2"  , "OffPt_L1Pt_2"  , 300,0,600, 300,0,600 );

  DeltaPt_1       = fs->make<TH1F>("L1_DeltaPt_1"   , "L1_DeltaPt_1"    , 100,-2,2);
  DeltaPt_2       = fs->make<TH1F>("DeltaPt_2"   , "DeltaPt_2"    , 100,-2,2 );

  DeltaPt_1_pt    = fs->make<TH2F>("DeltaPt_1_pt"   , "DeltaPt_1_pt"    , 100,-2,2, 300,0,600 );
  CurrDeltaPt_pt  = fs->make<TH2F>("CurrDeltaPt_pt" , "CurrDeltaPt_pt"  , 100,-2,2, 300,0,600 );
  DeltaPt_2_pt    = fs->make<TH2F>("DeltaPt_2_pt"   , "DeltaPt_2_pt"    , 100,-2,2, 300,0,600 );

  DeltaPt_1_nvtx  = fs->make<TH2F>("DeltaPt_1_nvtx"  , "DeltaPt_1_nvtx"  , 100,-2,2, 100,0,100 );
  CurrDeltaPt_nvtx= fs->make<TH2F>("CurrDeltaPt_nvtx", "CurrDeltaPt_nvtx", 100,-2,2, 100,0,100 );
  DeltaPt_2_nvtx  = fs->make<TH2F>("DeltaPt_2_nvtx"  , "DeltaPt_2_nvtx"  , 100,-2,2, 100,0,100 );

  DeltaPt_1_eta   = fs->make<TH2F>("DeltaPt_1_eta"  , "DeltaPt_1_eta"  , 100,-2,2, 100,-5.0,5.0 );
  CurrDeltaPt_eta = fs->make<TH2F>("CurrDeltaPt_eta", "CurrDeltaPt_eta", 100,-2,2, 100,-5.0,5.0 );
  DeltaPt_2_eta   = fs->make<TH2F>("DeltaPt_2_eta"  , "DeltaPt_2_eta"  , 100,-2,2, 100,-5.0,5.0 );







  OfflinePt3DIC40 = fs->make<TH3F>("OfflinePt3DIC40 ","OfflinePt3DIC40 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC45 = fs->make<TH3F>("OfflinePt3DIC45 ","OfflinePt3DIC45 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC50 = fs->make<TH3F>("OfflinePt3DIC50 ","OfflinePt3DIC50 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC60 = fs->make<TH3F>("OfflinePt3DIC60 ","OfflinePt3DIC60 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC80 = fs->make<TH3F>("OfflinePt3DIC80 ","OfflinePt3DIC80 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC100= fs->make<TH3F>("OfflinePt3DIC100","OfflinePt3DIC100", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC120= fs->make<TH3F>("OfflinePt3DIC120","OfflinePt3DIC120", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC140= fs->make<TH3F>("OfflinePt3DIC140","OfflinePt3DIC140", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DICDen= fs->make<TH3F>("OfflinePt3DICDen","OfflinePt3DICDen", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);

                                                                                                    


   OfflinePt3DIC40_1  = fs->make<TH3F>("OfflinePt3DIC40_1 ","OfflinePt3DIC40_1 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
   OfflinePt3DIC45_1  = fs->make<TH3F>("OfflinePt3DIC45_1 ","OfflinePt3DIC45_1 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
   OfflinePt3DIC50_1  = fs->make<TH3F>("OfflinePt3DIC50_1 ","OfflinePt3DIC50_1 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
   OfflinePt3DIC60_1  = fs->make<TH3F>("OfflinePt3DIC60_1 ","OfflinePt3DIC60_1 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
   OfflinePt3DIC80_1  = fs->make<TH3F>("OfflinePt3DIC80_1 ","OfflinePt3DIC80_1 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
   OfflinePt3DIC100_1 = fs->make<TH3F>("OfflinePt3DIC100_1","OfflinePt3DIC100_1", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
   OfflinePt3DIC120_1 = fs->make<TH3F>("OfflinePt3DIC120_1","OfflinePt3DIC120_1", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
   OfflinePt3DIC140_1 = fs->make<TH3F>("OfflinePt3DIC140_1","OfflinePt3DIC140_1", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
   OfflinePt3DICDen_1 = fs->make<TH3F>("OfflinePt3DICDen_1","OfflinePt3DICDen _1", 80,-4.0,4.0, 80,-4.0,4.0,300,0,300);


  OfflinePt3DCurr40 = fs->make<TH3F>("OfflinePt3DCurr40 ","OfflinePt3DCurr40 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr45 = fs->make<TH3F>("OfflinePt3DCurr45 ","OfflinePt3DCurr45 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr50 = fs->make<TH3F>("OfflinePt3DCurr50 ","OfflinePt3DCurr50 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr60 = fs->make<TH3F>("OfflinePt3DCurr60 ","OfflinePt3DCurr60 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr80 = fs->make<TH3F>("OfflinePt3DCurr80 ","OfflinePt3DCurr80 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr100= fs->make<TH3F>("OfflinePt3DCurr100","OfflinePt3DCurr100", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr120= fs->make<TH3F>("OfflinePt3DCurr120","OfflinePt3DCurr120", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr140= fs->make<TH3F>("OfflinePt3DCurr140","OfflinePt3DCurr140", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurrDen= fs->make<TH3F>("OfflinePt3DCurrDen","OfflinePt3DCurrDen", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);

  pt_value     = fs->make<TH3F>("pt_value"   ,"pt_value"   , 500,0,500, 80,-4.0,4.0, 1000,0,1000);
  large_value  = fs->make<TH3F>("large_value","large_value", 500,0,500, 80,-4.0,4.0, 1000,0,1000);
  pt_value_2D = fs->make<TH2F>("pt_value_2D"   ,"pt_value_2D"   , 500,0,500, 80,-4.0,4.0);



  ICjets_hist[0] = fs->make <TH1D>("ICjets_1_distro", "ICjets_1_distro", 255, 0., 255.);
  ICjets_hist[1] = fs->make <TH1D>("ICjets_2_distro", "ICjets_2_distro", 255, 0., 255.);
  ICjets_hist[2] = fs->make <TH1D>("ICjets_3_distro", "ICjets_3_distro", 255, 0., 255.);
  ICjets_hist[3] = fs->make <TH1D>("ICjets_4_distro", "ICjets_4_distro", 255, 0., 255.);
  ICjets_rate[0] = fs->make <TH1D>("ICSingleJet_rate", "Single Jet", 255, 0., 255.);
  ICjets_rate[1] = fs->make <TH1D>("ICDoubleJet_rate", "Double Jet", 255, 0., 255.);
  ICjets_rate[2] = fs->make <TH1D>("ICTripleJet_rate", "Triple Jet", 255, 0., 255.);
  ICjets_rate[3] = fs->make <TH1D>("ICQuadJet_rate", "Quad Jet", 255, 0., 255.);

  CurrL1jets_hist[0] = fs->make <TH1D>("CurrL1jets_1_distro", "CurrL1jets_1_distro", 255, 0., 255.);
  CurrL1jets_hist[1] = fs->make <TH1D>("CurrL1jets_2_distro", "CurrL1jets_2_distro", 255, 0., 255.);
  CurrL1jets_hist[2] = fs->make <TH1D>("CurrL1jets_3_distro", "CurrL1jets_3_distro", 255, 0., 255.);
  CurrL1jets_hist[3] = fs->make <TH1D>("CurrL1jets_4_distro", "CurrL1jets_4_distro", 255, 0., 255.);
  CurrL1jets_rate[0] = fs->make <TH1D>("CurrL1SingleJet_rate", "Single Jet", 255, 0., 255.);
  CurrL1jets_rate[1] = fs->make <TH1D>("CurrL1DoubleJet_rate", "Double Jet", 255, 0., 255.);
  CurrL1jets_rate[2] = fs->make <TH1D>("CurrL1TripleJet_rate", "Triple Jet", 255, 0., 255.);
  CurrL1jets_rate[3] = fs->make <TH1D>("CurrL1QuadJet_rate", "Quad Jet", 255, 0., 255.);



  IC_HT_40  = fs->make<TH1F>("IC_HT_40 " ,"IC_HT_40 " , 600,0,600);
  IC_HT_50  = fs->make<TH1F>("IC_HT_50 " ,"IC_HT_50 " , 600,0,600);
  IC_HT_60  = fs->make<TH1F>("IC_HT_60 " ,"IC_HT_60 " , 600,0,600);
  IC_HT_75  = fs->make<TH1F>("IC_HT_75 " ,"IC_HT_75 " , 600,0,600);
  IC_HT_100 = fs->make<TH1F>("IC_HT_100" ,"IC_HT_100" , 600,0,600);
  IC_HT_120 = fs->make<TH1F>("IC_HT_120" ,"IC_HT_120" , 600,0,600);
  IC_HT_150 = fs->make<TH1F>("IC_HT_150" ,"IC_HT_150" , 600,0,600);
  IC_HT_175 = fs->make<TH1F>("IC_HT_175" ,"IC_HT_175" , 600,0,600);
  IC_HT_200 = fs->make<TH1F>("IC_HT_200" ,"IC_HT_200" , 600,0,600);
  IC_HT_300 = fs->make<TH1F>("IC_HT_300" ,"IC_HT_300" , 600,0,600);
  HT_IC_all = fs->make<TH1F>("HT_IC_all" ,"HT_IC_all" , 600,0,600);
  
  //HT turn on curves
  Curr_HT_40   = fs->make<TH1F>("Curr_HT_40 " ,"Curr_HT_40 " , 600,0,600); 
  Curr_HT_50   = fs->make<TH1F>("Curr_HT_50 " ,"Curr_HT_50 " , 600,0,600);
  Curr_HT_60   = fs->make<TH1F>("Curr_HT_60 " ,"Curr_HT_60 " , 600,0,600);
  Curr_HT_75   = fs->make<TH1F>("Curr_HT_75 " ,"Curr_HT_75 " , 600,0,600);
  Curr_HT_100  = fs->make<TH1F>("Curr_HT_100" ,"Curr_HT_100" , 600,0,600);
  Curr_HT_120  = fs->make<TH1F>("Curr_HT_120" ,"Curr_HT_120" , 600,0,600);
  Curr_HT_150  = fs->make<TH1F>("Curr_HT_150" ,"Curr_HT_150" , 600,0,600);
  Curr_HT_175  = fs->make<TH1F>("Curr_HT_175" ,"Curr_HT_175" , 600,0,600);
  Curr_HT_200  = fs->make<TH1F>("Curr_HT_200" ,"Curr_HT_200" , 600,0,600);
  Curr_HT_300  = fs->make<TH1F>("Curr_HT_300" ,"Curr_HT_300" , 600,0,600);
  Curr_HT_all  = fs->make<TH1F>("Curr_HT_all" ,"Curr_HT_all" , 600,0,600);


  Offline_HT  = fs->make<TH1F>("Offline_HT"  ,"Offline_HT"    , 800,0.0,800.0);
  upgrade_HT  = fs->make<TH1F>("upgrade_HT"  ,"upgrade_HT"    , 800,0.0,800.0);
  current_HT  = fs->make<TH1F>("current_HT"  ,"current_HT"    , 800,0.0,800.0);

  Offline_HT_vs_upgrade_HT = fs->make<TH2F>("Offline_HT_vs_upgrade_HT","Offline_HT_vs_upgrade_HT", 800,0,800, 800,0,800);
  Offline_HT_vs_current_HT = fs->make<TH2F>("Offline_HT_vs_current_HT","Offline_HT_vs_current_HT", 800,0,800, 800,0,800);













}


CalibrateJets::~CalibrateJets()
{
 
	// do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
std::pair< double, int> CalibrateJets::Leading_Match(vector<TLorentzVector> offlineJets,vector<TLorentzVector> L1Jets )
{
  std::pair< double, int> dist_ref;
  dist_ref.first = 999; //this is the distance between jets
  dist_ref.second = 999; //dummy initialization of jet reference

  oPt_=offlineJets[0].Pt();
  oEta_=offlineJets[0].Eta();
  oPhi_=offlineJets[0].Phi();
  

  for(unsigned int i = 0 ; i<L1Jets.size(); ++i){
    if(offlineJets[0].Pt()<=0 || L1Jets[i].Pt()<=0) continue;

  //cout<<" ic eta: "<< upgradeJets[i].Eta()<<" off eta:  "<<offlineJets[0].Eta()<<endl;
  //cout<<" ic pt: "<< upgradeJets[i].Pt()<<" off pt:  "<<offlineJets[0].Pt()<<endl;
  //cout<<" calulcating DR: ";
  //cout<<" DR = "<<upgradeJets[i].DeltaR(offlineJets[0])<<endl;

    if(L1Jets[i].DeltaR(offlineJets[0])<dist_ref.first) {
      dist_ref.first = L1Jets[i].DeltaR(offlineJets[0]) ;
      dist_ref.second = i;
    }
  }

  return dist_ref;
}

float CalibrateJets::get_rho(double L1_rho)
{

  //get the rho multiplication factor:
  if(L1_rho<=40.5)   return rho_cal_vec[L1_rho*2].second;
  else return 1.44576*L1_rho;

}
float
CalibrateJets::val_pt_cal(float l1pt, float l1eta)
{
  l1Eta=l1eta;
  l1Pt=l1pt;

  Float_t val = (reader->EvaluateRegression(TString("BDT method") ))[0];

  //cout<<"l1 pt: " << l1pt <<" corr_pt_1 "<<corr_pt_1<<endl;  


  return val;


}

double refine_val( float x )
{

  double cal_fac(0);

  double mul (1);

  if( x >10 && x< 250 ){
    double expo[] = {1.58085e+00 , -6.16652e-02};
    double pol[] = {-3.55742e+00 , 1.48633e-01 , -2.58165e-03 , 2.35495e-05 , -1.18144e-07 , 3.07489e-10, -3.23791e-13}; 

    cal_fac+=exp( expo[0]+expo[1]*x);
 
    for(unsigned int i=0; i<7; i++){
      cal_fac+=pol[i]*pow(x,i);
    }

  mul =  1/(1- cal_fac) ;

  }


  return mul;

}


double CalibrateJets::Median( vector<double> aVec){
  sort( aVec.begin(), aVec.end() );
  double median(0);
  int size = aVec.size();
  if(size ==0){
    median = 0;
  }
  else if(size==1){
    median = aVec[size-1];
  }
  else if( size%2 == 0 ){
    median = ( aVec[ (size/2)-1  ] + aVec[ (size /2) ] )/2;
  }else{
    median = aVec [size/2+1 ];
  }
  return median;
}


void CalibrateJets::DoRateCalc(TH1D* h1, TH1D *h2, int preScale, int nLumis)
{
//short function to calculate rate plots

    Int_t nbins = h1->GetNbinsX();
    
    //This is a constant factor, in s, defining the time per lumi section
    Double_t tLumi = 23.3570304;

    //average run lumi in units of 10^31 cm^-2 s^-1
    //Double_t runAvgLumi = 4.7; //2011 HPU Run
    Double_t runAvgLumi = 2; //2012 66PU Run - rough estimate!
    //Double_t runAvgLumi = 1.5; //2012 45PU Run - rough estimate!

    Double_t lumiAdjust = 2000./runAvgLumi; //translation factor for 2e34 inst lumi
    //Double_t lumiAdjust = 1; //no luminosity scaling
     
    Double_t normalization = (1.*preScale*lumiAdjust)/(nLumis*tLumi); //check!! add in tLumi (~23.45secs)

    h2->GetXaxis()->SetTitle("Threshold (GeV)");
    h2->GetYaxis()->SetTitle("Rate (Hz)");
    //h2->Sumw2();
    for(Int_t i=1; i < nbins; ++i)
    {
      float rate = h1->Integral(i, nbins+1);
      h2->SetBinContent(i, rate*normalization);
    }
    
}
// ------------ method called for each event  ------------
void
CalibrateJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


  bool evValid =true;

  double DR = 0.5;
  float Pt_min=30.;

  //Handle for offline AK5 calo jets
  edm::Handle<reco::CaloJetCollection> Calojets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CaloJets"), Calojets);// cout<<" Calojets: "<< Calojets.isValid() ;
  if(!Calojets.isValid()){evValid=false;}

//   edm::Handle<double> rhoCALOjets;  
//   iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RhoCaloJets"), rhoCALOjets); cout<<" rhoCALOjets.: "<< rhoCALOjets.isValid() ;
//   if(!rhoCALOjets.isValid()){evValid=false;}

  edm::Handle<float> rhoTOWjets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RhoTowerJets"), rhoTOWjets); //cout<<" rhoTOWjets.: "<< rhoTOWjets.isValid() ;

	//Handle for uncalibrated central jets
  edm::Handle<L1TowerJetCollection > UnCalib;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle8"), UnCalib);//cout<<" UnCalib.: "<< UnCalib.isValid() ;


	//Handle to get fwd jet collection: uncalibrated, non pu subtracted
  edm::Handle<L1TowerJetCollection > UnCalibFwd;
  iEvent.getByLabel("L1TowerFwdJetFilter2D", UnCalibFwd);

	//Handle for vertices
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel("offlinePrimaryVertices", vtx); //cout<<"vtx: "<< vtx.isValid() ;
  if(!vtx.isValid()){evValid=false;     }


  vector<InputTag>  l1extraparticlesjet= conf_.getParameter< vector < InputTag > >("extrajet");
  vector<InputTag>  l1extraparticlesem= conf_.getParameter< vector < InputTag > >("extraem");

   edm::Handle< edm::ValueMap<reco::JetID> > jetID;
   iEvent.getByLabel( conf_.getParameter<edm::InputTag>("jetIDHelperConfig"), jetID );


 // handle to the jet ID variables
   edm::Handle<reco::JetIDValueMap> hJetIDMap;
   iEvent.getByLabel( conf_.getParameter<edm::InputTag>("jetIDHelperConfig") , hJetIDMap );


   edm::Handle<edm::View< reco::CaloJet > > hJets; // uncorrected jets!
   iEvent.getByLabel("ak5CaloJets", hJets );

   edm::Handle<l1extra::L1EtMissParticleCollection> L1MHT;
   iEvent.getByLabel("l1extraParticles","MHT",L1MHT);

  if( !evValid ) {/*cout<<"ERROR: invalid event....................."<<endl;*/}
  else{

//cout<<"valid event....................."<<endl;


    int PVSIZE=vtx->size();
    NVTX->Fill(PVSIZE);
    if (PVSIZE>5){

      ///////////////////////////////////////////////////////////////////////////////////////
      // L1 RHO
      ///////////////////////////////////////////////////////////////////////////////////////

      int count(0);
      vector<double> Jet2Energies;
      for (L1TowerJetCollection::const_iterator il1 = UnCalib->begin();
        il1!= UnCalib->end() ; ++il1 ){
        if( abs(il1->p4().eta() )>3) continue;
        if(count>1) {
          Jet2Energies.push_back(il1->p4().Pt());   
        }
        count++;
      }
      double areaPerJet = 52 * (0.087 * 0.087) ;

      float raw_rho2 = ( Median( Jet2Energies ) / areaPerJet );

      double cal_rhoL1 = raw_rho2 * get_rho(raw_rho2);

      //cout<<"Jet energies.size : "<<Jet2Energies.size() <<" raw rho : "<< raw_rho2 <<" calibrated rho: "<<cal_rhoL1<<endl;
     

      TLorentzVector off,curr,IC;
      vector<TLorentzVector> offlineJets;

      ///////////////////////////////////////////////////////////////////////////////////////
      //Offline Jets
      ///////////////////////////////////////////////////////////////////////////////////////

      vector<double> fHPD;
      vector<double> n90hits;
      vector<double> emfrac;

      unsigned int idx;
    
      //L2L3 PU corrected jet collection
      reco::CaloJetCollection::const_iterator it = Calojets->begin() ;
    
      //uncorrected jet collection for jet ID reference
      for ( edm::View<reco::CaloJet>::const_iterator ibegin = hJets->begin(),
              iend = hJets->end(), ijet = ibegin;
            ijet != iend; ++ijet , ++it) {
    
        idx = ijet - ibegin;
        edm::RefToBase<reco::CaloJet> jetRef = hJets->refAt(idx);  
        reco::JetID const & jetId = (*hJetIDMap)[ jetRef ];
    
        //set the TLorentz vector with corrected energy
        off.SetPtEtaPhiM(it->p4().Pt(),it->p4().eta(),it->p4().phi(),it->p4().M());
    
        //loose jet selections
    
        if (it->emEnergyFraction()<0.01) continue;
        if (jetId.fHPD>0.98  ) continue;
        if (jetId.n90Hits<=1) continue;
        if (off.Pt()<10) continue;

        //Fill vector with loose ID's calojets
        offlineJets.push_back(off);
        //out<<" offline (pt,eta,phi) =("<<off.Pt()<<","<<off.Eta()<<","<<off.Phi()<<")"<<endl;
         
      }
       
      if(offlineJets.size()>0)  sort (offlineJets.begin(), offlineJets.end(),sortTLorentz);



      ///////////////////////////////////////////////////////////////////////////////////////
      //Upgrade central Jets
      ///////////////////////////////////////////////////////////////////////////////////////
      
      vector<TLorentzVector> upgradeJets;

      for (L1TowerJetCollection::const_iterator il1 = UnCalib->begin(); 
          il1!= UnCalib->end() ; ++il1 ){
   
        //PU subtract:
        double l1PU_pt = il1->p4().Pt() - areaPerJet* cal_rhoL1;
        
        //weighted eta and phi are slightly wrong: adjust
        double w_eta = il1->WeightedEta();
        if( il1->WeightedEta() >= 0.1 ) w_eta += 0.110481;
        if( il1->WeightedEta() < 0 ) w_eta += -0.0167664;
        double w_phi = il1->WeightedPhi() + 8.38659e-02;
        //Keep jet if it has a positive energy
			  if(  l1PU_pt < 0.5 ) continue;

	//			if( corr_pt < 5 ) continue;
        IC.SetPtEtaPhiM(l1PU_pt,w_eta,w_phi,il1->p4().M()); 
        //cout<<" upgrade (pt,eta,phi) =("<<IC.Pt()<<","<<IC.Eta()<<","<<IC.Phi()<<")"<<endl;
        upgradeJets.push_back(IC); 
      }

      if(upgradeJets.size()>0)  sort (upgradeJets.begin(), upgradeJets.end(),sortTLorentz);




      ///////////////////////////////////////////////////////////////////////////////////////
      //Current L1 Jets
      ///////////////////////////////////////////////////////////////////////////////////////
      
      vector<TLorentzVector> currentJets;

      for (uint i=0; i< l1extraparticlesjet.size();++i){
        Handle<l1extra::L1JetParticleCollection> l1Jets;
        iEvent.getByLabel(l1extraparticlesjet[i],l1Jets);
 
        for(l1extra::L1JetParticleCollection::const_iterator itr = l1Jets->begin();  itr != l1Jets->end(); ++itr ) {

          curr.SetPtEtaPhiM(itr->p4().Pt(),itr->p4().eta(),itr->p4().phi(),itr->p4().M());
          currentJets.push_back(  curr );

        }
      }
  
      if(currentJets.size()>0)  sort (currentJets.begin(), currentJets.end(),sortTLorentz);


      //Match jets
  
      bool leading_jet_only = true;
      //try calibrating to the leading jet only
      if(offlineJets.size()>0 && leading_jet_only){  //For the turn on curves
  


      ///////////////////////////////////////////////////////////////////////////////////////
      //Upgrade L1 Jets: Fill histograms
      ///////////////////////////////////////////////////////////////////////////////////////


        std::pair< double, int> dist_ref_up = Leading_Match(offlineJets, upgradeJets);

        //now we have a matched jet, fill tree with values
        if( dist_ref_up.first<DR){
          
          //cout<<"Matched jet. "<<endl;
          //offline jet: offlineJets[0] , upgrade jet: upgradeJets[i]
          
          oPt_ =offlineJets[0].Pt();
          oEta_=offlineJets[0].Eta();
          oPhi_=offlineJets[0].Phi();
          
          l1Pt_ =upgradeJets[dist_ref_up.second].Pt();
          l1Eta_=upgradeJets[dist_ref_up.second].Eta();
          l1Phi_=upgradeJets[dist_ref_up.second].Phi();

          ratio_=oPt_/l1Pt_;

          ///Fill the tree with L1 and Off values
          tree_->Fill();

          //Fill a few histograms to check
          float _dpt4=1-((l1Pt_ )/oPt_);

          DeltaEta       ->Fill(oEta_-l1Eta_) ;
          DeltaPhi       ->Fill(oPhi_-l1Phi_) ;
          DeltaPt        ->Fill(_dpt4) ;

          //offline pt vs l1 pt
          OffPt_L1Pt -> Fill(oPt_ , l1Pt_);

          //cout<<"calibrated pt: "<<val_pt_cal(l1Pt_, l1Eta_, l1Phi_)<<" -> "<< val_pt_cal(l1Pt_, l1Eta_, l1Phi_)*l1Pt_<<endl; 
          double corr_pt = val_pt_cal(l1Pt_, l1Eta_)*l1Pt_;
          _dpt4=1-((corr_pt )/oPt_);
          
        //  if(fabs(_dpt4)> 2) cout<<"calibrated pt: "<<val_pt_cal(l1Pt_, l1Eta_)<<" * "<<l1Pt_<<" gives pt res: "<<_dpt4<<endl;
          
          DeltaPt_1	    ->Fill(_dpt4);
          DeltaPt_1_pt      ->Fill(_dpt4 , corr_pt);
          DeltaPt_1_nvtx    ->Fill(_dpt4 , PVSIZE);
          DeltaPt_1_eta     ->Fill(_dpt4 , l1Eta_);
    
          //Fill a few turn on curves
          
          if(corr_pt>40)     OfflinePt3DIC40->Fill( oEta_ , oPhi_ , oPt_ );
          if(corr_pt>45)     OfflinePt3DIC45->Fill( oEta_ , oPhi_ , oPt_ );
          if(corr_pt>50)     OfflinePt3DIC50->Fill( oEta_ , oPhi_ , oPt_ );
          if(corr_pt>60)     OfflinePt3DIC60->Fill( oEta_ , oPhi_ , oPt_ );
          if(corr_pt>80)     OfflinePt3DIC80->Fill( oEta_ , oPhi_ , oPt_ );
          if(corr_pt>100)    OfflinePt3DIC100->Fill( oEta_ , oPhi_ , oPt_ );
          if(corr_pt>120)    OfflinePt3DIC120->Fill( oEta_ , oPhi_ , oPt_ );
          if(corr_pt>140)    OfflinePt3DIC140->Fill( oEta_ , oPhi_ , oPt_ );

          OfflinePt3DICDen->Fill( oEta_ , oPhi_ , oPt_ );
          OffPt_L1Pt_1 -> Fill(oPt_ , corr_pt);
 
          //L1 Pt vs offline Pt

//          cout<<"multi factor using tmva: "<<val_pt_cal(l1Pt_, l1Eta_)<<endl;
//          cout<<"multi factor using lookup: "<<get_cal_pt(l1Pt_, l1Eta_)<<endl;

          double l1Pt_1 = corr_pt* refine_val(corr_pt);
          _dpt4=1-((l1Pt_1 )/oPt_);
          DeltaPt_2	  ->Fill(_dpt4);
          DeltaPt_2_pt    ->Fill(_dpt4 , l1Pt_1);
          DeltaPt_2_nvtx  ->Fill(_dpt4 , PVSIZE);
          DeltaPt_2_eta   ->Fill(_dpt4 , l1Eta_);
         
          if(l1Pt_1>40)        OfflinePt3DIC40_1  ->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_1>45)        OfflinePt3DIC45_1  ->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_1>50)        OfflinePt3DIC50_1  ->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_1>60)        OfflinePt3DIC60_1  ->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_1>80)        OfflinePt3DIC80_1  ->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_1>100)       OfflinePt3DIC100_1 ->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_1>120)       OfflinePt3DIC120_1 ->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_1>140)       OfflinePt3DIC140_1 ->Fill( oEta_ , oPhi_ , oPt_ );
          OfflinePt3DICDen_1 ->Fill( oEta_ , oPhi_ , oPt_ );

          OffPt_L1Pt_2 -> Fill(oPt_ , l1Pt_1);

        }//DR matching


      ///////////////////////////////////////////////////////////////////////////////////////
      //Current L1 Jets: Fill histograms
      ///////////////////////////////////////////////////////////////////////////////////////

        std::pair< double, int> dist_ref_curr = Leading_Match(offlineJets, currentJets);

        //now we have a matched jet, fill tree with values
        if( dist_ref_curr.first<DR){
          
          //cout<<"Matched jet. "<<endl;
          //offline jet: offlineJets[0] , upgrade jet: upgradeJets[i]
          
          oPt_ =offlineJets[0].Pt();
          oEta_=offlineJets[0].Eta();
          oPhi_=offlineJets[0].Phi();
          
          l1Pt_ =currentJets[dist_ref_curr.second].Pt();
          l1Eta_=currentJets[dist_ref_curr.second].Eta();
          l1Phi_=currentJets[dist_ref_curr.second].Phi();

          ratio_=oPt_/l1Pt_;

          ///Fill the tree with L1 and Off values

          //Fill a few histograms to check
          float _dpt4=1-((l1Pt_ )/oPt_);

          CurrDeltaEta       ->Fill(oEta_-l1Eta_) ;
          CurrDeltaPhi       ->Fill(oPhi_-l1Phi_) ;
          CurrDeltaPt        ->Fill(_dpt4) ;

          //cout<<"calibrated pt: "<<val_pt_cal(l1Pt_, l1Eta_, l1Phi_)<<" -> "<< val_pt_cal(l1Pt_, l1Eta_, l1Phi_)*l1Pt_<<endl; 
          
          CurrDeltaPt_pt ->Fill(_dpt4 , l1Pt_);
          CurrDeltaPt_nvtx   ->Fill(_dpt4, PVSIZE) ;
          CurrDeltaPt_eta    ->Fill(_dpt4, l1Eta_ ) ;


          //Fill a few turn on curves
          
          if(l1Pt_>40)     OfflinePt3DCurr40->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_>45)     OfflinePt3DCurr45->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_>50)     OfflinePt3DCurr50->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_>60)     OfflinePt3DCurr60->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_>80)     OfflinePt3DCurr80->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_>100)    OfflinePt3DCurr100->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_>120)    OfflinePt3DCurr120->Fill( oEta_ , oPhi_ , oPt_ );
          if(l1Pt_>140)    OfflinePt3DCurr140->Fill( oEta_ , oPhi_ , oPt_ );

          OfflinePt3DCurrDen->Fill( oEta_ , oPhi_ , oPt_ );


        }//DR matching




      } //offline jets size > 0 && leading jet


//////////////////////////////////////////////////////////////////////////////////////////////////////
// RATE PLOTS 
//////////////////////////////////////////////////////////////////////////////////////////////////////


      
      vector<double> Cal_upgradeJets_pt;

      for (unsigned int i = 0; i< upgradeJets.size(); i++)
      {
        float corr_pt = val_pt_cal(l1Pt_, l1Eta_)*upgradeJets[i].Pt();
          
   
        Cal_upgradeJets_pt.push_back(corr_pt); 
      }

      if(Cal_upgradeJets_pt.size()>0)  sort (Cal_upgradeJets_pt.begin(), Cal_upgradeJets_pt.end());
 

      int nLumis = lumis.size();
  
   //THIS IS THE ZERO BIAS PRESCALE: MAY BE WRONG
      //Upgrade Jets
   
       for(unsigned int i=0; i<4; i++){
         if (Cal_upgradeJets_pt.size() > i) {
           ICjets_hist[i]->Fill(Cal_upgradeJets_pt.at(i));
           }
       }
   
       for(int i=0;i<4;i++){
         DoRateCalc(ICjets_hist[i], ICjets_rate[i], preScale, nLumis);
       }
   

      //Current L1 Jets

       for(unsigned int i=0; i<4; i++){
         if (currentJets.size() > i) {
           CurrL1jets_hist[i]->Fill(currentJets.at(i).Pt());
         }
       }
 
       //At the moment doing this for all lumi sections: can change
       for(int i=0;i<4;i++){
         DoRateCalc(CurrL1jets_hist[i], CurrL1jets_rate[i], preScale, nLumis);
       }
     
 
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HT PLOTS 
//////////////////////////////////////////////////////////////////////////////////////////////////////


      //uncorrected jet collection for jet ID reference
      double Off_HT(0);
      for (unsigned int i = 0; i< offlineJets.size(); i++)
      {
        if(offlineJets[i].Pt()<40) continue;
        if(fabs(offlineJets[i].Eta())>3) continue;
        Off_HT+=offlineJets[i].Pt();
      }
      if(   Off_HT >= 40 ){

        //calculate upgrade HT
        vector<double> Cal_upgradeJets_ht;
        double up_HT(0);
        for (unsigned int i = 0; i< upgradeJets.size(); i++)
        {
          float corr_pt = val_pt_cal(l1Pt_, l1Eta_)*upgradeJets[i].Pt();
          if( corr_pt > 40 )  up_HT+=corr_pt;
        }
  
        float Curr_HT = L1MHT->begin()->etTotal();
  
        Offline_HT->Fill(Off_HT);
        upgrade_HT->Fill(up_HT);
        current_HT->Fill(Curr_HT);
  
        Offline_HT_vs_upgrade_HT->Fill(Off_HT ,up_HT );
        Offline_HT_vs_current_HT->Fill(Off_HT , Curr_HT);



       if(up_HT>40)  IC_HT_40->Fill(Off_HT);
       if(up_HT>50)  IC_HT_50->Fill(Off_HT);
       if(up_HT>60)  IC_HT_60->Fill(Off_HT);
       if(up_HT>75)  IC_HT_75->Fill(Off_HT);
       if(up_HT>100) IC_HT_100->Fill(Off_HT);
       if(up_HT>120) IC_HT_120->Fill(Off_HT);
       if(up_HT>150) IC_HT_150->Fill(Off_HT);
       if(up_HT>175) IC_HT_175->Fill(Off_HT);
       if(up_HT>200) IC_HT_200->Fill(Off_HT);
       if(up_HT>300) IC_HT_300->Fill(Off_HT);
       HT_IC_all->Fill(Off_HT);

      // cout<<"Current HT "<< Curr_HT<< endl;
       //HT turn on curves
       if(Curr_HT>40)  Curr_HT_40->Fill(Off_HT);
       if(Curr_HT>50)  Curr_HT_50->Fill(Off_HT);
       if(Curr_HT>60)  Curr_HT_60->Fill(Off_HT);
       if(Curr_HT>75)  Curr_HT_75->Fill(Off_HT);
       if(Curr_HT>100) Curr_HT_100->Fill(Off_HT);
       if(Curr_HT>120) Curr_HT_120->Fill(Off_HT);
       if(Curr_HT>150) Curr_HT_150->Fill(Off_HT);
       if(Curr_HT>175) Curr_HT_175->Fill(Off_HT);
       if(Curr_HT>200) Curr_HT_200->Fill(Off_HT);
       if(Curr_HT>300) Curr_HT_300->Fill(Off_HT);
       Curr_HT_all->Fill(Off_HT);
      }



////////////////////////////////////////////////////////////////////////////////////////////
//HF Jets
////////////////////////////////////////////////////////////////////////////////////////////


      vector<TLorentzVector> upgradeFwdJets;
      TLorentzVector ICfwd;
      for (L1TowerJetCollection::const_iterator il1 = UnCalibFwd->begin(); 
          il1!= UnCalibFwd->end() ; ++il1 ){
        areaPerJet = 0.087*8*1; 
        //PU subtract:
        double l1PU_pt = il1->p4().Pt() - areaPerJet* cal_rhoL1;
        
        //weighted eta and phi are slightly wrong: adjust

        if(  l1PU_pt < 0.5 ) continue;

	//			if( corr_pt < 5 ) continue;
        ICfwd.SetPtEtaPhiM(l1PU_pt,il1->p4().Eta(),il1->p4().Phi(),il1->p4().M()); 
        //cout<<" upgrade (pt,eta,phi) =("<<IC.Pt()<<","<<IC.Eta()<<","<<IC.Phi()<<")"<<endl;
        upgradeFwdJets.push_back(IC); 
      }

      if(upgradeFwdJets.size()>0)  sort (upgradeFwdJets.begin(), upgradeFwdJets.end(),sortTLorentz);












    } //PU>5
  } //Valid event

}


// ------------ method called once each job just before starting event loop  ------------
void 
CalibrateJets::beginJob()
{
  std::cout<<"Begin job: declare TTree."<<endl;
 
  tree_=new TTree("t1","t1");
  tree_->Branch("offPt",&oPt_,"offPt/F");
  tree_->Branch("offEta",&oEta_,"offEta/F");
  tree_->Branch("offPhi",&oPhi_,"offPhi/F");
  tree_->Branch("l1Pt",&l1Pt_,"l1Pt/F");
  tree_->Branch("l1Eta",&l1Eta_,"l1Eta/F");
  tree_->Branch("l1Phi",&l1Phi_,"l1Phi/F");
  tree_->Branch("ratio",&ratio_,"ratio/F");
  tree_->Branch("RhoRaw",&rhoraw_,"RhoRaw/F");
  tree_->Branch("MVA_Rho",&mva_rho_,"MVA_Rho/F");
  tree_->Branch("RhoRatio",&rhoratio_,"RhoRatio/F");
  tree_->Branch("RhoOff",&orho_,"RhoOff/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
CalibrateJets::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
CalibrateJets::beginRun(edm::Run const&, edm::EventSetup const&)
{
  inrhodata.open(inRhoData_edm.fullPath().c_str());

  if(!inrhodata){
    cerr<<" Unable to open rho lookup file. "<<endl;
    exit(1);
  }
  //read into a vector
  pair<double, double> rho_cal;
  double L1rho_(9999), calFac_(9999);
  while ( !inrhodata.eof() ) { // keep reading until end-of-file
    // sets EOF flag if no value found
    inrhodata >> L1rho_ >> calFac_ ;

    //bin of rho: bin = L1rho/2   

    rho_cal.first = L1rho_;
    rho_cal.second= calFac_;
    rho_cal_vec.push_back(rho_cal);
  }
  inrhodata.close();

  cout<<" Read in rho lookup table"<<endl;


  TMVA_calibration();
  

  //Convert to look up table in log pt
  //look_up_table();
  read_lookup_pt();


}

void CalibrateJets::TMVA_calibration()
{
  cout<<"Getting lookup from MVA"<<endl;

  reader = new TMVA::Reader("!Color:Silent");
  reader->AddVariable( "l1Pt", &l1Pt);
  reader->AddVariable( "l1Eta", &l1Eta);
//  reader->AddVariable( "l1Phi", &l1Phi);
//  reader->AddVariable( "MVA_Rho", &MVA_Rho );

   cout<<"Booking MVA reader"<<endl;

  reader->BookMVA("BDT method",inMVAweights_edm.fullPath().c_str());

  cout<<"MVA reader booked: start processing events."<<endl;

}

void CalibrateJets::look_up_table()
{

//  float corr_pt = val_pt_cal(l1Pt_, l1Eta_)*l1Pt_;
  int count(0);
  for(float pt=0.1;pt<400; ){
    if(abs(pt)<0.01 ) pt = 0.0;      
    
    for(float eta=-3;eta<3.2;eta+=0.1){          
      if(fabs(eta)<0.01 )        eta = 0;

      double val = val_pt_cal(pt, eta) ;

      pt_value->Fill(pt, eta, val);
      pt_value_2D->Fill(pt, eta, val);
      if( val_pt_cal(pt, eta) > 150 ){
        double large = val_pt_cal(pt, eta);
        for( double et=-3; et<3; et+=0.1){
           if( val_pt_cal(pt, eta+et) < large ) large =  val_pt_cal(pt, eta+et);
           if( large < 150 ){ break; }
        }
        cout<<"   "<<pt<<"\t"<<eta<<"\t\t\t"<<large<<endl;
        continue;
      }

      cout<<" "<<pt<<"\t"<<eta<<"\t\t\t"<<val_pt_cal(pt, eta)<<endl;

    }

    if(pt<15){
      pt = exp( log(pt) + 0.2 ); 
// 	//cout<<"  15 pt: "<<pt<<" "<<count<<endl;
    } else if(pt<105){ 
      pt = exp( log(pt) + 0.1 ) ; 
// 	//cout<<"  105 pt: "<<pt<<" "<<count<<endl; 
    } else { 
      pt = exp( log(pt) + 0.05 ) ;            
// 	//cout<<"  else pt: "<<pt<<" "<<count<<" log "<<pt<<"- log(109.663) ="<< log(pt)-log(109.663)<<" /0.05 = "<<(log(pt)- log(109.663))/0.05 <<endl;
    }
  }

}   

void CalibrateJets::read_lookup_pt()
{

  indata.open(inPtData_edm.fullPath().c_str());

  if(!indata){
    cerr<<" Unable to open pt lookup file. "<<endl;
    exit(1);
  }
  indata.seekg(0);

  double prevpt(0), preveta(0);
  int ptcount(0), etacount(0);
         
  while ( !indata.eof() ) { // keep reading until end-of-file

    // sets EOF flag if no value found
    indata >> pt_ >> eta_ >>  cal_pt_;
 
    if(pt_!=prevpt)  ptcount++;
    if(eta_!=preveta) etacount++;   

    prevpt = pt_;
    preveta = eta_ ;

    if(eta_/(-3) <1.001 && eta_/(-3) >0.9999) etacount=0;

    

    corrFactors[ptcount][etacount] = cal_pt_;
  }

   indata.close();

   cout<<" pt lookup table read in. "<<endl;

}

float CalibrateJets::get_cal_pt( float pt , float eta )
{


  vector<double> pt_values;
  //find the closest pt value in file


  for(float p=0.1;p<400; ){

    pt_values.push_back(p);
    
    if(p<15){
      p = exp( log(p) + 0.2 ); 

    } else if(pt<105){ 
      p = exp( log(p) + 0.1 ) ; 

    } else { 
      p = exp( log(p) + 0.05 ) ;  
    }

  }

  //get pt bin
  pair<double, int> diff;
  diff.first=999;
  diff.second=100;
  for(unsigned int i =0; i<pt_values.size(); i++)
  {
    if(  fabs(pt_values[i] - pt) < diff.first) {
      diff.first = fabs( pt_values[i] - pt) ;
      diff.second = i;
    } 
  }
  int pt_bin = diff.second;
  int eta_bin = (10*eta + 30); 

//  cout<<" eta: "<<eta<<" bin "<<eta_bin<<"  pt: "<< pt <<" bin: "<<pt_bin <<" corrFactor: "<<      corrFactors[pt_bin][eta_bin]   <<endl;
  return corrFactors[pt_bin][eta_bin];
}
    
    
// ------------ method called when ending the processing of a run  ------------
void 
CalibrateJets::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CalibrateJets::beginLuminosityBlock(edm::LuminosityBlock const& iLumiBlock, edm::EventSetup const& iSetup)
{  
  int thisLumi = iLumiBlock.id().luminosityBlock();
  cout<<"luminosity block: "<<    iLumiBlock.id().luminosityBlock()    <<endl;  
  bool newLumi = true;

  for(unsigned int i=0; i<lumis.size(); ++i){
    if(lumis.at(i)==thisLumi) newLumi = false;
  }
  if(newLumi) lumis.push_back(thisLumi);
  if(thisLumi<4) preScale = 9;
  if(thisLumi<103) preScale = 2;
  if(thisLumi<766) preScale = 3;
  else preScale = 4;

}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CalibrateJets::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CalibrateJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CalibrateJets);
