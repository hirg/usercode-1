// -*- C++ -*-
//
// Package:    NtupleProducer
// Class:      NtupleProducer
// 
/**\class NtupleProducer NtupleProducer.cc JetSLHC/NtupleProducer/src/NtupleProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Robyn Lucas
//         Created:  Mon Jul  9 12:08:30 BST 2012
// $Id$
//
//

//#include "../interface/NtupleProducer.h"
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include "TTree.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//root
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TProfile.h"

//File references
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "SimCalorimetry/EcalTrigPrimProducers/plugins/EcalTrigPrimProducer.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveSample.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"  

//SLHC stuff
#include "SimDataFormats/SLHC/interface/L1CaloTriggerSetup.h"
#include "SimDataFormats/SLHC/interface/L1CaloCluster.h"
#include "SimDataFormats/SLHC/interface/L1CaloClusterFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloJet.h"
#include "SimDataFormats/SLHC/interface/L1CaloJetFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloRegion.h"
#include "SimDataFormats/SLHC/interface/L1CaloRegionFwd.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "SimDataFormats/SLHC/interface/L1CaloTowerFwd.h"    
#include "SimDataFormats/SLHC/interface/EtaPhiContainer.h"
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/L1CaloAlgoBase.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#define DEBUG std::cout<< __FILE__ << " : " << __LINE__ << std::endl;

#include "../interface/ClusterTree.h"

double median( vector<double> aVec)
{
//sort the members of the vector
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
    median = aVec [ double( (size/2) ) +0.5 ];
  }
  return median;
}

double rho_scale_op2_L1towBF(double rho, int j)
{ 
  double p0(9999),p1(9999),p2(9999),p3(9999),p4(9999),p5(9999),p6(9999),p7(9999),p8(9999),p9(9999); 
  double a(9999), mean(9999), sigma(9999);   
  if(j==0){ 
    p0  =      1.97713;
    p1  =      1.57091;
    p2  =     0.290832;
    p3  =    -0.115215;
    p4  =    0.0120728;
    p5  = -0.000614244;
    p6  =  1.68188e-05;
    p7  = -2.38646e-07;
    p8  =  1.38035e-09;
    p9  =  0;    
  }if(j==1){
    p0     =      2.03717;
    p1     =       1.5874;
    p2     =     0.503526;
    p3     =    -0.194814;
    p4     =    0.0228755;
    p5     =  -0.00135892;
    p6     =  4.58219e-05;
    p7     = -8.86491e-07;
    p8     =  9.13297e-09;
    p9     =  -3.8559e-11;

  }if(j==2){
    p0  =      2.04799;
    p1  =      1.51676;
    p2  =     0.702718;
    p3  =    -0.265936;
    p4  =      0.03302;
    p5  =  -0.00210803;
    p6  =   7.7203e-05;
    p7  = -1.64081e-06;
    p8  =  1.88477e-08;
    p9  = -9.06008e-11;

  }if(j==3){  //12x12 circ
    p0   =       2.4159;
    p1   =      2.02434;
    p2   =       0.5504;
    p3   =    -0.260468;
    p4   =    0.0348197;
    p5   =  -0.00234878;
    p6   =  9.05368e-05;
    p7   = -2.02555e-06;
    p8   =  2.45321e-08;
    p9   = -1.24604e-10;

  }if(j==4){//8x8 sq 
    p0  =      1.94963;
    p1  =     0.799762;
    p2  =      1.05884;
    p3  =    -0.356832;
    p4  =    0.0467175;
    p5  =  -0.00327205;
    p6  =  0.000133561;
    p7  = -3.18654e-06;
    p8  =  4.12263e-08;
    p9  = -2.23528e-10;

  }if(j==5){//9x9 sq
    p0  =      1.85128;
    p1  =      1.29029;
    p2  =      1.04034;
    p3  =    -0.406005;
    p4  =    0.0568464;
    p5  =  -0.00415508;
    p6  =  0.000174695;
    p7  = -4.25775e-06;
    p8  =  5.59643e-08;
    p9  = -3.07123e-10;

  }if(j==6){ //10x10     sq
    p0  =       1.6354; 
    p1  =      2.88123; 
    p2  =     0.214806; 
    p3  =    -0.233812; 
    p4  =    0.0370448; 
    p5  =   -0.0027886; 
    p6  =  0.000117193; 
    p7  = -2.81934e-06; 
    p8  =  3.63607e-08; 
    p9  = -1.95213e-10;

  }if(j==7)   {   //12x12 sq 
    p0  =      1.38482; 
    p1  =      4.27415; 
    p2  =     -0.44515; 
    p3  =    -0.113427; 
    p4  =    0.0253515; 
    p5  =  -0.00211578; 
    p6  =  9.35611e-05; 
    p7  = -2.32261e-06; 
    p8  =  3.06279e-08; 
    p9  = -1.67308e-10; 

  }
  double scaled_rho=9999;
  if(a<9999){
    scaled_rho = a*exp(- ((pow( (rho-mean),2 ) )/2*sigma*sigma) );
  }else if(p0<9999) 
  { 
    scaled_rho = p0+rho*p1+pow(rho,2)*p2+pow(rho,3)*p3+
            pow(rho,4)*p4+pow(rho,5)*p5+pow(rho,6)*p6+
            pow(rho,7)*p7+pow(rho,8)*p8+pow(rho,9)*p9;
  }
  return scaled_rho; 

}

double rho_Rescale_op2_L1towBF(double rho, int j)
{
  double p0(0),p1(0),p2(0);
  
  if(j==0){ //8x8 circle
    p0=0.00628503;  p1= 0.786142;  p2=0.0263271;
  }if(j==1){
    p0 = -0.452651; p1 = 0.892006;
  }if(j==2){
    p0=0.0450745;   p1=0.755794;   p2=0.022053;
  }if(j==3){  //12x12 circ    
    p0=-0.0164094;  p1=0.75261;    p2=0.0209709;
  }  if(j==4){//8x8 sq  
    p0=0.335104;    p1=0.69551;    p2=0.0264787;
  }if(j==5){//9x9 sq
    p0= 0.274005;   p1= 0.67706;   p2=0.0285173;
  }if(j==6){ //10x10     sq  
    p0=0.0631925;   p1=0.665901;   p2=0.0334327;
  }if(j==7)   {   //12x12 sq
    p0=0.301525;    p1=0.595759;   p2=0.037671;
  }
  return p0+rho*p1+pow(rho,2)*p2;
}

//
// class declaration
//

class NtupleProducer : public edm::EDAnalyzer {
   public:
      explicit NtupleProducer(const edm::ParameterSet&);
      ~NtupleProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      virtual void tow_ntuple_pt(int j, float goodPt, int no_jets);
      virtual void tow_ntuple_eta(int j, float goodEta, int no_jets);
      virtual void tow_ntuple_phi(int j, float goodPhi, int no_jets);
      virtual void tow_ntuple_multi( int j, int no_jets);
      virtual void tow_ntuple_rho( int j, float rho);

      virtual void Off_matched_tow_ntuple_pt( int j, float scaled_pt, int no_jets);
      virtual void Off_matched_tow_ntuple_eta( int j, float eta, int no_jets);
      virtual void Off_matched_tow_ntuple_phi( int j, float phi, int no_jets);
      // ----------member data ---------------------------
      edm::InputTag CaloJets_;
      edm::InputTag CaloJets_rho_;
      edm::InputTag CaloJets_corrections_;
      //edm::InputTag PFJets_;
      //edm::InputTag PFJets_rho_;
      //edm::InputTag PFJets_corrections_;
      edm::InputTag L1ak5Jets_;
      vector < edm::InputTag > L1extras_;
      edm::InputTag Circle8_; 
      edm::InputTag Circle9_; 
      edm::InputTag Circle10_;
      edm::InputTag Circle12_;
      edm::InputTag Square8_; 
      edm::InputTag Square9_; 
      edm::InputTag Square10_;
      edm::InputTag Square12_;
      
      std::string fileName_;

      TFile       *clusterFile_;
      ClusterTree *clusterTree_;
      edm::ParameterSet conf_;

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
NtupleProducer::NtupleProducer(const edm::ParameterSet& iConfig):
conf_(iConfig)
{
   //now do what ever initialization is needed

  //Get in the names from the config file
  
  //calo jets
  CaloJets_ = iConfig.getParameter<edm::InputTag>("CaloJets");
  CaloJets_rho_ = iConfig.getParameter<edm::InputTag>("RhoCaloJets");
//CaloJets_corrections_ = iConfig.getParameter<string>("JetCorrectorCalo");
  //PF jets
  //PFJets_ = iConfig.getParameter<edm::InputTag>("PFJets");
  //PFJets_rho_ = iConfig.getParameter<edm::InputTag>("RhoPFJets");
  //PFJets_corrections_ = iConfig.getParameter<string>("JetCorrectorPF");
  //L1 AK5 jets
  L1ak5Jets_ = iConfig.getParameter<edm::InputTag>("ReRecoPfJets");
  //L1 extra particles: current L1
  L1extras_ = iConfig.getParameter< vector <edm::InputTag> >("extra"); 
  //L1 tower jets
  Circle8_ = iConfig.getParameter<edm::InputTag>("FilteredCircle8");	
  Circle9_ = iConfig.getParameter<edm::InputTag>("FilteredCircle9");	
  Circle10_ = iConfig.getParameter<edm::InputTag>("FilteredCircle10"); 
  Circle12_ = iConfig.getParameter<edm::InputTag>("FilteredCircle12");
  Square8_ = iConfig.getParameter<edm::InputTag>("FilteredSquare8"); 
  Square9_ = iConfig.getParameter<edm::InputTag>("FilteredSquare9");
  Square10_ = iConfig.getParameter<edm::InputTag>("FilteredSquare10");
  Square12_ = iConfig.getParameter<edm::InputTag>("FilteredSquare12");



  edm::Service<TFileService> fs;
  //Filename stuff
  fileName_ = iConfig.getParameter<std::string>("fileName");
  const char * c = fileName_.c_str();
  clusterFile_ = TFile::Open(c,"RECREATE");
  clusterTree_ = new ClusterTree();
  clusterTree_->CreateTree();

 

}


NtupleProducer::~NtupleProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  clusterFile_->cd();
  clusterTree_->tree_->Write();
  clusterFile_->Close();
}



//
// member functions
//

// ------------ method called for each event  ------------
void
NtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVec_;
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVec_PU_;
  bool evValid = true;

  //number of vertices
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel("offlinePrimaryVertices", vtx);
  int PVSIZE=vtx->size();


  //Get the OFFLINE CALO JETS
  edm::Handle<reco::CaloJetCollection> Calojets;
  iEvent.getByLabel(CaloJets_, Calojets);
  if(!Calojets.isValid()){evValid=false;}
  
  //Corrections for offline calo jets
  //const JetCorrectorCalo* jetCorrectorCalo = JetCorrector::getJetCorrector(CaloJets_corrections_,iSetup);
  const JetCorrector* jetCorrectorCalo = JetCorrector::getJetCorrector(conf_.getParameter<string>("JetCorrectorCalo"),iSetup);

  //rho of offline calo jets
  edm::Handle<double> rhoCALOjets;  
  iEvent.getByLabel(CaloJets_rho_, rhoCALOjets);
  double rhoCalo = *rhoCALOjets;
  if(!rhoCALOjets.isValid()){evValid=false;}
  //rho of offline PF jets

  //edm::Handle<double> rhoPFjets;  
  //iEvent.getByLabel(PFJets_rho_, rhoPFjets);
  //double rhoPF = *rhoPFjets;
  //if(!rhoPFjets.isValid()){evValid=false;}

  //
  // set tree variables
  //
  
  clusterTree_->run_ = iEvent.id().run();
  clusterTree_->event_ = iEvent.id().event();
  clusterTree_->lumi_ = iEvent.luminosityBlock();
  clusterTree_->nvtx_ = PVSIZE;
  clusterTree_->rho_Calo_ = rhoCalo;
  //clusterTree_->rho_PF_ = rhoPF;
  clusterTree_->primaryVertex_ = vtx->begin()->position();

  //OFFILNE CALO JETS
  int CaloJet=0;
  
  for(reco::CaloJetCollection::const_iterator it = Calojets->begin(); it != Calojets->end() ; ++it) {
    if(abs(it->eta())<3){
      double scale = jetCorrectorCalo->correction(*it,iEvent,iSetup);
      double scaled_pt = it->pt()*scale;
      //cout<<"Calo jet pt: "<<it->pt()<<" scale factor: "<<scale<<" scaled pt: "<<scaled_pt;
      //cout<<" (eta, phi) = ("<<it->eta()<<","<<it->phi()<<")"<<endl;
      LorentzVec_PU_.SetCoordinates( scaled_pt,it->eta(),it->phi(),0);
      LorentzVec_.SetCoordinates( it->pt(),it->eta(),it->phi(),0);
      clusterTree_->Calojet_indx_->push_back( CaloJet );
      clusterTree_->Calojet4_mom_->push_back( LorentzVec_ );
      clusterTree_->Calojet4_mom_PU_->push_back( LorentzVec_PU_ );
      clusterTree_->Calojet_area_->push_back( it->jetArea() );
     
      
 //     clusterTree_-> Calojet_jet_pt_[CaloJet] = scaled_pt ;      
 //     clusterTree_-> Calojet_jet_eta_[CaloJet] = it->eta() ;
 //     clusterTree_-> Calojet_jet_phi_[CaloJet] = it->phi() ;
       ++CaloJet;
    }
  }



///////////////////////////////////////////////////////////////////////////////
/// Fill the ntuple with tower jets matched to offline calo jets only  

//  edm::Handle<reco::CaloJetCollection> Calojets;
//  iEvent.getByLabel(CaloJets_, Calojets);
//  if(!Calojets.isValid()){evValid=false;}
  
  //Corrections for offline calo jets
  //const JetCorrectorCalo* jetCorrectorCalo = JetCorrector::getJetCorrector(CaloJets_corrections_,iSetup);
 
  //-------------------------------------------------------------------------------------------//
  //Jet collection already filtered

  std::vector< const  EtaPhiContainer<l1slhc::L1TowerJet>* > TowerJets;
  edm::Handle< EtaPhiContainer <l1slhc::L1TowerJet> > hk;
  
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle8"), hk);  if(!hk.isValid()){evValid=false;}	
  TowerJets.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle9"), hk);       if(!hk.isValid()){evValid=false;}	
  TowerJets.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle10"), hk); if(!hk.isValid()){evValid=false;}	
  TowerJets.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle12"), hk);      if(!hk.isValid()){evValid=false;} 
  TowerJets.push_back( &(*hk) );
  
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredSquare8"), hk);  if(!hk.isValid()){evValid=false;}	
  TowerJets.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredSquare9"), hk);       if(!hk.isValid()){evValid=false;}	
  TowerJets.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredSquare10"), hk); if(!hk.isValid()){evValid=false;}	
  TowerJets.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredSquare12"), hk);	if(!hk.isValid()){evValid=false;}	
  TowerJets.push_back( &(*hk) );
  

  if(evValid){}
  else{cout<<"-------------------------> ERROR: Event content corrupted"<<endl;}
  

  double no_tows[]={52,69,80,112,64,81,100,144};
  int Typecount (0);
  TLorentzVector PF, TRPF;
  //Jet collection contains bitonic sorted & filtered 8x8 circ, 9x9 circ, etc; no mean or median yet
  for(std::vector< const  EtaPhiContainer<l1slhc::L1TowerJet>*>::const_iterator Type= TowerJets.begin() ;
      Type!= TowerJets.end() ; ++Type, Typecount++  )
  {
    ///--- Rho calculations---//
    vector<double> JetEnergies;
    for( EtaPhiContainer<l1slhc::L1TowerJet>::const_iterator itr = (**Type).begin(); itr != (**Type).end() ;  ++ itr )
    {
      if( abs(itr->p4().eta() )<3){
        //Evaluate rhos from non scaled L1 jet energies
        JetEnergies.push_back(itr->p4().pt());
      }
    }
    double areaPerJet = no_tows[Typecount] * (0.087 * 0.087) ;
    //option 2
    double median_rho_noscaled = ( median( JetEnergies ) / areaPerJet );
    //scale rho
    double median_BF_rho_noscaled_corr = rho_scale_op2_L1towBF( median_rho_noscaled , Typecount);
    float median_BF_rho_noscaled_Recorr = rho_Rescale_op2_L1towBF( median_BF_rho_noscaled_corr , Typecount);

    //Fill ntuple with scaled rho values
    tow_ntuple_rho( Typecount, median_BF_rho_noscaled_Recorr);

    ///---Match tower jets to calo jets and fill ntuple
    //TLorentzVector PF, TRPF;
    int no_jets=0;
    for(reco::CaloJetCollection::const_iterator it = Calojets->begin(); it != Calojets->end() ; ++it) {
      double scale = jetCorrectorCalo->correction(*it,iEvent,iSetup);
      double scaled_pt = it->pt()*scale;
      if( abs(it->eta())<3){
        PF.SetPtEtaPhiM( scaled_pt,it->eta(),it->phi(),0);
        float maxDist=999;
        float goodPt=-100;
        float goodEta=-100;
        float goodPhi=-100;
        for( EtaPhiContainer<l1slhc::L1TowerJet>::const_iterator itr = (**Type).begin() ;
            itr != (**Type).end() ; ++ itr 
          ) 
        {
          double tow_scaled_pt = itr->p4().pt() - median_BF_rho_noscaled_Recorr*areaPerJet;
          TRPF.SetPtEtaPhiM(tow_scaled_pt,itr->p4().eta(),itr->p4().phi(),0);
          float dr=PF.DeltaR(TRPF);
          if (dr>0.5) continue;
          if (dr<maxDist){//Performs an iteration. maxDist only gets a new value 
            maxDist=dr;
            goodEta=itr->p4().eta(); 
            goodPhi=itr->p4().phi();
            goodPt=tow_scaled_pt;
          }
        }
        //matched jets, fill ntuple

        if (maxDist<998){
          ++no_jets;
          //cout<<"Matched Calo jet: "<<scaled_pt<<","<<it->eta()<<","<<it->phi()<<" no jets: "<<no_jets<<endl;
          //cout<<"Matched tower jet: "<<goodPt<<","<<goodEta<<","<<goodPhi<<endl;

          ///functions to fill ntuple branches
          tow_ntuple_pt( Typecount, goodPt, no_jets);
          Off_matched_tow_ntuple_pt( Typecount, scaled_pt, no_jets);

          tow_ntuple_eta( Typecount, goodEta, no_jets);
          Off_matched_tow_ntuple_phi( Typecount, it->eta(), no_jets);

          tow_ntuple_phi( Typecount, goodPhi, no_jets);
          Off_matched_tow_ntuple_eta( Typecount, it->phi(), no_jets);

        }
      }//eta < 3
    }//calo jets

    tow_ntuple_multi( Typecount, no_jets );

  }//tower jet type
  
  if(evValid){}
  else{cout<<"-------------------------> ERROR: Event content corrupted"<<endl;}

  //Have now read in all the jet collections.

  clusterTree_->tree_->Fill();
  clusterTree_->clearVectors();

}

void NtupleProducer::tow_ntuple_rho( int j, float rho)
{
  if(j==0){     //8x8 circle
    clusterTree_->L1_tow_8circ_scaledrho_ = rho;
  }if(j==1){    //9x9 circle
     clusterTree_->L1_tow_9circ_scaledrho_ = rho;
  }if(j==2){    //10x10 circle
    clusterTree_->L1_tow_10circ_scaledrho_ = rho;
  }if(j==3){    //12x12 circ    
    clusterTree_->L1_tow_12circ_scaledrho_ = rho;
  }  if(j==4){  //8x8 sq  
    clusterTree_->L1_tow_8sq_scaledrho_ = rho;
  }if(j==5){    //9x9 sq
    clusterTree_->L1_tow_9sq_scaledrho_ = rho;
  }if(j==6){    //10x10 sq  
    clusterTree_->L1_tow_10sq_scaledrho_ = rho;
  }if(j==7){    //12x12 sq
    clusterTree_->L1_tow_12sq_scaledrho_ = rho;
  }
}

void NtupleProducer::tow_ntuple_pt( int j, float goodPt, int no_jets)
{
  if(j==0){     //8x8 circle
    clusterTree_->L1_tow_8circ_jet_pt_ [no_jets] = goodPt;
  }if(j==1){    //9x9 circle
    clusterTree_->L1_tow_9circ_jet_pt_ [no_jets] = goodPt;
  }if(j==2){    //10x10 circle
    clusterTree_->L1_tow_10circ_jet_pt_ [no_jets] = goodPt;
  }if(j==3){    //12x12 circ    
    clusterTree_->L1_tow_12circ_jet_pt_ [no_jets] = goodPt;
  }  if(j==4){  //8x8 sq  
    clusterTree_->L1_tow_8sq_jet_pt_ [no_jets] = goodPt;
  }if(j==5){    //9x9 sq
    clusterTree_->L1_tow_9sq_jet_pt_ [no_jets] = goodPt;
  }if(j==6){    //10x10 sq  
    clusterTree_->L1_tow_10sq_jet_pt_ [no_jets] = goodPt;
  }if(j==7){    //12x12 sq
    clusterTree_->L1_tow_12sq_jet_pt_ [no_jets] = goodPt;
  }
}

void NtupleProducer::Off_matched_tow_ntuple_pt( int j, float scaled_pt, int no_jets)
{
  if(j==0){     //8x8 circle
    clusterTree_->Calo_match8circ_jet_pt_ [no_jets] = scaled_pt;
  }if(j==1){    //9x9 circle
    clusterTree_->Calo_match9circ_jet_pt_ [no_jets] = scaled_pt;
  }if(j==2){    //10x10 circle
    clusterTree_->Calo_match10circ_jet_pt_ [no_jets] = scaled_pt;
  }if(j==3){    //12x12 circ    
    clusterTree_->Calo_match12circ_jet_pt_ [no_jets] = scaled_pt;
  }  if(j==4){  //8x8 sq  
    clusterTree_->Calo_match8sq_jet_pt_ [no_jets] = scaled_pt;
  }if(j==5){    //9x9 sq
    clusterTree_->Calo_match9sq_jet_pt_ [no_jets] = scaled_pt;
  }if(j==6){    //10x10 sq  
    clusterTree_->Calo_match10sq_jet_pt_ [no_jets] = scaled_pt;
  }if(j==7){    //12x12 sq
    clusterTree_->Calo_match12sq_jet_pt_ [no_jets] = scaled_pt;
  }
}

void NtupleProducer::tow_ntuple_eta(int j, float goodEta, int no_jets)
{
  if(j==0){     //8x8 circle
    clusterTree_->L1_tow_8circ_jet_eta_ [no_jets] = goodEta;
  }if(j==1){    //9x9 circle
    clusterTree_->L1_tow_9circ_jet_eta_ [no_jets] = goodEta;
  }if(j==2){    //10x10 circle
    clusterTree_->L1_tow_10circ_jet_eta_ [no_jets] = goodEta;
  }if(j==3){    //12x12 circ    
    clusterTree_->L1_tow_12circ_jet_eta_ [no_jets] = goodEta;
  }  if(j==4){  //8x8 sq  
    clusterTree_->L1_tow_8sq_jet_eta_ [no_jets] = goodEta;
  }if(j==5){    //9x9 sq
    clusterTree_->L1_tow_9sq_jet_eta_ [no_jets] = goodEta;
  }if(j==6){    //10x10 sq  
    clusterTree_->L1_tow_10sq_jet_eta_ [no_jets] = goodEta;
  }if(j==7){    //12x12 sq
    clusterTree_->L1_tow_12sq_jet_eta_ [no_jets] = goodEta;
  }
}

void NtupleProducer::Off_matched_tow_ntuple_eta( int j, float eta, int no_jets)
{
  if(j==0){     //8x8 circle
    clusterTree_->Calo_match8circ_jet_eta_ [no_jets] = eta;
  }if(j==1){    //9x9 circle
    clusterTree_->Calo_match9circ_jet_eta_ [no_jets] = eta;
  }if(j==2){    //10x10 circle
    clusterTree_->Calo_match10circ_jet_eta_ [no_jets] = eta;
  }if(j==3){    //12x12 circ    
    clusterTree_->Calo_match12circ_jet_eta_ [no_jets] = eta;
  }  if(j==4){  //8x8 sq  
    clusterTree_->Calo_match8sq_jet_eta_ [no_jets] = eta;
  }if(j==5){    //9x9 sq
    clusterTree_->Calo_match9sq_jet_eta_ [no_jets] = eta;
  }if(j==6){    //10x10 sq  
    clusterTree_->Calo_match10sq_jet_eta_ [no_jets] = eta;
  }if(j==7){    //12x12 sq
    clusterTree_->Calo_match12sq_jet_eta_ [no_jets] = eta;
  }
}

void NtupleProducer::tow_ntuple_phi( int j, float goodPhi, int no_jets)
{
  if(j==0){     //8x8 circle
    clusterTree_->L1_tow_8circ_jet_phi_ [no_jets] = goodPhi;
  }if(j==1){    //9x9 circle
    clusterTree_->L1_tow_9circ_jet_phi_ [no_jets] = goodPhi;
  }if(j==2){    //10x10 circle
    clusterTree_->L1_tow_10circ_jet_phi_ [no_jets] = goodPhi;
  }if(j==3){    //12x12 circ    
    clusterTree_->L1_tow_12circ_jet_phi_ [no_jets] = goodPhi;
  }  if(j==4){  //8x8 sq  
    clusterTree_->L1_tow_8sq_jet_phi_ [no_jets] = goodPhi;
  }if(j==5){    //9x9 sq
    clusterTree_->L1_tow_9sq_jet_phi_ [no_jets] = goodPhi;
  }if(j==6){    //10x10 sq  
    clusterTree_->L1_tow_10sq_jet_phi_ [no_jets] = goodPhi;
  }if(j==7){    //12x12 sq
    clusterTree_->L1_tow_12sq_jet_phi_ [no_jets] = goodPhi;
  }
}

void NtupleProducer::Off_matched_tow_ntuple_phi( int j, float phi, int no_jets)
{
  if(j==0){     //8x8 circle
    clusterTree_->Calo_match8circ_jet_phi_ [no_jets] = phi;
  }if(j==1){    //9x9 circle
    clusterTree_->Calo_match9circ_jet_phi_ [no_jets] = phi;
  }if(j==2){    //10x10 circle
    clusterTree_->Calo_match10circ_jet_phi_ [no_jets] = phi;
  }if(j==3){    //12x12 circ    
    clusterTree_->Calo_match12circ_jet_phi_ [no_jets] = phi;
  }  if(j==4){  //8x8 sq  
    clusterTree_->Calo_match8sq_jet_phi_ [no_jets] = phi;
  }if(j==5){    //9x9 sq
    clusterTree_->Calo_match9sq_jet_phi_ [no_jets] = phi;
  }if(j==6){    //10x10 sq  
    clusterTree_->Calo_match10sq_jet_phi_ [no_jets] = phi;
  }if(j==7){    //12x12 sq
    clusterTree_->Calo_match12sq_jet_phi_ [no_jets] = phi;
  }
}


void NtupleProducer::tow_ntuple_multi( int j, int no_jets)
{
  if(j==0){     //8x8 circle
    clusterTree_->L1_tow_8circ_jet_multi_= no_jets;
    clusterTree_->Calo_match8circ_jet_multi_= no_jets;
  }if(j==1){    //9x9 circle
    clusterTree_->L1_tow_9circ_jet_multi_= no_jets;
    clusterTree_->Calo_match9circ_jet_multi_= no_jets;
  }if(j==2){    //10x10 circle
    clusterTree_->L1_tow_10circ_jet_multi_= no_jets;
    clusterTree_->Calo_match10circ_jet_multi_= no_jets;
  }if(j==3){    //12x12 circ    
    clusterTree_->L1_tow_12circ_jet_multi_= no_jets;
    clusterTree_->Calo_match12circ_jet_multi_= no_jets;
  }  if(j==4){  //8x8 sq  
    clusterTree_->L1_tow_8sq_jet_multi_= no_jets;    
    clusterTree_->Calo_match8sq_jet_multi_= no_jets;
  }if(j==5){    //9x9 sq
    clusterTree_->L1_tow_9sq_jet_multi_= no_jets;    
    clusterTree_->Calo_match9sq_jet_multi_= no_jets;
  }if(j==6){    //10x10 sq  
    clusterTree_->L1_tow_10sq_jet_multi_= no_jets;    
    clusterTree_->Calo_match10sq_jet_multi_= no_jets;
  }if(j==7){    //12x12 sq
    clusterTree_->L1_tow_12sq_jet_multi_= no_jets;
    clusterTree_->Calo_match12sq_jet_multi_= no_jets;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
NtupleProducer::beginJob()
{
 
} 
  
// ------------ method called once each job just after ending the event loop  ------------
void 
NtupleProducer::endJob() 
{ 
} 
  
// ------------ method called when starting to processes a run  ------------
void 
NtupleProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
NtupleProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
NtupleProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
NtupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleProducer);
