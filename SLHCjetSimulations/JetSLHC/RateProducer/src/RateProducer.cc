// -*- C++ -*-
//
// Package:    RateProducer
// Class:      RateProducer
// 
/**\class RateProducer RateProducer.cc JetSLHC/RateProducer/src/RateProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michele Pioppi,510 1-017,+41227671792,
//         Created:  Fri Jul 27 15:20:08 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloJet.h"
#include "SimDataFormats/SLHC/interface/L1CaloJetFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
bool myfunction (float i,float j) { return (i>j); }
//
// class declaration
//
using namespace l1slhc;
using namespace std;
using namespace edm;
using namespace reco;
class RateProducer : public edm::EDAnalyzer {
   public:
      explicit RateProducer(const edm::ParameterSet&);
      ~RateProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
 ParameterSet conf_;
 TH2F *ICSingleJet;
 TH2F *ICDoubleJet;
 TH2F *ICTripleJet;
 TH2F *ICQuadJet;
TH2F *UWSingleJet;
 TH2F *UWDoubleJet;
 TH2F *UWTripleJet;
 TH2F *UWQuadJet;
TH2F *CUSingleJet;
 TH2F *CUDoubleJet;
 TH2F *CUTripleJet;
 TH2F *CUQuadJet;
 TH2F *ICDeltaEta;
TH2F *ICInvMass;
 TH2F *CUDeltaEta; 
TH2F *CUInvMass; 
 TH2F *UWDeltaEta; 
TH2F *UWInvMass; 
      // ----------member data ---------------------------
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
RateProducer::RateProducer(const edm::ParameterSet& iConfig):
  conf_(iConfig)

{
   //now do what ever initialization is needed
  Service<TFileService> fs;
  ICDeltaEta   = fs->make<TH2F>("ICDeltaEta", 
			      "ICDeltaEta", 51,-0.1,5.,20,0,100);
  ICInvMass = fs->make<TH2F>("ICInvMass","ICInvMass",501,-1,500,20,0,100);
  ICSingleJet = fs->make<TH2F>("ICSingleJet","ICSingleJet",101,-10,1000,20,0,100);
  ICDoubleJet = fs->make<TH2F>("ICDoubleJet","ICDoubleJet",101,-10,1000,20,0,100);
  ICTripleJet = fs->make<TH2F>("ICTripleJet","ICTripleJet",101,-10,1000,20,0,100);
  ICQuadJet = fs->make<TH2F>("ICQuadJet","ICQuadJet",101,-10,1000,20,0,100);

  CUDeltaEta   = fs->make<TH2F>("CUDeltaEta",
                              "CUDeltaEta", 51,-0.1,5.,20,0,100);
  CUInvMass = fs->make<TH2F>("CUInvMass","CUInvMass",501,-1,500,20,0,100);
  CUSingleJet = fs->make<TH2F>("CUSingleJet","CUSingleJet",101,-10,1000,20,0,100);
  CUDoubleJet = fs->make<TH2F>("CUDoubleJet","CUDoubleJet",101,-10,1000,20,0,100);
  CUTripleJet = fs->make<TH2F>("CUTripleJet","CUTripleJet",101,-10,1000,20,0,100);
  CUQuadJet = fs->make<TH2F>("CUQuadJet","CUQuadJet",101,-10,1000,20,0,100);	 


  UWDeltaEta   = fs->make<TH2F>("UWDeltaEta",
                              "UWDeltaEta", 51,-0.1,5.,20,0,100);
  UWInvMass = fs->make<TH2F>("UWInvMass","UWInvMass",501,-1,500,20,0,100);
  UWSingleJet = fs->make<TH2F>("UWSingleJet","UWSingleJet",101,-10,1000,20,0,100);
  UWDoubleJet = fs->make<TH2F>("UWDoubleJet","UWDoubleJet",101,-10,1000,20,0,100);
  UWTripleJet = fs->make<TH2F>("UWTripleJet","UWTripleJet",101,-10,1000,20,0,100);
  UWQuadJet = fs->make<TH2F>("UWQuadJet","UWQuadJet",101,-10,1000,20,0,100);

}


RateProducer::~RateProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RateProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel("offlinePrimaryVertices", vtx);
  float nv=float(vtx->size());



    //ASSUME HERE THAT JETS ARE ALREADY CALIBRATED AND ORDERED IN PT
    Handle<L1TowerJetCollection > hk;
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle8"), hk);
    if (hk->size()<1) ICSingleJet->Fill(-5.,nv);
    if (hk->size()<2) ICDoubleJet->Fill(-5.,nv);
    if (hk->size()<3) ICTripleJet->Fill(-5.,nv);
    if (hk->size()<4) ICQuadJet->Fill(-5.,nv);
    int ijet=0;
    for (L1TowerJetCollection::const_iterator il1 =hk->begin();
	 il1!=hk->end();++il1){
         ijet++;
         if (ijet==1) ICSingleJet->Fill(il1->p4().Pt(),nv);
        if (ijet==2) ICDoubleJet->Fill(il1->p4().Pt(),nv);
        if (ijet==3) ICTripleJet->Fill(il1->p4().Pt(),nv);
        if (ijet==4) ICQuadJet->Fill(il1->p4().Pt(),nv);
         }


    Handle<L1CaloJetCollection > uw;
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UWJets"), uw);
    if (uw->size()<1) UWSingleJet->Fill(-5.,nv);
    if (uw->size()<2) UWDoubleJet->Fill(-5.,nv);
    if (uw->size()<3) UWTripleJet->Fill(-5.,nv);
    if (uw->size()<4) UWQuadJet->Fill(-5.,nv);
    ijet=0;
    for (L1CaloJetCollection::const_iterator il1 =uw->begin();
	 il1!=uw->end();++il1){
         ijet++;
         if (ijet==1) UWSingleJet->Fill(il1->p4().Pt(),nv);
        if (ijet==2) UWDoubleJet->Fill(il1->p4().Pt(),nv);
        if (ijet==3) UWTripleJet->Fill(il1->p4().Pt(),nv);
        if (ijet==4) UWQuadJet->Fill(il1->p4().Pt(),nv);
         }

vector<float> ptcurrl1;
    vector<InputTag>  l1extraparticles= conf_.getParameter< vector < InputTag > >("extra"); 
    for (uint i=0; i< l1extraparticles.size();++i){
      Handle<l1extra::L1JetParticleCollection> l1Jets;
      iEvent.getByLabel(l1extraparticles[i],l1Jets);
 for(l1extra::L1JetParticleCollection::const_iterator itr = l1Jets->begin();  itr != l1Jets->end(); ++itr ) {
   ptcurrl1.push_back(itr->p4().Pt());
}}
sort(ptcurrl1.begin(),ptcurrl1.end(),myfunction);
    if (ptcurrl1.size()<1) CUSingleJet->Fill(-5.,nv); else CUSingleJet->Fill(ptcurrl1[0],nv);
    if (ptcurrl1.size()<2) CUDoubleJet->Fill(-5.,nv); else CUDoubleJet->Fill(ptcurrl1[1],nv);
    if (ptcurrl1.size()<3) CUTripleJet->Fill(-5.,nv); else CUTripleJet->Fill(ptcurrl1[2],nv);
    if (ptcurrl1.size()<4) CUQuadJet->Fill(-5.,nv); else CUQuadJet->Fill(ptcurrl1[3],nv);
}


// ------------ method called once each job just before starting event loop  ------------
void 
RateProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RateProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
RateProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
RateProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
RateProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
RateProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RateProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RateProducer);
