// -*- C++ -*-
//
// Package:    VBFHBBEffAnalyzer
// Class:      VBFHBBEffAnalyzer
// 
/**\class VBFHBBEffAnalyzer VBFHBBEffAnalyzer.cc JetSLHC/VBFHBBEffAnalyzer/src/VBFHBBEffAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michele Pioppi,510 1-017,+41227671792,
//         Created:  Wed Jul 25 12:23:09 CEST 2012
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloJet.h"
#include "SimDataFormats/SLHC/interface/L1CaloJetFwd.h"
//
// class declaration
//
using namespace l1slhc;
using namespace std;
using namespace edm;
using namespace reco;
class VBFHBBEffAnalyzer : public EDAnalyzer {
   public:
      explicit VBFHBBEffAnalyzer(const ParameterSet&);
      ~VBFHBBEffAnalyzer();

      static void fillDescriptions(ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const Event&, const EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(Run const&, EventSetup const&);
      virtual void endRun(Run const&, EventSetup const&);
      virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&);
      virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&);
  ParameterSet conf_;
  TH1F* DeltaEta;
  TH1F* CurrDeltaEta;
  TH1F* UWDeltaEta;

  TH1F* InvMassL1;
  TH1F* InvMassCurrL1;
  TH1F* InvMassUWL1;
  bool bjets,vbf;
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
VBFHBBEffAnalyzer::VBFHBBEffAnalyzer(const ParameterSet& iConfig):
  conf_(iConfig)

{   //now do what ever initialization is needed
  Service<TFileService> fs;
  DeltaEta   = fs->make<TH1F>("L1_DeltaEta", 
			      "L1_DeltaEta", 51,-0.1,5.);
  CurrDeltaEta   = fs->make<TH1F>("CurrL1_DeltaEta", 
			      "CurrL1_DeltaEta", 51,-0.1,5.);
  UWDeltaEta   = fs->make<TH1F>("UWL1_DeltaEta", 
			      "UWL1_DeltaEta", 51,-0.1,5.);
  InvMassL1 = fs->make<TH1F>("InvMassL1","InvMassL1",501,-1,500);
  InvMassCurrL1 = fs->make<TH1F>("InvMassCurrL1","InvMassCurrL1",501,-1,500);
  InvMassUWL1 = fs->make<TH1F>("InvMassUWL1","InvMassUWL1",501,-1,500);
  bjets=conf_.getParameter<bool>("JetsFromHiggsDecay");
  vbf=conf_.getParameter<bool>("JetsFromVBF");
}


VBFHBBEffAnalyzer::~VBFHBBEffAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  

}


//
// member functions
//

// ------------ method called for each event  ------------
void
VBFHBBEffAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
{

  Handle<GenParticleCollection> GEN;
  iEvent.getByLabel("genParticles",GEN);
  //  cout<<"NUOVO EVENTO "<<endl;
  bool isb1ok=false;
  bool isb2ok=false;
  TLorentzVector genb1,genb2;
  
  bool isv1ok=false;
  bool isv2ok=false;
  TLorentzVector genv1,genv2;
  int ig=0;
  for(GenParticleCollection::const_iterator it = GEN->begin(); it != GEN->end(); ++it) {
    ig++;
    //BJET SELECTION FROM H-->bb
    if (bjets){
      if (fabs(it->pdgId())!=5) continue;
      if (it->motherRefVector().size()==0) continue;
      if ((*it->motherRefVector()[0]).pdgId()!=25) continue;
      if (fabs(it->momentum().eta())>3.) continue;
      if (it->pdgId()==5){
	isb1ok=true;
	//	cout<<"GENB1 "<<it->momentum().eta()<<" "<<it->momentum().phi()<<endl;
	genb1.SetPxPyPzE(it->momentum().x(),it->momentum().y(),
			 it->momentum().z(),it->energy());
      }
      if (it->pdgId()==-5){
	isb2ok=true;
	//	cout<<"GENB2 "<<it->momentum().eta()<<" "<<it->momentum().phi()<<endl;
	genb2.SetPxPyPzE(it->momentum().x(),it->momentum().y(),
			 it->momentum().z(),it->energy());
      } 
    }
    if (vbf && ig<30){
      bool isjetok=false;
      if (fabs(it->pdgId())>4) continue;
      if (it->motherRefVector().size()==0) continue;
      if ((*it->motherRefVector()[0]).daughterRefVector().size()==0) continue;
      cout<<"PDG "<<it->pdgId()<<" "<<it->motherRefVector().size()<<" "<<(*it->motherRefVector()[0]).daughterRefVector().size()<<endl;
      for (uint ida=0;ida< (*it->motherRefVector()[0]).daughterRefVector().size();ida++){
	if (isjetok)continue;
	if ((*it->motherRefVector()[0]).daughterRefVector()[ida]->pdgId()==25)
	  isjetok=true;
      }
      cout<<" "<<isjetok<<"pt "<<it->pt()<<endl;	
      if (isjetok && !isv1ok && !isv2ok){
	cout<<"DENTRO 1 "<<endl;
	isv1ok=true;
	genv1.SetPxPyPzE(it->momentum().x(),it->momentum().y(),
			 it->momentum().z(),it->energy());
      }else if (isjetok && isv1ok && !isv2ok){
	cout<<"DENTRO 2 "<<endl;
	isv2ok=true;
	genv2.SetPxPyPzE(it->momentum().x(),it->momentum().y(),
			 it->momentum().z(),it->energy());
	
      }
    }
  }

  if (bjets &&isb1ok &&isb2ok){
    cout<<"HIGGS MASS "<<(genb1+genb2).M()<<endl;
    TLorentzVector L1B,L2B;
    bool isl1b1ok=false;
    bool isl1b2ok=false;
    //JETS HERE MUST BE ALREADY CALIBRATED AND PU SUBTRACTED
    Handle<L1TowerJetCollection > hk;
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle8"), hk);
    for (L1TowerJetCollection::const_iterator il1 =hk->begin();
	 il1!=hk->end();++il1){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(il1->p4().Pt(),il1->p4().eta(),il1->p4().phi(),il1->p4().M());
      cout<<"ETA,phi "<<il1->p4().eta()<<" "<<il1->p4().phi()<<" "<<tmp.DeltaR(genb1)<<" "<<tmp.DeltaR(genb2)<<endl;
      if (!isl1b1ok && tmp.DeltaR(genb1)<0.4){
	isl1b1ok=true;
	L1B=tmp;
      }
      if (!isl1b2ok && tmp.DeltaR(genb2)<0.4){
	isl1b2ok=true;
	L2B=tmp;
      }  
    }
    cout<<"IC MMM"<<isl1b1ok<<isl1b2ok<<endl;
    if (isl1b1ok && isl1b2ok) InvMassL1->Fill((L1B+L2B).M());
    else InvMassL1->Fill(-0.5);
      

    // DO THE SAME FOR CURRENT L1 JETS

    TLorentzVector CurrL1B,CurrL2B;
    bool iscurrl1b1ok=false;
    bool iscurrl1b2ok=false;
    vector<InputTag>  l1extraparticles= conf_.getParameter< vector < InputTag > >("extra"); 
    for (uint i=0; i< l1extraparticles.size();++i){
      Handle<l1extra::L1JetParticleCollection> l1Jets;
      iEvent.getByLabel(l1extraparticles[i],l1Jets);
      cout<<"SS "<<l1extraparticles[i]<<" "<<l1Jets->size()<<endl;
      for(l1extra::L1JetParticleCollection::const_iterator itr = l1Jets->begin();  itr != l1Jets->end(); ++itr ) {
	TLorentzVector tmp;
	tmp.SetPtEtaPhiM(itr->p4().Pt(),itr->p4().eta(),itr->p4().phi(),itr->p4().M());
	cout<<"ETA,phi "<<itr->p4().eta()<<" "<<itr->p4().phi()<<" "<<tmp.DeltaR(genb1)<<" "<<tmp.DeltaR(genb2)<<endl;
	if (!iscurrl1b1ok && tmp.DeltaR(genb1)<0.4){
	  iscurrl1b1ok=true;
	  CurrL1B=tmp;
	}
	if (!iscurrl1b2ok && tmp.DeltaR(genb2)<0.4){
	  iscurrl1b2ok=true;
	  CurrL2B=tmp;
	}  
      }
    }
    cout<<"CU MMM"<<iscurrl1b1ok<<iscurrl1b2ok<<endl;
    if (iscurrl1b1ok && iscurrl1b2ok) InvMassCurrL1->Fill((CurrL1B+CurrL2B).M());
    else InvMassCurrL1->Fill(-0.5);
  
    // DO THE SAME FOR  UW JETS
    TLorentzVector UWL1B,UWL2B;
    bool isuwl1b1ok=false;
    bool isuwl1b2ok=false;
    //JETS HERE MUST BE ALREADY CALIBRATED AND PU SUBTRACTED
    Handle<L1CaloJetCollection > uw;
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UWJets"), uw);
    for (L1CaloJetCollection::const_iterator il1 =uw->begin();
	 il1!=uw->end();++il1){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(il1->p4().Pt(),il1->p4().eta(),il1->p4().phi(),il1->p4().M());
      cout<<"ETA,phi "<<il1->p4().eta()<<" "<<il1->p4().phi()<<" "<<tmp.DeltaR(genb1)<<" "<<tmp.DeltaR(genb2)<<endl;
      if (!isuwl1b1ok && tmp.DeltaR(genb1)<0.4){
	isuwl1b1ok=true;
	UWL1B=tmp;
      }
      if (!isuwl1b2ok && tmp.DeltaR(genb2)<0.4){
	isuwl1b2ok=true;
	UWL2B=tmp;
      }  
    }
    cout<<"UW MMM"<<isuwl1b1ok<<isuwl1b2ok<<endl;
    if (isuwl1b1ok && isuwl1b2ok) InvMassUWL1->Fill((UWL1B+UWL2B).M());
    else InvMassUWL1->Fill(-0.5);
  }
  //VBF
  if (vbf &&isv1ok &&isv2ok){
    cout<<"Delta ETA "<<fabs(genv1.Eta()  - genv2.Eta())<<endl;
    TLorentzVector L1B,L2B;
    bool isl1b1ok=false;
    bool isl1b2ok=false;
    //JETS HERE MUST BE ALREADY CALIBRATED AND PU SUBTRACTED
    Handle<L1TowerJetCollection > hk;
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle8"), hk);
    for (L1TowerJetCollection::const_iterator il1 =hk->begin();
	 il1!=hk->end();++il1){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(il1->p4().Pt(),il1->p4().eta(),il1->p4().phi(),il1->p4().M());
      cout<<"ETA,phi "<<il1->p4().eta()<<" "<<il1->p4().phi()<<" "<<tmp.DeltaR(genv1)<<" "<<tmp.DeltaR(genv2)<<endl;
      if (!isl1b1ok && tmp.DeltaR(genv1)<0.4){
	isl1b1ok=true;
	L1B=tmp;
      }
      if (!isl1b2ok && tmp.DeltaR(genv2)<0.4){
	isl1b2ok=true;
	L2B=tmp;
      }  
    }
    cout<<"IC ETA"<<isl1b1ok<<isl1b2ok<<endl;
    if (isl1b1ok && isl1b2ok)  DeltaEta->Fill(fabs(L1B.Eta()-L2B.Eta()));
    else DeltaEta->Fill(-0.05); 
      

    // DO THE SAME FOR CURRENT L1 JETS

    TLorentzVector CurrL1B,CurrL2B;
    bool iscurrl1b1ok=false;
    bool iscurrl1b2ok=false;
    vector<InputTag>  l1extraparticles= conf_.getParameter< vector < InputTag > >("extra"); 
    for (uint i=0; i< l1extraparticles.size();++i){
      Handle<l1extra::L1JetParticleCollection> l1Jets;
      iEvent.getByLabel(l1extraparticles[i],l1Jets);
      cout<<"SS "<<l1extraparticles[i]<<" "<<l1Jets->size()<<endl;
      for(l1extra::L1JetParticleCollection::const_iterator itr = l1Jets->begin();  itr != l1Jets->end(); ++itr ) {
	TLorentzVector tmp;
	tmp.SetPtEtaPhiM(itr->p4().Pt(),itr->p4().eta(),itr->p4().phi(),itr->p4().M());
	cout<<"ETA,phi "<<itr->p4().eta()<<" "<<itr->p4().phi()<<" "<<tmp.DeltaR(genv1)<<" "<<tmp.DeltaR(genv2)<<endl;
	if (!iscurrl1b1ok && tmp.DeltaR(genv1)<0.4){
	  iscurrl1b1ok=true;
	  CurrL1B=tmp;
	}
	if (!iscurrl1b2ok && tmp.DeltaR(genv2)<0.4){
	  iscurrl1b2ok=true;
	  CurrL2B=tmp;
	}  
      }
    }
    cout<<"CU ETA"<<iscurrl1b1ok<<iscurrl1b2ok<<endl;
    if (iscurrl1b1ok && iscurrl1b2ok)  CurrDeltaEta->Fill(fabs(CurrL1B.Eta()-CurrL2B.Eta()));
    else CurrDeltaEta->Fill(-0.05); 
    // DO THE SAME FOR  UW JETS
    TLorentzVector UWL1B,UWL2B;
    bool isuwl1b1ok=false;
    bool isuwl1b2ok=false;
    //JETS HERE MUST BE ALREADY CALIBRATED AND PU SUBTRACTED
    Handle<L1CaloJetCollection > uw;
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("UWJets"), uw);
    for (L1CaloJetCollection::const_iterator il1 =uw->begin();
	 il1!=uw->end();++il1){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(il1->p4().Pt(),il1->p4().eta(),il1->p4().phi(),il1->p4().M());
      cout<<"ETA,phi "<<il1->p4().eta()<<" "<<il1->p4().phi()<<" "<<tmp.DeltaR(genv1)<<" "<<tmp.DeltaR(genv2)<<endl;
      if (!isuwl1b1ok && tmp.DeltaR(genv1)<0.4){
	isuwl1b1ok=true;
	UWL1B=tmp;
      }
      if (!isuwl1b2ok && tmp.DeltaR(genv2)<0.4){
	isuwl1b2ok=true;
	UWL2B=tmp;
      }  
    }
    cout<<"UW ETA"<<isuwl1b1ok<<isuwl1b2ok<<endl;
    if (isuwl1b1ok && isuwl1b2ok)  UWDeltaEta->Fill(fabs(UWL1B.Eta()-UWL2B.Eta()));
    else UWDeltaEta->Fill(-0.05);
  }


}

// ------------ method called once each job just before starting event loop  ------------
void 
VBFHBBEffAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHBBEffAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
VBFHBBEffAnalyzer::beginRun(Run const&, EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
VBFHBBEffAnalyzer::endRun(Run const&, EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
VBFHBBEffAnalyzer::beginLuminosityBlock(LuminosityBlock const&, EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
VBFHBBEffAnalyzer::endLuminosityBlock(LuminosityBlock const&, EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VBFHBBEffAnalyzer::fillDescriptions(ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VBFHBBEffAnalyzer);
