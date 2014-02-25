#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "DataFormats/JetReco/interface/JetID.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "DataFormats/Common/interface/ValueMap.h"

//
// class declaration
//

using namespace l1slhc;
using namespace edm;
using namespace std;
using namespace reco;

bool myfunction (float i,float j) { return (i>j); }

bool sortTLorentz (TLorentzVector i,TLorentzVector j) { return ( i.Pt()>j.Pt() ); }


class AnalyzeJets : public edm::EDAnalyzer {
   public:
      explicit AnalyzeJets(const edm::ParameterSet&);
      ~AnalyzeJets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      //------------member functions---------------------
      std::pair< double, int> Leading_Match(vector<TLorentzVector> offlineJets,vector<TLorentzVector> L1Jets );

      // ----------member data ---------------------------
      ParameterSet conf_;
    edm::InputTag m_jetCollection;
    edm::InputTag m_jetIDMap;
			float oPt_;
			float oPhi_;
			float oEta_;
			float l1Pt_;
			float l1Phi_;
			float l1Eta_;

      //Jet Rates
      TH1D *upgradeJets_hist[4];
      TH1D *upgradeJets_rate[4];

      TH1D *CurrL1jets_hist[4];
      TH1D *CurrL1jets_rate[4];

      //HT, MHT Rates
      TH1D *UpgradeHT_hist;
      TH1D *UpgradeHT_rate;
      TH1D *UpgradeMHT_hist;
      TH1D *UpgradeMHT_rate;

      TH1D *CurrL1HT_hist;
      TH1D *CurrL1HT_rate;
      TH1D *CurrL1MHT_hist;
      TH1D *CurrL1MHT_rate;

			//Turn on Curves
      TH3F*    OfflinePt3DIC40  ;
      TH3F*    OfflinePt3DIC45  ;
      TH3F*    OfflinePt3DIC50  ;
      TH3F*    OfflinePt3DIC60  ;
      TH3F*    OfflinePt3DIC80  ;
      TH3F*    OfflinePt3DIC100 ;
      TH3F*    OfflinePt3DIC120 ;
      TH3F*    OfflinePt3DIC140 ;
      TH3F*    OfflinePt3DICDen ;


      TH3F*    OfflinePt3DCurr40  ;
      TH3F*    OfflinePt3DCurr45  ;
      TH3F*    OfflinePt3DCurr50  ;
      TH3F*    OfflinePt3DCurr60  ;
      TH3F*    OfflinePt3DCurr80  ;
      TH3F*    OfflinePt3DCurr100 ;
      TH3F*    OfflinePt3DCurr120 ;
      TH3F*    OfflinePt3DCurr140 ;
      TH3F*    OfflinePt3DCurrDen ;


			//Resolutions
			TH1F*			DeltaEta          ; 
			TH1F*			CurrDeltaEta      ;
			TH1F*			DeltaPhi          ;
			TH1F*			CurrDeltaPhi      ;
			TH1F*			DeltaPt           ;
			TH1F*			CurrDeltaPt       ;
			TH2F*			OffPt_upL1Pt      ;
			TH2F*			OffPt_currL1Pt    ;
			TH2F*			DeltaPt_pt        ;
			TH2F*			CurrDeltaPt_pt    ; 
			TH2F*			DeltaPt_nvtx      ;
			TH2F*			CurrDeltaPt_nvtx  ;
			TH2F*			DeltaPt_eta       ;
			TH2F*			CurrDeltaPt_eta   ;

			TH1F*    NVTX;

};


AnalyzeJets::AnalyzeJets(const edm::ParameterSet& iConfig):
	conf_(iConfig)
{
  Service<TFileService> fs;

	//Rates: Jets
  upgradeJets_hist[0] = fs->make <TH1D>("upgradeJets_1_distro" , "upgradeJets_1_distro", 400, 0., 400.);
  upgradeJets_hist[1] = fs->make <TH1D>("upgradeJets_2_distro" , "upgradeJets_2_distro", 400, 0., 400.);
  upgradeJets_hist[2] = fs->make <TH1D>("upgradeJets_3_distro" , "upgradeJets_3_distro", 400, 0., 400.);
  upgradeJets_hist[3] = fs->make <TH1D>("upgradeJets_4_distro" , "upgradeJets_4_distro", 400, 0., 400.);
  upgradeJets_rate[0] = fs->make <TH1D>("UpgradeSingleJet_rate", "Single Jet"          , 400, 0., 400.);
  upgradeJets_rate[1] = fs->make <TH1D>("UpgradeDoubleJet_rate", "Double Jet"          , 400, 0., 400.);
  upgradeJets_rate[2] = fs->make <TH1D>("UpgradeTripleJet_rate", "Triple Jet"          , 400, 0., 400.);
  upgradeJets_rate[3] = fs->make <TH1D>("UpgradeQuadJet_rate"  , "Quad Jet"            , 400, 0., 400.);

  CurrL1jets_hist[0] = fs->make <TH1D>("CurrL1jets_1_distro"   , "CurrL1jets_1_distro", 400, 0., 400.);
  CurrL1jets_hist[1] = fs->make <TH1D>("CurrL1jets_2_distro"   , "CurrL1jets_2_distro", 400, 0., 400.);
  CurrL1jets_hist[2] = fs->make <TH1D>("CurrL1jets_3_distro"   , "CurrL1jets_3_distro", 400, 0., 400.);
  CurrL1jets_hist[3] = fs->make <TH1D>("CurrL1jets_4_distro"   , "CurrL1jets_4_distro", 400, 0., 400.);
  CurrL1jets_rate[0] = fs->make <TH1D>("CurrL1SingleJet_rate"  , "Single Jet"         , 400, 0., 400.);
  CurrL1jets_rate[1] = fs->make <TH1D>("CurrL1DoubleJet_rate"  , "Double Jet"         , 400, 0., 400.);
  CurrL1jets_rate[2] = fs->make <TH1D>("CurrL1TripleJet_rate"  , "Triple Jet"         , 400, 0., 400.);
  CurrL1jets_rate[3] = fs->make <TH1D>("CurrL1QuadJet_rate"    , "Quad Jet"           , 400, 0., 400.);


	//Rates: MHT, HT
  UpgradeHT_hist = fs->make <TH1D>("UpgradeHT_distro"  , "Upgrade HT distro"  , 400, 0., 400.);
  UpgradeHT_rate = fs->make <TH1D>("UpgradeHT_rate"    , "Upgrade HT rate"    , 400, 0., 400.);
  UpgradeMHT_hist= fs->make <TH1D>("UpgradeMHT_distro" , "Upgrade MHT distro" , 400, 0., 400.);
  UpgradeMHT_rate= fs->make <TH1D>("UpgradeMHT_rate"   , "Upgrade MHT rate"   , 400, 0., 400.);

  CurrL1HT_hist = fs->make <TH1D>("CurrL1HT_distro" , "Curr HT distro"  ,400, 0., 400. );
  CurrL1HT_rate = fs->make <TH1D>("CurrL1HT_rate"   , "Curr HT rate"    ,400, 0., 400. );
  CurrL1MHT_hist= fs->make <TH1D>("CurrL1MHT_distro", "Curr MHT distro" ,400, 0., 400. );
  CurrL1MHT_rate= fs->make <TH1D>("CurrL1MHT_rate"  , "Curr MHT rate"   ,400, 0., 400. );


	//Turn on curves

  OfflinePt3DIC40 = fs->make<TH3F>("OfflinePt3DIC40 ","OfflinePt3DIC40 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC45 = fs->make<TH3F>("OfflinePt3DIC45 ","OfflinePt3DIC45 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC50 = fs->make<TH3F>("OfflinePt3DIC50 ","OfflinePt3DIC50 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC60 = fs->make<TH3F>("OfflinePt3DIC60 ","OfflinePt3DIC60 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC80 = fs->make<TH3F>("OfflinePt3DIC80 ","OfflinePt3DIC80 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC100= fs->make<TH3F>("OfflinePt3DIC100","OfflinePt3DIC100", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC120= fs->make<TH3F>("OfflinePt3DIC120","OfflinePt3DIC120", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DIC140= fs->make<TH3F>("OfflinePt3DIC140","OfflinePt3DIC140", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DICDen= fs->make<TH3F>("OfflinePt3DICDen","OfflinePt3DICDen", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);


  OfflinePt3DCurr40 = fs->make<TH3F>("OfflinePt3DCurr40 ","OfflinePt3DCurr40 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr45 = fs->make<TH3F>("OfflinePt3DCurr45 ","OfflinePt3DCurr45 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr50 = fs->make<TH3F>("OfflinePt3DCurr50 ","OfflinePt3DCurr50 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr60 = fs->make<TH3F>("OfflinePt3DCurr60 ","OfflinePt3DCurr60 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr80 = fs->make<TH3F>("OfflinePt3DCurr80 ","OfflinePt3DCurr80 ", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr100= fs->make<TH3F>("OfflinePt3DCurr100","OfflinePt3DCurr100", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr120= fs->make<TH3F>("OfflinePt3DCurr120","OfflinePt3DCurr120", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurr140= fs->make<TH3F>("OfflinePt3DCurr140","OfflinePt3DCurr140", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);
  OfflinePt3DCurrDen= fs->make<TH3F>("OfflinePt3DCurrDen","OfflinePt3DCurrDen", 80,-4.0,4.0, 80,-4.0,4.0, 300,0,300);


  //resolutions
  DeltaEta        = fs->make<TH1F>("L1_DeltaEta"    ,  "L1_DeltaEta"    , 100,-0.5,0.5);
  CurrDeltaEta    = fs->make<TH1F>("CurrL1_DeltaEta", "CurrL1_DeltaEta" , 100,-0.5,0.5);
  DeltaPhi        = fs->make<TH1F>("L1_DeltaPhi"    ,  "L1_DeltaPhi"    , 100,-0.5,0.5);
  CurrDeltaPhi    = fs->make<TH1F>("CurrL1_DeltaPhi", "CurrL1_DeltaPhi" , 100,-0.5,0.5);
  DeltaPt         = fs->make<TH1F>("L1_DeltaPt"     , "L1_DeltaPt"      , 100,-2,2);
  CurrDeltaPt     = fs->make<TH1F>("CurrL1_DeltaPt" ,  "CurrL1_DeltaPt" , 100,-2,2);

  OffPt_upL1Pt    = fs->make<TH2F>("OffPt_upL1Pt"    , "OffPt_upL1Pt"    , 300,0,600, 300,0,600 );
  OffPt_currL1Pt  = fs->make<TH2F>("OffPt_currL1Pt"  , "OffPt_currL1Pt"  , 300,0,600, 300,0,600 );

  DeltaPt_pt      = fs->make<TH2F>("DeltaPt_pt"   , "DeltaPt_1_pt"    , 100,-2,2, 300,0,600 );
  CurrDeltaPt_pt  = fs->make<TH2F>("CurrDeltaPt_pt" , "CurrDeltaPt_pt"  , 100,-2,2, 300,0,600 );

  DeltaPt_nvtx    = fs->make<TH2F>("DeltaPt_nvtx"  , "DeltaPt_1_nvtx"  , 100,-2,2, 100,0,100 );
  CurrDeltaPt_nvtx= fs->make<TH2F>("CurrDeltaPt_nvtx", "CurrDeltaPt_nvtx", 100,-2,2, 100,0,100 );

  DeltaPt_eta     = fs->make<TH2F>("DeltaPt_eta"  , "DeltaPt_1_eta"  , 100,-2,2, 100,-5.0,5.0 );
  CurrDeltaPt_eta = fs->make<TH2F>("CurrDeltaPt_eta", "CurrDeltaPt_eta", 100,-2,2, 100,-5.0,5.0 );

	//PU distribution
  NVTX						= fs->make<TH1F>("NVTX","Number of Primary Vertices",100,0,100);



}

AnalyzeJets::~AnalyzeJets()
{

}


void
AnalyzeJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool evValid =true;

	//////////////////////////////////
	// Handles to L1 information
	/////////////////////////////////


  //calibrated tower jet collection: L1TowerJetCollection
  edm::Handle<L1TowerJetCollection> CalibJets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibCircle8"), CalibJets);
  if(!CalibJets.isValid()){
		evValid=false;
		edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibCircle8") << std::endl;
  }	
  //calibrated tower MHT/HT collection: L1EtMissParticleCollection
  edm::Handle<l1extra::L1EtMissParticleCollection> L1MHT_up;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibCircle8MHT"),L1MHT_up);
  if(!L1MHT_up.isValid()){
    evValid=false;
		edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibCircle8MHT") << std::endl;
	}
  //calibrated tower jet collection: l1extra:: L1JetParticleCollection
  edm::Handle<l1extra::L1JetParticleCollection> CalibJets_L1extra;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalibCircle8l1extra"), CalibJets_L1extra);
  if(!CalibJets_L1extra.isValid()){
		evValid=false;
		edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalibCircle8l1extra") << std::endl;
	}


  //Current L1 extra jets
  vector<InputTag>  l1extraparticles= conf_.getParameter< vector < InputTag > >("extrajet"); 

  for (uint i=0; i< l1extraparticles.size();++i){
    Handle<l1extra::L1JetParticleCollection> currl1Jets;
    iEvent.getByLabel(l1extraparticles[i], currl1Jets );
    if(!currl1Jets.isValid()){evValid=false;}
  }

  //Current HT&MHT
  edm::Handle<l1extra::L1EtMissParticleCollection> L1MHT_curr;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("L1extraMHT"),L1MHT_curr);
  if(!L1MHT_curr.isValid()){
		evValid=false;
		edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("L1extraMHT") << std::endl;

	}


	//////////////////////////////////
	// Handles to offline information
	//////////////////////////////////

	//Need this for information about PU
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RecoVertices"), vtx); 
  if(!vtx.isValid()){
		edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RecoVertices") << std::endl;
		evValid=false;
	}

	//PU subtracted AK5 calo jets-must be in root file read in
  edm::Handle<reco::CaloJetCollection> Calojets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PUsubCaloJets"), Calojets);
  if(!Calojets.isValid()){
		edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("PUsubCaloJets") << std::endl;
		evValid=false;
	}

	//information about offline calo rho-must be in root file read in
  edm::Handle<double> rhoCALOjets;  
//  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RhoCaloJets"), rhoCALOjets);
	iEvent.getByLabel("ak5CaloJets","rho", rhoCALOjets);
  if(!rhoCALOjets.isValid()){
		edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("RhoCaloJets") << std::endl;
		evValid=false;
	}


  // handle to the jet ID variables
  edm::Handle<reco::JetIDValueMap> hJetIDMap;
  iEvent.getByLabel( conf_.getParameter<edm::InputTag>("jetIDHelperConfig") , hJetIDMap );
  if(!hJetIDMap.isValid()){
		edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("jetIDHelperConfig")<< std::endl;
		evValid=false;
	}

   edm::Handle<edm::View< reco::CaloJet > > hJets; // uncorrected jets!
   iEvent.getByLabel("ak5CaloJets", hJets );
   if(!hJets.isValid()){
		edm::LogWarning("MissingProduct") << "ak5CaloJets"<< std::endl;
		evValid=false;
	}





	/////////////////////////////////////////////////////////////
  // Start filling histograms etc if event is valid
	//     ie input root file contains all necessary collections
	/////////////////////////////////////////////////////////////
 
  if( !evValid ) {

	cout<<"invalid event "<<endl;	
		
 
  }else{

   	int PVSIZE=vtx->size();
    NVTX->Fill(PVSIZE);

	  float min_l1_pt(10);
    float min_off_pt(20);

    //Fill a vector of TLorentzVectors with upgrade calibrated L1 jets

    vector<TLorentzVector> upgradeJets;
    TLorentzVector up_l1;
    for (l1extra::L1JetParticleCollection::const_iterator il1 = CalibJets_L1extra->begin(); il1!= CalibJets_L1extra->end() ; ++il1 ){

      up_l1.SetPtEtaPhiM(il1->p4().Pt(),il1->p4().Eta(),il1->p4().Phi(),il1->p4().M());
      //cout<<" upgrade (pt,eta,phi) =("<<up_l1.Pt()<<","<<up_l1.Eta()<<","<<up_l1.Phi()<<")"<<endl;
			if(up_l1.Pt() < min_l1_pt) continue;
			if(up_l1.Eta() > 3.0) continue;

      upgradeJets.push_back(up_l1);

    }








    //Fill a vector of TLorentzVecotrs with current L1 jets

    vector<TLorentzVector> currentJets;
    TLorentzVector curr_l1;

    for (uint i=0; i< l1extraparticles.size();++i){
			Handle<l1extra::L1JetParticleCollection> currl1Jets;
      iEvent.getByLabel(l1extraparticles[i],currl1Jets); 
      for(l1extra::L1JetParticleCollection::const_iterator itr = currl1Jets->begin();  itr != currl1Jets->end(); ++itr ) {

        curr_l1.SetPtEtaPhiM(itr->p4().Pt(),itr->p4().eta(),itr->p4().phi(),itr->p4().M());
        //cout<<"Current (pt,eta,phi): "<<curr_l1.Pt()<<" "<<curr_l1.Eta()<<" "<<curr_l1.Phi()<<endl;
			  if(curr_l1.Pt() < min_l1_pt) continue;
        if(curr_l1.Eta() > 3.0) continue;

				currentJets.push_back(curr_l1);

      }
    }

    //Fill a vector of TLorentzVecotrs with offline PU corrected AK5 calo jets
    vector<TLorentzVector> offlineJets;
    TLorentzVector off;

    unsigned int idx;

    //L2L3 PU corrected jet collection
    reco::CaloJetCollection::const_iterator it = Calojets->begin() ;

    //uncorrected jet collection for jet ID reference
    for ( edm::View<reco::CaloJet>::const_iterator ibegin = hJets->begin(),
          iend = hJets->end(), ijet = ibegin;
          ijet != iend; ++ijet , ++it) {

      idx = ijet - ibegin; 

//      edm::RefToBase<reco::CaloJet> jetRef = hJets->refAt(idx);
//      reco::JetID const & jetId = (*hJetIDMap)[ jetRef ];
      
			//set the TLorentz vector with corrected energy
      off.SetPtEtaPhiM(it->p4().Pt(),it->p4().eta(),it->p4().phi(),it->p4().M());

      //loose jet selections
      if (it->emEnergyFraction()<0.01) continue;
      //if (jetId.fHPD>0.98  ) continue;
      //if (jetId.n90Hits<=1) continue;
      
      if( off.Pt() < min_off_pt ) continue;
			if( off.Eta() > 3.0 ) continue;

      //Fill vector with loose ID's calojets
      offlineJets.push_back(off);
      //out<<" offline (pt,eta,phi) =("<<off.Pt()<<","<<off.Eta()<<","<<off.Phi()<<")"<<endl;

    }

    if(offlineJets.size()>0)  sort (offlineJets.begin(), offlineJets.end(),sortTLorentz);


 		///////////////////////////////////////////
		//Rate Plots: require L1 information only
		///////////////////////////////////////////





		//upgrade jets
    if(upgradeJets.size()>0)  sort (upgradeJets.begin(), upgradeJets.end(),sortTLorentz);

		//Jet triggers: 1, 2, 3, 4 jets
    for(unsigned int i=0; i<4; i++){
      if (upgradeJets.size() > i) {
        upgradeJets_hist[i]->Fill(upgradeJets.at(i).Pt());
      }else{
        upgradeJets_hist[i]->Fill(0);
      }
    }
    //Upgrade L1 HT & MHT
    float Up_HT = L1MHT_up->begin()->etTotal();
    float Up_MHT= L1MHT_up->begin()->pt();

    UpgradeHT_rate->Fill(Up_HT);
    UpgradeMHT_hist->Fill(Up_MHT);


    //Current L1 Jets
    if(currentJets.size()>0)  sort (currentJets.begin(), currentJets.end(),sortTLorentz);

		//Jet triggers: 1, 2, 3, 4 jets
    for(unsigned int i=0; i<4; i++){
      if (currentJets.size() > i) {
        CurrL1jets_hist[i]->Fill(currentJets.at(i).Pt());
      }else{
        CurrL1jets_hist[i]->Fill(0);
      }
    }

    //Current L1 HT & MHT
    float Curr_HT = L1MHT_curr->begin()->etTotal();
    float Curr_MHT= L1MHT_curr->begin()->pt();

    CurrL1HT_hist->Fill(Curr_HT);
    CurrL1MHT_hist->Fill(Curr_MHT);



		//////////////////////////////////////////////////////
		// Match L1 jet to offline jet: turn ons, resolutions
		//////////////////////////////////////////////////////

		double DR = 0.5;

		//match to leading jet only
    bool leading_jet_only = true;
    if(offlineJets.size()>0 && leading_jet_only){  //For the turn on curves

			//upgrade jets

      std::pair< double, int> dist_ref_up = Leading_Match(offlineJets, upgradeJets);
			//matched jet: distance between jets is dist_ref_up.first and ref of matched L1 jet is dist_ref_up.second

      if( dist_ref_up.first<DR){
				//leading offline jet
        oPt_ =offlineJets[0].Pt();
        oEta_=offlineJets[0].Eta();
        oPhi_=offlineJets[0].Phi();

        l1Pt_ =upgradeJets[dist_ref_up.second].Pt();
        l1Eta_=upgradeJets[dist_ref_up.second].Eta();
        l1Phi_=upgradeJets[dist_ref_up.second].Phi();

			  float _dpt4=1-((l1Pt_ )/oPt_);

        DeltaEta       ->Fill(oEta_-l1Eta_) ;
        DeltaPhi       ->Fill(oPhi_-l1Phi_) ;
        DeltaPt        ->Fill(_dpt4) ; 

				//offline pt vs l1 pt, res Pt vs various things
				OffPt_upL1Pt      -> Fill(oPt_  , l1Pt_ );
        DeltaPt_pt    -> Fill(_dpt4 , l1Pt_ );
        DeltaPt_nvtx  -> Fill(_dpt4 , PVSIZE);
        DeltaPt_eta   -> Fill(_dpt4 , l1Eta_);


				//Fill a few turn on curves
				if(l1Pt_>40)     OfflinePt3DIC40->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>45)     OfflinePt3DIC45->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>50)     OfflinePt3DIC50->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>60)     OfflinePt3DIC60->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>80)     OfflinePt3DIC80->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>100)    OfflinePt3DIC100->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>120)    OfflinePt3DIC120->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>140)    OfflinePt3DIC140->Fill( oEta_ , oPhi_ , oPt_ );

				OfflinePt3DICDen->Fill( oEta_ , oPhi_ , oPt_ );
			}



			//current L1 jets

      dist_ref_up = Leading_Match(offlineJets, currentJets);
			//matched jet: distance between jets is dist_ref_up.first and ref of matched L1 jet is dist_ref_up.second

      if( dist_ref_up.first<DR){
				//leading offline jet
        oPt_ =offlineJets[0].Pt();
        oEta_=offlineJets[0].Eta();
        oPhi_=offlineJets[0].Phi();

        l1Pt_ =currentJets[dist_ref_up.second].Pt();
        l1Eta_=currentJets[dist_ref_up.second].Eta();
        l1Phi_=currentJets[dist_ref_up.second].Phi();

			  float _dpt4=1-((l1Pt_ )/oPt_);

        CurrDeltaEta       ->Fill(oEta_-l1Eta_) ;
        CurrDeltaPhi       ->Fill(oPhi_-l1Phi_) ;
        CurrDeltaPt        ->Fill(_dpt4) ; 

				//offline pt vs l1 pt, res Pt vs various things
				OffPt_currL1Pt  -> Fill(oPt_  , l1Pt_ );
        CurrDeltaPt_pt    -> Fill(_dpt4 , l1Pt_ );
        CurrDeltaPt_nvtx  -> Fill(_dpt4 , PVSIZE);
        CurrDeltaPt_eta   -> Fill(_dpt4 , l1Eta_);


				//Fill a few turn on curves
				if(l1Pt_>40)     OfflinePt3DCurr40->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>45)     OfflinePt3DCurr45->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>50)     OfflinePt3DCurr50->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>60)     OfflinePt3DCurr60->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>80)     OfflinePt3DCurr80->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>100)    OfflinePt3DCurr100->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>120)    OfflinePt3DCurr120->Fill( oEta_ , oPhi_ , oPt_ );
				if(l1Pt_>140)    OfflinePt3DCurr140->Fill( oEta_ , oPhi_ , oPt_ );

				OfflinePt3DICDen->Fill( oEta_ , oPhi_ , oPt_ );
			}


		}



  }//valid event


}

//////////////  member functions  ///////////////

std::pair< double, int> AnalyzeJets::Leading_Match(vector<TLorentzVector> offlineJets,vector<TLorentzVector> L1Jets )
{
  std::pair< double, int> dist_ref;
  dist_ref.first = 999; //this is the distance between jets
  dist_ref.second = 999; //dummy initialization of jet reference

  oPt_=offlineJets[0].Pt();
  oEta_=offlineJets[0].Eta();
  oPhi_=offlineJets[0].Phi();
  

  for(unsigned int i = 0 ; i<L1Jets.size(); ++i){
    if(offlineJets[0].Pt()<=0 || L1Jets[i].Pt()<=0) continue;

    if(L1Jets[i].DeltaR(offlineJets[0])<dist_ref.first) {
      dist_ref.first = L1Jets[i].DeltaR(offlineJets[0]) ;
      dist_ref.second = i;
    }
  }

  return dist_ref;
}

// ------------ method called once each job just before starting event loop  ------------
void 
AnalyzeJets::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnalyzeJets::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
AnalyzeJets::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
AnalyzeJets::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AnalyzeJets::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AnalyzeJets::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnalyzeJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzeJets);
