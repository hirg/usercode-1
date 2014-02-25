#ifndef ClusterTree_H
#define ClusterTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Math/LorentzVector.h"
#include "Math/PositionVector3D.h"

#include <cmath>
#include <vector>
#include <string>
#include "assert.h"

using namespace std;
//
// Ntuple structure:
//

// Ntuple content:

class ClusterTree {

  public:
    /// float doesn't have dictionary by default, so use double
    //typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > Point;

    /// event identification
    unsigned int   event_;
    unsigned int   run_;
    unsigned int   lumi_;

    // event properties

    unsigned int   nvtx_;
    float          rho_Calo_;
    float          rho_PF_;
    Point          primaryVertex_;


    // jet properties
    vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > *Calojet4_mom_;    vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > *Calojet4_mom_PU_;
    vector<int>                                                          *Calojet_indx_;
    unsigned int                                                          Calojet_multi_;
    vector<float>                                                        *Calojet_area_;
    
    //recast as an array of floats:
    float Calojet_jet_pt_[200];
    float Calojet_jet_eta_[200];
    float Calojet_jet_phi_[200];

    float Calo_match8circ_jet_pt_[200];
    float Calo_match9circ_jet_pt_[200];
    float Calo_match10circ_jet_pt_[200];
    float Calo_match12circ_jet_pt_[200];
    float Calo_match8sq_jet_pt_[200];
    float Calo_match9sq_jet_pt_[200];
    float Calo_match10sq_jet_pt_[200];
    float Calo_match12sq_jet_pt_[200];
      
    float Calo_match8circ_jet_eta_[200];
    float Calo_match9circ_jet_eta_[200];
    float Calo_match10circ_jet_eta_[200];
    float Calo_match12circ_jet_eta_[200];
    float Calo_match8sq_jet_eta_[200];
    float Calo_match9sq_jet_eta_[200];
    float Calo_match10sq_jet_eta_[200];
    float Calo_match12sq_jet_eta_[200];

    float Calo_match8circ_jet_phi_[200];
    float Calo_match9circ_jet_phi_[200];
    float Calo_match10circ_jet_phi_[200];
    float Calo_match12circ_jet_phi_[200];
    float Calo_match8sq_jet_phi_[200];
    float Calo_match9sq_jet_phi_[200];
    float Calo_match10sq_jet_phi_[200];
    float Calo_match12sq_jet_phi_[200];

    //jet multiplicity
    unsigned int Calo_match8circ_jet_multi_;
    unsigned int Calo_match9circ_jet_multi_;      
    unsigned int Calo_match10circ_jet_multi_;
    unsigned int Calo_match12circ_jet_multi_;       
            
    unsigned int Calo_match8sq_jet_multi_;       
    unsigned int Calo_match9sq_jet_multi_;
    unsigned int Calo_match10sq_jet_multi_;       
    unsigned int Calo_match12sq_jet_multi_;       

    //tower jets
    vector<int>                                                          * L1_tow_8circ_jet_indx_;
    unsigned int                                                           L1_tow_8circ_jet_multi_;
    float                                                                  L1_tow_8circ_scaledrho_;
    float L1_tow_8circ_jet_pt_[200];
    float L1_tow_8circ_jet_eta_[200];
    float L1_tow_8circ_jet_phi_[200];


    vector<int>                                                          * L1_tow_9circ_jet_indx_;
    unsigned int                                                           L1_tow_9circ_jet_multi_;
    float                                                                  L1_tow_9circ_scaledrho_;
    float L1_tow_9circ_jet_pt_[200];
    float L1_tow_9circ_jet_eta_[200];
    float L1_tow_9circ_jet_phi_[200];

    vector<int>                                                          * L1_tow_10circ_jet_indx_;
    unsigned int                                                           L1_tow_10circ_jet_multi_;
    float                                                                  L1_tow_10circ_scaledrho_;
    float L1_tow_10circ_jet_pt_[200];
    float L1_tow_10circ_jet_eta_[200];
    float L1_tow_10circ_jet_phi_[200];


    vector<int>                                                          * L1_tow_12circ_jet_indx_;
    unsigned int                                                           L1_tow_12circ_jet_multi_;
    float                                                                  L1_tow_12circ_scaledrho_;
    float L1_tow_12circ_jet_pt_[200];
    float L1_tow_12circ_jet_eta_[200];
    float L1_tow_12circ_jet_phi_[200];

    vector<int>                                                          * L1_tow_8sq_jet_indx_;
    unsigned int                                                           L1_tow_8sq_jet_multi_;
    float                                                                  L1_tow_8sq_scaledrho_;
    float L1_tow_8sq_jet_pt_[200];
    float L1_tow_8sq_jet_eta_[200];
    float L1_tow_8sq_jet_phi_[200];

    vector<int>                                                          * L1_tow_9sq_jet_indx_;
    unsigned int                                                           L1_tow_9sq_jet_multi_;
    float                                                                  L1_tow_9sq_scaledrho_;
    float L1_tow_9sq_jet_pt_[200];
    float L1_tow_9sq_jet_eta_[200];
    float L1_tow_9sq_jet_phi_[200];

    vector<int>                                                          * L1_tow_10sq_jet_indx_;
    unsigned int                                                           L1_tow_10sq_jet_multi_;
    float                                                                  L1_tow_10sq_scaledrho_;
    float L1_tow_10sq_jet_pt_[200];
    float L1_tow_10sq_jet_eta_[200];
    float L1_tow_10sq_jet_phi_[200];

    vector<int>                                                          * L1_tow_12sq_jet_indx_;
    unsigned int                                                           L1_tow_12sq_jet_multi_;
    float                                                                  L1_tow_12sq_scaledrho_;
    float L1_tow_12sq_jet_pt_[200];
    float L1_tow_12sq_jet_eta_[200];
    float L1_tow_12sq_jet_phi_[200];

  public:
    /// this is the main element
    TTree *tree_;
    TFile *f_;

    /// hold the names of variables to facilitate things (filled during Init)
    std::vector<std::string> variables_;

    /// default constructor  

    /// default destructor
    ~ClusterTree(){ 

      delete Calojet4_mom_;
      delete Calojet4_mom_PU_;
      delete Calojet_indx_;

      delete L1_tow_8circ_jet_indx_;
      delete L1_tow_9circ_jet_indx_;
      delete L1_tow_10circ_jet_indx_;
      delete L1_tow_12circ_jet_indx_;

      delete L1_tow_8sq_jet_indx_;
      delete L1_tow_9sq_jet_indx_;
      delete L1_tow_10sq_jet_indx_;
      delete L1_tow_12sq_jet_indx_;

      if (f_) f_->Close();
 
    };

    /// initialize varibles and fill list of available variables
    void InitVariables();
    void clearVectors();
    void resizeVectors(int);

    /// load a ClusterTree
    void LoadTree(const char* file){
      f_ = TFile::Open(file);
      assert(f_);
      tree_ = dynamic_cast<TTree*>(f_->Get("tree"));
      assert(tree_);
    }

    /// create a ClusterTree
    void CreateTree(){
      edm::Service<TFileService> fs;
	tree_ = fs->make <TTree>("Jet_ntuple","Jet_ntuple");
    //  tree_ = new TTree("Jet_ntuple","Jet_ntuple");
      f_ = 0;
      InitVariables();

      //book the branches
      tree_->Branch("event"         , &event_            ,   "event/i");
      tree_->Branch("run"           , &run_              ,   "run/i");
      tree_->Branch("lumi"          , &lumi_             ,   "lumi/i");

      tree_->Branch("nvtx"          , &nvtx_             ,   "nvtx/i");
      tree_->Branch("rho_Calo"      , &rho_Calo_         ,   "rho_Calo/F");
      tree_->Branch("rho_PF"        , &rho_PF_           ,   "rho_PF/F");
      tree_->Branch("primaryVertex" , "Point"            , &primaryVertex_);
      
      //Calo jets
      tree_->Branch("Offline_Calojet"        , "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >", &Calojet4_mom_);
      tree_->Branch("Offline_Calojet_PU"        , "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >", &Calojet4_mom_PU_);
      tree_->Branch("Offline_Calojet_index"         , "std::vector<int>"         , &Calojet_indx_ );
      tree_->Branch("Offline_Calojet_multiplicity"  , &Calojet_multi_            , "Offline_Calojet_multiplicity/i");
      tree_->Branch("Offline_Calojet_area"          , "std::vector<float>"       , &Calojet_area_ );

      tree_->Branch("Offline_Calojet_pt"     , &Calojet_jet_pt_ , "Offline_Calojet_pt[Offline_Calojet_multiplicity]/f"  );
      tree_->Branch("Offline_Calojet_eta"    , &Calojet_jet_eta_, "Offline_Calojet_eta[Offline_Calojet_multiplicity]/f" );
      tree_->Branch("Offline_Calojet_phi"    , &Calojet_jet_phi_, "Offline_Calojet_phi[Offline_Calojet_multiplicity]/f" );

      tree_->Branch("Calojet_match8circ_multi"  , &Calo_match8circ_jet_multi_,   "Calojet_match8circ_multi/i"   );
      tree_->Branch("Calojet_match9circ_multi"  , &Calo_match9circ_jet_multi_,   "Calojet_match9circ_multi/i"   );
      tree_->Branch("Calojet_match10circ_multi" , &Calo_match10circ_jet_multi_,  "Calojet_match10circ_multi/i"  );
      tree_->Branch("Calojet_match12circ_multi" , &Calo_match12circ_jet_multi_,  "Calojet_match12circ_multi/i"  );
      tree_->Branch("Calojet_match8sq_multi"    , &Calo_match8sq_jet_multi_,     "Calojet_match8sq_multi/i"     );
      tree_->Branch("Calojet_match9sq_multi"    , &Calo_match9sq_jet_multi_,     "Calojet_match9sq_multi/i"     );
      tree_->Branch("Calojet_match10sq_multi"   , &Calo_match10sq_jet_multi_,    "Calojet_match10sq_multi/i"    );
      tree_->Branch("Calojet_match12sq_multi"   , &Calo_match12sq_jet_multi_,    "Calojet_match12sq_multi/i"    );

      //pt for matched jets
      tree_->Branch("Calojet_match8circ_pt"     , &Calo_match8circ_jet_pt_,
             "Calojet_match8circ_pt[Calojet_match8circ_multi]/f");
      tree_->Branch("Calojet_match9circ_pt"     , &Calo_match9circ_jet_pt_,
             "Calojet_match9circ_pt[Calojet_match9circ_multi]/f");
      tree_->Branch("Calojet_match10circ_pt"    , &Calo_match10circ_jet_pt_,
             "Calojet_match10circ_pt[Calojet_match10circ_multi]/f");
      tree_->Branch("Calojet_match12circ_pt"    , &Calo_match12circ_jet_pt_,
             "Calojet_match12circ_pt[Calojet_match12circ_multi]/f");
      tree_->Branch("Calojet_match8sq_pt"       , &Calo_match8sq_jet_pt_,
             "Calojet_match8sq_pt[Calojet_match8sq_multi]/f");
      tree_->Branch("Calojet_match9sq_pt"       , &Calo_match9sq_jet_pt_,
             "Calojet_match9sq_pt[Calojet_match9sq_multi]/f");
      tree_->Branch("Calojet_match10sq_pt"      , &Calo_match10sq_jet_pt_,
             "Calojet_match10sq_pt[Calojet_match10sq_multi]/f");
      tree_->Branch("Calojet_match12sq_pt"      , &Calo_match12sq_jet_pt_,
             "Calojet_match12sq_pt[Calojet_match12sq_multi]/f");
      //eta for matched jets
      tree_->Branch("Calojet_match8circ_eta"    , &Calo_match8circ_jet_eta_,
             "Calojet_match8circ_eta[Calojet_match8circ_multi]/f");
      tree_->Branch("Calojet_match9circ_eta"    , &Calo_match9circ_jet_eta_,
             "Calojet_match9circ_eta[Calojet_match9circ_multi]/f");
      tree_->Branch("Calojet_match10circ_eta"   , &Calo_match10circ_jet_eta_,
             "Calojet_match10circ_eta[Calojet_match10circ_multi]/f");
      tree_->Branch("Calojet_match12circ_eta"   , &Calo_match12circ_jet_eta_,
             "Calojet_match12circ_eta[Calojet_match12circ_multi]/f");
      tree_->Branch("Calojet_match8sq_eta"      , &Calo_match8sq_jet_eta_,
             "Calojet_match8sq_eta[Calojet_match8sq_multi]/f");
      tree_->Branch("Calojet_match9sq_eta"      , &Calo_match9sq_jet_eta_,
             "Calojet_match9sq_eta[Calojet_match9sq_multi]/f");
      tree_->Branch("Calojet_match10sq_eta"     , &Calo_match10sq_jet_eta_,
             "Calojet_match10sq_eta[Calojet_match10sq_multi]/f");
      tree_->Branch("Calojet_match12sq_eta"     , &Calo_match12sq_jet_eta_,
             "Calojet_match12sq_eta[Calojet_match12sq_multi]/f");
      //phi for matched jets
      tree_->Branch("Calojet_match8circ_phi"    , &Calo_match8circ_jet_phi_,
             "Calojet_match8circ_phi[Calojet_match8circ_multi]/f");
      tree_->Branch("Calojet_match9circ_phi"    , &Calo_match9circ_jet_phi_,
             "Calojet_match9circ_phi[Calojet_match9circ_multi]/f");
      tree_->Branch("Calojet_match10circ_phi"   , &Calo_match10circ_jet_phi_,
             "Calojet_match10circ_phi[Calojet_match10circ_multi]/f");
      tree_->Branch("Calojet_match12circ_phi"   , &Calo_match12circ_jet_phi_,
             "Calojet_match12circ_phi[Calojet_match12circ_multi]/f");
      tree_->Branch("Calojet_match8sq_phi"      , &Calo_match8sq_jet_phi_,
             "Calojet_match8sq_phi[Calojet_match8sq_multi]/f");
      tree_->Branch("Calojet_match9sq_phi"      , &Calo_match9sq_jet_phi_,
             "Calojet_match9sq_phi[Calojet_match9sq_multi]/f");
      tree_->Branch("Calojet_match10sq_phi"     , &Calo_match10sq_jet_phi_,
             "Calojet_match10sq_phi[Calojet_match10sq_multi]/f");
      tree_->Branch("Calojet_match12sq_phi"     , &Calo_match12sq_jet_phi_,
             "Calojet_match12sq_phi[Calojet_match12sq_multi]/f");



      //Tower jets: circular
      tree_->Branch("L1_tow_8circ_jet_index"       , "std::vector<int>"       , &L1_tow_8circ_jet_indx_ );
      tree_->Branch("L1_tow_8circ_jet_multiplicity", &L1_tow_8circ_jet_multi_ , "L1_tow_8circ_jet_multiplicity/i");
      tree_->Branch("L1_tow_8circ_rho"             , &L1_tow_8circ_scaledrho_ , "L1_tow_8circ_rho_/f");

      tree_->Branch("L1_tow_8circ_pt"     , &L1_tow_8circ_jet_pt_ , "L1_tow_8circ_pt[Calojet_match8circ_multi]/f"  );
      tree_->Branch("L1_tow_8circ_eta"    , &L1_tow_8circ_jet_eta_, "L1_tow_8circ_eta[Calojet_match8circ_multi]/f" );
      tree_->Branch("L1_tow_8circ_phi"    , &L1_tow_8circ_jet_phi_, "L1_tow_8circ_phi[Calojet_match8circ_multi]/f" );      

      tree_->Branch("L1_tow_9circ_jet_index"       , "std::vector<int>"     , &L1_tow_9circ_jet_indx_ );
      tree_->Branch("L1_tow_9circ_jet_multiplicity", &L1_tow_9circ_jet_multi_ , "L1_tow_9circ_jet_multiplicity/i");
      tree_->Branch("L1_tow_9circ_rho"             , &L1_tow_9circ_scaledrho_ , "L1_tow_9circ_rho_/f");

      tree_->Branch("L1_tow_9circ_pt"     , &L1_tow_9circ_jet_pt_ , "L1_tow_9circ_pt[Offline_Calojet_multiplicity]/f"  );
      tree_->Branch("L1_tow_9circ_eta"    , &L1_tow_9circ_jet_eta_, "L1_tow_9circ_eta[Offline_Calojet_multiplicity]/f" );
      tree_->Branch("L1_tow_9circ_phi"    , &L1_tow_9circ_jet_phi_, "L1_tow_9circ_phi[Offline_Calojet_multiplicity]/f" );

      tree_->Branch("L1_tow_10circ_jet_index"       , "std::vector<int>"        , &L1_tow_10circ_jet_indx_ );
      tree_->Branch("L1_tow_10circ_jet_multiplicity", &L1_tow_10circ_jet_multi_ , "L1_tow_10circ_jet_multiplicity/i");
      tree_->Branch("L1_tow_10circ_rho"             , &L1_tow_10circ_scaledrho_ , "L1_tow_10circ_rho_/f");

      tree_->Branch("L1_tow_10circ_pt"   , &L1_tow_10circ_jet_pt_ , "L1_tow_10circ_pt[Offline_Calojet_multiplicity]/f"  );
      tree_->Branch("L1_tow_10circ_eta"  , &L1_tow_10circ_jet_eta_, "L1_tow_10circ_eta[Offline_Calojet_multiplicity]/f" );
      tree_->Branch("L1_tow_10circ_phi"  , &L1_tow_10circ_jet_phi_, "L1_tow_10circ_phi[Offline_Calojet_multiplicity]/f" );

      tree_->Branch("L1_tow_12circ_jet_index"       , "std::vector<int>"     , &L1_tow_12circ_jet_indx_ );
      tree_->Branch("L1_tow_12circ_jet_multiplicity", &L1_tow_12circ_jet_multi_ , "L1_tow_12circ_jet_multiplicity/i");
      tree_->Branch("L1_tow_12circ_rho"             , &L1_tow_12circ_scaledrho_ , "L1_tow_12circ_rho_/f");

      tree_->Branch("L1_tow_12circ_pt"   , &L1_tow_12circ_jet_pt_ , "L1_tow_12circ_pt[Offline_Calojet_multiplicity]/f"  );
      tree_->Branch("L1_tow_12circ_eta"  , &L1_tow_12circ_jet_eta_, "L1_tow_12circ_eta[Offline_Calojet_multiplicity]/f" );
      tree_->Branch("L1_tow_12circ_phi"  , &L1_tow_12circ_jet_phi_, "L1_tow_12circ_phi[Offline_Calojet_multiplicity]/f" );

      //Tower jets: square
      tree_->Branch("L1_tow_8sq_jet_index"       , "std::vector<int>"       , &L1_tow_8sq_jet_indx_ );
      tree_->Branch("L1_tow_8sq_jet_multiplicity", &L1_tow_8sq_jet_multi_ , "L1_tow_8sq_jet_multiplicity/i");
      tree_->Branch("L1_tow_8sq_rho"             , &L1_tow_8sq_scaledrho_ , "L1_tow_8sq_rho_/f");

      tree_->Branch("L1_tow_8sq_pt"     , &L1_tow_8sq_jet_pt_ , "L1_tow_8sq_pt[Offline_Calojet_multiplicity]/f"  );
      tree_->Branch("L1_tow_8sq_eta"    , &L1_tow_8sq_jet_eta_, "L1_tow_8sq_eta[Offline_Calojet_multiplicity]/f" );
      tree_->Branch("L1_tow_8sq_phi"    , &L1_tow_8sq_jet_phi_, "L1_tow_8sq_phi[Offline_Calojet_multiplicity]/f" );

      
      tree_->Branch("L1_tow_9sq_jet_index"       , "std::vector<int>"     , &L1_tow_9sq_jet_indx_ );
      tree_->Branch("L1_tow_9sq_jet_multiplicity", &L1_tow_9sq_jet_multi_ , "L1_tow_9sq_jet_multiplicity/i");
      tree_->Branch("L1_tow_9sq_rho"             , &L1_tow_9sq_scaledrho_ , "L1_tow_9sq_rho_/f");

      tree_->Branch("L1_tow_9sq_pt"     , &L1_tow_9sq_jet_pt_ , "L1_tow_9sq_pt[Offline_Calojet_multiplicity]/f"  );
      tree_->Branch("L1_tow_9sq_eta"    , &L1_tow_9sq_jet_eta_, "L1_tow_9sq_eta[Offline_Calojet_multiplicity]/f" );
      tree_->Branch("L1_tow_9sq_phi"    , &L1_tow_9sq_jet_phi_, "L1_tow_9sq_phi[Offline_Calojet_multiplicity]/f" );


      tree_->Branch("L1_tow_10sq_jet_index"       , "std::vector<int>"        , &L1_tow_10sq_jet_indx_ );
      tree_->Branch("L1_tow_10sq_jet_multiplicity", &L1_tow_10sq_jet_multi_ , "L1_tow_10sq_jet_multiplicity/i");
      tree_->Branch("L1_tow_10sq_rho"             , &L1_tow_10sq_scaledrho_ , "L1_tow_10sq_rho_/f");

      tree_->Branch("L1_tow_10sq_pt"     , &L1_tow_10sq_jet_pt_ , "L1_tow_10sq_pt[Offline_Calojet_multiplicity]/f"  );
      tree_->Branch("L1_tow_10sq_eta"    , &L1_tow_10sq_jet_eta_, "L1_tow_10sq_eta[Offline_Calojet_multiplicity]/f" );
      tree_->Branch("L1_tow_10sq_phi"    , &L1_tow_10sq_jet_phi_, "L1_tow_10sq_phi[Offline_Calojet_multiplicity]/f" );


      tree_->Branch("L1_tow_12sq_jet_index"       , "std::vector<int>"     , &L1_tow_12sq_jet_indx_ );
      tree_->Branch("L1_tow_12sq_jet_multiplicity", &L1_tow_12sq_jet_multi_ , "L1_tow_12sq_jet_multiplicity/i");
      tree_->Branch("L1_tow_12sq_rho"             , &L1_tow_12sq_scaledrho_ , "L1_tow_12sq_rho_/f");

      tree_->Branch("L1_tow_12sq_pt"     , &L1_tow_12sq_jet_pt_ , "L1_tow_12sq_pt[Offline_Calojet_multiplicity]/f"  );
      tree_->Branch("L1_tow_12sq_eta"    , &L1_tow_12sq_jet_eta_, "L1_tow_12sq_eta[Offline_Calojet_multiplicity]/f" );
      tree_->Branch("L1_tow_12sq_phi"    , &L1_tow_12sq_jet_phi_, "L1_tow_12sq_phi[Offline_Calojet_multiplicity]/f" );


    }

    // initialze a ClusterTree
    void InitTree(){
      assert(tree_);
      // don't forget to set pointers to zero before you set address
      // or you will fully appreciate that "ROOT sucks" :)
      InitVariables();
      //Set branch address
      Int_t currentState = gErrorIgnoreLevel;
      // gErrorIgnoreLevel = kError;
      gErrorIgnoreLevel = kBreak;
      tree_->SetBranchAddress("event",                  &event_);
      tree_->SetBranchAddress("run",                    &run_);
      tree_->SetBranchAddress("lumi",                   &lumi_);

      tree_->SetBranchAddress("nvtx",                   &nvtx_);  
      tree_->SetBranchAddress("rho_Calo",               &rho_Calo_);
      tree_->SetBranchAddress("rho_PF",                 &rho_PF_);
      tree_->SetBranchAddress("primaryVertex",          &primaryVertex_);

      tree_->SetBranchAddress("Offline_Calojet"             , &Calojet4_mom_);
      tree_->SetBranchAddress("Offline_Calojet_PU"          , &Calojet4_mom_PU_);
      tree_->SetBranchAddress("Offline_Calojet_index"       , &Calojet_indx_);
      tree_->SetBranchAddress("Offline_Calojet_multiplicity", &Calojet_multi_);
      tree_->SetBranchAddress("Offline_Calojet_area"        , &Calojet_area_);
      
      //tree_->SetBranchAddress("Offline_Calojet_pt"             , &Calojet_jet_pt_  );
      //tree_->SetBranchAddress("Offline_Calojet_eta"            , &Calojet_jet_eta_ );
      //tree_->SetBranchAddress("Offline_Calojet_phi"            , &Calojet_jet_phi_ );
      


    
      
      //Tower jets: circle

      tree_->SetBranchAddress("L1_tow_8circ_jet_index"        , &L1_tow_8circ_jet_indx_ );
      tree_->SetBranchAddress("L1_tow_8circ_jet_multiplicity" , &L1_tow_8circ_jet_multi_);
      tree_->SetBranchAddress("L1_tow_8circ_rho"              , &L1_tow_8circ_scaledrho_);
      
      //pt,eta,phi
      //tree_->SetBranchAddress("L1_tow_8circ_pt"               , &L1_tow_8circ_jet_pt_ );
      //tree_->SetBranchAddress("L1_tow_8circ_eta"              , &L1_tow_8circ_jet_eta_);
      //tree_->SetBranchAddress("L1_tow_8circ_phi"              , &L1_tow_8circ_jet_phi_);
      
 
      tree_->SetBranchAddress("L1_tow_9circ_jet_index"        , &L1_tow_9circ_jet_indx_ );
      tree_->SetBranchAddress("L1_tow_9circ_jet_multiplicity" , &L1_tow_9circ_jet_multi_);
      tree_->SetBranchAddress("L1_tow_9circ_rho"              , &L1_tow_9circ_scaledrho_);
      
      tree_->SetBranchAddress("L1_tow_10circ_jet_index"       , &L1_tow_10circ_jet_indx_ );
      tree_->SetBranchAddress("L1_tow_10circ_jet_multiplicity", &L1_tow_10circ_jet_multi_);
      tree_->SetBranchAddress("L1_tow_10circ_rho"             , &L1_tow_10circ_scaledrho_);
      
      tree_->SetBranchAddress("L1_tow_12circ_jet_index"       , &L1_tow_12circ_jet_indx_ );
      tree_->SetBranchAddress("L1_tow_12circ_jet_multiplicity", &L1_tow_12circ_jet_multi_);
      tree_->SetBranchAddress("L1_tow_12circ_rho"             , &L1_tow_12circ_scaledrho_);
      
      //Tower jets: square
      tree_->SetBranchAddress("L1_tow_8sq_jet_index"        , &L1_tow_8sq_jet_indx_ );
      tree_->SetBranchAddress("L1_tow_8sq_jet_multiplicity" , &L1_tow_8sq_jet_multi_);
      tree_->SetBranchAddress("L1_tow_8sq_rho"              , &L1_tow_8sq_scaledrho_);
      
      tree_->SetBranchAddress("L1_tow_9sq_jet_index"        , &L1_tow_9sq_jet_indx_ );
      tree_->SetBranchAddress("L1_tow_9sq_jet_multiplicity" , &L1_tow_9sq_jet_multi_);
      tree_->SetBranchAddress("L1_tow_9sq_rho"              , &L1_tow_9sq_scaledrho_);
      
      tree_->SetBranchAddress("L1_tow_10sq_jet_index"       , &L1_tow_10sq_jet_indx_ );
      tree_->SetBranchAddress("L1_tow_10sq_jet_multiplicity", &L1_tow_10sq_jet_multi_);
      tree_->SetBranchAddress("L1_tow_10sq_rho"             , &L1_tow_10sq_scaledrho_);
      
      tree_->SetBranchAddress("L1_tow_12sq_jet_index"       , &L1_tow_12sq_jet_indx_ );
      tree_->SetBranchAddress("L1_tow_12sq_jet_multiplicity", &L1_tow_12sq_jet_multi_);
      tree_->SetBranchAddress("L1_tow_12sq_rho"             , &L1_tow_12sq_scaledrho_);




      gErrorIgnoreLevel = currentState;
    }

      /// get a built in type variable by name
      double Get(std::string value);
      /// compare two ClusterTrees for a given event on a given level of precision; 
      /// returns the variables that failed the comparison 
      std::vector<std::string> Compare(ClusterTree* value, double prec=0.005);

    private:



}; 

inline void ClusterTree::InitVariables(){
  // create list of available variables
  if(variables_.empty()){
    //make sure that this is only done once
    variables_.push_back(std::string("event"                ));
    variables_.push_back(std::string("run"                  ));
    variables_.push_back(std::string("lumi"                 ));

    variables_.push_back(std::string("rho_Calo"             ));
    variables_.push_back(std::string("rho_PF"               ));
    variables_.push_back(std::string("nvtx"                 ));
    variables_.push_back(std::string("primaryVertex"        ));

    variables_.push_back(std::string("Offline_Calojet"             ));
    variables_.push_back(std::string("Offline_Calojet_PU"          ));
    variables_.push_back(std::string("Offline_Calojet_index"       ));
    variables_.push_back(std::string("Offline_Calojet_multiplicity"));
    variables_.push_back(std::string("Offline_Calojet_area"        ));
    
    variables_.push_back(std::string("Offline_Calojet_pt"  ));  
    variables_.push_back(std::string("Offline_Calojet_eta" ));  
    variables_.push_back(std::string("Offline_Calojet_phi" ));  
    
    //pt for matched jets
    variables_.push_back("Calojet_match8circ_pt" );
    variables_.push_back("Calojet_match9circ_pt" );
    variables_.push_back("Calojet_match10circ_pt");
    variables_.push_back("Calojet_match12circ_pt");
    variables_.push_back("Calojet_match8sq_pt"   );
    variables_.push_back("Calojet_match9sq_pt"   );
    variables_.push_back("Calojet_match10sq_pt"  );
    variables_.push_back("Calojet_match12sq_pt"  );
    //eta for matched jets
    variables_.push_back("Calojet_match8circ_eta"   );
    variables_.push_back("Calojet_match9circ_eta"   );
    variables_.push_back("Calojet_match10circ_eta"  );
    variables_.push_back("Calojet_match12circ_eta"  );
    variables_.push_back("Calojet_match8sq_eta"     );
    variables_.push_back("Calojet_match9sq_eta"     );
    variables_.push_back("Calojet_match10sq_eta"    );
    variables_.push_back("Calojet_match12sq_eta"    );
    //phi for matched jets
    variables_.push_back("Calojet_match8circ_phi"   );
    variables_.push_back("Calojet_match9circ_phi"   );
    variables_.push_back("Calojet_match10circ_phi"  );
    variables_.push_back("Calojet_match12circ_phi"  );
    variables_.push_back("Calojet_match8sq_phi"     );
    variables_.push_back("Calojet_match9sq_phi"     );
    variables_.push_back("Calojet_match10sq_phi"    );
    variables_.push_back("Calojet_match12sq_phi"    );
    //multiplicity
    variables_.push_back("Calojet_match8circ_multi"   ); 
    variables_.push_back("Calojet_match9circ_multi"   );
    variables_.push_back("Calojet_match10circ_multi"  );
    variables_.push_back("Calojet_match12circ_multi"  );
    variables_.push_back("Calojet_match8sq_multi"     );
    variables_.push_back("Calojet_match9sq_multi"     );
    variables_.push_back("Calojet_match10sq_multi"    );
    variables_.push_back("Calojet_match12sq_multi"    );


    variables_.push_back("L1_tow_8circ_jet_index"   );
    variables_.push_back("L1_tow_8circ_jet_multiplicity");  
    variables_.push_back("L1_tow_8circ_rho_");          
    
    //pt,eta,phi
    variables_.push_back(std::string("L1_tow_8circ_pt"  ));
    variables_.push_back(std::string("L1_tow_8circ_eta" ));
    variables_.push_back(std::string("L1_tow_8circ_phi" ));
    
    
    variables_.push_back("L1_tow_9circ_jet_index"   );
    variables_.push_back("L1_tow_9circ_jet_multiplicity");  
    variables_.push_back("L1_tow_9circ_rho_");  
    
    variables_.push_back(std::string("L1_tow_9circ_pt"  ));
    variables_.push_back(std::string("L1_tow_9circ_eta" ));
    variables_.push_back(std::string("L1_tow_9circ_phi" ));

    variables_.push_back("L1_tow_10circ_jet_index"  );
    variables_.push_back("L1_tow_10circ_jet_multiplicity");  
    variables_.push_back("L1_tow_10circ_rho_");

    variables_.push_back(std::string("L1_tow_10circ_pt"  ));
    variables_.push_back(std::string("L1_tow_10circ_eta" ));
    variables_.push_back(std::string("L1_tow_10circ_phi" ));
    
    variables_.push_back("L1_tow_12circ_jet_index"  );
    variables_.push_back("L1_tow_12circ_jet_multiplicity");  
    variables_.push_back("L1_tow_12circ_rho_");    

    variables_.push_back(std::string("L1_tow_12circ_pt"  ));
    variables_.push_back(std::string("L1_tow_12circ_eta" ));
    variables_.push_back(std::string("L1_tow_12circ_phi" ));     

    //Tower jets: square
    variables_.push_back(std::string("L1_tow_8sq_jet_index"   ));
    variables_.push_back(std::string("L1_tow_8sq_jet_multiplicity"));  
    variables_.push_back(std::string("L1_tow_8sq_rho_"));
    //pt,eta,phi
    variables_.push_back(std::string("L1_tow_8sq_pt"  ));
    variables_.push_back(std::string("L1_tow_8sq_eta" ));
    variables_.push_back(std::string("L1_tow_8sq_phi" ));

    variables_.push_back(std::string("L1_tow_9sq_jet_index"   ));
    variables_.push_back(std::string("L1_tow_9sq_jet_multiplicity"));  
    variables_.push_back(std::string("L1_tow_9sq_rho_"));  
    variables_.push_back(std::string("L1_tow_9sq_pt"  ));
    variables_.push_back(std::string("L1_tow_9sq_eta" ));
    variables_.push_back(std::string("L1_tow_9sq_phi" ));

    variables_.push_back(std::string("L1_tow_10sq_jet_index"  ));
    variables_.push_back(std::string("L1_tow_10sq_jet_multiplicity"));  
    variables_.push_back(std::string("L1_tow_10sq_rho_"));
    variables_.push_back(std::string("L1_tow_10sq_pt"  ));
    variables_.push_back(std::string("L1_tow_10sq_eta" ));
    variables_.push_back(std::string("L1_tow_10sq_phi" ));

    variables_.push_back(std::string("L1_tow_12sq_jet_index"  ));
    variables_.push_back(std::string("L1_tow_12sq_jet_multiplicity"));  
    variables_.push_back(std::string("L1_tow_12sq_rho_"));
    variables_.push_back(std::string("L1_tow_12sq_pt"  ));
    variables_.push_back(std::string("L1_tow_12sq_eta" ));
    variables_.push_back(std::string("L1_tow_12sq_phi" ));


  }   
    // inizialize variables
    event_                = 0;
    run_                  = 0;
    lumi_                 = 0;
      
    rho_Calo_             = -999;
    rho_PF_               = -999;
    nvtx_                 = -999;
    primaryVertex_        = Point();
      
    Calojet4_mom_         = new std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >();
    Calojet4_mom_PU_      = new std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >();   
    Calojet_indx_         = new std::vector<int>();
    Calojet_multi_        = 0;
    Calojet_area_         = new std::vector<float>();

    std::fill_n(Calojet_jet_pt_, 200, 0);
    std::fill_n(Calojet_jet_eta_, 200, 0);
    std::fill_n(Calojet_jet_phi_, 200, 0);

    std::fill_n(Calo_match8circ_jet_pt_, 200, 0);
    std::fill_n(Calo_match9circ_jet_pt_, 200, 0);
    std::fill_n(Calo_match10circ_jet_pt_, 200, 0);
    std::fill_n(Calo_match12circ_jet_pt_, 200, 0);
    std::fill_n(Calo_match8sq_jet_pt_, 200, 0);
    std::fill_n(Calo_match9sq_jet_pt_, 200, 0);
    std::fill_n(Calo_match10sq_jet_pt_, 200, 0);
    std::fill_n(Calo_match12sq_jet_pt_, 200, 0);
      
    std::fill_n(Calo_match8circ_jet_eta_, 200, 0);
    std::fill_n(Calo_match9circ_jet_eta_, 200, 0);
    std::fill_n(Calo_match10circ_jet_eta_, 200, 0);
    std::fill_n(Calo_match12circ_jet_eta_, 200, 0);
    std::fill_n(Calo_match8sq_jet_eta_, 200, 0);
    std::fill_n(Calo_match9sq_jet_eta_, 200, 0);
    std::fill_n(Calo_match10sq_jet_eta_, 200, 0);
    std::fill_n(Calo_match12sq_jet_eta_, 200, 0);

    std::fill_n(Calo_match8circ_jet_phi_, 200, 0);
    std::fill_n(Calo_match9circ_jet_phi_, 200, 0);
    std::fill_n(Calo_match10circ_jet_phi_, 200, 0);
    std::fill_n(Calo_match12circ_jet_phi_, 200, 0);
    std::fill_n(Calo_match8sq_jet_phi_, 200, 0);
    std::fill_n(Calo_match9sq_jet_phi_, 200, 0);
    std::fill_n(Calo_match10sq_jet_phi_, 200, 0);
    std::fill_n(Calo_match12sq_jet_phi_, 200, 0);

    Calo_match8circ_jet_multi_= 0;
    Calo_match9circ_jet_multi_= 0; 
    Calo_match10circ_jet_multi_= 0;
    Calo_match12circ_jet_multi_= 0;
    Calo_match8sq_jet_multi_= 0; 
    Calo_match9sq_jet_multi_= 0;
    Calo_match10sq_jet_multi_= 0; 
    Calo_match12sq_jet_multi_= 0; 


    L1_tow_8circ_jet_indx_    = new std::vector<int>();
    L1_tow_8circ_jet_multi_  = 0;
    L1_tow_8circ_scaledrho_  = 0;

    std::fill_n(L1_tow_8circ_jet_pt_, 200, 0);
    std::fill_n(L1_tow_8circ_jet_phi_, 200, 0);
    std::fill_n(L1_tow_8circ_jet_eta_, 200, 0);

    L1_tow_9circ_jet_indx_    = new std::vector<int>();
    L1_tow_9circ_jet_multi_   = 0;
    L1_tow_9circ_scaledrho_   = 0;

    std::fill_n(L1_tow_9circ_jet_pt_, 200, 0);
    std::fill_n(L1_tow_9circ_jet_phi_, 200, 0);
    std::fill_n(L1_tow_9circ_jet_eta_, 200, 0);
       
    L1_tow_10circ_jet_indx_    = new std::vector<int>();
    L1_tow_10circ_jet_multi_   = 0;
    L1_tow_10circ_scaledrho_   = 0;

    std::fill_n(L1_tow_10circ_jet_pt_, 200, 0);
    std::fill_n(L1_tow_10circ_jet_phi_, 200, 0);
    std::fill_n(L1_tow_10circ_jet_eta_, 200, 0);
  
    L1_tow_12circ_jet_indx_ = new std::vector<int>();
    L1_tow_12circ_jet_multi_ = 0;
    L1_tow_12circ_scaledrho_ = 0;
    
    std::fill_n(L1_tow_12circ_jet_pt_, 200, 0);
    std::fill_n(L1_tow_12circ_jet_phi_, 200, 0);
    std::fill_n(L1_tow_12circ_jet_eta_, 200, 0);

    L1_tow_8sq_jet_indx_  = new std::vector<int>();
    L1_tow_8sq_jet_multi_ = 0;
    L1_tow_8sq_scaledrho_ = 0;
    std::fill_n(L1_tow_8sq_jet_pt_, 200, 0);
    std::fill_n(L1_tow_8sq_jet_phi_, 200, 0); 
    std::fill_n(L1_tow_8sq_jet_eta_, 200, 0); 
    
    L1_tow_9sq_jet_indx_  = new std::vector<int>();
    L1_tow_9sq_jet_multi_ = 0;
    L1_tow_9sq_scaledrho_ = 0;
    std::fill_n(L1_tow_9sq_jet_pt_, 200, 0);
    std::fill_n(L1_tow_9sq_jet_phi_, 200, 0); 
    std::fill_n(L1_tow_9sq_jet_eta_, 200, 0); 
                       
    L1_tow_10sq_jet_indx_ = new std::vector<int>();
    L1_tow_10sq_jet_multi_ = 0;
    L1_tow_10sq_scaledrho_ = 0;
    std::fill_n(L1_tow_10sq_jet_pt_, 200, 0);
    std::fill_n(L1_tow_10sq_jet_phi_, 200, 0);
    std::fill_n(L1_tow_10sq_jet_eta_, 200, 0);  
    
    L1_tow_12sq_jet_indx_ = new std::vector<int>();
    L1_tow_12sq_jet_multi_ = 0;
    L1_tow_12sq_scaledrho_ = 0;
    std::fill_n(L1_tow_12sq_jet_pt_, 200, 0);
    std::fill_n(L1_tow_12sq_jet_phi_, 200, 0);
    std::fill_n(L1_tow_12sq_jet_eta_, 200, 0);    
}

inline void ClusterTree::clearVectors(){

  Calojet4_mom_  ->clear();
  Calojet4_mom_PU_  ->clear();
  Calojet_indx_  ->clear();
  Calojet_area_  ->clear();
  /*
  L1_tow_8circ_jet_indx_ ->clear();
  L1_tow_9circ_jet_indx_ ->clear();
  L1_tow_10circ_jet_indx_ ->clear();
  L1_tow_12circ_jet_indx_ ->clear();
  L1_tow_8sq_jet_indx_  ->clear();
  L1_tow_9sq_jet_indx_    ->clear();
  L1_tow_10sq_jet_indx_ ->clear();
  L1_tow_12sq_jet_indx_ ->clear();
  */ 
}  
inline void ClusterTree::resizeVectors(int size){
   
  Calojet4_mom_     ->resize(size);
  Calojet_indx_     ->resize(size);
/*
  L1_tow_8circ_jet_indx_ ->resize(size);
  L1_tow_9circ_jet_indx_ ->resize(size);
  L1_tow_10circ_jet_indx_->resize(size);
  L1_tow_12circ_jet_indx_->resize(size);
  
  L1_tow_8sq_jet_indx_   ->resize(size);
  L1_tow_9sq_jet_indx_   ->resize(size);
  L1_tow_10sq_jet_indx_  ->resize(size);
  L1_tow_12sq_jet_indx_  ->resize(size);
*/
} 
  
inline double ClusterTree::Get(std::string value)
{ 
    if(value=="event"            ) { return this->event_;	   }
    if(value=="run"              ) { return this->run_;	           }
    if(value=="lumi"             ) { return this->lumi_;	   }
  
    if(value=="rho_Calo"         ) { return this->rho_Calo_;        }
    if(value=="rho_PF"           ) { return this->rho_PF_;          }
    if(value=="nvtx"             ) { return this->nvtx_;	   }
  
    return -9999.; 
} 

inline std::vector<std::string> 
ClusterTree::Compare(ClusterTree* value, double prec){
    std::vector<std::string> fails;
    // this should alway fit with ultimate precision
    if( this->event_ != value->event_ ){ fails.push_back( "event" ); }
    if( this->run_   != value->run_   ){ fails.push_back( "run"   ); }
    if( this->lumi_  != value->lumi_  ){ fails.push_back( "lumi"  ); }
  
    // check within  (relative) precision
    for(std::vector<std::string>::const_iterator var=variables_.begin(); var!=variables_.end(); ++var){
        if( fabs(Get(*var)-value->Get(*var))/Get(*var)>prec ) fails.push_back(*var);
    }
    return fails;
} 


#endif
