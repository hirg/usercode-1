
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 26 17:44:53 2012 by ROOT version 5.34/00
// from TTree t1/t1
// found on file: L1Tree.root
//////////////////////////////////////////////////////////

#ifndef REG_h
#define REG_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class REG {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *fChain_rho;
   Int_t           fCurrent; //!current Tree number in a TChain

   Float_t eta_;
   float phi_;
   float pt_;
   float rho_;
   float ratio_;
   // Declaration of leaf types
   Float_t         RhoRaw;
   Float_t         RhoOff;
   Float_t         mva_rho;
   Float_t         offPt;
   Float_t         offEta;
   Float_t         offPhi;
   Float_t         l1Pt;
   Float_t         l1Eta;
   Float_t         l1Phi;

   // List of branches
   TBranch        *b_RhoRaw;   //!
   TBranch        *b_mva_rho;
   TBranch        *b_RhoOff;   //!
   TBranch        *b_offPt;   //!
   TBranch        *b_offEta;   //!
   TBranch        *b_offPhi;   //!
   TBranch        *b_l1Pt;   //!
   TBranch        *b_l1Eta;   //!
   TBranch        *b_l1Phi;   //!

   REG(TTree *tree=0 , TTree *tree_rho=0);
   virtual ~REG();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init_(TTree *tree);
   virtual void     Init_rho_(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef REG_cxx
REG::REG(TTree *tree_ , TTree *tree_rho_ ) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.


   if (tree_ == 0 ) {

      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("L1Tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("L1Tree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("L1Tree.root:/demo");
      dir->GetObject("t1",tree_);
  }
  Init_(tree_);

//    if (tree_rho_ == 0 ) {
//       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("L1Tree.root");
//       if (!f || !f->IsOpen()) {
//          f = new TFile("L1Tree.root");
//       }
//       TDirectory * dir = (TDirectory*)f->Get("L1Tree.root:/demo");
//       dir->GetObject("t2",tree_rho_);
// 
//    }
// 
//    Init_rho_(tree_rho_);

}

REG::~REG()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t REG::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t REG::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void REG::Init_(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
//   fChain->SetBranchAddress("RhoL1", &RhoL1, &b_RhoL1);
//    fChain->SetBranchAddress("RhoRaw", &RhoRaw, &b_RhoRaw);
//    fChain->SetBranchAddress("RhoOff", &RhoOff, &b_RhoOff);
//    fChain->SetBranchAddress("MVA_Rho", &mva_rho, &b_mva_rho);
//    fChain->SetBranchAddress("offPt", &offPt, &b_offPt);
//    fChain->SetBranchAddress("offEta", &offEta, &b_offEta);
//    fChain->SetBranchAddress("offPhi", &offPhi, &b_offPhi);
   fChain->SetBranchAddress("l1Pt", &l1Pt, &b_l1Pt);
//   fChain->SetBranchAddress("l1Eta", &l1Eta, &b_l1Eta);
//   fChain->SetBranchAddress("l1Phi", &l1Phi, &b_l1Phi);
   Notify();


  cout<<"initialized"<<endl;
}

// void REG::Init_rho_(TTree *tree)
// {
//    // The Init() function is called when the selector needs to initialize
//    // a new tree or chain. Typically here the branch addresses and branch
//    // pointers of the tree will be set.
//    // It is normally not necessary to make changes to the generated
//    // code, but the routine can be extended by the user if needed.
//    // Init() will be called many times when running on PROOF
//    // (once per file to be processed).
// 
//    // Set branch addresses and branch pointers
//    if (!tree) return;
//    fChain = tree;
//    fCurrent = -1;
//    fChain->SetMakeClass(1);
// 
//    fChain->SetBranchAddress("RhoRaw", &RhoRaw, &b_RhoRaw);
//    fChain->SetBranchAddress("RhoOff", &RhoOff, &b_RhoOff);
//    Notify();
//    cout<<"rho"<<endl;
// }

Bool_t REG::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void REG::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t REG::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef REG_cxx

