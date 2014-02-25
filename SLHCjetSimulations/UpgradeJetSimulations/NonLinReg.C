#define REG_cxx 
#include "REG.h" 
#include <TH2.h> 
#include <TStyle.h> 
#include <TCanvas.h> 
#include <TTree.h> 
void REG::Loop() {

  if (fChain == 0) return;


  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  TMVA::Tools::Instance();
  std::cout << std::endl;
  std::cout << "==> Start TMVARegression" << std::endl;
  TString outfileName( "TMVAReg.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile, 
					      "!V:!Silent:Color:DrawProgressBar" );

  factory->AddVariable( "l1Pt", "Variable 1", "units", 'F' );
  factory->AddVariable( "l1Eta", "Variable 2", "units", 'F' );
  factory->AddVariable( "l1Phi", "Variable 1", "units", 'F' );
//  factory->AddVariable( "RhoL1", "Variable 2", "units", 'F' );

  
  factory->AddTarget( "ratio" ); 
  factory->AddRegressionTree(fChain, 1.0 );
  TCut mycut = "";
  factory->PrepareTrainingAndTestTree( mycut, 
				       "nTrain_Regression=10000:nTest_Regression=0:SplitMode=Block:NormMode=NumEvents:!V" );

  factory->BookMethod( TMVA::Types::kBDT, "BDT",
		       "!H:!V:NTrees=100:nEventsMin=5:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );
  factory->TrainAllMethods();
 // ---- Evaluate all MVAs using the set of test events
  //   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
  //  factory->EvaluateAllMethods();    

   // --------------------------------------------------------------
   
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;      

   delete factory;
}
