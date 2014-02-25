// -*- C++ -*-
//
// Package:    CalibTowerJetProducer
// Class:      CalibTowerJetProducer
// 
/**\class CalibTowerJetProducer CalibTowerJetProducer.cc JetSLHC/CalibTowerJetProducer/src/CalibTowerJetProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michele Pioppi,510 1-017,+41227671792,
//         Created:  Mon Aug 20 13:22:58 CEST 2012
// $Id: CalibTowerJetProducer.cc,v 1.4 2012/08/22 11:19:04 pioppi Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "TLorentzVector.h"
using namespace std;
using namespace edm;
using namespace l1slhc;
using namespace l1extra;
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


///Function to assign the absolute L2L3 energy corrections in bins of eta 
void AssignParameters(vector<double> &eta_,
                      double goodPt,
                      double goodEta) 
{

  if(  abs(goodEta) > 0 && abs(goodEta) < 0.348 &&   goodPt > 4.61818 &&    goodPt < 438.647 ){

    eta_.push_back(1.03803);   eta_.push_back(5.52273);   eta_.push_back(3.12868);    eta_.push_back(1.09865);

  }else if (abs(goodEta) > 0.348 && abs(goodEta) < 0.695 && goodPt > 4.33423 &&    goodPt < 411.04 ){

    eta_.push_back(1.01282);   eta_.push_back(5.66508);   eta_.push_back(3.03777);    eta_.push_back(1.02177);

  }else if (abs(goodEta) > 0.695 && abs(goodEta) < 1.044 && goodPt > 3.94817 &&    goodPt < 416.864 ){

    eta_.push_back(1.04077);   eta_.push_back(6.14283);   eta_.push_back(3.09797);    eta_.push_back(1.03392);

  }else if (abs(goodEta) > 1.044 && abs(goodEta) < 1.392 && goodPt > 3.54529 &&    goodPt < 379.031 ){

    eta_.push_back(1.08491);   eta_.push_back(7.04365);   eta_.push_back(3.16871);    eta_.push_back(1.10764);

  }else if (abs(goodEta) > 1.392 && abs(goodEta) < 1.74 && goodPt > 4.9 &&    goodPt < 426.055 ){

    eta_.push_back(-0.755097); eta_.push_back(0.0643375); eta_.push_back(0.0187837);  eta_.push_back(-0.9887);

  }else if (abs(goodEta) > 1.74 && abs(goodEta) < 2.172 && goodPt > 4.40558 &&    goodPt < 287.101 ){

    eta_.push_back(1.45345);   eta_.push_back(8.69944);   eta_.push_back(4.13872);    eta_.push_back(1.93416);

  }else if (abs(goodEta) > 2.172 && abs(goodEta) < 3 &&   goodPt > 5.15565 &&    goodPt < 168.777 ){

    eta_.push_back(1.48724);   eta_.push_back(8.41225);   eta_.push_back(4.3892);     eta_.push_back(2.4617);

  }

  return;
}






//
// class declaration
//

class CalibTowerJetProducer : public edm::EDProducer {
   public:
      explicit CalibTowerJetProducer(const edm::ParameterSet&);
      ~CalibTowerJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  ParameterSet conf_;
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
CalibTowerJetProducer::CalibTowerJetProducer(const edm::ParameterSet& iConfig):
  conf_(iConfig)
{


  produces<L1TowerJetCollection>("CalibJets");
  produces<float>("Rho");
  produces< L1JetParticleCollection >( "Tower" ) ;
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


CalibTowerJetProducer::~CalibTowerJetProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CalibTowerJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  auto_ptr< L1TowerJetCollection > outputColl(new L1TowerJetCollection());
  auto_ptr< L1JetParticleCollection > outputExtra(new L1JetParticleCollection());
  
  Handle<L1TowerJetCollection > fc8;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle8"), fc8);

  float outrho;

////////////////////////////////
// Calibrate rho. Use the median jet energy / area and apply 3 rounds of calibrations:

    vector<double> JetEnergies;
    for (L1TowerJetCollection::const_iterator il1 = fc8->begin(); il1!= fc8->end() ; ++il1 ){
      if( abs(il1->p4().eta() )<3){
        JetEnergies.push_back(il1->p4().Pt());
      } 
    }
    double areaPerJet = 52 * (0.087 * 0.087) ;
    float raw_rho = ( median( JetEnergies ) / areaPerJet );
    //cout<<"median rho producer: "<<raw_rho<<endl;
    //scale rho: 3 rounds of calibs
    double corr_rho(0);
    double corr_rho_1(0);
    double corr_rho_2(0);
    if( raw_rho<15 ){ //calibrate values for raw_rho<15: when above this is usually 0 PU event so rho=0;

      ///1st round of calibs: take_l1raw_l1corr_
      corr_rho = 0;
      vector<double> L1rawrho_L1_corrected_rho;
      if(raw_rho<3) corr_rho = raw_rho; 
      else if(raw_rho >= 3 && raw_rho<11 ){
        //pol6(0)+[7]*pow(x,[8])
        L1rawrho_L1_corrected_rho.push_back(4.73071);      L1rawrho_L1_corrected_rho.push_back(0.249325); 
        L1rawrho_L1_corrected_rho.push_back(-0.00768474);  L1rawrho_L1_corrected_rho.push_back(0.0438207); 
        L1rawrho_L1_corrected_rho.push_back(-0.00489053);  L1rawrho_L1_corrected_rho.push_back(0.000293992); 
        L1rawrho_L1_corrected_rho.push_back(-8.41092e-06); L1rawrho_L1_corrected_rho.push_back(-0.0434629); 
        L1rawrho_L1_corrected_rho.push_back(2.34468); 
      
        for(int a=0; a<7; ++a){
          corr_rho+=L1rawrho_L1_corrected_rho[a]*pow(raw_rho,a);
        }
        corr_rho+=L1rawrho_L1_corrected_rho[7]*pow(raw_rho,L1rawrho_L1_corrected_rho[8]);

      }else if(raw_rho>=11){
        //pol3(0)
        L1rawrho_L1_corrected_rho.push_back(10.1793);    L1rawrho_L1_corrected_rho.push_back(0.744746); 
        L1rawrho_L1_corrected_rho.push_back(-0.0476675); L1rawrho_L1_corrected_rho.push_back(0.000904294); 

        for(int a=0; a<4; ++a){
          corr_rho+=L1rawrho_L1_corrected_rho[a]*pow(raw_rho,a);
        }
      }else corr_rho = 0;
      //cout<<"in producer: rho1 = "<<corr_rho<<endl;
      ///second round of rho calibs:      L1corrected_rho_again = Offline_correctedL1rho(L1corrected_rho, 0);
      vector <double> L1_corrected_rho_Calorho;
      corr_rho_1 = 0;
  
      if(corr_rho<3){
        L1_corrected_rho_Calorho.push_back(-0.642985); L1_corrected_rho_Calorho.push_back(1.50127); 
        L1_corrected_rho_Calorho.push_back(0.308268);
      }else if(corr_rho >= 3 && corr_rho < 15){
        //pol5(0)
        L1_corrected_rho_Calorho.push_back(-6.04152);   L1_corrected_rho_Calorho.push_back(9.80707); 
        L1_corrected_rho_Calorho.push_back(-2.82378);   L1_corrected_rho_Calorho.push_back(0.382044);
        L1_corrected_rho_Calorho.push_back(-0.0233116); L1_corrected_rho_Calorho.push_back(0.000522253); 
      }else if (corr_rho>=15){
        //pol4(0)
        L1_corrected_rho_Calorho.push_back(19.0038);    L1_corrected_rho_Calorho.push_back(-1.37124); 
        L1_corrected_rho_Calorho.push_back(0.00248121); L1_corrected_rho_Calorho.push_back(0.00197638); 
        L1_corrected_rho_Calorho.push_back(-3.83044e-05); 
      }
    
      if(L1_corrected_rho_Calorho.size()>0){  //if raw rho in the calibrated range, assign value to corr_rho_1
        for(unsigned int a=0; a<L1_corrected_rho_Calorho.size(); a++){
          corr_rho_1 += L1_corrected_rho_Calorho[a]*pow(corr_rho,a);
        }
      }else corr_rho_1 = corr_rho;            //else keep original value
      //cout<<"in producer: rho2 = "<<corr_rho_1<<endl;
      
      ///third round of calibs
      vector<double> L1_corrected_rho_Calorho_2;
      corr_rho_2 = 0;
      //pol9(0)+[10]*pow(x,[11])
      if(corr_rho_1>=5 && corr_rho_1<13){
        L1_corrected_rho_Calorho_2.push_back(24.2307); L1_corrected_rho_Calorho_2.push_back(-8.89551); 
        L1_corrected_rho_Calorho_2.push_back(1.52808); L1_corrected_rho_Calorho_2.push_back(-0.103572);  
        L1_corrected_rho_Calorho_2.push_back(0.00255142); 
      }else if(corr_rho_1<5){
        L1_corrected_rho_Calorho_2.push_back(2.02139); L1_corrected_rho_Calorho_2.push_back(0.509741); 
        L1_corrected_rho_Calorho_2.push_back(0.131989); L1_corrected_rho_Calorho_2.push_back(-0.199705);
        L1_corrected_rho_Calorho_2.push_back(0.0376892); 
      }
  
      if( L1_corrected_rho_Calorho_2.size() > 0){ //if raw rho in the calibrated range, assign value to corr_rho_2
        for(unsigned int a=0; a<L1_corrected_rho_Calorho_2.size(); a++){
          corr_rho_2 += L1_corrected_rho_Calorho_2[a]*pow(corr_rho_1,a);
        }
      }else corr_rho_2 = corr_rho_1;            //else keep original value
    
    } //raw rho > 15: if not, corr_2_rho = 0;

///////////////////////////////////////////////////
// SET VALUE OF RHO USING THESE CALIBRATIONS
  //cout<<"in producer: rho3 = "<<corr_rho_2<<endl;
  outrho=corr_rho_2;

///////////////////////////////////////////////////
// Now PU subtract and apply calibration to the Pt
///////////////////////////////////////////////////

  TLorentzVector PF;
  for (L1TowerJetCollection::const_iterator il =fc8->begin();
       il!=fc8->end();++il){
    L1TowerJet h=(*il);

    //DEFINE CORRECT ENERGY, AND ALL THE OTHER VARIABLES NEEDED 
    //TO CREATE A TOWERJET 

    //PU subtract using calibrated rho
    double areaPerJet = 52 * (0.087 * 0.087); //8x8 circular jet made of 52 towers
    double PU_sub_pt = il->p4().Pt() - corr_rho_2*areaPerJet;

    //Calculate correction factor using official jet calibration package; fill vector above: 
    vector<double> correction_parameters;  
    AssignParameters( correction_parameters , PU_sub_pt , il->p4().eta() ); 

    double correction_factor(999);

    if(correction_parameters.size()>0){
      correction_factor = correction_parameters[0]+
        correction_parameters[1]/(pow(log10(PU_sub_pt),correction_parameters[2])+correction_parameters[3]);
    }

    double pt_corrected(PU_sub_pt);
    double pt_corrected_again(PU_sub_pt);
   // cout<<"PU sub pt:"<<PU_sub_pt<<endl;
   // cout<<"correction factor: "<<correction_factor<<endl;
    if(correction_factor<999){ //good jet and can correct
      //calculate the calibrated pt 
      pt_corrected = PU_sub_pt * correction_factor ;
    //  cout<<"corrected pt: "<<pt_corrected<<endl;
      //apply a small correction to ensure flat response across pt
       if(pt_corrected<20) { pt_corrected_again = -1; }
       else if(pt_corrected>20 && pt_corrected<200){
         vector<double> calibs; 
         pt_corrected_again = 0;
        //pol8(0)
        calibs.push_back(78.016);      calibs.push_back(-3.88863);     calibs.push_back(0.116187); 
        calibs.push_back(-0.00143401); calibs.push_back(1.01852e-05);  calibs.push_back(-4.32462e-08);
        calibs.push_back(1.08326e-10); calibs.push_back(-1.47494e-13); calibs.push_back(8.4071e-17); 

        for(unsigned int a=0; a<calibs.size(); ++a){ 
          pt_corrected_again+=calibs[a]*pow(pt_corrected,a);
        }

      }else pt_corrected_again = pt_corrected;
   
    }//else pt_corrected_again = PU_sub_pt with no corrections as is out of range  

    //cout<<"raw pt: "<<il->p4().Pt()<<" PU sub pt: "<<PU_sub_pt<<" calibrated pt: "<<pt_corrected_again<<endl;

    
    //////////////////////////////////////////////////////
    // SET VALUE OF LORENTZ VECTOR USING THESE CALIBRATIONS
    math::PtEtaPhiMLorentzVector p4;

    if(pt_corrected_again > 0 ){
      p4.SetCoordinates(pt_corrected_again ,il->p4().eta(),il->p4().phi(),il->p4().M());// corrTowerJet;
      h.setP4(p4); 
      outputColl->insert( il->iEta(), il->iPhi(), h ); 
      outputExtra->push_back(L1JetParticle(h.p4()));
    }
  }

  auto_ptr<float> outRho(new float(outrho));

  iEvent.put(outputColl,"CalibJets");
  iEvent.put(outRho,"Rho");
  iEvent.put(outputExtra,"Tower");
}

// ------------ method called once each job just before starting event loop  ------------
void 
CalibTowerJetProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CalibTowerJetProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
CalibTowerJetProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CalibTowerJetProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CalibTowerJetProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CalibTowerJetProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CalibTowerJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CalibTowerJetProducer);
