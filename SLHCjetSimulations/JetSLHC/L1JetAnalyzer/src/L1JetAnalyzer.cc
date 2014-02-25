//Calo jets. change all offline PF jets to offline calo jets

#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <string>
#include <stdio.h>
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

using namespace std;
using namespace edm;


// Container for L1 calibration contants
struct L1TowerCalibrationContants{
  double a,b,c;
};
struct L1TowerCalibrationContants2{
  double a,b,c,d,e;
};
struct L1AK5CalibrationContants{
  double a,b,c,d;
};
struct L1AK5CalibrationContantsRescale{
  double a,b,c,d,e,f,g,h,i,j;
};
struct L1AK5CalibrationContants2{
  double a,b,c,d,e,f,g,h,i,j;
};

//function to find mean from vector of doubles
double mean(vector<double> aVec)
{
  double sum(0);
  for(unsigned int i = 0; i< aVec.size(); ++i){
    sum+=aVec[i];
  } 
  double mean = sum / aVec.size();
  return mean;
}

//function to find median from vector of doubles
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

//Function to sort std pair by first element
struct  pair_sort: public std::binary_function<pair <double,double>, pair <double,double>, bool> {
  bool operator()(const pair <double,double>& x, const pair <double,double>& y) {
    return ( x.first > y.first ) ;
  }
};

//function to find median given a std pair
pair <double,double> median_area( vector< pair <double,double> > aVec)
{
//sort the members of the vector: the first is the energy
  sort(aVec.begin(),aVec.end(),pair_sort());
  double median(0), area(0);
  int size = aVec.size();
  if(size ==0){
    median = 0;
    area = 1;
  }
  else if(size==1){
    median = aVec[size-1].first;
    area = aVec[size-1].second;
  }
  else if( size%2 == 0 ){
    median = ( aVec[ (size/2)-1  ].first + aVec[ (size /2) ].first )/2;
    area = ( aVec[ (size/2)-1  ].second + aVec[ (size /2) ].second )/2;
  }else{
    median = aVec[ double( (size/2) ) +0.5 ].first ;
    area = aVec[ double( (size/2) ) +0.5 ].second ;
  }
  pair <double,double> MedianArea;
  MedianArea.first = median;
  MedianArea.second = area;
  return MedianArea;
}

//Function to return the eta reference
int eta_ref(double eta)
{

  double i(0), j(0);
  if(eta < -3. )  j=i; 
  ++i;

  for(double eta_low = -3. ; eta_low < 3. ; eta_low+=0.5 , ++i){
    double eta_high = eta_low + 0.5;
    //eta in range of this bin
    if( eta >= eta_low && eta < eta_high ) {
      j=i;
      break;
    } 
  }
  if(eta >= 3.0 ) j=13;

  return j;
}

//Function to return the pt reference
int pt_ref(double pt)
{
  int k=0, indx=0;
  double min=0, max=10;

  if(pt < max && pt >= min) indx=k; // 0>pt>30
  ++k;
  for(double min = 10 ; min < 250 ; min=min+10 , ++k){
    double max = min+10;
    //pt in range of this bin
    if( pt< max && pt >= min){ 
      indx=k; 
      break;
    }
  }
  if(pt > 250) indx=25;
  return indx;

}


//Function to sort Lorentz vector by Pt
struct  Lorentz_sort: public std::binary_function<TLorentzVector, TLorentzVector, bool> {
  bool operator()(const TLorentzVector& x, const TLorentzVector& y) {
    return ( x.Pt() > y.Pt() ) ;
  }
};

//Function to filter and sort tower jets given a vector containing the entire unsorted collection
//and the number of towers it contains
vector<TLorentzVector> sort_filter ( vector <TLorentzVector> jetcoll, double towers )
{

  vector <TLorentzVector> sortedjetcoll;

  if(jetcoll.size()>0){

    sort(jetcoll.begin(),jetcoll.end(),Lorentz_sort());

    sortedjetcoll.push_back( jetcoll[0] );

    //at this point the jet collection is sorted into energies.
    //need to take the first element, greatest pt, and put it into filtered jet collection.
    //take the next greatest energy. if this is overlapping the preceeding elements don't add it 

    double overlap = towers*pow(0.087,2); 
    unsigned int maxCollSize=40;
    for(unsigned k=1; k<jetcoll.size(); k++){

      bool noOverlap = true;
      for(unsigned p=0; p<sortedjetcoll.size(); p++){
        noOverlap = noOverlap && (sortedjetcoll[p].DeltaR( jetcoll[k] ) > overlap );
      } 

      // then none of the jets overlap
      if(noOverlap){sortedjetcoll.push_back( jetcoll[k] );} // so this is a good jet 
      if(sortedjetcoll.size() == maxCollSize){break;}
    }
  }
  return sortedjetcoll;
}

double rho_scale_op2_L1ak5(double rho)
{ 
  double p0(0),p1(0),p2(0),p3(0),p4(0),p5(0),p6(0),p7(0),p8(0),p9(0);    
  p0 =     0.344009;   p1 =      0.996484;  p2 =     -0.338334;  p3 =      0.124934;
  p4 =    -0.0169141;  p5 =    0.00110314;  p6 =  -3.93872e-05;  p7 =   7.89364e-07;
  p8 =   -8.3169e-09;  p9 =   3.56686e-11;

 double scaled_rho = p0+rho*p1+pow(rho,2)*p2+pow(rho,3)*p3+pow(rho,4)*p4+pow(rho,5)*p5+pow(rho,6)*p6+pow(rho,7)*p7+pow(rho,8)*p8+pow(rho,9)*p9;
 
 return scaled_rho;
}

double rho_op2_again( double rho)
{
  double p0,p1;
  p0  =     0.352535; 
  p1  =     0.963717;
  double scaled_rho = p0 + p1*(rho);
  return scaled_rho;   
}


double rho_scale_op1_L1ak5(double rho)
{ 
  double p0(0),p1(0),p2(0),p3(0),p4(0),p5(0),p6(0),p7(0),p8(0),p9(0);    
  p0      =     6.24945  ;  p1      =     -60.3155 ;  p2      =      23.0541 ;  p3      =     -3.65462 ;  
  p4      =      0.31131 ;  p5      =    -0.015488 ;  p6      =  0.000459161 ;  p7      = -7.88241e-06 ;
  p8      =  7.08122e-08 ;  p9      = -2.47611e-10 ;

  double scaled_rho = p0+rho*p1+pow(rho,2)*p2+pow(rho,3)*p3+pow(rho,4)*p4+pow(rho,5)*p5+pow(rho,6)*p6+pow(rho,7)*p7+pow(rho,8)*p8+pow(rho,9)*p9;
  return scaled_rho;
}

double rho_Rescale_op1_L1towBF(double rho, int j)
{
  double p0(0),p1(0),p2(0),p3(0),p4(0),p5(0),p6(0); 
  
  if(j==0||j==1||j==2){ //8x8 circle
    p0 = 0.386014;   p1 = 0.560608 ;  p2 = 0.04073;
  }if(j==3||j==4||j==5){
    p0 = 1.37294;    p1 = 0.428241 ;  p2 = 0.0429987;
  }if(j==6||j==7||j==8){
    p0 = 6.1046;     p1 = 8.17838;    p2 = -7.37965;       p3 = 2.76408;
    p4 = -0.495444;  p5 = 0.0423746;  p6 = -0.00139121; 
  }if(j==9||j==10||j==11){  //12x12 circ
    p0 = 6.7248;     p1 = 0.390695;   
  }  if(j==12||j==13||j==14){//8x8 sq  
    p0 = 0.0642713;  p1 = 0.745763;   p2 = 0.0244281;
  }if(j==15||j==16||j==17){//9x9 sq
    p0 = 0.786133;   p1 = 0.394278;   p2 = 0.055773;
  }if(j==18||j==19||j==20){ //10x10     sq  
    p0 = 1.00966;    p1 = 0.33125;    p2 = 0.0600329;
  }if(j==21||j==22||j==23)   {   //12x12 sq
    p0 = -0.105544;  p1 = 0.357505;   p2 = 0.0707033;
  }
  return p0+rho*p1+pow(rho,2)*p2+pow(rho,3)*p3+pow(rho,4)*p4+pow(rho,5)*p5+pow(rho,6);

}
double rho_Rescale_op2_L1towBF(double rho, int j)
{
  double p0(0),p1(0),p2(0);
  
  if(j==0||j==1||j==2){ //8x8 circle
    p0=0.00628503;  p1= 0.786142;  p2=0.0263271;
  }if(j==3||j==4||j==5){
    p0 = -0.452651; p1 = 0.892006;
  }if(j==6||j==7||j==8){
    p0=0.0450745;   p1=0.755794;   p2=0.022053;
  }if(j==9||j==10||j==11){  //12x12 circ    
    p0=-0.0164094;  p1=0.75261;    p2=0.0209709;
  }  if(j==12||j==13||j==14){//8x8 sq  
    p0=0.335104;    p1=0.69551;    p2=0.0264787;
  }if(j==15||j==16||j==17){//9x9 sq
    p0= 0.274005;   p1= 0.67706;   p2=0.0285173;
  }if(j==18||j==19||j==20){ //10x10     sq  
    p0=0.0631925;   p1=0.665901;   p2=0.0334327;
  }if(j==21||j==22||j==23)   {   //12x12 sq
    p0=0.301525;    p1=0.595759;   p2=0.037671;
  }
  return p0+rho*p1+pow(rho,2)*p2;
}


double rho_scale_op1_L1towBF(double rho, int j)
{
  double p0(0),p1(0),p2(0),p3(0),p4(0),p5(0),p6(0),p7(0),p8(0),p9(0); 
  double a(0), mean(0), sigma(0); 
    
  if(j==0||j==1||j==2){ //8x8 circle
    p0  =      4.65023 ;    p1  =     -28.2042 ;    p2  =      11.5195 ;    p3  =     -2.07126 ;
    p4  =       0.2076 ;    p5  =   -0.0124822 ;    p6  =  0.000459418 ;    p7  = -1.01333e-05 ;
    p8  =  1.22991e-07 ;    p9  = -6.31718e-10 ;
  }if(j==3||j==4||j==5){
    p0  =     -2646.68;     p1  =      1253.88;     p2  =     -256.808;     p3  =      29.7449; 
    p4  =     -2.13998;     p5  =    0.0991147;     p6  =  -0.00295938;     p7  =  5.50565e-05; 
    p8  = -5.80621e-07;     p9  =   2.6513e-09; 
  }if(j==6||j==7||j==8){
    a     = 9.35932e+00;       mean  = 1.26950e+01;        sigma = 3.32765e+00;  
  }if(j==9||j==10||j==11){  //12x12 circ    
    a     = 9.94674e+00;       mean  = 9.54629e+00;        sigma = 2.59992e+00;    
  }  if(j==12||j==13||j==14){//8x8 sq  
    p0   =      2.19715;    p1   =     -2.41229;    p2   =     -1.29498;    p3   =     0.427967;
    p4   =   -0.0460086;    p5   =   0.00243636;    p6   = -6.93362e-05;    p7   =  1.01877e-06;
    p8   = -6.07893e-09;
  }if(j==15||j==16||j==17){//9x9 sq
    p0   =      3.73926;     p1   =     -18.3593;     p2   =       4.0438;     p3   =    -0.163512; 
    p4   =    -0.025613;     p5   =   0.00329448;     p6   = -0.000165975;     p7   =  4.35768e-06; 
    p8   = -5.90098e-08;     p9   =  3.26298e-10; 
  }if(j==18||j==19||j==20){ //10x10     sq  
    p0   =      6.48213;     p1   =     -47.1878;     p2   =      18.3842;    p3   =     -2.94168;
    p4   =     0.255431;     p5   =   -0.0133051;     p6   =  0.000428748;    p7   =  -8.3886e-06;
    p8   =  9.14491e-08;     p9   = -4.26541e-10;
  }if(j==21||j==22||j==23)   {   //12x12 sq
    p0    =      4.24352;    p1    =     -24.4068;    p2    =      12.1743;    p3    =       -2.347;
    p4    =     0.236346;    p5    =   -0.0139171;    p6    =  0.000497757;    p7    =   -1.066e-05;
    p8    =  1.25827e-07;    p9    = -6.29896e-10;    
  } 
  
  if(a!=0){
    return a*exp(- ((pow( (rho-mean),2 ) )/2*sigma*sigma) );
  }else if(p0!=0) 
  {  return     p0+rho*p1+pow(rho,2)*p2+pow(rho,3)*p3+pow(rho,4)*p4+pow(rho,5)*p5+pow(rho,6)*p6+pow(rho,7)*p7+pow(rho,8)*p8+pow(rho,9)*p9;
  }else return rho; 
  
}

double rho_scale_op2_L1towBF(double rho, int j)
{
  double p0(0),p1(0),p2(0),p3(0),p4(0),p5(0),p6(0),p7(0),p8(0),p9(0); 
  double a(0), mean(0), sigma(0);   
  if(j==0||j==1||j==2){ 
    p0  =      1.97713;    p1  =      1.57091;     p2  =     0.290832;       p3  =    -0.115215;
    p4  =    0.0120728;    p5  = -0.000614244;     p6  =  1.68188e-05;       p7  = -2.38646e-07;
    p8  =  1.38035e-09;
  }if(j==3||j==4||j==5){
    p0     =      2.03717; p1     =       1.5874;  p2     =     0.503526;    p3     =    -0.194814;
    p4     =    0.0228755; p5     =  -0.00135892;  p6     =  4.58219e-05;    p7     = -8.86491e-07;
    p8     =  9.13297e-09; p9     =  -3.8559e-11;
  }if(j==6||j==7||j==8){
    p0  =      2.04799;    p1  =      1.51676;     p2  =     0.702718;       p3  =    -0.265936;
    p4  =      0.03302;    p5  =  -0.00210803;     p6  =   7.7203e-05;       p7  = -1.64081e-06;
    p8  =  1.88477e-08;    p9  = -9.06008e-11;
  }if(j==9||j==10||j==11){  //12x12 circ
    p0   =       2.4159;   p1   =      2.02434;    p2   =       0.5504;      p3   =    -0.260468;
    p4   =    0.0348197;   p5   =  -0.00234878;    p6   =  9.05368e-05;      p7   = -2.02555e-06;
    p8   =  2.45321e-08;   p9   = -1.24604e-10;
  }if(j==12||j==13||j==14){//8x8 sq 
    p0  =      1.94963;    p1  =     0.799762;     p2  =      1.05884;      p3  =    -0.356832;
    p4  =    0.0467175;    p5  =  -0.00327205;     p6  =  0.000133561;      p7  = -3.18654e-06;
    p8  =  4.12263e-08;    p9  = -2.23528e-10;
  }if(j==15||j==16||j==17){//9x9 sq
    p0  =      1.85128;    p1  =      1.29029;     p2  =      1.04034;      p3  =    -0.406005;
    p4  =    0.0568464;    p5  =  -0.00415508;     p6  =  0.000174695;      p7  = -4.25775e-06;
    p8  =  5.59643e-08;    p9  = -3.07123e-10;
  }if(j==18||j==19||j==20){ //10x10     sq
    p0  =       1.6354;    p1  =      2.88123;     p2  =     0.214806;      p3  =    -0.233812; 
    p4  =    0.0370448;    p5  =   -0.0027886;     p6  =  0.000117193;      p7  = -2.81934e-06; 
    p8  =  3.63607e-08;    p9  = -1.95213e-10;
  }if(j==21||j==22||j==23)   {   //12x12 sq 
    p0  =      1.38482;    p1  =      4.27415;     p2  =     -0.44515;      p3  =    -0.113427; 
    p4  =    0.0253515;    p5  =  -0.00211578;     p6  =  9.35611e-05;      p7  = -2.32261e-06; 
    p8  =  3.06279e-08;    p9  = -1.67308e-10; 

  }
  if(a!=0){
    return a*exp(- ((pow( (rho-mean),2 ) )/2*sigma*sigma) );
  }else if(p0!=0) 
  {  return p0+rho*p1+pow(rho,2)*p2+pow(rho,3)*p3+
      pow(rho,4)*p4+pow(rho,5)*p5+pow(rho,6)*p6+
      pow(rho,7)*p7+pow(rho,8)*p8+pow(rho,9)*p9;
  }else return rho; 

}   



//
// class declaration
//
class L1JetAnalyzer : public edm::EDAnalyzer {
public:
  explicit L1JetAnalyzer(const edm::ParameterSet&);
  ~L1JetAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  
private:
   virtual void beginJob() ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
  
   virtual void beginRun(edm::Run const&, edm::EventSetup const&);
   virtual void endRun(edm::Run const&, edm::EventSetup const&);
   virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
   virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
   L1TowerCalibrationContants getL1TowerCalibrationContants(int l1JetColId, int l1JetEtaRange);
   L1TowerCalibrationContants2 getL1TowerCalibrationContants2(int l1JetColId, int l1JetEtaRange);
   L1AK5CalibrationContants getL1AK5CalibrationContants(int l1JetColId, int l1JetEtaRange);
   L1AK5CalibrationContantsRescale getL1AK5CalibrationContantsRescale(int l1JetColId, int l1JetEtaRange, double l1JetPt);
   L1AK5CalibrationContants2 getL1AK5CalibrationContants2(int l1JetColId, int l1JetEtaRange, double l1JetPt);

  // ----------member data ---------------------------

  ParameterSet conf_;

  //Offline calo Jets compared to Current L1 jets 
  TH1F* Offline_L1_eta;
  TH1F* Offline_L1_phi;
  TH1F* Offline_L1_pt;
  TH1F* Offline_L1_pt_PU;
  TH2F* Pt_res_vertices_op1_CurrentL1;

  //Offline calo Jets compared to L1 AK5 jets 
  TH1F* Offline_L1ak5_eta;
  TH1F* Offline_L1ak5_phi;
  TH1F* Offline_L1ak5_pt;
  TH1F* Offline_L1ak5_pt_PU;
  TH1F* Offline_L1ak5_pt_scale;
  TH1F* Offline_L1ak5_pt_PU_scale; 

  TH1F* Offline_L1ak5_Op1_pt_res;
  TH1F* Offline_L1ak5_Op2_pt_res;
  TH1F* Offline_L1ak5_Op1Rescale_pt_res;  

  //scale: L1ak5 op1, op2, rescaleop1
  vector <TProfile*> L1ak5_Pt_scale_eta;
  vector <TProfile*> L1ak5_Pt_scale_eta_op2; 
  vector <TProfile*> L1ak5_Pt_scale_op1_rescale;

  //Rho resolutions
  TH1F* Rho_L1AK5_Offline_scaled;
  TH1F* Rho_L1AK5_Offline_nonscaled;
  TH1F* Rho_L1AK5_Offline_Option1;
  TH1F* Rho_L1AK5_Offline_Option2;
  TH1F* Rho_L1AK5_Offline_Option2Again;  
    //PU sub of these rhos:
  TH1F* Offline_L1ak5_pt_PU_scale_median_scaled;
  TH1F* Offline_L1ak5_pt_PU_scale_median_nonscaled; 
  //Pt res as a function of PU
  TH2F* Pt_res_vertices_op1;
  TH2F* Pt_res_vertices_op2;

  vector<string> TowerJet;
  //Offline calo Jets compared to bitonic filtered and sorted Tower level L1 jets 
  vector <TH1F*> L1tow_BitFilt_eta ;        
  vector <TH1F*> L1tow_BitFilt_phi;
  vector <TH1F*> L1tow_BitFilt_pt;
  vector <TH1F*> L1tow_BitFilt_pt_PU;
  vector <TH1F*> L1tow_BitFilt_pt_scale;
  vector <TH1F*> L1tow_BitFilt_pt_PU_scale; 
  vector <TH1F*> Offline_L1tow_Op1_pt_res;
  vector <TH1F*> Offline_L1tow_Op2_pt_res;

  //binned res in terms of pt
  //binned in terms of PU
  vector <TH2F*> Off_BFL1tow_Pt_PU_op1_binpt;
  vector <TH2F*> Off_BFL1tow_Pt_PU_op2_binpt;

  vector <TH1F*> L1Tow_BitFilt_rho_res_Option1;
  vector <TH1F*> L1Tow_BitFilt_rho_res_Option2;

  //scale
  vector <TProfile*> L1towBF_Pt_scale_eta_op1;
  vector <TProfile*> L1towBF_Pt_scale_eta_op2;
  vector <TProfile*> L1towBF_Pt_scale_eta_op1Rescale;

  //rho
  vector <TH1F*> L1Tower_rho_res_BitFilt_scaled;
  vector <TH1F*> L1Tower_rho_res_BitFilt_noscaled;
  vector <TH1F*> L1Tower_rho;

  //Trigger turn on curves for L1 jets
  vector <TH1F*> L1Tower_num;
  vector <TH1F*> L1Tower_den;
  TH1F* L1AK5_num;
  TH1F* L1AK5_den;
  TH1F* CurrL1_num;
  TH1F* CurrL1_den; 

  //Pt, eta, PU binning:
  vector <string> etast, ptst, PU;

  //rho scaling:
  TProfile* Rho_Scaling_Option1;
  TProfile* Rho_Scaling_Option2;
  TProfile* Rho_Scaling_Option2_Again;
  
  vector<TProfile*> Rho_Scaling_Option1_BitFilt;
  vector<TProfile*> Rho_Scaling_Option2_BitFilt;
  vector<TProfile*> Rho_ReScaling_Option1_BitFilt;
  vector<TProfile*> Rho_ReScaling_Option2_BitFilt;
  vector<TProfile*> Rho_ReReScaling_Option2_BitFilt;

  //pt res in bins of pt, PU
  TH2F* Off_L1ak5_binnedPtres_op1;
  TH2F* Off_L1ak5_binnedPtres_op1Rescale;
  TH2F* Off_L1ak5_binnedPtres_op2;
  TH2F* Offline_L1_binnedPtres;
  vector <TH2F*> Off_BFL1tow_binnedPtres_op2;
  
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
L1JetAnalyzer::L1JetAnalyzer(const edm::ParameterSet& iConfig):
  conf_(iConfig)
{
  Service<TFileService> fs;
//------------------------------------------------------------------------------------------//
   string name, title;

//create vectors for eta and pt binning

  //Eta
  for(double k=-3; k<3.6; k+=0.5){
    std::ostringstream ref;
    ref << k;
    std::string refst = ref.str();
    etast.push_back( refst) ;
  }
  //Momentum
  for(int k=0; k<251; k+=10){
    std::ostringstream ref;
    ref << k;
    std::string refst = ref.str();
    ptst.push_back( refst) ;
  }
  //PU
  for(int k=0; k<55; k+=5){
    std::ostringstream ref;
    ref << k;
    std::string refst = ref.str();
    PU.push_back( refst) ;
  }

//-------------------------------------------------------------------------------------------//
  //Offline calo Jets compared to Current L1 jets 
  
  TFileDirectory iDir0 = fs->mkdir( "CurrentL1" );
  
  Offline_L1_eta   = iDir0.make<TH1F>("Offline_L1_eta", 
                                      "#eta resolution of Offline Calo jets vs Current L1 jets; #eta res; ", 150,-1.5,1.5);
  Offline_L1_phi   = iDir0.make<TH1F>("Offline_L1_phi", 
                                      "#phi resolution of Offline Calo jets vs Current L1 jets; #phi res; ", 150,-1.5,1.5);
  Offline_L1_pt    = iDir0.make<TH1F>("Offline_L1_pt", 
                                      "P_{T} resolution of Offline Calo jets vs Current L1 jets; P_{T} res; ", 150,-1.5,1.5);
  Offline_L1_pt_PU = iDir0.make<TH1F>("Offline_L1_pt_PU", 
                                      "P_{T} resolution of PU sub Offline Calo jets vs Current L1 jets; P_{T} res; ", 150,-1.5,1.5);  
  
  Pt_res_vertices_op1_CurrentL1= iDir0.make<TH2F>("Offline_L1_pt_res_PrimVert",
      "P_{T} resolution of L1 jets vs Offline jets as a function of PU; No. of vertices; P_{T} res", 
      40,0,40, 150, -1.5, 1.5); 
  

  Offline_L1_binnedPtres = iDir0.make<TH2F>("Offline_L1_binnedPtres",
                                            "P_{T} res of Offline Calo jets vs Current L1 jet in binned P_{T}; Current L1 P_{T}; P_{T} resolution",
                                            25,0,250,150,-1.5,1.5);

  
  //Offline calo jets compared to L1 AK5 jets 
  TFileDirectory iDir1 = fs->mkdir( "L1AK5" );

  Offline_L1ak5_eta   = iDir1.make<TH1F>("Offline_L1ak5_eta", 
                                         "#eta resolution of Offline Calo jets vs L1ak5 jets; #eta res; ", 150,-1.5,1.5);
  Offline_L1ak5_phi   = iDir1.make<TH1F>("Offline_L1ak5_phi", 
                                         "#phi resolution of Offline Calo jets vs L1ak5 jets; #phi res; ", 150,-1.5,1.5);
  Offline_L1ak5_pt    = iDir1.make<TH1F>("Offline_L1ak5_pt", 
                                         "P_{T} resolution of Offline Calo jets vs L1ak5 jets; P_{T} res; ", 150,-1.5,1.5);
  Offline_L1ak5_pt_PU = iDir1.make<TH1F>("Offline_L1ak5_pt_PU", 
                                         "P_{T} resolution of PU sub Offline Calo jets vs L1ak5 jets; P_{T} res; ", 150,-1.5,1.5);    
  Offline_L1ak5_pt_scale = iDir1.make<TH1F>("Offline_L1ak5_pt_scale",
                                            "P_{T} resolution of Offline Calo jets vs scaled L1ak5 jets; P_{T} res; ", 150,-1.5,1.5);
  Offline_L1ak5_pt_PU_scale = iDir1.make<TH1F>("Offline_L1ak5_pt_PU_scale",
                                               "P_{T} resolution of PU sub Offline Calo jets vs scaled L1ak5 jets; P_{T} res; ", 150,-1.5,1.5);    

  //non scaled pu calculation
  //non scaled rho
  Offline_L1ak5_pt_PU_scale_median_nonscaled = iDir1.make<TH1F>("Offline_L1ak5_pt_PU_scale_median",
      "P_{T} resolution of PU sub Offline Calo jets vs PU sub median scaled L1ak5 jets; P_{T} res; ", 
      150,-1.5,1.5);    
  
  //scaled PU calculation
  //non scaled rho
  Offline_L1ak5_pt_PU_scale_median_scaled = iDir1.make<TH1F>("Offline_L1ak5_pt_PU_scale_median_scaled",
      "P_{T} resolution of PU sub Offline Calo jets vs PU sub median scaled L1ak5 jets; P_{T} res; ", 
      150,-1.5,1.5);    

  Offline_L1ak5_Op1_pt_res = iDir1.make<TH1F>("Offline_L1ak5_Op1_pt_res",
                                             "P_{T} resolution of L1ak5 jets vs Offline jets op1; P_{T} res; ", 
                                             150,-1.5,1.5);    
  Offline_L1ak5_Op2_pt_res = iDir1.make<TH1F>("Offline_L1ak5_Op2_pt_res",
                                              "P_{T} resolution of L1ak5 jets vs Offline jets op2; P_{T} res; ", 
                                             150,-1.5,1.5);    
  Offline_L1ak5_Op1Rescale_pt_res= iDir1.make<TH1F>("Offline_L1ak5_Op1Rescale_pt_res",
      "P_{T} resolution of L1ak5 jets vs Offline jets op1Rescale; P_{T} res; ", 
      150,-1.5,1.5);    
  
  //To scale in bins of eta
  for( unsigned int i =0; i<etast.size(); ++i){
    name  = "L1ak5_scale_uncorrOff_"+etast[i];
    title = "Scale of Non PU sub Offline Calo jets vs L1AK5 jet in P_{T} "+etast[i]+" jets; L1 AK5 P_{T}; Offline P_{T}  ";
    L1ak5_Pt_scale_eta.push_back( iDir1.make<TProfile>( name.c_str() , title.c_str()  ,200,-0,200) );
    
    name  = "L1ak5_scale_corrOff_"+etast[i];
    title = "Scale of PU sub Offline Calo jets vs L1AK5 jet in P_{T} "+etast[i]+" jets; L1 AK5 P_{T}; Offline P_{T}  ";
    L1ak5_Pt_scale_eta_op2.push_back( iDir1.make<TProfile>( name.c_str() , title.c_str()  ,200,-0,200) );

    name = "L1ak5_scale_op1Rescale_"+etast[i];
    title = "Rescale PU sub L1ak5 jets to offline PU jets, op1 "+etast[i]+" jets; L1AK5 P_{T}; Offline P_{T}";
    L1ak5_Pt_scale_op1_rescale.push_back( iDir1.make<TProfile>( name.c_str() , title.c_str()  ,200,-0,200) );

  }
  
  Off_L1ak5_binnedPtres_op1 = iDir1.make<TH2F>("Offline_L1ak5_binnedPtres_op1",
                  "P_{T} res of Offline Calo jets vs L1 AK5 jet in binned P_{T}; L1 AK5 P_{T}; P_{T} resolution",
                  25,0,250, 150,-1.5,1.5);
    
  Off_L1ak5_binnedPtres_op1Rescale = iDir1.make<TH2F>("Offline_L1ak5_binnedPtres_op1Rescale",
                  "P_{T} res of Offline Calo jets vs L1 AK5 jet in binned P_{T}; L1 AK5 P_{T}; P_{T} resolution",
                  25,0,250, 150,-1.5,1.5);
  
  Off_L1ak5_binnedPtres_op2 = iDir1.make<TH2F>("Offline_L1ak5_binnedPtres_op2",
                  "P_{T} res of Offline Calo jets vs L1 AK5 jet in binned P_{T}; L1 AK5 P_{T}; P_{T} resolution",
                   25,0,250, 150,-1.5,1.5);
  
  //L1 ak5 rho histograms
  
  TFileDirectory iDir2 = fs->mkdir( "Rho" );
    
  //using scaled L1AK5 jets
  Rho_L1AK5_Offline_scaled = iDir2.make<TH1F>("Rho_L1AK5_Offline_scaled",
      "#rho resolution using median of rhos of scaled L1ak5 jets; P_{T} res; ", 150,-1.5,1.5);
  //using non scaled L1AK5 jets
  Rho_L1AK5_Offline_nonscaled = iDir2.make<TH1F>("Rho_L1AK5_Offline_nonscaled",
      "#rho resolution using median of rhos of non scaled L1ak5 jets; P_{T} res; ", 150,-1.5,1.5);

  Rho_Scaling_Option1 = iDir2.make<TProfile>("Rho_Scaling_Option1" ,
      "Scale of L1 #rho (from scaled L1 jets) to offline rho: option1 ; L1 #rho; Offline #rho",
      200,-0,40);
  Rho_Scaling_Option2 = iDir2.make<TProfile>("Rho_Scaling_Option2" , 
       "Scale of L1 #rho (from unscaled L1 jets) to offline rho: option2 ; L1 #rho; Offline #rho",
       200,-0,40);

  Rho_Scaling_Option2_Again = iDir2.make<TProfile>("Rho_Scaling_Option2_again" , 
      "Scale of L1 #rho (from unscaled L1 jets) to offline rho again: option2 ; L1 #rho;Offline #rho",
      200,-0,40);
  
  
  //option 1 rho resolution
  Rho_L1AK5_Offline_Option1= iDir2.make<TH1F>("Rho_L1AK5_Offline_Option1",
      "#rho resolution using option 1; #rho res; ", 150,-1.5,1.5);
  //option 2 rho resolution
  Rho_L1AK5_Offline_Option2= iDir2.make<TH1F>("Rho_L1AK5_Offline_Option2",
      "#rho resolution using option 2; #rho res; ", 150,-1.5,1.5);
  Rho_L1AK5_Offline_Option2Again= iDir2.make<TH1F>("Rho_L1AK5_Offline_Option2Again",
      "#rho resolution using option 2 rescaled; #rho res; ", 150,-1.5,1.5);
  
  //Pt resolution as a function of PU
  Pt_res_vertices_op1  = iDir1.make<TH2F>("Offline_L1ak5_Op1_pt_res_PrimVert",
      "P_{T} resolution of L1ak5 jets vs Offline jets as a function of PU for op1; No. of vertices; P_{T} res", 
       40,-0.5,39.5, 150, -1.5, 1.5);       
  Pt_res_vertices_op2= iDir1.make<TH2F>("Offline_L1ak5_Op2_pt_res_PrimVert",
      "P_{T} resolution of L1ak5 jets vs Offline jets as a function of PU for op2; No. of vertices; P_{T} res", 
      40,-0.5,39.5, 150, -1.5, 1.5); 


  //Trigger turn on curves for Currnet L1 and L1 AK5
  L1AK5_num = iDir1.make<TH1F>( "MatchedOffline40", "Offline P_{T} of matched jets L1 AK5 > 40GeV",
                                  250,0,250);
  L1AK5_den = iDir1.make<TH1F>( "MatchedOfflineAll", "Offline P_{T} of all matched L1 jets",
                                  250,0,250); 
  //for current L1 jets
  CurrL1_num = iDir0.make<TH1F>( "MatchedOffline40", "Offline P_{T} of matched L1 jets > 40GeV",
                                   250,0,250);
  CurrL1_den = iDir0.make<TH1F>( "MatchedOfflineAll", "Offline P_{T} of all matched current L1 jets",
                                   250,0,250);
  
      
  //Bitonic filtered and sorted jets compared to current L1 jets
  TowerJet.push_back("8x8_circle");
  TowerJet.push_back("8x8_circle_mean");
  TowerJet.push_back("8x8_circle_median");
  TowerJet.push_back("9x9_circle");
  TowerJet.push_back("9x9_circle_mean");
  TowerJet.push_back("9x9_circle_median");
  TowerJet.push_back("10x10_circle");
  TowerJet.push_back("10x10_circle_mean");
  TowerJet.push_back("10x10_circle_median");
  TowerJet.push_back("12x12_circle");
  TowerJet.push_back("12x12_circle_mean");
  TowerJet.push_back("12x12_circle_median");
  TowerJet.push_back("8x8_square");
  TowerJet.push_back("8x8_square_mean");
  TowerJet.push_back("8x8_square_median");
  TowerJet.push_back("9x9_square");
  TowerJet.push_back("9x9_square_mean");
  TowerJet.push_back("9x9_square_median");
  TowerJet.push_back("10x10_square");
  TowerJet.push_back("10x10_square_mean");
  TowerJet.push_back("10x10_square_median");
  TowerJet.push_back("12x12_square");
  TowerJet.push_back("12x12_square_mean");
  TowerJet.push_back("12x12_square_median");


  for(std::vector<std::string>::iterator lIt = TowerJet.begin(); lIt != TowerJet.end(); ++lIt){
    TFileDirectory iDir = fs->mkdir( (*lIt) );

    //Bitonic sorted jets: eta, phi, pt resolution vs offline calo jets
    name = "L1tow_BitFilt_eta"+*lIt;
    title = "#eta resolution of Offline Calo jets vs L1 tower "+(*lIt)+" jets; #eta Resolution; Count";
    L1tow_BitFilt_eta.push_back(iDir.make<TH1F>( name.c_str() , title.c_str()  ,150,-1.5,1.5) );

    name = "L1tow_BitFilt_phi"+*lIt;
    title = "#phi resolution of Offline Calo jets vs L1 tower "+(*lIt)+" jets; #phi Resolution; Count";
    L1tow_BitFilt_phi.push_back(iDir.make<TH1F>( name.c_str() , title.c_str()  ,150,-1.5,1.5) );

    name = "L1tow_BitFilt_pt"+*lIt;
    title = "P_{T} resolution of Offline Calo jets vs L1 tower "+(*lIt)+" jets; P_{T} Resolution; Count";
    L1tow_BitFilt_pt.push_back(iDir.make<TH1F>( name.c_str() , title.c_str()  ,150,-1.5,1.5) );

    name = "L1tow_BitFilt_pt_PU_"+*lIt;
    title = "P_{T} resolution of PU sub Offline Calo jets vs filtered L1 tower "+(*lIt)+" jets; P_{T} Resolution; Count";
    L1tow_BitFilt_pt_PU.push_back(iDir.make<TH1F>( name.c_str() , title.c_str()  ,150,-1.5,1.5) );

    name = "L1tow_BitFilt_pt_scale_"+*lIt;
    title = "P_{T} resolution of Offline Calo jets vs L1 tower "+(*lIt)+" jets; P_{T} Resolution; Count";
    L1tow_BitFilt_pt_scale.push_back(iDir.make<TH1F>( name.c_str() , title.c_str()  ,150,-1.5,1.5) );

    name = "L1tow_BitFilt_pt_PU_scale_"+*lIt;
    title = "P_{T} resolution of PU sub Offline Calo jets vs filtered L1 tower "+(*lIt)+" jets; P_{T} Resolution; Count";
    L1tow_BitFilt_pt_PU_scale.push_back(iDir.make<TH1F>( name.c_str() , title.c_str()  ,150,-1.5,1.5) );
    
    name = "Offline_L1towBF_Op1_pt_res_"+*lIt;
    title = "P_{T} resolution using option 1 for "+(*lIt)+" jets; P_{T} Resolution; ";
    Offline_L1tow_Op1_pt_res.push_back(iDir.make<TH1F>( name.c_str() , title.c_str()  ,150,-1.5,1.5) );
    
    name = "Offline_L1towBF_Op2_pt_res_"+*lIt;    
    title = "P_{T} resolution using option 2 for "+(*lIt)+" jets; P_{T} Resolution; ";
    Offline_L1tow_Op2_pt_res.push_back(iDir.make<TH1F>( name.c_str() , title.c_str()  ,150,-1.5,1.5) );

    //scales in bins of eta
    for( unsigned int i =0; i<etast.size(); ++i){
      name  = "L1towBF_scale_pt_uncorrOff_"+etast[i];
      title = "Scale of Offline Calo jets vs L1tow jet in P_{T} "+etast[i]+" jets; L1 AK5 P_{T}; Offline P_{T}  ";
      L1towBF_Pt_scale_eta_op1.push_back( iDir.make<TProfile>( 
          name.c_str() , title.c_str()  ,200,-0,200) );

      name  = "L1towBF_scale_pt_corrOff_"+etast[i];
      title = "Scale of Offline Calo jets vs L1tow jet in P_{T} "+etast[i]+" jets;L1 AK5 P_{T};Offline P_{T}";
      L1towBF_Pt_scale_eta_op2.push_back( iDir.make<TProfile>( 
          name.c_str() , title.c_str()  ,200,-0,200) );

      name  = "L1towBF_scale_pt_RescaleOp1_"+etast[i];
      title = "Scale of Offline Calo jets vs L1tow jet in P_{T} "+etast[i]+" jets rescale op1; L1 AK5 P_{T}; Offline P_{T}  ";
      L1towBF_Pt_scale_eta_op1Rescale.push_back( iDir.make<TProfile>( 
          name.c_str() , title.c_str()  ,200,-0,200) );
    }  
    
    TFileDirectory iDirRho = iDir.mkdir( "Rho" );
    //TProfiles for scaling
    name = "Rho_Scaling_Option1_BitFilt_"+*lIt;
    title = "Scale of L1 #rho (from scaled L1 jets) to offline rho: BitFilt option1 ; L1 #rho; Offline #rho ";
    Rho_Scaling_Option1_BitFilt.push_back(iDirRho.make<TProfile>( name.c_str() , title.c_str()  ,200,-0,40));
        
    name = "Rho_Scaling_Option2_BitFilt_"+*lIt;
    title = "Scale of L1 #rho (from unscaled L1 jets) to offline rho: BitFilt option2 ; L1 #rho; Offline #rho ";
    Rho_Scaling_Option2_BitFilt.push_back(iDirRho.make<TProfile>( name.c_str() , title.c_str()  ,200,-0,40));
    
    name = "Rho_ReScaling_Option1_BitFilt_"+*lIt;
    title = "Scale of L1 #rho (from scaled L1 jets) to offline rho: BitFilt option1 rescale; L1 #rho; Offline #rho ";
    Rho_ReScaling_Option1_BitFilt.push_back(iDirRho.make<TProfile>( name.c_str() , title.c_str()  ,200,-0,40));
    
    name = "Rho_ReScaling_Option2_BitFilt_"+*lIt;
    title= "Scale of L1 #rho (from scaled L1 jets) to offline rho: BitFilt option2 rescale; L1 #rho; Offline #rho ";
    Rho_ReScaling_Option2_BitFilt.push_back(iDirRho.make<TProfile>( name.c_str() , title.c_str()  ,200,-0,40));
    
    name = "Rho_ReReScaling_Option2_BitFilt_"+*lIt;
    title = "ReReScale of L1 #rho (from scaled L1 jets) to offline rho: BitFilt option2 rescale; L1 #rho; Offline #rho ";
    Rho_ReReScaling_Option2_BitFilt.push_back(iDirRho.make<TProfile>( name.c_str() , title.c_str()  ,200,-0,40));
    
    //rho resolution
    name = "L1tow_BitFilt_rho_res_scaled"+*lIt;
    title = "#rho resolution of Offline Calo jets vs scaled L1 tower rho (scaled jets) by median method, "+(*lIt)+" jets; #rho Resolution; Count";
    L1Tower_rho_res_BitFilt_scaled.push_back(iDirRho.make<TH1F>( 
        name.c_str() , title.c_str()  ,150,-1.5,1.5) );
    
    name = "L1tow_BitFilt_rho_res_noscaled"+*lIt;
    title = "#rho resolution of Offline Calo jets vs scaled L1 tower rho (non scaled jets) by median method, "+(*lIt)+" jets; #rho Resolution; Count";   
    L1Tower_rho_res_BitFilt_noscaled.push_back(iDirRho.make<TH1F>( 
        name.c_str() , title.c_str()  ,150,-1.5,1.5) );
    
    name = "L1tow_BitFilt_rho_"+*lIt;
    title = "#rho of L1 tower rho jets by median (scaled jets) method, "+(*lIt)+" jets; #rho ; Count";
    L1Tower_rho.push_back(iDirRho.make<TH1F>( 
        name.c_str() , title.c_str()  ,200,0,40) );

    name = "L1tow_BitFilt_rho_res_Option1_"+*lIt;
    title = "#rho res of L1 tower jets by option1 method, "+(*lIt)+" jets; #rho res; ";
    L1Tow_BitFilt_rho_res_Option1.push_back( iDirRho.make<TH1F>( 
        name.c_str() , title.c_str()  ,150,-1.5,1.5) );
  //option 2 rho resolution
    name = "L1tow_BitFilt_rho_res_Option2_"+*lIt;
    title = "#rho res of L1 tower jets by option2 method, "+(*lIt)+" jets; #rho res ; ";
    L1Tow_BitFilt_rho_res_Option2.push_back( iDirRho.make<TH1F>( 
        name.c_str() , title.c_str()  ,150,-1.5,1.5) );

    TFileDirectory iDirBin = iDir.mkdir( "Binned" );
    //binned in Pt
    name ="Off_BFL1tow_binnedPtres_op2";
    title="P_{T} res of Offline Calo jets vs L1 tower jets in binned P_{T}; L1 tower jet P_{T}; P_{T} res";
    Off_BFL1tow_binnedPtres_op2.push_back( iDirBin.make<TH2F>(
        name.c_str() , title.c_str() , 25,0,250, 150,-1.5,1.5));

    //binned in terms of PU      
    name ="Off_BFL1tow_Pt_PU_op1_binpt"+(*lIt);
    title="P_{T} res of Offline Calo jets vs Current L1 jet in binned PU op1"+(*lIt)+";no. PU;P_{T} res ";
    Off_BFL1tow_Pt_PU_op1_binpt.push_back( iDirBin.make<TH2F>( 
        name.c_str() , title.c_str(),  50,-0.5,49.5, 150,-1.5,1.5) );

    name  = "Off_BFL1tow_Pt_PU_op2_binpt"+(*lIt);
    title = "P_{T} res of Offline Calo jets vs Current L1 jet in binned PU op2"+(*lIt)+";no. PU;P_{T} res";
    Off_BFL1tow_Pt_PU_op2_binpt.push_back( iDirBin.make<TH2F>( 
        name.c_str() , title.c_str(),  50,-0.5,49.5, 150,-1.5,1.5) );

    //TH1F's for the trigger turn on curves
    name="MatchedOffline40_"+(*lIt) ;
    title="Offline jet P_{T} matched to L1 tower jet above 40GeV;Offline Jet P_{T};Count ";
    L1Tower_num.push_back( iDir.make<TH1F>( name.c_str() , title.c_str(),  250,0,250) );

    name="MatchedOfflineAll_"+(*lIt) ;
    title="Offline jet P_{T} matched to L1 tower jets;Offline Jet P_{T};Count ";
    L1Tower_den.push_back( iDir.make<TH1F>( name.c_str() , title.c_str(),  250,0,250) );

  }


}


L1JetAnalyzer::~L1JetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}
//
// member functions
//

//_________________________________________________________________________________________________
// Function: analyze
//_________________________________________________________________________________________________
void L1JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool evValid = true;
  
  const JetCorrector* jetCorrector = JetCorrector::getJetCorrector(conf_.getParameter<string>("JetCorrector"),iSetup);

  //OFFLINE CALO JETS
  edm::Handle<reco::CaloJetCollection> Calojets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CaloJets"), Calojets);
  if(!Calojets.isValid()){evValid=false;}
  
  //RHO OF OFFLINE CALO JETS
  edm::Handle<double> rhoCALOjets;  
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("RhoCaloJets"), rhoCALOjets);
  double rhoCalo = *rhoCALOjets;
  if(!rhoCALOjets.isValid()){evValid=false;}

    //PF COLLECTION FROM ECAL AND HCAL DIGI PRIMITIVES: L1 AK5 JETS
  edm::Handle<reco::PFJetCollection> L1ak5jets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("ReRecoPfJets"), L1ak5jets);
  if(!L1ak5jets.isValid()){evValid=false;}   
    
  //RHO OF L1 AK5 JETS
/*  edm::Handle<double> rhoL1AK5jets;  
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>(" RhoL1PFJets"), rhoL1AK5jets);
  double rhoL1ak5 = *rhoL1AK5jets;
  if(!rhoL1AK5jets.isValid()){evValid=false;}
*/
  //L1EXTRA PARTICLES 
  vector<InputTag>  l1extraparticles= conf_.getParameter< vector < InputTag > >("extra"); 
    //This vector called l1extraparticles references "extra" in the config file: contains Central, Forward, Tau jets

  //--------------------------------Vertex Multiplicity--------------------------------------//

  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel("offlinePrimaryVertices", vtx);
  int PVSIZE=vtx->size();
  // cout<<"PVSIZE: "<<PVSIZE<<endl;

  //-------------------------------------------------------------------------------------------------//
  
  double DR    = 0.5;
  double PTlim = 30;

  //-------------------------------------------------------------------------------------------//
  //Jet collection already filtered
  
  std::vector< const  EtaPhiContainer<l1slhc::L1TowerJet>* > FilteredBitonicSortJetColl;
  edm::Handle< EtaPhiContainer <l1slhc::L1TowerJet> > hk;
  
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle8"), hk);  if(!hk.isValid()){evValid=false;}	
  FilteredBitonicSortJetColl.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle9"), hk);       if(!hk.isValid()){evValid=false;}	
  FilteredBitonicSortJetColl.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle10"), hk); if(!hk.isValid()){evValid=false;}	
  FilteredBitonicSortJetColl.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredCircle12"), hk);      if(!hk.isValid()){evValid=false;} 
  FilteredBitonicSortJetColl.push_back( &(*hk) );
  
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredSquare8"), hk);  if(!hk.isValid()){evValid=false;}	
  FilteredBitonicSortJetColl.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredSquare9"), hk);       if(!hk.isValid()){evValid=false;}	
  FilteredBitonicSortJetColl.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredSquare10"), hk); if(!hk.isValid()){evValid=false;}	
  FilteredBitonicSortJetColl.push_back( &(*hk) );
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("FilteredSquare12"), hk);	if(!hk.isValid()){evValid=false;}	
  FilteredBitonicSortJetColl.push_back( &(*hk) );
  

  if(evValid){}
  else{cout<<"-------------------------> ERROR: Event content corrupted"<<endl;}
  
//-------------------------------------------------------------------------------------------//

  float maxDist;
  float goodPt(0), goodEta(0), goodPhi(0), goodArea(0);
  double no_tows[]={52,52,52,69,69,69,80,80,80,112,112,112,64,64,64,81,81,81,100,100,100,144,144,144};
  
  int Typecount ( 0 );
  TLorentzVector PF,TRPF;

  //-----PU bins-----//
  int m=0, PVindx=0;
  int PVmin=0, PVmax=5;
  if(PVSIZE < PVmax && PVSIZE >= PVmin){
    PVindx=m; // 5>PU>0
  }
  ++m;
  for( PVmin = 5 ; PVmin < 50 ; PVmin=PVmin+5 , ++m){    
    PVmax = PVmin+5;
    if(PVSIZE < PVmax && PVSIZE >= PVmin){     PVindx=m;        }
  }
  if(PVSIZE >= 50)   PVindx=m;

//--------------------------------------------------------------------------------------------------//
//  Median rho of L1 ak5 jets
//--------------------------------------------------------------------------------------------------//

  vector <double> JetRhosScaled;
  vector <double> JetRhosNoScaled;
  
  for(reco::PFJetCollection::const_iterator itr = L1ak5jets->begin(); itr != L1ak5jets->end() ; ++itr) {
    //only take jets with a eta of less than |3|
    if( abs( itr->eta() ) < 3 ){ 

      //Option 1:
      //Scale L1 jets to offline non PU sub jets DONE
      //Evaluate L1 rho from scaled L1 jets
      //scale rho L1 (scaled jets) to offline rho
      //Apply scaling to L1 rho
      //Correct L1 jets for L1 rho
      int j = eta_ref( itr->eta() );
      L1AK5CalibrationContants jetCalibration = getL1AK5CalibrationContants(Typecount,j);
      double scaled_pt_ak5 = jetCalibration.a
            +jetCalibration.b*(itr->pt())
            +jetCalibration.c*(itr->pt())*(itr->pt());
      JetRhosScaled.push_back( scaled_pt_ak5/itr->jetArea() );

      //Option 2:
      //Evaluate rho from unscaled L1 jets
      //scale L1 rho (unscaled jets) to offline rho
      //apply rho correction to L1 jets
      //scale jets to offline PU subtracted jets
      //apply scaling to L1 jets
      JetRhosNoScaled.push_back( itr->pt()/itr->jetArea() );
    }
  }

  //Evaluate L1 rho from scaled L1 jets:
  double median_rho_scaled = median(JetRhosScaled);
  Rho_L1AK5_Offline_scaled->Fill( (rhoCalo - median_rho_scaled)/rhoCalo );

  //Evaluate rho from unscaled L1 jets
  double median_rho_noscaled = median(JetRhosNoScaled);
  Rho_L1AK5_Offline_nonscaled->Fill( (rhoCalo - median_rho_noscaled)/rhoCalo );

  //Option 1: scale rho L1 (scaled jets) to offline rho
   Rho_Scaling_Option1->Fill( median_rho_scaled , rhoCalo );   
  //Option 2: Scale L1 rho(unscaled jets) to Offline rho
  Rho_Scaling_Option2->Fill( median_rho_noscaled , rhoCalo );

  //Apply rho correction to the medain_rho_scaled
  double median_rho_scaled_corr = rho_scale_op1_L1ak5( median_rho_scaled ); 
  //Apply rho correction to median_rho_noscaled
  double median_rho_noscaled_corr = rho_scale_op2_L1ak5( median_rho_noscaled );

  //Fill histograms of the rho resolution
  Rho_L1AK5_Offline_Option1 -> Fill( ( median_rho_scaled_corr - rhoCalo )/rhoCalo );
  if(rhoCalo>5) Rho_L1AK5_Offline_Option2->Fill( ( median_rho_noscaled_corr - rhoCalo )/rhoCalo ); 

  //Rescale op2 rhos: its not peaked at 0
  Rho_Scaling_Option2_Again -> Fill (median_rho_noscaled_corr , rhoCalo);
  double rescaled_op2Rho =   rho_op2_again (  median_rho_scaled_corr );
  Rho_L1AK5_Offline_Option2Again ->Fill( ( rescaled_op2Rho - rhoCalo)/rhoCalo  );

//-------------------------------------------------------------------------------------------------//
//  Scale of L1ak5 jets to non PU subtracted offline jets
//  Scale of L1ak5 jets to PU subtracted offline jets
//-------------------------------------------------------------------------------------------------//

  for(reco::CaloJetCollection::const_iterator it = Calojets->begin(); it != Calojets->end() ; ++it) {
    double scale = jetCorrector->correction(*it,iEvent,iSetup);
    double scaled_pt = it->pt()*scale;
    PF.SetPtEtaPhiM( scaled_pt,it->eta(),it->phi(),0);
    maxDist=999;
    goodPt=-100;
    goodEta=-100;
    goodPhi=-100;
    goodArea = 1000;
    //eta reference for scale
    int j = eta_ref( it->eta() );

    for(reco::PFJetCollection::const_iterator itr = L1ak5jets->begin(); itr != L1ak5jets->end() ; ++itr) {
      TRPF.SetPtEtaPhiM(itr-> pt(),itr->eta(),itr->phi(),0);
      float dr=PF.DeltaR(TRPF);
      if (dr>DR) continue;
      if (dr<maxDist){
        maxDist=dr;
        goodEta=itr->eta();
        goodPhi=itr->phi();
        goodPt=itr->pt();
        goodArea=itr->jetArea();
      }
    }

    if(maxDist<998){
      //cout<<"Offline scaled pt: "<<scaled_pt<<" eta: "<<it->eta()<<" phi:"<< it->phi() << endl;
      //cout<<"Macthed OPTION 2: pt: "<<goodPt<<"PU sub pt: "<<(goodPt-median_rho_noscaled_corr*goodArea);
      //cout <<" eta: "<<goodEta<<" phi:"<< goodPhi << endl;

      if(goodPt>0){
        //TProfiles for the scale: scale L1 jets to non PU subtracted offline jets OPTION 1
        L1ak5_Pt_scale_eta[j]->Fill( goodPt, it->pt() );
      }

      double corrected_pt = (goodPt - median_rho_noscaled_corr*goodArea) ;
      if(corrected_pt>0){
        //TProfiles for the scale: scale L1 jets to non PU subtracted offline jets OPTION 2 
        L1ak5_Pt_scale_eta_op2[j]->Fill( corrected_pt , scaled_pt );
      }

      //Going to scale op1 jets again because they peak at -0.5: weird.
      L1AK5CalibrationContants jetCalibration = getL1AK5CalibrationContants(Typecount,j);
      double scaled_pt_ak5 = jetCalibration.a
            +jetCalibration.b*(goodPt)
            +jetCalibration.c*(goodPt)*(goodPt);

      double l1ak5_Corrpt_op1 = (scaled_pt_ak5 - median_rho_scaled_corr*goodArea );
      if(l1ak5_Corrpt_op1>0){
        L1ak5_Pt_scale_op1_rescale[j]->Fill( l1ak5_Corrpt_op1 , scaled_pt );
      }

      //Trigger turn on curves
      double scaled_pt_ak5_op2 (0); 
      double l1ak5_Corrpt_op2 = (goodPt - median_rho_noscaled_corr*goodArea);
      L1AK5CalibrationContants2 jetCalibration2 = getL1AK5CalibrationContants2(Typecount,j,l1ak5_Corrpt_op2);

      scaled_pt_ak5_op2 = jetCalibration2.a
          +jetCalibration2.b*(l1ak5_Corrpt_op2)
          +jetCalibration2.c*pow(l1ak5_Corrpt_op2,2)
          +jetCalibration2.d*pow(l1ak5_Corrpt_op2,3)
          +jetCalibration2.e*pow(l1ak5_Corrpt_op2,4)
          +jetCalibration2.f*pow(l1ak5_Corrpt_op2,5)
          +jetCalibration2.g*pow(l1ak5_Corrpt_op2,6)
          +jetCalibration2.h*pow(l1ak5_Corrpt_op2,7)
          +jetCalibration2.i*pow(l1ak5_Corrpt_op2,8)
          +jetCalibration2.j*pow(l1ak5_Corrpt_op2,9); 

      if(scaled_pt_ak5_op2 > 40) L1AK5_num->Fill(scaled_pt);
      L1AK5_den->Fill(scaled_pt);

    }
  }

//-----------------------------------------------------------------------------------------//
//  L1 AK5 jets vs offline jets
//-----------------------------------------------------------------------------------------//

  for(reco::CaloJetCollection::const_iterator it = Calojets->begin(); it != Calojets->end() ; ++it) {
    double scale = jetCorrector->correction(*it,iEvent,iSetup);
    double scaled_pt = it->pt()*scale;
    if( scaled_pt > PTlim && abs(it->eta() ) < 3 ) {
      PF.SetPtEtaPhiM( scaled_pt,it->eta(),it->phi(),0);
      maxDist=999;
      goodPt=-100;
      goodEta=-100;
      goodPhi=-100;
      goodArea=100;

      int j = eta_ref( it->eta() );

      for(reco::PFJetCollection::const_iterator itr = L1ak5jets->begin(); itr != L1ak5jets->end() ; ++itr) {
        TRPF.SetPtEtaPhiM(itr-> pt(),itr->eta(),itr->phi(),0);
        float dr=PF.DeltaR(TRPF);
        if (dr>DR) continue;
        if (dr<maxDist){
          maxDist=dr;
          goodEta=itr->eta();
          goodPhi=itr->phi();
          goodPt=itr->pt();
          goodArea=itr->jetArea();
        }
      }

      if(maxDist<998 && goodPt>0){
       
     //   std::cout <<"Calo Jet Pt = " << it->pt() << ", eta = "<<it->eta()<<", phi = "<< it->phi();
     //   std::cout <<", Correction factor = "<< scale <<", ptCorr = "<<scale * it->pt() <<std::endl;
     //   std::cout <<"Matched L1ak5 raw Jet Pt = " << goodPt << ", eta = "<< goodEta;
     //   std::cout <<", phi = "<<goodPhi <<", jet area: "<<goodArea<<std::endl;

        //If we have found the closest jet, then fill eta, phi and pt resolutions
        Offline_L1ak5_eta   -> Fill( it->eta() - goodEta );
        Offline_L1ak5_phi   -> Fill( it->phi() - goodPhi );
        Offline_L1ak5_pt    -> Fill(  (it->pt() - goodPt ) / it->pt() );
        Offline_L1ak5_pt_PU -> Fill( (scaled_pt - goodPt) / scaled_pt );

        L1AK5CalibrationContants jetCalibration = getL1AK5CalibrationContants(Typecount,j);
        double scaled_pt_ak5 = jetCalibration.a
              +jetCalibration.b*(goodPt)
              +jetCalibration.c*(goodPt)*(goodPt);

        Offline_L1ak5_pt_scale->Fill( ( scaled_pt_ak5 - it->pt())/ it->pt() );
        Offline_L1ak5_pt_PU_scale->Fill( ( scaled_pt_ak5 - scaled_pt)/ scaled_pt );
        Offline_L1ak5_pt_PU_scale_median_scaled->Fill( ((scaled_pt_ak5 - median_rho_scaled*goodArea) - scaled_pt)/scaled_pt);

        double l1ak5_Corrpt_op1 = (scaled_pt_ak5 - median_rho_scaled_corr*goodArea );
        //don't want the bin at zero:
        if(l1ak5_Corrpt_op1>0){
          Offline_L1ak5_Op1_pt_res->Fill( (l1ak5_Corrpt_op1 - scaled_pt)/scaled_pt); 
        }

        //rescale l1 pu sub jets to get peak at 0
        L1AK5CalibrationContantsRescale jetCalibrationRescale = getL1AK5CalibrationContantsRescale(Typecount,j,l1ak5_Corrpt_op1);

        double scaled_pt_ak5_op1Rescale = jetCalibrationRescale.a
            +jetCalibrationRescale.b*(l1ak5_Corrpt_op1)
            +jetCalibrationRescale.c*(l1ak5_Corrpt_op1)*(l1ak5_Corrpt_op1)
            +jetCalibrationRescale.d*pow(l1ak5_Corrpt_op1,3)
            +jetCalibrationRescale.e*pow(l1ak5_Corrpt_op1,4)
            +jetCalibrationRescale.f*pow(l1ak5_Corrpt_op1,5)
            +jetCalibrationRescale.g*pow(l1ak5_Corrpt_op1,6)
            +jetCalibrationRescale.h*pow(l1ak5_Corrpt_op1,7)
            +jetCalibrationRescale.i*pow(l1ak5_Corrpt_op1,8)
            +jetCalibrationRescale.j*pow(l1ak5_Corrpt_op1,9);     

        if(scaled_pt_ak5_op1Rescale>0){
          Offline_L1ak5_Op1Rescale_pt_res->Fill( (scaled_pt_ak5_op1Rescale - scaled_pt)/scaled_pt); 
        }

        double l1ak5_Corrpt_op2 = (goodPt - median_rho_noscaled_corr*goodArea);
        double scaled_pt_ak5_op2 (0); 
        //scale the op2 jets
        L1AK5CalibrationContants2 jetCalibration2 = getL1AK5CalibrationContants2(Typecount,j,l1ak5_Corrpt_op2);
        scaled_pt_ak5_op2 = jetCalibration2.a
            +jetCalibration2.b*(l1ak5_Corrpt_op2)
            +jetCalibration2.c*pow(l1ak5_Corrpt_op2,2)
            +jetCalibration2.d*pow(l1ak5_Corrpt_op2,3)
            +jetCalibration2.e*pow(l1ak5_Corrpt_op2,4)
            +jetCalibration2.f*pow(l1ak5_Corrpt_op2,5)
            +jetCalibration2.g*pow(l1ak5_Corrpt_op2,6)
            +jetCalibration2.h*pow(l1ak5_Corrpt_op2,7)
            +jetCalibration2.i*pow(l1ak5_Corrpt_op2,8)
            +jetCalibration2.j*pow(l1ak5_Corrpt_op2,9);          

        if(scaled_pt_ak5_op2>1 && abs(goodEta)<3){
          Offline_L1ak5_Op2_pt_res->Fill( (scaled_pt_ak5_op2 - scaled_pt)/ scaled_pt );
        }
        
        //pt res in bins of pt
        Off_L1ak5_binnedPtres_op1->Fill( l1ak5_Corrpt_op1 , ((l1ak5_Corrpt_op1 - scaled_pt)/ scaled_pt ) );
        Off_L1ak5_binnedPtres_op1Rescale->Fill( scaled_pt_ak5_op1Rescale, ((scaled_pt_ak5_op1Rescale - scaled_pt)/scaled_pt ));
        Off_L1ak5_binnedPtres_op2->Fill( scaled_pt_ak5_op2 , ((scaled_pt_ak5_op2 - scaled_pt) / scaled_pt) );

        //pt resolution vs number of vertices
        //option 1
        Pt_res_vertices_op1 -> Fill( PVSIZE  , ( (scaled_pt_ak5_op1Rescale - scaled_pt)/ scaled_pt)  );
        Pt_res_vertices_op2 -> Fill( PVSIZE  , ( (scaled_pt_ak5_op2 - scaled_pt) / scaled_pt) ); 

      }
    } //jet > 30GeV and eta < +- 3
  }

//-----------------------------------------------------------------------------------------//
//  L1 extra jets vs offline jets
//-----------------------------------------------------------------------------------------//
  for(reco::CaloJetCollection::const_iterator it = Calojets->begin(); it != Calojets->end() ; ++it) {
    double scale = jetCorrector->correction(*it,iEvent,iSetup);
    double scaled_pt = it->pt()*scale;  
    if( scaled_pt > PTlim && abs(it->eta()) < 3){
      PF.SetPtEtaPhiM( scaled_pt,it->eta(),it->phi(),0);
      maxDist=999;
      goodPt=-100;
      goodEta=-100;
      goodPhi=-100;

      //Get Central, Forward, tau jet collections
      for (uint i=0; i< l1extraparticles.size();++i){
        Handle<l1extra::L1JetParticleCollection> l1Jets;
        iEvent.getByLabel(l1extraparticles[i],l1Jets);

        for(l1extra::L1JetParticleCollection::const_iterator itr = l1Jets->begin();  itr != l1Jets->end(); ++itr ) {
          TRPF.SetPtEtaPhiM(itr-> pt(),itr->eta(),itr->phi(),0);  
          float dr=PF.DeltaR(TRPF);                
          if (dr>DR) continue;
          if (dr<maxDist){      
            maxDist=dr;
            goodEta=itr->eta();
            goodPhi=itr->phi();
            goodPt=itr->pt();
          }
        }
      }

      if(maxDist<998){

        //cout<<"Calo Jet Pt = " << it->pt() << ", eta = "<<it->eta()<<", phi = "<< it->phi();
        //cout<<", Correction factor = "<< scale <<", ptCorr = "<<scale * it->pt() <<std::endl;
        //cout<<"Matched L1 extra Jet Pt = " <<goodPt<< ", eta = "<<goodEta<<", phi = "<<goodPhi<<endl;

        //If we have found the closest jet, then fill eta, phi and pt resolutions
        Offline_L1_eta   -> Fill( it->eta() - goodEta );
        Offline_L1_phi   -> Fill( it->phi() - goodPhi );
        Offline_L1_pt    -> Fill(  (it->pt() - goodPt ) / it->pt() );  
        Offline_L1_pt_PU -> Fill( (scaled_pt - goodPt) / scaled_pt );

        //Pt resolution as a function of PU
        Pt_res_vertices_op1_CurrentL1->Fill( PVSIZE, (( goodPt - scaled_pt ) / scaled_pt ) );
        //Pt resolution as a function of Pt
        Offline_L1_binnedPtres->Fill( goodPt , (( goodPt - scaled_pt ) / scaled_pt ) );
      }
    } //jet > 30GeV  

    //no 30GeV limit for turn on curves: rematch
    PF.SetPtEtaPhiM( scaled_pt,it->eta(),it->phi(),0);
    maxDist=999;
    goodPt=-100;
    goodEta=-100;
    goodPhi=-100;

    //Get Central, Forward, tau jet collections
    for (uint i=0; i< l1extraparticles.size();++i){
      Handle<l1extra::L1JetParticleCollection> l1Jets;
      iEvent.getByLabel(l1extraparticles[i],l1Jets);

      for(l1extra::L1JetParticleCollection::const_iterator itr = l1Jets->begin();  itr != l1Jets->end(); ++itr ) {
        TRPF.SetPtEtaPhiM(itr-> pt(),itr->eta(),itr->phi(),0);  
        float dr=PF.DeltaR(TRPF);                
        if (dr>DR) continue;
        if (dr<maxDist){      
          maxDist=dr;
          goodEta=itr->eta();
          goodPhi=itr->phi();
          goodPt=itr->pt();
        }
      }
    }

    if(maxDist<998){
      if(goodPt>40) { CurrL1_num -> Fill( scaled_pt ) ; }//numerator
      CurrL1_den -> Fill(scaled_pt) ;//denominator
    }
  }

//------------------------------------------------------------------------------------------------//
// Tower level jets: bitonic sort
//-----------------------------------------------------------------------------------------------------//


  Typecount=0;    //relates to the jet collection

  //Jet collection contains bitonic sorted & filtered 8x8 circ, 9x9 circ, etc; no mean or median yet
  for(std::vector< const  EtaPhiContainer<l1slhc::L1TowerJet>*>::const_iterator PreFilteredJetType= FilteredBitonicSortJetColl.begin() ;
      PreFilteredJetType!= FilteredBitonicSortJetColl.end() ; 
      ++PreFilteredJetType,       Typecount+=3 
     )
  {

    //------------- Rho calculations: loop through tower jets and fill vectors ----------------//
    vector<double> TowerEnergiesScaled;
    vector<double> TowerEnergiesUnScaled;
    for( EtaPhiContainer<l1slhc::L1TowerJet>::const_iterator PreFilteredTowJet = (**PreFilteredJetType).begin() ;
        PreFilteredTowJet != (**PreFilteredJetType).end() ; 
        ++ PreFilteredTowJet )
    {
      if( abs(PreFilteredTowJet->p4().eta() )<3){
        //Option 1:
        //Scale L1 jets to offline non PU sub jets
        int j = eta_ref( PreFilteredTowJet->p4().eta() );
        L1TowerCalibrationContants jetCalibration = getL1TowerCalibrationContants(Typecount,j);
        double scaled_pt_tow = jetCalibration.a
              +jetCalibration.b*(PreFilteredTowJet->p4().pt())
              +jetCalibration.c*(PreFilteredTowJet->p4().pt())*(PreFilteredTowJet->p4().pt());

        //Option 2:
        //Evaluate rhos from scaled and non scaled L1 jet energies
        TowerEnergiesScaled.push_back(scaled_pt_tow);
        TowerEnergiesUnScaled.push_back(PreFilteredTowJet->p4().pt());

      }
    }
    double towers = no_tows[Typecount];
    double areaPerJet = towers* (0.087 * 0.087) ;
    //option 1
    double median_rho_scaled = ( median( TowerEnergiesScaled ) / areaPerJet );
    //option 2
    double median_rho_noscaled = ( median( TowerEnergiesUnScaled ) / areaPerJet );
 
    //Fill TProfiles for rho scaling
    Rho_Scaling_Option1_BitFilt[Typecount]->Fill( median_rho_scaled , rhoCalo );
    Rho_Scaling_Option2_BitFilt[Typecount]->Fill( median_rho_noscaled , rhoCalo );

    //Apply rho correction to the medain_rho_scaled
    double median_BF_rho_scaled_corr = rho_scale_op1_L1towBF( median_rho_scaled , Typecount ); 
    //Apply rho correction to median_rho_noscaled
    double median_BF_rho_noscaled_corr = rho_scale_op2_L1towBF( median_rho_noscaled , Typecount);

    L1Tower_rho[Typecount]->Fill( median_rho_scaled );
    L1Tower_rho_res_BitFilt_scaled[Typecount]->Fill( (rhoCalo - median_BF_rho_scaled_corr) / rhoCalo );
    L1Tower_rho_res_BitFilt_noscaled[Typecount]->Fill( (rhoCalo - median_BF_rho_noscaled_corr) / rhoCalo );

    Rho_ReScaling_Option1_BitFilt[Typecount]->Fill( median_BF_rho_scaled_corr , rhoCalo );
    Rho_ReScaling_Option2_BitFilt[Typecount]->Fill( median_BF_rho_noscaled_corr , rhoCalo );

    double median_BF_rho_scaled_Recorr = rho_Rescale_op1_L1towBF( median_BF_rho_scaled_corr , Typecount ); 
    double median_BF_rho_noscaled_Recorr = rho_Rescale_op2_L1towBF( median_BF_rho_noscaled_corr , Typecount);

    Rho_ReReScaling_Option2_BitFilt[Typecount]->Fill( median_BF_rho_noscaled_Recorr , rhoCalo );

    if(median_BF_rho_scaled_corr>0.1){
      L1Tow_BitFilt_rho_res_Option1[Typecount]->Fill((median_BF_rho_noscaled_Recorr-rhoCalo)/ rhoCalo);
    }
    if(median_BF_rho_noscaled_corr>0.1){
      L1Tow_BitFilt_rho_res_Option2[Typecount]->Fill((median_BF_rho_noscaled_Recorr-rhoCalo)/ rhoCalo);
    }
    
    //--------------------------------------------------------------------------------------------//
    //  Tower level jets: Pt scale for options 1 and 2
    //--------------------------------------------------------------------------------------------------//

    for(reco::CaloJetCollection::const_iterator it = Calojets->begin(); it != Calojets->end() ; ++it) {
  
      double scale = jetCorrector->correction(*it,iEvent,iSetup);
      double scaled_pt = it->pt()*scale; 
      //eta ref
      int j = eta_ref( it->eta() );
      PF.SetPtEtaPhiM( scaled_pt,it->eta(),it->phi(),0);
      float maxDist=999;
      float goodPt=-100;
      float goodEta=-100;
      float goodPhi=-100;
      for( EtaPhiContainer<l1slhc::L1TowerJet>::const_iterator PreFilteredTowJet = (**PreFilteredJetType).begin() ;
           PreFilteredTowJet != (**PreFilteredJetType).end() ; ++ PreFilteredTowJet 
         ) 
      {
        TRPF.SetPtEtaPhiM(PreFilteredTowJet->p4().pt(),PreFilteredTowJet->p4().eta(),PreFilteredTowJet->p4().phi(),0);
        float dr=PF.DeltaR(TRPF);
        if (dr>DR) continue;
        if (dr<maxDist){//Performs an iteration. maxDist only gets a new value 
          maxDist=dr;
          goodEta=PreFilteredTowJet->p4().eta(); 
          goodPhi=PreFilteredTowJet->p4().phi();
          goodPt=PreFilteredTowJet->p4().pt();
        }
      }
      //matched jets, fill TProfiles for the scaling
      if (maxDist<998){
          
        L1towBF_Pt_scale_eta_op1[Typecount*etast.size() + j]->Fill( goodPt , it->pt());
        //op2:
        L1towBF_Pt_scale_eta_op2[Typecount*etast.size() + j]->Fill( (goodPt -median_BF_rho_noscaled_Recorr) , scaled_pt);

        //Rescale op1 jets perhaps to get peak at 0
        L1TowerCalibrationContants jetCalibration = getL1TowerCalibrationContants(Typecount,j);
        double scaled_pt_tow = jetCalibration.a
              +jetCalibration.b*(goodPt)
              +jetCalibration.c*(goodPt)*(goodPt);
        double l1tow_Corrpt_op1 = (scaled_pt_tow - median_BF_rho_scaled_Recorr*areaPerJet );

        L1towBF_Pt_scale_eta_op1Rescale[Typecount*etast.size() + j]->Fill(l1tow_Corrpt_op1, scaled_pt);

        double l1tow_Corrpt_op2 = (goodPt - median_BF_rho_noscaled_Recorr*areaPerJet);
        L1TowerCalibrationContants2 jetCalibration2 = getL1TowerCalibrationContants2(Typecount,j);
        double scaled_pt_tow_op2 = jetCalibration2.a
              +jetCalibration2.b*(l1tow_Corrpt_op2)
              +jetCalibration2.c*pow(l1tow_Corrpt_op2, 2)
              +jetCalibration2.d*pow(l1tow_Corrpt_op2 , 3)
              +jetCalibration2.e*pow(l1tow_Corrpt_op2 , 4);

        //Trigger turn on curves: not within 30GeV offline jet cut
        if(scaled_pt_tow_op2>40) L1Tower_num[Typecount]->Fill(scaled_pt);    //numerator: L1 jet>40GeV
        L1Tower_den[Typecount]->Fill(scaled_pt);                          //demonimator: all matched jets

      }

    //----------Match the L1 tower jet to offline jet-------------------------------//

      if(scaled_pt>PTlim && abs(it->eta()) < 3 ){
        PF.SetPtEtaPhiM( scaled_pt,it->eta(),it->phi(),0);
        float maxDist=999;
        float goodPt=-100;
        float goodEta=-100;
        float goodPhi=-100;
        for( EtaPhiContainer<l1slhc::L1TowerJet>::const_iterator PreFilteredTowJet = (**PreFilteredJetType).begin() ;
            PreFilteredTowJet != (**PreFilteredJetType).end() ; ++ PreFilteredTowJet 
          ) 
        {
          TRPF.SetPtEtaPhiM(PreFilteredTowJet->p4().pt(),PreFilteredTowJet->p4().eta(),PreFilteredTowJet->p4().phi(),0);
          float dr=PF.DeltaR(TRPF);
          if (dr>DR) continue;
          if (dr<maxDist){ //Performs an iteration to match tower jet to closest offline jet 
            maxDist=dr;
            goodEta=PreFilteredTowJet->p4().eta(); 
            goodPhi=PreFilteredTowJet->p4().phi();
            goodPt=PreFilteredTowJet->p4().pt();
          }
        }

        if (maxDist<998){ //if jet is matched
          L1tow_BitFilt_eta[Typecount]->Fill(it->eta() - goodEta);
          L1tow_BitFilt_phi[Typecount]->Fill(it->phi() - goodPhi);
          L1tow_BitFilt_pt[Typecount]->Fill( (it->pt() - goodPt)/ it->pt() );
          L1tow_BitFilt_pt_PU[Typecount]->Fill( ( scaled_pt-(goodPt - median_rho_scaled*areaPerJet) )/ scaled_pt );

          L1TowerCalibrationContants jetCalibration = getL1TowerCalibrationContants(Typecount,j);
          double scaled_pt_tow = jetCalibration.a
               +jetCalibration.b*(goodPt)
               +jetCalibration.c*(goodPt)*(goodPt);
          double scaled_pt_tow_PU = (scaled_pt_tow - median_rho_scaled*areaPerJet);

          L1tow_BitFilt_pt_scale[Typecount]->Fill( (scaled_pt_tow - it->pt())/ it->pt() );
          L1tow_BitFilt_pt_PU_scale[Typecount]->Fill( (scaled_pt_tow_PU - scaled_pt)/ scaled_pt );
          //option 1
          double l1tow_Corrpt_op1 = (scaled_pt_tow - median_BF_rho_scaled_Recorr*areaPerJet );
          if( l1tow_Corrpt_op1>0 ){
            Offline_L1tow_Op1_pt_res[Typecount]->Fill( (l1tow_Corrpt_op1 - scaled_pt)/scaled_pt); 
          }

          //option 2
          double l1tow_Corrpt_op2 = (goodPt - median_BF_rho_noscaled_Recorr*areaPerJet); 

          L1TowerCalibrationContants2 jetCalibration2 = getL1TowerCalibrationContants2(Typecount,j);
          double scaled_pt_tow_op2 = jetCalibration2.a
                  +jetCalibration2.b*(l1tow_Corrpt_op2)
                  +jetCalibration2.c*(l1tow_Corrpt_op2)*(l1tow_Corrpt_op2)
                  +jetCalibration2.d*pow(l1tow_Corrpt_op2 , 3)
                  +jetCalibration2.e*pow(l1tow_Corrpt_op2 , 4);

          //put this into a pt res histogram
          if(scaled_pt_tow_op2>0){ //get rid of -1 bin
            Offline_L1tow_Op2_pt_res[Typecount]->Fill( (scaled_pt_tow_op2 - scaled_pt)/ scaled_pt ); 
          }

          //option 2:in bins of Pt
          Off_BFL1tow_binnedPtres_op2[Typecount]->Fill( scaled_pt_tow_op2 , ((scaled_pt_tow_op2 - scaled_pt)/ scaled_pt)  );
          //in bins of PU
          Off_BFL1tow_Pt_PU_op2_binpt[Typecount]->Fill(  PVSIZE, ((scaled_pt_tow_op2 - scaled_pt)/ scaled_pt)  ); 

        }  
      }//>30GeV
    }//calo jets
  }
}

//=================================================================================================
// Functions
//=================================================================================================

//option 1 tower calibration constants
L1TowerCalibrationContants L1JetAnalyzer::getL1TowerCalibrationContants(int l1JetColId, int l1JetEtaRange){

  L1TowerCalibrationContants out;
  out.a = 0; out.b = 0; out.c = 0; 

   if(l1JetColId==0 || l1JetColId==1 || l1JetColId==2 ){   //8x8 out.cirout.cle
     if( l1JetEtaRange==1){  out.a = 13.2258;  out.b= 2.3969;   out.c= -0.0151855;   }
     if( l1JetEtaRange==2){  out.a = 7.25494;  out.b= 1.44515;  out.c= 0.000267118;  }
     if( l1JetEtaRange==3){  out.a = 4.71518;  out.b= 1.30174;  out.c= 0.000131428;  }
     if( l1JetEtaRange==4){  out.a = 4.14074;  out.b= 1.31278;  out.c= -0.001055;    }
     if( l1JetEtaRange==5){  out.a = 4.23914;  out.b= 1.30269;  out.c= -0.00118351;  }
     if( l1JetEtaRange==6){  out.a = 4.27434;  out.b= 1.315;    out.c= -0.00109691;  }
     if( l1JetEtaRange==7){  out.a = 4.28461;  out.b= 1.31528;  out.c= -0.000867265; }
     if( l1JetEtaRange==8){  out.a = 4.17597;  out.b= 1.29612;  out.c= -0.000821865; }
     if( l1JetEtaRange==9){  out.a = 3.99956;  out.b= 1.30952;  out.c= -0.00150007;  }
     if( l1JetEtaRange==10){ out.a = 4.54639;  out.b= 1.31673;  out.c= 0.000278103;  }
     if( l1JetEtaRange==11){ out.a = 6.81068;  out.b= 1.40814;  out.c= 3.99861e-05;  }
     if( l1JetEtaRange==12){ out.a = 11.6192;  out.b= 2.27774;  out.c= -0.0103032;   }
   }if(l1JetColId==3 || l1JetColId==4 || l1JetColId==5 ){   //9x9 out.cirout.cle
     if( l1JetEtaRange==1){ out.a = 13.7046;   out.b= 1.72282;    out.c= -0.00631639;      }
     if( l1JetEtaRange==2){ out.a = 6.77388;   out.b= 1.42266;    out.c= 0.000810111;      }
     if( l1JetEtaRange==3){ out.a = 4.80397;   out.b= 1.168;      out.c= -0.00104256;      }
     if( l1JetEtaRange==4){ out.a = 4.28681;   out.b= 1.14205;    out.c= 0.000524676;      }
     if( l1JetEtaRange==5){ out.a = 4.46379;   out.b= 1.11792;    out.c= -0.00411943;      }
     if( l1JetEtaRange==6){ out.a = 4.44028;   out.b= 1.16678;    out.c= -0.000334375;     }
     if( l1JetEtaRange==7){ out.a = 4.47806;   out.b= 1.15379;    out.c= 0.000174234;      }
     if( l1JetEtaRange==8){ out.a = 4.30516;   out.b= 1.15518;    out.c= 0.000142683;      }
     if( l1JetEtaRange==9){ out.a = 3.9672;    out.b= 1.26836;    out.c= -0.00724705;      }
     if( l1JetEtaRange==10){out.a = 4.70676;   out.b= 1.15819;    out.c= 0.00216844;       }
     if( l1JetEtaRange==11){out.a = 7.08029;   out.b= 1.16898;    out.c= -0.0088703;       }
     if( l1JetEtaRange==12){out.a = 10.8914;   out.b= 2.0546;     out.c= -0.00264239;      }
   }if(l1JetColId==6 || l1JetColId==7 || l1JetColId==8){   //10x10 out.cirout.cle
     if( l1JetEtaRange==2){  out.a = 6.96831;  out.b= 1.24766;  out.c= 0.00238192;      }
     if( l1JetEtaRange==3){  out.a = 4.87064;  out.b= 1.10341;  out.c= 0.000974906;     }
     if( l1JetEtaRange==4){  out.a = 4.42627;  out.b= 1.04538;  out.c= 0.00157314;      }
     if( l1JetEtaRange==5){  out.a = 4.44696;  out.b= 1.06347;  out.c= 0.000232833;     }
     if( l1JetEtaRange==6){  out.a = 4.50476;  out.b= 1.05951;  out.c= 0.000229663;     }
     if( l1JetEtaRange==7){  out.a = 4.50495;  out.b= 1.0644;   out.c= 0.000572767;     }
     if( l1JetEtaRange==8){  out.a = 4.29545;  out.b= 1.11928;  out.c= -0.00396351;     }
     if( l1JetEtaRange==9){  out.a = 4.24735;  out.b= 1.06513;  out.c= 0.000249408;     }
     if( l1JetEtaRange==10){ out.a = 4.61083;  out.b= 1.18202;  out.c= -0.00518618;     }
     if( l1JetEtaRange==11){ out.a = 6.12731;  out.b= 1.55299;  out.c= -0.0116564;      }
     if( l1JetEtaRange==12){ out.a = 9.92609;  out.b= 2.40657;  out.c= -0.0133076;      }
   }if(l1JetColId==9 || l1JetColId==10 || l1JetColId==11){   //12x12 out.cirout.cle
     if( l1JetEtaRange==2){ out.a = 6.13973;  out.b= 1.30568;    out.c= -0.00769519;      }
     if( l1JetEtaRange==3){ out.a = 4.99152;  out.b= 0.931042;   out.c= 0.00318799;       }
     if( l1JetEtaRange==4){ out.a = 4.57265;  out.b= 0.876347;   out.c= 0.000459323;      }
     if( l1JetEtaRange==5){ out.a = 4.67392;  out.b= 0.844406;   out.c= 0.00172336;       }
     if( l1JetEtaRange==6){ out.a = 4.51301;  out.b= 0.944927;   out.c= -0.00507978;      }
     if( l1JetEtaRange==7){ out.a = 4.53566;  out.b= 0.94017;    out.c= -0.00464729;      }
     if( l1JetEtaRange==8){ out.a = 4.59243;  out.b= 0.862914;   out.c= -6.89403e-05;     }
     if( l1JetEtaRange==9){ out.a = 4.37393;  out.b= 0.916107;   out.c= -0.00497283;      }
     if( l1JetEtaRange==10){out.a = 4.75126;  out.b= 1.02541;    out.c= -0.00601781;      }
     if( l1JetEtaRange==11){out.a = 6.45743;  out.b= 1.1434;     out.c= -0.00814567;      }
   }if(l1JetColId==12 || l1JetColId==13 || l1JetColId==14){  //8x8 squout.are
     if( l1JetEtaRange==1){  out.a = 14.0907;   out.b= 2.20944;  out.c= -0.00887922;    }
     if( l1JetEtaRange==2){  out.a = 7.13183;   out.b= 1.25564;  out.c= 0.00228916;     }
     if( l1JetEtaRange==3){  out.a = 4.78378;   out.b= 1.19668;  out.c= 0.00171253;     }
     if( l1JetEtaRange==4){  out.a = 4.52582;   out.b= 1.10433;  out.c= 0.00115559;     }
     if( l1JetEtaRange==5){  out.a = 4.61356;   out.b= 1.10011;  out.c= 0.000433501;    }
     if( l1JetEtaRange==6){  out.a = 4.64989;   out.b= 1.11497;  out.c= 0.00016716;     }
     if( l1JetEtaRange==7){  out.a = 4.6114;    out.b= 1.12664;  out.c= -3.59491e-05;   }
     if( l1JetEtaRange==8){  out.a = 4.47007;   out.b= 1.11238;  out.c= -3.31398e-05;   }
     if( l1JetEtaRange==9){  out.a = 4.18184;   out.b= 1.17559;  out.c= -0.00223842;    }
     if( l1JetEtaRange==10){ out.a = 4.49961;   out.b= 1.33933;  out.c= -0.0116934;     }
     if( l1JetEtaRange==11){ out.a = 4.8388;    out.b= 2.35682;  out.c= -0.0137188;     }
     if( l1JetEtaRange==12){ out.a = 11.7885;   out.b= 2.12557;  out.c= -0.0070071;     }
   }if(l1JetColId==15 || l1JetColId==16 || l1JetColId==17){  //9x9 squout.are
     if( l1JetEtaRange==1){  out.a = 12.5656;   out.b= 2.06762;   out.c= -0.00739286;    }
     if( l1JetEtaRange==2){  out.a = 24.2954;   out.b= -5.72677;  out.c= 0.0639203;      }
     if( l1JetEtaRange==3){  out.a = 4.8236;    out.b= 1.11877;   out.c= 0.00189857;     }
     if( l1JetEtaRange==4){  out.a = 4.61783;   out.b= 0.986284;  out.c= 0.00184376;     }
     if( l1JetEtaRange==5){  out.a = 4.73558;   out.b= 0.973591;  out.c= -0.000656711;    }
     if( l1JetEtaRange==6){  out.a = 4.74088;   out.b= 1.0175;    out.c= 0.000445015;    }
     if( l1JetEtaRange==7){  out.a = 4.7089;    out.b= 1.01932;   out.c= 0.00102063;     }
     if( l1JetEtaRange==8){  out.a = 4.47217;   out.b= 1.04528;   out.c= -0.00201983;    }
     if( l1JetEtaRange==9){  out.a = 4.40425;   out.b= 0.993732;  out.c= 0.00143311;     }
     if( l1JetEtaRange==10){ out.a = 4.86523;   out.b= 1.05992;   out.c= 0.00245171;     }
     if( l1JetEtaRange==11){ out.a = 4.75617;   out.b= 2.24675;   out.c= -0.0123394;     }
     if( l1JetEtaRange==12){ out.a = 10.4637;   out.b= 2.16339;   out.c= -0.00556654;    }
   }if(l1JetColId==18 || l1JetColId==19 || l1JetColId==20){  //10x10 squout.are
     if( l1JetEtaRange==2){   out.a = 10.0298;  out.b= -0.157889; out.c= 0.00145047;    }
     if( l1JetEtaRange==3){   out.a = 5.50012;  out.b= 0.819637;  out.c= 0.000935719;   }
     if( l1JetEtaRange==4){   out.a = 4.70199;  out.b= 0.921028;  out.c= 0.00126585;    }
     if( l1JetEtaRange==5){   out.a = 4.7418;   out.b= 0.922145;  out.c= -0.000561485;  } 
     if( l1JetEtaRange==6){   out.a = 4.74681;  out.b= 0.9218;    out.c= 0.0010017;     }
     if( l1JetEtaRange==7){   out.a = 4.78948;  out.b= 0.892598;  out.c= 0.00159056;    }
     if( l1JetEtaRange==8){   out.a = 4.45831;  out.b= 0.985886;  out.c= -0.00492871;   }
     if( l1JetEtaRange==9){   out.a = 4.38918;  out.b= 0.922127;  out.c= -3.2851e-05;   }
     if( l1JetEtaRange==10){  out.a = 4.84705;  out.b= 1.03928;   out.c= -0.00763978;   }
     if( l1JetEtaRange==11){  out.a = 6.24422;  out.b= 1.40707;   out.c= -0.0109499;    }
     if( l1JetEtaRange==12){  out.a = 9.60314;  out.b= 2.30307;   out.c= -0.0110794;    }
   }if(l1JetColId==21 || l1JetColId==22 || l1JetColId==23){  //12x12 squout.are
     if( l1JetEtaRange==2){  out.a = 6.9798;   out.b= 0.766655;  out.c= -0.00739624;     }
     if( l1JetEtaRange==3){  out.a = 5.06892;  out.b= 0.865894;  out.c= 0.00315562;      }
     if( l1JetEtaRange==4){  out.a = 4.85475;  out.b= 0.720402;  out.c= 0.00375579;      }
     if( l1JetEtaRange==5){  out.a = 4.86363;  out.b= 0.73327;   out.c= 0.00124869;      }
     if( l1JetEtaRange==6){  out.a = 4.92773;  out.b= 0.718131;  out.c= 0.00213193;      }
     if( l1JetEtaRange==7){  out.a = 4.51582;  out.b= 0.86787;   out.c= -0.0046014;      }
     if( l1JetEtaRange==8){  out.a = 4.50126;  out.b= 0.817494;  out.c= -0.00228342;     }
     if( l1JetEtaRange==9){  out.a = 4.29426;  out.b= 0.833016;  out.c= -0.00423125;     }
     if( l1JetEtaRange==10){ out.a = 4.75989;  out.b= 0.967867;  out.c= -0.00651742;     }
     if( l1JetEtaRange==11){ out.a = 6.17994;  out.b= 1.15279;   out.c= -0.00957535;     }
   }
  return out;

}

//option 2 tower calibration constants
L1TowerCalibrationContants2 L1JetAnalyzer::getL1TowerCalibrationContants2(int l1JetColId, int l1JetEtaRange){
               
  // Initialization of output
  L1TowerCalibrationContants2 out;
  out.a = 0; out.b = 0; out.c = 0; out.d = 0; out.e = 0;

   //Fits of the profiles: from the local file, run without 30GeV cutoff
  if(l1JetColId==0 || l1JetColId==1 || l1JetColId==2 ){   //8x8 out.cirout.cle
    if( l1JetEtaRange==1){ out.a = 23.7558;  out.b= 1.93145;  out.c= -0.00512296;   }
    if( l1JetEtaRange==2){ out.a = 19.2431;  out.b= 1.51314;  out.c= -0.000105093;  }
    if( l1JetEtaRange==3){ out.a = 22.2783;  out.b= 1.72528;  out.c= -0.00145616;   }
    if( l1JetEtaRange==4){ out.a = 24.2513;  out.b= 1.77265;  out.c= -0.00247916;   }
    if( l1JetEtaRange==5){ out.a = 20.9402;  out.b= 1.66363;  out.c= -0.00178097;   }
    if( l1JetEtaRange==6){ out.a = 18.9574;  out.b= 1.66627;  out.c= -0.00211727;   }
    if( l1JetEtaRange==7){ out.a = 19.2048;  out.b= 1.6858;   out.c= -0.00228489;   }
    if( l1JetEtaRange==8){ out.a = 20.7764;  out.b= 1.67469;  out.c= -0.00199804;   }
    if( l1JetEtaRange==9){ out.a = 21.6223;  out.b= 2.12562;  out.c= -0.00427456;   }
    if( l1JetEtaRange==10){out.a = 18.378;   out.b= 1.54955;  out.c= -0.000972697;  }
    if( l1JetEtaRange==11){out.a = 18.9834;  out.b= 2.01681;  out.c= -0.00684155;   }
    if( l1JetEtaRange==12){ out.a = 11.6192;  out.b= 2.27774;  out.c= -0.0103032;   }
  }if(l1JetColId==3 || l1JetColId==4 || l1JetColId==5 ){   //9x9 out.cirout.cle
    if( l1JetEtaRange==1){ out.a = 13.7046;   out.b= 1.72282;    out.c= -0.00631639;      }
    if( l1JetEtaRange==2){ out.a = 6.77388;   out.b= 1.42266;    out.c= 0.000810111;      }
    if( l1JetEtaRange==3){ out.a = 4.80397;   out.b= 1.168;      out.c= -0.00104256;      }
    if( l1JetEtaRange==4){ out.a = 4.28681;   out.b= 1.14205;    out.c= 0.000524676;      }
    if( l1JetEtaRange==5){ out.a = 4.46379;   out.b= 1.11792;    out.c= -0.00411943;      }
    if( l1JetEtaRange==6){ out.a = 4.44028;   out.b= 1.16678;    out.c= -0.000334375;     }
    if( l1JetEtaRange==7){ out.a = 4.47806;   out.b= 1.15379;    out.c= 0.000174234;      }
    if( l1JetEtaRange==8){ out.a = 4.30516;   out.b= 1.15518;    out.c= 0.000142683;      }
    if( l1JetEtaRange==9){ out.a = 3.9672;    out.b= 1.26836;    out.c= -0.00724705;      }
    if( l1JetEtaRange==10){out.a = 4.70676;   out.b= 1.15819;    out.c= 0.00216844;       }
    if( l1JetEtaRange==11){out.a = 7.08029;   out.b= 1.16898;    out.c= -0.0088703;       }
    if( l1JetEtaRange==12){out.a = 10.8914;   out.b= 2.0546;     out.c= -0.00264239;      }
  }if(l1JetColId==6 || l1JetColId==7 || l1JetColId==8){   //10x10 out.cirout.cle
    if( l1JetEtaRange==1){ out.a = 18.1737;  out.b= 1.35362;  out.c= 0.0011629;   }
    if( l1JetEtaRange==2){ out.a = 19.421;  out.b= 1.51568;  out.c= -8.07696e-05; }
    if( l1JetEtaRange==3){ out.a = 19.5876;  out.b= 1.69937;  out.c= -0.00284367; }
    if( l1JetEtaRange==4){ out.a = 17.2096;  out.b= 1.53229;  out.c= -0.00112623; }
    if( l1JetEtaRange==5){ out.a = 15.5081;  out.b= 1.50451;  out.c= -0.00116657; }
    if( l1JetEtaRange==6){ out.a = 15.6505;  out.b= 1.54713;  out.c= -0.00151503; }
    if( l1JetEtaRange==7){ out.a = 16.1767;  out.b= 1.79387;  out.c= -0.00634596; }
    if( l1JetEtaRange==8){ out.a = 20.3783;  out.b= 1.58301;  out.c= -0.000980894;}
    if( l1JetEtaRange==9){ out.a = 19.4159;  out.b= 1.62985;  out.c= -0.00125462; }
    if( l1JetEtaRange==10){out.a = 16.7135;  out.b= 1.48148;  out.c= 0.000515484; }
    if( l1JetEtaRange==11){out.a = 19.2585;  out.b= 2.04719;  out.c= -0.00708857; }
    if( l1JetEtaRange==12){out.a = 19.4072;  out.b= 1.99785;  out.c= -0.00675134; }
  }if(l1JetColId==9 || l1JetColId==10 || l1JetColId==11){   //12x12 out.cirout.cle
    if( l1JetEtaRange==2){ out.a = 14.8886;  out.b= 1.34894;  out.c= 0.000875287;                       }
    if( l1JetEtaRange==3){ out.a = 15.0335;  out.b= 1.73809;  out.c= -0.00954032;                       }
    if( l1JetEtaRange==4){ out.a = 14.6043;  out.b= 1.92568;  out.c= -0.0195457;   out.d=0.000112846;  }
    if( l1JetEtaRange==5){ out.a = 13.4065;  out.b= 1.44792;  out.c= -0.00292822;                       }
    if( l1JetEtaRange==6){ out.a = 12.297;   out.b= 1.37412;  out.c= -0.00323287;                       }
    if( l1JetEtaRange==7){ out.a = 11.5848;  out.b= 1.70966;  out.c= -0.0130503;   out.d = 4.51248e-05; }
    if( l1JetEtaRange==8){ out.a = 12.0733;  out.b= 1.91448;  out.c= -0.0206547;   out.d = 0.000149631;     out.e =  -2.97192e-07; }
    if( l1JetEtaRange==9){ out.a = 14.9504;  out.b=1.48058;  out.c= -0.00284384; }
    if( l1JetEtaRange==10){ out.a = 16.0877; out.b= 1.41414; out.c= -0.00213214; }
    if( l1JetEtaRange==11){ out.a = 14.644;  out.b= 1.34298; out.c= 0.000797422; }
  }if(l1JetColId==12 || l1JetColId==13 || l1JetColId==14){  //8x8 squout.are
    if( l1JetEtaRange==1){  out.a = 21.8934;  out.b= 1.90989;  out.c= -0.003622;    }
    if( l1JetEtaRange==2){  out.a = 16.962;   out.b= 1.39447;  out.c= 0.000624719;  }
    if( l1JetEtaRange==3){  out.a = 18.8847;  out.b= 1.68679;  out.c= -0.001322;    }
    if( l1JetEtaRange==4){  out.a = 18.8055;  out.b= 1.74495;  out.c= -0.00198856;  }
    if( l1JetEtaRange==5){  out.a = 16.5087;  out.b= 1.66905;  out.c= -0.002084;    }
    if( l1JetEtaRange==6){  out.a = 15.4088;  out.b= 1.63109;  out.c= -0.00206412;  }
    if( l1JetEtaRange==7){  out.a = 15.4536;  out.b= 1.65226;  out.c= -0.0022026;   }
    if( l1JetEtaRange==8){  out.a = 16.1638;  out.b= 1.6986;  out.c= -0.00260055;   }
    if( l1JetEtaRange==9){  out.a = 18.1566;  out.b= 1.80647;  out.c= -0.0027926;   }
    if( l1JetEtaRange==10){ out.a = 18.7645;  out.b= 1.7464;  out.c= -0.0019839;    }
    if( l1JetEtaRange==11){ out.a = 14.8315;  out.b= 1.5122;  out.c= -0.00032842;   }
    if( l1JetEtaRange==12){ out.a = 15.1368;  out.b= 2.05485;  out.c= -0.00750572;  }
  }if(l1JetColId==15 || l1JetColId==16 || l1JetColId==17){  //9x9 squout.are
    if( l1JetEtaRange==1){  out.a = 12.5656;   out.b= 2.06762;   out.c= -0.00739286;    }
    if( l1JetEtaRange==2){  out.a = 24.2954;   out.b= -5.72677;  out.c= 0.0639203;      }
    if( l1JetEtaRange==3){  out.a = 4.8236;    out.b= 1.11877;   out.c= 0.00189857;     }
    if( l1JetEtaRange==4){  out.a = 4.61783;   out.b= 0.986284;  out.c= 0.00184376;     }
    if( l1JetEtaRange==5){  out.a = 4.73558;   out.b= 0.973591;  out.c= -0.000656711;   }
    if( l1JetEtaRange==6){  out.a = 4.74088;   out.b= 1.0175;    out.c= 0.000445015;    }
    if( l1JetEtaRange==7){  out.a = 4.7089;    out.b= 1.01932;   out.c= 0.00102063;     }
    if( l1JetEtaRange==8){  out.a = 4.47217;   out.b= 1.04528;   out.c= -0.00201983;    }
    if( l1JetEtaRange==9){  out.a = 4.40425;   out.b= 0.993732;  out.c= 0.00143311;     }
    if( l1JetEtaRange==10){ out.a = 4.86523;   out.b= 1.05992;   out.c= 0.00245171;     }
    if( l1JetEtaRange==11){ out.a = 4.75617;   out.b= 2.24675;   out.c= -0.0123394;     }
    if( l1JetEtaRange==12){ out.a = 10.4637;   out.b= 2.16339;   out.c= -0.00556654;    }
  }if(l1JetColId==18 || l1JetColId==19 || l1JetColId==20){  //10x10 squout.are
    if( l1JetEtaRange==2){   out.a = 10.0298;  out.b= -0.157889; out.c= 0.00145047;    }
    if( l1JetEtaRange==3){   out.a = 5.50012;  out.b= 0.819637;  out.c= 0.000935719;   }
    if( l1JetEtaRange==4){   out.a = 4.70199;  out.b= 0.921028;  out.c= 0.00126585;    }
    if( l1JetEtaRange==5){   out.a = 4.7418;   out.b= 0.922145;  out.c= -0.000561485;  } 
    if( l1JetEtaRange==6){   out.a = 4.74681;  out.b= 0.9218;    out.c= 0.0010017;     }
    if( l1JetEtaRange==7){   out.a = 4.78948;  out.b= 0.892598;  out.c= 0.00159056;    }
    if( l1JetEtaRange==8){   out.a = 4.45831;  out.b= 0.985886;  out.c= -0.00492871;   }
    if( l1JetEtaRange==9){   out.a = 4.38918;  out.b= 0.922127;  out.c= -3.2851e-05;   }
    if( l1JetEtaRange==10){  out.a = 4.84705;  out.b= 1.03928;   out.c= -0.00763978;   }
    if( l1JetEtaRange==11){  out.a = 6.24422;  out.b= 1.40707;   out.c= -0.0109499;    }
    if( l1JetEtaRange==12){  out.a = 9.60314;  out.b= 2.30307;   out.c= -0.0110794;    }
  }if(l1JetColId==21 || l1JetColId==22 || l1JetColId==23){  //12x12 squout.are
    if( l1JetEtaRange==2){  out.a = 6.24008 ;  out.b= 1.1747 ;  out.c= -0.00674199 ;}
    if( l1JetEtaRange==3){  out.a = 12.9523;  out.b= 1.49474;  out.c= -0.000462372; }
    if( l1JetEtaRange==4){  out.a = 12.4347;  out.b= 1.3471;  out.c= 9.66848e-05;   }
    if( l1JetEtaRange==5){  out.a = 10.7919;  out.b= 1.33801;  out.c= -0.000933652; }  
    if( l1JetEtaRange==6){  out.a = 9.85858;  out.b= 1.30005;  out.c= -0.00247947;  }
    if( l1JetEtaRange==7){  out.a = 10.0347;  out.b= 1.23912;  out.c= -0.000183111; }
    if( l1JetEtaRange==8){  out.a = 10.4139;  out.b= 1.36981;  out.c= -0.00146992;  }
    if( l1JetEtaRange==9){  out.a = 11.4242;  out.b= 1.35691;  out.c= -0.000165442; }
    if( l1JetEtaRange==10){ out.a = 12.696;  out.b= 1.44643;  out.c=-0.00201867;    }
    if( l1JetEtaRange==11){ out.a = 8.57317;  out.b= 2.06295;  out.c= -0.033153; out.d = 0.000391576; out.e =  -1.26876e-06;  }
  }
  return out;

}

//L1 AK5 calibration constants op1
L1AK5CalibrationContants L1JetAnalyzer::getL1AK5CalibrationContants(int l1JetColId, int l1JetEtaRange){

  // Initialization of output
  L1AK5CalibrationContants out;
  out.a = 0; out.b = 0; out.c = 0; out.d = 0;

  //using results from jobs  1,3,22,57,68,89-90
  if(l1JetEtaRange==0) { out.a = 4.51201;  out.b= 0.335363;  out.c= 0.00459238;    }
  if(l1JetEtaRange==1) { out.a = 9.07023;  out.b= 0.0195037; out.c= 0.0111327;     }
  if(l1JetEtaRange==2) { out.a = 2.52005;  out.b= 1.63931;   out.c= -0.00327911;   }
  if(l1JetEtaRange==3) { out.a = 2.02806;  out.b= 1.32588;   out.c= -0.000637534;  }   
  if(l1JetEtaRange==4) { out.a = 2.20294;  out.b= 1.22756;   out.c= -0.000735539;  }   
  if(l1JetEtaRange==5) { out.a = 2.9497;   out.b= 1.18908;   out.c= -0.000553623;  }   
  if(l1JetEtaRange==6) { out.a = 3.60766;  out.b= 1.07561;   out.c= -0.00105512;   }   
  if(l1JetEtaRange==7) { out.a = 2.88895;  out.b= 1.2178;    out.c= -0.000885785;  }   
  if(l1JetEtaRange==8) { out.a = 2.82073;  out.b= 1.20378;   out.c= -0.000796205;  }   
  if(l1JetEtaRange==9) { out.a = 1.47378;  out.b= 1.40259;   out.c= -0.00561334;   }   
  if(l1JetEtaRange==10){ out.a = 1.97746;  out.b= 1.3261;    out.c= -0.000945595;  }   
  if(l1JetEtaRange==11){ out.a = 2.7193;   out.b= 1.6431;    out.c= -0.00278475;   }    
  if(l1JetEtaRange==12){ out.a = 8.89084;  out.b= 0.051971;  out.c= 0.00844144;    }   
  if(l1JetEtaRange==13){ out.a = 4.74277;  out.b= 0.278456;  out.c= 0.00628881;    }    

  return out;
}

//L1 AK5 calibration constants op1 rescale
L1AK5CalibrationContantsRescale L1JetAnalyzer::getL1AK5CalibrationContantsRescale(int l1JetColId, int l1JetEtaRange, double l1JetPt){
  L1AK5CalibrationContantsRescale out;
  out.a = 0;  out.b = 0;  out.c = 0;  out.d = 0;  out.e = 0;  
  out.f = 0;  out.g = 0;  out.h = 0;  out.i = 0;  out.j = 0;
  //CALIBRATION FOR OPTION 2: CALIBRATING RHO CORRECTED L1 JETS TO PU SUB OFFLINE JETS

  if(l1JetEtaRange==0) { //-3
    if(l1JetPt < 60){
      out.a = 0.682091;    out.b = -0.205164;    out.c = 0.112949;     out.d = -0.0192092;
      out.e = 0.00208274;  out.f = -9.59003e-05; out.g = 1.79527e-06;  out.h = -2.79567e-09;
      out.i = -3.05133e-10;out.j =  2.65215e-12;
    }else{
      out.a = -14.5523;    out.b = 1.52599;  
    }
  }if(l1JetEtaRange==1) {     
    if(l1JetPt < 60){
      out.a = 0.393037;     out.b = -0.372476;    out.c = 0.217104;     out.d = -0.0524845;
      out.e = 0.00625194;   out.f = -0.00038275;  out.g = 1.29208e-05;  out.h = -2.42716e-07;
      out.i = 2.37433e-09;  out.j = -9.42285e-12;
    }else{ //>100
      out.a = 16.9244;      out.b = 0.361365;
    }
  }if(l1JetEtaRange==2) { //-2
    if(l1JetPt < 120){
      out.a = 0.456893;     out.b = -0.324723;      out.c = 0.0986209;   out.d = -0.0021623; 
      out.e = 1.92619e-05;  out.f = -6.04222e-08; 
    }else{
      out.a = -12.4889;     out.b = 1.28213;
    } 
  }if(l1JetEtaRange==3) { 
    if(l1JetPt < 150){
      out.a = 1.31978;      out.b = -0.0952089;     out.c = 0.350186;     out.d = -0.0219245;
      out.e = 0.000666506;  out.f = -1.14448e-05;   out.g = 1.1624e-07;   out.h = -6.91227e-10;
      out.i = 2.22195e-12;  out.j = -2.97864e-15;
    }else{  
      out.a = 26.8031;      out.b = 1.05623;
    }
  }if(l1JetEtaRange==4) { //-1
    out.a = 2.20987 ;        out.b = 0.833414 ;        out.c = 0.267548 ;       out.d = -0.0186549 ;
    out.e = 0.000579946 ;    out.f = -9.77537e-06 ;    out.g = 9.50211e-08 ;    out.h = -5.31444e-10 ;
    out.i = 1.58667e-12 ;    out.j = -1.95722e-15 ;
  }if(l1JetEtaRange==5) { 
      out.a = 1.04226;       out.b = 2.0093;    out.c = -0.0118415;     out.d = 4.78419e-05 ;
  }if(l1JetEtaRange==6) { //0
      out.a = 1.0659 ;    out.b = 1.99774;    out.c = -0.0111557;   out.d = 3.99304e-05;
  }if(l1JetEtaRange==7) { 
      out.a = 0.92738;    out.b = 1.87435;    out.c = -0.00898217;  out.d = 3.05085e-05;
  }if(l1JetEtaRange==8) { //1
      out.a = 1.07099;    out.b= 1.98441;     out.c= -0.0106089;    out.d= 3.78115e-05;
  }if(l1JetEtaRange==9) { 
      out.a = 0.425475;   out.b= 2.20663;     out.c= -0.0166944;    out.d= 0.000134281;
  }if(l1JetEtaRange==10){ //2
      out.a = -0.868285;  out.b= 1.95554;     out.c= -0.00933064;   out.d= 3.37264e-05;
  }if(l1JetEtaRange==11){ 
    out.a = 0.386037;    out.b = -0.297345;    out.c = 0.093929;     out.d = -0.00206772;
    out.e = 2.00733e-05; out.f = -8.88054e-08; out.g = 1.4618e-10; 
  }if(l1JetEtaRange==12){ //3
    if(l1JetPt < 60){
      out.a = -0.0544346;     out.b = 0.231386;    out.c = -0.0636845;     out.d = 0.00668052;
      out.e = -0.00016757;    out.f = 1.27613e-06;
    }
    else{
      out.a = -116.709 ;      out.b = 3.69675;    out.c = -0.012712;
    }
  }if(l1JetEtaRange==13){ 
    if( l1JetPt < 50){
      out.a = 0.644164;      out.b = -0.162112;      out.c = 0.0990824;      out.d = -0.0168025; 
      out.e = 0.00207218;    out.f = -0.000122797;   out.g = 3.85766e-06;    out.h = -6.63742e-08; 
      out.i = 5.89888e-10;   out.j = -2.11311e-12; 
   }else{
     out.a =     19.6176;out.b =   0.676703 ;
   }
 }

 return out;
}  

L1AK5CalibrationContants2 L1JetAnalyzer::getL1AK5CalibrationContants2(int l1JetColId, int l1JetEtaRange, double l1JetPt){

  // Initialization of output
  L1AK5CalibrationContants2 out;
  out.a = 0; out.b = 0; out.c = 0; out.d = 0; out.e = 0; out.f = 0;  
  out.g = 0; out.h = 0; out.i = 0; out.j = 0;
  //CALIBRATION FOR OPTION 2: CALIBRATING RHO CORRECTED L1 JETS TO PU SUB OFFLINE JETS

  //Using 700 events run locally. Linear fits best I think
   if(l1JetEtaRange==0) { 
     out.a = 0.461685;      out.b= 0.000528405;   out.c = 0.0199043;    out.d = -0.00239633;
     out.e = 0.000131892;   out.f = -2.92814e-06; out.g=  3.21297e-08;  out.h = -1.87253e-10;
     out.i = 5.57671e-13;   out.j= -6.69448e-16;
   }if(l1JetEtaRange==1) { 
     if(l1JetPt < 100){
      out.a = 0.245306;     out.b = -0.0175123;   out.c = 0.000924904;  out.d = 0.000564157;
      out.e = -7.11798e-05; out.f = 3.60797e-06;  out.g = -7.45287e-08; out.h = 6.72995e-10;
      out.i = -2.22022e-12;
     }else{ //>100
      out.a   =     -39.5836 ; out.b   =      1.40207 ;
     }
  }if(l1JetEtaRange==2) { out.a = 0.661592;  out.b= 1.39546;   out.c= 0.0839384;   out.d= -0.00294603;
    out.e=3.84166e-05; out.f=-2.14728e-07; out.g=4.30121e-10;        }
   if(l1JetEtaRange==3) { out.a = 4.69352;   out.b= 3.43598;   out.c= -0.0574304;  out.d= 0.000736731;
     out.e=-4.71088e-06; out.f=1.39882e-08; out.g=-1.46108e-11;   }   
   if(l1JetEtaRange==4) { out.a = 8.73216;   out.b= 3.32547;   out.c= -0.0804369;  out.d= 0.00161277;
     out.e=-1.68876e-05; out.f=8.46531e-08; out.g=-1.59173e-10;   }   
   if(l1JetEtaRange==5) { out.a = 9.38982;   out.b= 2.80158;   out.c= -0.0589441;  out.d= 0.00119714;
     out.e=-1.31281e-05; out.f=6.84816e-08; out.g=-1.32038e-10;   }   
   if(l1JetEtaRange==6) { out.a = 8.88216;   out.b= 2.60116;   out.c= -0.048385;   out.d= 0.00101379;
     out.e=-1.21377e-05; out.f=6.82927e-08; out.g=-1.396e-10;     }   
   if(l1JetEtaRange==7) { out.a = 8.85892;   out.b= 2.64762;   out.c= -0.0502263;  out.d= 0.00104999;
     out.e=-1.25472e-05; out.f=7.07149e-08; out.g=-1.45041e-10;   }   
   if(l1JetEtaRange==8) { out.a = 9.2381;    out.b= 2.85245;   out.c= -0.0617271;  out.d= 0.00126547;
     out.e=-1.4001e-05; out.f=7.38509e-08; out.g=-1.44073e-10;    }   
   if(l1JetEtaRange==9) { out.a = 8.56558;   out.b= 3.361;     out.c= -0.0802383;  out.d= 0.00155765;
     out.e=-1.56934e-05; out.f=7.58546e-08; out.g=-1.37889e-10;   }   
   if(l1JetEtaRange==10){ out.a = 4.68284;   out.b= 3.50193;   out.c= -0.0643231;  out.d= 0.000942738;
     out.e=-7.25944e-06; out.f=2.7836e-08; out.g=-4.1973e-11;     }   
   if(l1JetEtaRange==11){ 
     if(l1JetPt < 135){
       out.a = 1.59046;    out.b = 0.0115098;     out.c = 0.455886;     out.d = -0.0332072;
       out.e = 0.0011856;  out.f = -2.40843e-05;  out.g = 2.90573e-07;  out.h = -2.05633e-09;
       out.i = 7.8714e-12; out.j = -1.25665e-14;
     }else{
       out.a = 24.9043;    out.b = 1.28753 ;
     } 
   }if(l1JetEtaRange==12){ 
     if(l1JetPt < 90){
      out.a = 0.246531;     out.b = -0.0587921;   out.c = 0.0196874;    out.d = -0.00300772;
      out.e = 0.000260743;  out.f = -1.28423e-05; out.g = 3.74493e-07;  out.h = -6.05221e-09;
      out.i = 4.94285e-11;  out.j = -1.58683e-13;
     }
     else{
       out.a = -205.6;      out.b = 4.37228;      out.c = -0.0133536;
     }
   }if(l1JetEtaRange==13){ 
     if( l1JetPt < 95){
      out.a =     0.448692;   out.b = 0.0158405;    out.c = 0.0145457;    out.d = -0.0015955;
      out.e =  7.94352e-05;   out.f = -1.44525e-06; out.g = 2.04844e-08;  out.h = -3.81301e-10;
      out.i =  4.30313e-12;   out.j = -1.7399e-14;
     }else{
       out.a = -538.74;       out.b = 12.3227;      out.c =   -0.0622738;
     }
   }

  return out;

}                                            
                                            
// ------------ method called once each job just before starting event loop  ------------
void L1JetAnalyzer::beginJob(){}                    
                                             
// ------------ method called once each job just after ending the event loop  ------------
void L1JetAnalyzer::endJob(){}                      
                                             
// ------------ method called when starting to processes a run  ------------
void L1JetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&){}
                                             
// ------------ method called when ending the processing of a run  ------------
void L1JetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&){}

// ------------ method called when starting to processes a luminosity block  ------------
void L1JetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

// ------------ method called when ending the processing of a luminosity block  ------------
void L1JetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1JetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1JetAnalyzer);
