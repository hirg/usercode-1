#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "CondFormats/DataRecord/interface/L1CaloEcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CalibCalorimetry/EcalTPGTools/interface/EcalTPGScale.h"
#include "JetSLHC/L1CaloCandidate/interface/TriggerTowerGeometry.h"
#include "TLorentzVector.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

using namespace std;
using namespace edm;
using namespace reco;
class L1CaloCandidate : public EDProducer {
   public:
      explicit L1CaloCandidate(const ParameterSet&);
      ~L1CaloCandidate();

      static void fillDescriptions(ConfigurationDescriptions& descriptions);

   private:
  virtual void beginJob() ;
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(Run&, EventSetup const&);
  virtual void endRun(Run&, EventSetup const&);
  virtual void beginLuminosityBlock(LuminosityBlock&, EventSetup const&);
  virtual void endLuminosityBlock(LuminosityBlock&, EventSetup const&);
  void addEcal( const int &, const int &, const int &, const bool & );
  
  ParameterSet conf_;
  const L1CaloEcalScale *mEcalScale;
  const L1CaloHcalScale *mHcalScale;
  TriggerTowerGeometry mGeometry;
  // ----------member data ---------------------------
};


L1CaloCandidate::L1CaloCandidate(const ParameterSet& iConfig):
  conf_(iConfig)
{
  produces<reco::PFCandidateCollection>();
}


L1CaloCandidate::~L1CaloCandidate()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1CaloCandidate::produce(Event& iEvent, const EventSetup& iSetup)
{
  

  std::auto_ptr< reco::PFCandidateCollection > outputColl(new reco::PFCandidateCollection());

  edm::ESHandle < L1CaloEcalScale > lEcalScaleHandle;
  iSetup.get < L1CaloEcalScaleRcd > (  ).get( lEcalScaleHandle );
  mEcalScale = lEcalScaleHandle.product(  );

  

  edm::ESHandle < L1CaloHcalScale > lHcalScaleHandle;
  iSetup.get < L1CaloHcalScaleRcd > (  ).get( lHcalScaleHandle );
  mHcalScale = lHcalScaleHandle.product(  );
  

  EcalTPGScale ecalScale ;
  ecalScale.setEventSetup(iSetup) ;

  Handle< EcalTrigPrimDigiCollection > lEcalDigiHandle;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("srcEcal"),lEcalDigiHandle);

  ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  


  //  cout<<"A "<<lEcalDigiHandle->size()<< endl;

  for ( EcalTrigPrimDigiCollection::const_iterator lEcalTPItr = lEcalDigiHandle->begin(  ); lEcalTPItr != lEcalDigiHandle->end(  ); ++lEcalTPItr ){

    const EcalTrigTowerDetId TPtowid=lEcalTPItr->id();
    float ENE= ecalScale.getTPGInGeV( lEcalTPItr->compressedEt(  ),TPtowid); 
   
    if (lEcalTPItr->compressedEt()){
      if (ENE>0.){

	//	cout<<"RRR "<<ENE<<" "<<
	double phi=mGeometry.phi( TPtowid.iphi() );
	double eta=mGeometry.eta( TPtowid.ieta() );

	PFCandidate a(0,math::XYZTLorentzVectorD(ENE*cos(phi),ENE*sin(phi),ENE*sinh(eta),ENE*cosh(eta)),reco::PFCandidate::gamma);
	outputColl->push_back(a);
      }
    }
  }
  edm::Handle < HcalTrigPrimDigiCollection > lHcalDigiHandle;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("srcHcal"),  lHcalDigiHandle );
  
  for ( HcalTrigPrimDigiCollection::const_iterator lHcalTPItr = lHcalDigiHandle->begin(  ); lHcalTPItr != lHcalDigiHandle->end(  ); ++lHcalTPItr ){
    //    int isHF = data_tp->id().ietaAbs() >= 29 ? 1 : 0;
    float ENE=mHcalScale->et( lHcalTPItr->SOI_compressedEt(),
				abs( lHcalTPItr->id().ieta()),
				( lHcalTPItr->id().ieta() > 0 ? +1 : -1 ) );
    
    if (ENE>0.){
      //	cout<<"RRR "<<ENE<<" "<<
      double phi=mGeometry.phi( lHcalTPItr->id().iphi() );
      double eta=mGeometry.eta( lHcalTPItr->id().ieta() );
      PFCandidate a(0,math::XYZTLorentzVectorD(ENE*cos(phi),ENE*sin(phi),ENE*sinh(eta),ENE*cosh(eta)),reco::PFCandidate::h0);
      outputColl->push_back(a);
    }
    
  }
 
  
  iEvent.put(outputColl);
  //  cout<<"C"<<endl;
}

// ------------ method called once each job just before starting event loop  ------------
void 
L1CaloCandidate::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1CaloCandidate::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
L1CaloCandidate::beginRun(Run&, EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
L1CaloCandidate::endRun(Run&, EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
L1CaloCandidate::beginLuminosityBlock(LuminosityBlock&, EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
L1CaloCandidate::endLuminosityBlock(LuminosityBlock&, EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1CaloCandidate::fillDescriptions(ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void L1CaloCandidate::addEcal( const int &aCompressedEt, const int &aIeta, const int &aIphi, const bool & aFG )
{
 
  if ( aCompressedEt > 0 )
    {
      
      int lET = ( int )( 2 * mEcalScale->et( aCompressedEt,
					     abs( aIeta ),
					     ( aIeta > 0 ? +1 : -1 ) ) );
      //cout<<"SS "<<aCompressedEt<<" "<<lET<<" "<<aIeta<<" "<<aIphi<< endl;
      if( lET ){
	//			l1slhc::L1CaloTower lCaloTower( aIeta, aIphi );
	//			lCaloTower.setEcal( lET, aFG );
	//			mCaloTowers.insert( aIeta , aIphi , lCaloTower );
      }
      
    }
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1CaloCandidate);
