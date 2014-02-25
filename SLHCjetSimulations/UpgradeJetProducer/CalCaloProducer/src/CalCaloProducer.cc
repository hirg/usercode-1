

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

using namespace edm;
using namespace std;
using namespace reco;

class CalCaloProducer : public edm::EDProducer {
   public:
      explicit CalCaloProducer(const edm::ParameterSet&);
      ~CalCaloProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
     ParameterSet conf_;
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
CalCaloProducer::CalCaloProducer(const edm::ParameterSet& iConfig):
conf_(iConfig)
{
   //register your products
produces<CaloJetCollection>();
  
}


CalCaloProducer::~CalCaloProducer()
{
 

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CalCaloProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  const JetCorrector* jetCorrector = JetCorrector::getJetCorrector(conf_.getParameter<string>("JetCorrector"),iSetup);
  edm::Handle<reco::CaloJetCollection> Calojets;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CaloJets"), Calojets);

auto_ptr< CaloJetCollection > outputColl(new CaloJetCollection());

  for(reco::CaloJetCollection::const_iterator it = Calojets->begin(); it != Calojets->end() ; ++it) {
    double scale = jetCorrector->correction(*it,iEvent,iSetup);
    //    cout<<" off (pt,eta,phi) =("<<it->p4().pt()<<","<<it->p4().eta()<<","<<it->p4().phi()<<") Area= "<<it->jetArea()<<" "<<(it->p4().e()/it->jetArea())<<endl;
   // double scaled_pt = it->pt()*scale;
//math::XYZTLorentzVector ff;
//    ff.SetXYZT(it->px()*scale,it->py()*scale,it->pz(),it->t());
  //  cout<<"QQ "<<it->pt()<<" "<<scale<<endl;
    outputColl->push_back(CaloJet(scale*it->p4(),it->vertex(),it->getSpecific ()));
   }
 iEvent.put(outputColl);

 
}

// ------------ method called once each job just before starting event loop  ------------
void 
CalCaloProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CalCaloProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
CalCaloProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CalCaloProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CalCaloProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CalCaloProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CalCaloProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CalCaloProducer);
