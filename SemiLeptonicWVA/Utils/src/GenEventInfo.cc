// -*- C++ -*-
//
// Package:    GenEventInfo
// Class:      GenEventInfo
// 
/**\class GenEventInfo GenEventInfo.cc RA2Classic/GenEventInfo/src/GenEventInfo.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 *     [Notes on implementation]
 */
//
// Original Author:  Arne-Rasmus Draeger,68/111,4719,
//         Created:  Fri Apr 11 16:35:33 CEST 2014
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//
// class declaration
//

class GenEventInfo : public edm::EDProducer {
public:
  explicit GenEventInfo( const edm::ParameterSet& iConfig);
	~GenEventInfo();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  edm::EDGetTokenT< GenEventInfoProduct > geneventToken_; 
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
GenEventInfo::GenEventInfo(const edm::ParameterSet& iConfig):
  geneventToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("geneventToken")))
{

	produces<double>("genEventWeight");
	/* Examples
	 *   produces<ExampleData2>();
	 * 
	 *   //if do put with a label
	 *   produces<ExampleData2>("label");
	 * 
	 *   //if you want to put into the Run
	 *   produces<ExampleData2,InRun>();
	 */
	//now do what ever other initialization is needed
	
}


GenEventInfo::~GenEventInfo()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenEventInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;	

	edm::Handle< GenEventInfoProduct > geneventInfo_;
	iEvent.getByToken(geneventToken_, geneventInfo_); 

	double eventWeight =1.;
	if (geneventInfo_.isValid())
	  eventWeight = geneventInfo_->weight();

	std::auto_ptr<double> htp(new double(eventWeight));
	iEvent.put(htp,"genEventWeight");

}

// ------------ method called once each job just before starting event loop  ------------
void 
GenEventInfo::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenEventInfo::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
GenEventInfo::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenEventInfo::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenEventInfo::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenEventInfo::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenEventInfo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenEventInfo);
