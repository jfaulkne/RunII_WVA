// -*- C++ -*-
//
// Package:    GenJet
// Class:      GenJet
// 
/**\class GenJet GenJet.cc RA2Classic/GenJet/src/GenJet.cc
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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//
// class declaration
//

class GenJet : public edm::EDProducer {
public:
	explicit GenJet(const edm::ParameterSet&);
	~GenJet();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag GenJetCollTag_;

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
GenJet::GenJet(const edm::ParameterSet& iConfig)
{
	//register your produc
	GenJetCollTag_ 				= 	iConfig.getParameter<edm::InputTag >("GenJetCollTag");	
	const std::string string1("GenJet");
	produces<std::vector<reco::GenJet> > (string1).setBranchAlias(string1);
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


GenJet::~GenJet()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenJet::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;	
	std::auto_ptr< std::vector<reco::GenJet> > selectedJet(new std::vector<reco::GenJet>);

	Handle<edm::View<reco::GenJet> > jet;
	iEvent.getByLabel(GenJetCollTag_,jet);

	for(size_t i=0; i<jet->size();i++)
	{
	  selectedJet->push_back(reco::GenJet(jet->at(i)));
	  //	  std::cout<<jet->at(i).pt()<<std::endl;
	}
	const std::string string1("GenJet");

	iEvent.put(selectedJet,string1);
	
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenJet::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
GenJet::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenJet::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenJet::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenJet::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenJet);
