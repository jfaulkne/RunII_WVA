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

class GenJetAK8 : public edm::EDProducer {
public:
	explicit GenJetAK8(const edm::ParameterSet&);
	~GenJetAK8();
	
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
GenJetAK8::GenJetAK8(const edm::ParameterSet& iConfig)
{
	//register your produc
	GenJetCollTag_ 				= 	iConfig.getParameter<edm::InputTag >("GenJetCollTag");	

	const std::string string1("GenJetAK8");
	produces<std::vector<pat::Jet> > (string1).setBranchAlias(string1);
        const std::string string2("GenAK8prunedMass");
        produces<std::vector<double> > (string2).setBranchAlias(string2);
        const std::string string3("GenAK8softdropMass");
	produces<std::vector<double> > (string3).setBranchAlias(string3);

        const std::string string4("GenAK8tau1");
        produces<std::vector<double> > (string4).setBranchAlias(string4);
        const std::string string5("GenAK8tau2");
        produces<std::vector<double> > (string5).setBranchAlias(string5);
        const std::string string6("GenAK8tau3");
        produces<std::vector<double> > (string6).setBranchAlias(string6);
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


GenJetAK8::~GenJetAK8()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenJetAK8::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;	

	std::auto_ptr< std::vector<pat::Jet> > selectedJet(new std::vector<pat::Jet>);

	std::auto_ptr< std::vector<double> > GenAK8prunedMass(new std::vector<double>);
	std::auto_ptr< std::vector<double> > GenAK8softdropMass(new std::vector<double>);

	std::auto_ptr< std::vector<double> > GenAK8tau1(new std::vector<double>);
	std::auto_ptr< std::vector<double> > GenAK8tau2(new std::vector<double>);
	std::auto_ptr< std::vector<double> > GenAK8tau3(new std::vector<double>);

	Handle<edm::View<pat::Jet> > jet;
	iEvent.getByLabel(GenJetCollTag_,jet);

	for(size_t i=0; i<jet->size();i++)
	{
	  selectedJet->push_back(pat::Jet(jet->at(i)));

	  GenAK8prunedMass->push_back( jet->at(i).userFloat("ak8GenJetsPrunedMass"));
	  GenAK8softdropMass->push_back( jet->at(i).userFloat("ak8GenJetsSoftDropMass"));//	  std::cout<<jet->at(i).pt()<<std::endl;

	  GenAK8tau1->push_back( jet->at(i).userFloat("NjettinessGenAK8:tau1"));
	  GenAK8tau2->push_back( jet->at(i).userFloat("NjettinessGenAK8:tau2"));
	  GenAK8tau3->push_back( jet->at(i).userFloat("NjettinessGenAK8:tau3"));
	}

	const std::string string1("GenJetAK8");
	iEvent.put(selectedJet,string1);
        const std::string string2("GenAK8prunedMass");
        iEvent.put(GenAK8prunedMass,string2);
        const std::string string3("GenAK8softdropMass");
        iEvent.put(GenAK8softdropMass,string3);	

        const std::string string4("GenAK8tau1");
        iEvent.put(GenAK8tau1,string4);	
        const std::string string5("GenAK8tau2");
        iEvent.put(GenAK8tau2,string5);	
        const std::string string6("GenAK8tau3");
        iEvent.put(GenAK8tau3,string6);	
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenJetAK8::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenJetAK8::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
GenJetAK8::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenJetAK8::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenJetAK8::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenJetAK8::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenJetAK8::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenJetAK8);
