// -*- C++ -*-
//
// Package:    Photon
// Class:      Photon
// 
/**\class Photon Photon.cc RA2Classic/Photon/src/Photon.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 *     [Notes on implementation]
 */
//
// Original Author:  James Faulkner
//         Created:  Fri June 5 2015
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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Common/interface/ValueMap.h"

//
// class declaration
//

class Photon : public edm::EDProducer {
public:
	explicit Photon(const edm::ParameterSet&);
	~Photon();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag PhoTag_;

  // ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> > phoLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoTightIdMapToken_;

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
Photon::Photon(const edm::ParameterSet& iConfig):
  phoLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoLooseIdMap"))),
  phoMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoMediumIdMap"))),
  phoTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoTightIdMap")))
{
	PhoTag_ = iConfig.getParameter<edm::InputTag>("PhoTag");

	//register your products
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
	//register your products
        produces<std::vector<pat::Photon> >();
        const std::string string0("pt");
	produces<std::vector<double> > (string0).setBranchAlias(string0);
        const std::string string1("eta");
	produces<std::vector<double> > (string1).setBranchAlias(string1);
        const std::string string2("phi");
	produces<std::vector<double> > (string2).setBranchAlias(string2);
        const std::string string3("e");
	produces<std::vector<double> > (string3).setBranchAlias(string3);
        const std::string string4("isLoose");
	produces<std::vector<bool> > (string4).setBranchAlias(string4);
        const std::string string5("isMedium");
        produces<std::vector<bool> > (string5).setBranchAlias(string5);
        const std::string string6("isTight");
        produces<std::vector<bool> > (string6).setBranchAlias(string6);
}


Photon::~Photon()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
Photon::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	std::auto_ptr<std::vector<pat::Photon> > prodPho(new std::vector<pat::Photon>());

	std::auto_ptr< std::vector<double> > pt(new std::vector<double>);
	std::auto_ptr< std::vector<double> > eta(new std::vector<double>);
	std::auto_ptr< std::vector<double> > phi(new std::vector<double>);
	std::auto_ptr< std::vector<double> > e(new std::vector<double>);
	std::auto_ptr< std::vector<bool> > isLoose(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > isMedium(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > isTight(new std::vector<bool>);

	using namespace edm;
	using namespace reco;
	using namespace pat;

	edm::Handle< edm::View<pat::Photon> > Photons;
	iEvent.getByLabel(PhoTag_, Photons);

	// Get the electron ID data from the event stream.
	// Note: this implies that the VID ID modules have been run upstream.
	// If you need more info, check with the EGM group.
	edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
	edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
	edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
	iEvent.getByToken(phoLooseIdMapToken_ ,loose_id_decisions);
	iEvent.getByToken(phoMediumIdMapToken_,medium_id_decisions);
	iEvent.getByToken(phoTightIdMapToken_,tight_id_decisions);

	if( Photons.isValid() ) {
		for(unsigned int i=0; i<Photons->size();i++)
		{
		  const auto ph = Photons->ptrAt(i);
		  bool isPassLoose = (*loose_id_decisions)[ph];
		  bool isPassMedium = (*medium_id_decisions)[ph];
		  bool isPassTight = (*tight_id_decisions)[ph];

		  prodPho->push_back(pat::Photon(Photons->at(i)));
		  isLoose->push_back( isPassLoose );
                  isMedium->push_back( isPassMedium );
                  isTight->push_back( isPassTight );
		  e->push_back(Photons->at(i).energy());
		  eta->push_back(Photons->at(i).superCluster()->eta());
		  pt->push_back(Photons->at(i).pt());
		  phi->push_back(Photons->at(i).phi());
		}
	}

	const std::string string00("");
        iEvent.put(prodPho );

	const std::string string0("pt");
	iEvent.put(pt,string0);
	const std::string string1("eta");
	iEvent.put(eta,string1);
	const std::string string2("phi");
	iEvent.put(phi,string2);
	const std::string string3("e");
	iEvent.put(e,string3);
	const std::string string4("isLoose");
	iEvent.put(isLoose,string4);
        const std::string string5("isMedium");
        iEvent.put(isMedium,string5);
        const std::string string6("isTight");
        iEvent.put(isTight,string6);
}


// ------------ method called once each job just before starting event loop  ------------
void 
Photon::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Photon::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
Photon::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Photon::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Photon::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Photon::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Photon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Photon);
