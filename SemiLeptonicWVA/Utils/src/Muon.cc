// -*- C++ -*-
//
// Package:    Muon
// Class:      Muon
// 
/**\class Muon Muon.cc RA2Classic/Muon/src/Muon.cc
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//
// class declaration
//
class Muon : public edm::EDProducer {
public:
	explicit Muon(const edm::ParameterSet&);
	~Muon();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	edm::InputTag MuTag_;
        edm::InputTag VertexTag_;
	edm::InputTag RhoTag_;
  double MinPt_;
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
Muon::Muon(const edm::ParameterSet& iConfig)
{
	MuTag_ = iConfig.getParameter<edm::InputTag>("MuTag");
        MinPt_ = iConfig.getParameter <double> ("MinPt");
        VertexTag_ = iConfig.getParameter<edm::InputTag>("VertexTag");
	RhoTag_ = iConfig.getParameter<edm::InputTag>("RhoTag");
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
        produces<std::vector<pat::Muon> >();
        const std::string string0("pt");
	produces<std::vector<double> > (string0).setBranchAlias(string0);
        const std::string string1("eta");
	produces<std::vector<double> > (string1).setBranchAlias(string1);
        const std::string string2("phi");
	produces<std::vector<double> > (string2).setBranchAlias(string2);
        const std::string string3("e");
	produces<std::vector<double> > (string3).setBranchAlias(string3);
        const std::string string4("isHighPt");
	produces<std::vector<bool> > (string4).setBranchAlias(string4);
        const std::string string5("isLoose");
	produces<std::vector<bool> > (string5).setBranchAlias(string5);
        const std::string string10("isMedium");
        produces<std::vector<bool> > (string10).setBranchAlias(string10);
        const std::string string11("isTight");
        produces<std::vector<bool> > (string11).setBranchAlias(string11);
        const std::string string6("type");
	produces<std::vector<int> > (string6).setBranchAlias(string6);
        const std::string string7("charge");
	produces<std::vector<int> > (string7).setBranchAlias(string7);
        const std::string string8("mass");
	produces<std::vector<double> > (string8).setBranchAlias(string8);
        const std::string string9("pfDeltaCorrRelIso");
	produces<std::vector<double> > (string9).setBranchAlias(string9);
        const std::string string12("pfRelIso");
	produces<std::vector<double> > (string12).setBranchAlias(string12);
        const std::string string13("photonIso");
	produces<std::vector<double> > (string13).setBranchAlias(string13);
        const std::string string14("neutralHadIso");
	produces<std::vector<double> > (string14).setBranchAlias(string14);
        const std::string string15("chargedHadIso");
	produces<std::vector<double> > (string15).setBranchAlias(string15);
        const std::string string16("trackIso");
	produces<std::vector<double> > (string16).setBranchAlias(string16);
        const std::string string17("isPFMuon");
	produces<std::vector<bool> > (string17).setBranchAlias(string17);
}


Muon::~Muon()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
Muon::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	std::auto_ptr<std::vector<pat::Muon> > prodMu(new std::vector<pat::Muon>());

	std::auto_ptr< std::vector<double> > pt(new std::vector<double>);
	std::auto_ptr< std::vector<double> > eta(new std::vector<double>);
	std::auto_ptr< std::vector<double> > phi(new std::vector<double>);
	std::auto_ptr< std::vector<double> > e(new std::vector<double>);
	std::auto_ptr< std::vector<bool> > isHighPt(new std::vector<bool>);
	std::auto_ptr< std::vector<bool> > isLoose(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > isMedium(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > isTight(new std::vector<bool>);
	std::auto_ptr< std::vector<bool> > isPFMuon(new std::vector<bool>);
	std::auto_ptr< std::vector<int> > type(new std::vector<int>);
	std::auto_ptr< std::vector<int> > charge(new std::vector<int>);
	std::auto_ptr< std::vector<double> > mass(new std::vector<double>);
	std::auto_ptr< std::vector<double> > pfDeltaCorrRelIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > pfRelIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > photonIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > neutralHadIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > chargedHadIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > trackIso(new std::vector<double>);

	using namespace edm;
	using namespace reco;
	using namespace pat;

	edm::Handle< edm::View<pat::Muon> > Muons;
	iEvent.getByLabel(MuTag_,Muons);
	edm::Handle< edm::View<reco::Vertex> > Vertices;
        iEvent.getByLabel(VertexTag_,Vertices);
	edm::Handle<double> rho_ ;
	iEvent.getByLabel(RhoTag_, rho_);

	if( Muons.isValid() ) {
		for(unsigned int i=0; i<Muons->size();i++)
		{
		  if (Muons->at(i).pt()<MinPt_) continue;

		  bool isHighPtId = Muons->at(i).isHighPtMuon(Vertices->at(0));
		  bool isLooseId = Muons->at(i).isLooseMuon();
                  bool isMediumId = Muons->at(i).isMediumMuon();
                  bool isTightId = Muons->at(i).isTightMuon(Vertices->at(0));
		  bool isPF = Muons->at(i).isPFMuon();
		
		  prodMu->push_back(pat::Muon(Muons->at(i)));
		  isHighPt->push_back( isHighPtId );
		  isLoose->push_back( isLooseId );
                  isMedium->push_back( isMediumId );
                  isTight->push_back( isTightId );
		  isPFMuon->push_back( isPF);
		  type->push_back(Muons->at(i).pdgId());
		  charge->push_back(Muons->at(i).charge());
		  e->push_back(Muons->at(i).energy());
		  eta->push_back(Muons->at(i).eta());
		  mass->push_back(Muons->at(i).mass());
		  pt->push_back(Muons->at(i).pt());
		  phi->push_back(Muons->at(i).phi());

		  pfDeltaCorrRelIso->push_back((Muons->at(i).chargedHadronIso() + std::max(0., Muons->at(i).neutralHadronIso() + Muons->at(i).photonIso() - 0.5*Muons->at(i).puChargedHadronIso()))/Muons->at(i).pt());
		  pfRelIso->push_back((Muons->at(i).chargedHadronIso() + Muons->at(i).neutralHadronIso()+ Muons->at(i).photonIso())/Muons->at(i).pt());
		  photonIso->push_back(Muons->at(i).photonIso());
		  neutralHadIso->push_back(Muons->at(i).neutralHadronIso());
		  chargedHadIso->push_back(Muons->at(i).chargedHadronIso());
		  trackIso->push_back(Muons->at(i).trackIso());		  

		}
	}

	const std::string string00("");
        iEvent.put(prodMu );

	const std::string string0("pt");
	iEvent.put(pt,string0);
	const std::string string1("eta");
	iEvent.put(eta,string1);
	const std::string string2("phi");
	iEvent.put(phi,string2);
	const std::string string3("e");
	iEvent.put(e,string3);
	const std::string string4("isHighPt");
	iEvent.put(isHighPt,string4);
	const std::string string5("isLoose");
	iEvent.put(isLoose,string5);
        const std::string string10("isMedium");
        iEvent.put(isMedium,string10);
        const std::string string11("isTight");
        iEvent.put(isTight,string11);
	const std::string string6("type");
	iEvent.put(type,string6);
	const std::string string7("charge");
	iEvent.put(charge,string7);
	const std::string string8("mass");
	iEvent.put(mass,string8);
	const std::string string9("pfDeltaCorrRelIso");
	iEvent.put(pfDeltaCorrRelIso,string9);
	const std::string string12("pfRelIso");
	iEvent.put(pfRelIso,string12);
	const std::string string13("photonIso");
	iEvent.put(photonIso,string13);
	const std::string string14("neutralHadIso");
	iEvent.put(neutralHadIso,string14);
	const std::string string15("chargedHadIso");
	iEvent.put(chargedHadIso,string15);
	const std::string string16("trackIso");
	iEvent.put(trackIso,string16);
	const std::string string17("isPFMuon");
	iEvent.put(isPFMuon,string17);
}


// ------------ method called once each job just before starting event loop  ------------
void 
Muon::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Muon::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
Muon::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Muon::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Muon::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Muon::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Muon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Muon);
