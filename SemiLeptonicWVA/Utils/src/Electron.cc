// -*- C++ -*-
//
// Package:    Electron
// Class:      Electron
// 
/**\class Electron Electron.cc RA2Classic/Electron/src/Electron.cc
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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//
// class declaration
//

class Electron : public edm::EDProducer {
public:
	explicit Electron(const edm::ParameterSet&);
	~Electron();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        float dEtaInSeed( const pat::Electron &ele );
	edm::InputTag EleTag_;
	edm::InputTag VertexTag_;
	edm::InputTag RhoTag_;
  double MinPt_;
  edm::EDGetTokenT<reco::ConversionCollection> ConversionTag_;
  edm::EDGetTokenT<reco::BeamSpot> BeamSpotTag_;

  // ID decisions objects
  //  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
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
Electron::Electron(const edm::ParameterSet& iConfig):
  ConversionTag_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("ConversionTag"))),
  BeamSpotTag_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotTag"))),
  //  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap")))
{
	EleTag_ = iConfig.getParameter<edm::InputTag>("EleTag");
        MinPt_ = iConfig.getParameter <double> ("MinPt");
	VertexTag_ = iConfig.getParameter<edm::InputTag>("VertexTag");
	RhoTag_ = iConfig.getParameter<edm::InputTag>("RhoTag");
	//	ConversionTag_ = iConfig.getParameter<edm::InputTag>("ConversionTag");
	//	BeamSpotTag_ = iConfig.getParameter<edm::InputTag>("BeamSpotTag");

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
        produces<std::vector<pat::Electron> >();
        const std::string string0("pt");
	produces<std::vector<double> > (string0).setBranchAlias(string0);
        const std::string string1("eta");
	produces<std::vector<double> > (string1).setBranchAlias(string1);
        const std::string string2("phi");
	produces<std::vector<double> > (string2).setBranchAlias(string2);
        const std::string string3("e");
	produces<std::vector<double> > (string3).setBranchAlias(string3);
        const std::string string4("isHEEP");
	produces<std::vector<bool> > (string4).setBranchAlias(string4);
	//  const std::string string5("isHEEPv50");
	//	produces<std::vector<bool> > (string5).setBranchAlias(string5);
        const std::string string6("type");
	produces<std::vector<int> > (string6).setBranchAlias(string6);
        const std::string string7("charge");
	produces<std::vector<int> > (string7).setBranchAlias(string7);
        const std::string string8("mass");
	produces<std::vector<double> > (string8).setBranchAlias(string8);
        const std::string string9("pfDeltaCorrRelIso");
	produces<std::vector<double> > (string9).setBranchAlias(string9);
        const std::string string10("pfRhoCorrRelIso04");
	produces<std::vector<double> > (string10).setBranchAlias(string10);
        const std::string string11("pfRhoCorrRelIso03");
	produces<std::vector<double> > (string11).setBranchAlias(string11);
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
        const std::string string17("isLoose");
	produces<std::vector<bool> > (string17).setBranchAlias(string17);
        const std::string string23("isMedium");
        produces<std::vector<bool> > (string23).setBranchAlias(string23);
        const std::string string24("isTight");
        produces<std::vector<bool> > (string24).setBranchAlias(string24);

        const std::string string18("SCEnergy");
	produces<std::vector<double> > (string18).setBranchAlias(string18);
        const std::string string19("deltaEtaSCTracker");
	produces<std::vector<double> > (string19).setBranchAlias(string19);
        const std::string string20("deltaPhiSCTracker");
	produces<std::vector<double> > (string20).setBranchAlias(string20);
        const std::string string21("sigmaIetaIeta");
	produces<std::vector<double> > (string21).setBranchAlias(string21);
        const std::string string22("sigmaIphiIphi");
	produces<std::vector<double> > (string22).setBranchAlias(string22);
}


Electron::~Electron()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

float Electron::dEtaInSeed( const pat::Electron &ele ){
  return ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() ?
    ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}

// ------------ method called to produce the data  ------------
void
Electron::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	std::auto_ptr<std::vector<pat::Electron> > prodEle(new std::vector<pat::Electron>());

	std::auto_ptr< std::vector<double> > pt(new std::vector<double>);
	std::auto_ptr< std::vector<double> > eta(new std::vector<double>);
	std::auto_ptr< std::vector<double> > phi(new std::vector<double>);
	std::auto_ptr< std::vector<double> > e(new std::vector<double>);
	std::auto_ptr< std::vector<bool> > isHEEP(new std::vector<bool>);
	//std::auto_ptr< std::vector<bool> > isHEEPv50(new std::vector<bool>);
	std::auto_ptr< std::vector<int> > type(new std::vector<int>);
	std::auto_ptr< std::vector<int> > charge(new std::vector<int>);
	std::auto_ptr< std::vector<double> > mass(new std::vector<double>);
	std::auto_ptr< std::vector<double> > pfDeltaCorrRelIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > pfRhoCorrRelIso04(new std::vector<double>);
	std::auto_ptr< std::vector<double> > pfRhoCorrRelIso03(new std::vector<double>);
	std::auto_ptr< std::vector<double> > pfRelIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > photonIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > neutralHadIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > chargedHadIso(new std::vector<double>);
	std::auto_ptr< std::vector<double> > trackIso(new std::vector<double>);
	std::auto_ptr< std::vector<bool> > isLoose(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > isMedium(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > isTight(new std::vector<bool>);

	std::auto_ptr< std::vector<double> > SCEnergy(new std::vector<double>);
	std::auto_ptr< std::vector<double> > deltaEtaSCTracker(new std::vector<double>);
	std::auto_ptr< std::vector<double> > deltaPhiSCTracker(new std::vector<double>);
	std::auto_ptr< std::vector<double> > sigmaIetaIeta(new std::vector<double>);
	std::auto_ptr< std::vector<double> > sigmaIphiIphi(new std::vector<double>);

	using namespace edm;
	using namespace reco;
	using namespace pat;

	edm::Handle< edm::View<pat::Electron> > Electrons;
	iEvent.getByLabel(EleTag_,Electrons);
	edm::Handle< edm::View<reco::Vertex> > Vertices;
	iEvent.getByLabel(VertexTag_,Vertices);
	edm::Handle<double> rho_ ;
	iEvent.getByLabel(RhoTag_, rho_);
	edm::Handle<reco::ConversionCollection> ConversionsHandle_;
	iEvent.getByToken(ConversionTag_, ConversionsHandle_);
	edm::Handle<reco::BeamSpot> BeamSpotHandle_;	
	iEvent.getByToken(BeamSpotTag_, BeamSpotHandle_);


	// Get the electron ID data from the event stream.
	// Note: this implies that the VID ID modules have been run upstream.
	// If you need more info, check with the EGM group.
	//	edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
	edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
	edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
	edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
	edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
	//iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
	iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
	iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
	iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
	iEvent.getByToken(eleHEEPIdMapToken_ ,heep_id_decisions);

	if( Electrons.isValid() ) {
	  //	  std::cout<<"Electrons: "<<Electrons->size()<<std::endl;
		for(unsigned int i=0; i<Electrons->size();i++)
		{
		  if (Electrons->at(i).pt()<MinPt_) continue;

		  const auto el = Electrons->ptrAt(i);
		  double rho = *(rho_.product());
		  //		  bool isPassVeto = (*veto_id_decisions)[el];
		  bool isPassLoose = (*loose_id_decisions)[el];
		  bool isPassMedium = (*medium_id_decisions)[el];
		  bool isPassTight = (*tight_id_decisions)[el];
		  bool isPassHEEP = (*heep_id_decisions)[el];

		  prodEle->push_back(pat::Electron(Electrons->at(i)));
		  isHEEP->push_back( isPassHEEP );
		  //		  isHEEPv50->push_back( isPassLooseOLD );
		  isLoose->push_back( isPassLoose );
                  isMedium->push_back( isPassMedium );
                  isTight->push_back( isPassTight );
		  type->push_back(Electrons->at(i).pdgId());
		  charge->push_back(Electrons->at(i).charge());
		  e->push_back(Electrons->at(i).energy());
		  eta->push_back(Electrons->at(i).superCluster()->eta());
		  mass->push_back(Electrons->at(i).mass());
		  pt->push_back(Electrons->at(i).pt());
		  phi->push_back(Electrons->at(i).phi());

		  double Aeff04 = 0.5;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) < 1.0 ) Aeff04 = 0.208;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 1.0 && fabs(Electrons->at(i).superCluster()->eta()) < 1.479 ) Aeff04 = 0.209;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 1.479 && fabs(Electrons->at(i).superCluster()->eta()) < 2.0 ) Aeff04 = 0.115;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 2.0 && fabs(Electrons->at(i).superCluster()->eta()) < 2.2 ) Aeff04 = 0.143;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 2.2 && fabs(Electrons->at(i).superCluster()->eta()) < 2.3 ) Aeff04 = 0.183;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 2.3 && fabs(Electrons->at(i).superCluster()->eta()) < 2.4 ) Aeff04 = 0.194;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 2.4 ) Aeff04 = 0.261;

		  double Aeff03 = 0.5;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) < 1.0 ) Aeff03 = 0.13;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 1.0 && fabs(Electrons->at(i).superCluster()->eta()) < 1.479 ) Aeff03 = 0.14;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 1.479 && fabs(Electrons->at(i).superCluster()->eta()) < 2.0 ) Aeff03 = 0.07;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 2.0 && fabs(Electrons->at(i).superCluster()->eta()) < 2.2 ) Aeff03 = 0.09;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 2.2 && fabs(Electrons->at(i).superCluster()->eta()) < 2.3 ) Aeff03 = 0.11;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 2.3 && fabs(Electrons->at(i).superCluster()->eta()) < 2.4 ) Aeff03 = 0.11;
		  if ( fabs(Electrons->at(i).superCluster()->eta()) > 2.4 ) Aeff03 = 0.14;

		  float DeltaCorrectedIso = (Electrons->at(i).chargedHadronIso() + std::max(0., Electrons->at(i).neutralHadronIso() + Electrons->at(i).photonIso() - 0.5*Electrons->at(i).puChargedHadronIso()))/Electrons->at(i).pt();
		  float RhoCorrectedIso04 = Electrons->at(i).chargedHadronIso() + std::max(Electrons->at(i).neutralHadronIso() + Electrons->at(i).photonIso() - rho*Aeff04, 0.);
		  float RhoCorrectedIso03 = Electrons->at(i).chargedHadronIso() + std::max(Electrons->at(i).neutralHadronIso() + Electrons->at(i).photonIso() - rho*Aeff03, 0.);

		  pfDeltaCorrRelIso->push_back(DeltaCorrectedIso);
		  pfRhoCorrRelIso04->push_back(RhoCorrectedIso04);
		  pfRhoCorrRelIso03->push_back(RhoCorrectedIso03);
		  pfRelIso->push_back((Electrons->at(i).chargedHadronIso() + Electrons->at(i).neutralHadronIso()+ Electrons->at(i).photonIso())/Electrons->at(i).pt());
		  photonIso->push_back(Electrons->at(i).photonIso());
		  neutralHadIso->push_back(Electrons->at(i).neutralHadronIso());
		  chargedHadIso->push_back(Electrons->at(i).chargedHadronIso());
		  trackIso->push_back(Electrons->at(i).trackIso());		  

		  //		  std::cout<<Electrons->at(i).superCluster()->energy()<<std::endl;
		  SCEnergy->push_back(Electrons->at(i).superCluster()->energy());
		  deltaEtaSCTracker->push_back(Electrons->at(i).deltaEtaSuperClusterTrackAtVtx());
		  deltaPhiSCTracker->push_back(Electrons->at(i).deltaPhiSuperClusterTrackAtVtx());
		  sigmaIetaIeta->push_back(Electrons->at(i).sigmaIetaIeta());
		  sigmaIphiIphi->push_back(Electrons->at(i).sigmaIphiIphi());
		}
	}

	const std::string string00("");
        iEvent.put(prodEle );

	const std::string string0("pt");
	iEvent.put(pt,string0);
	const std::string string1("eta");
	iEvent.put(eta,string1);
	const std::string string2("phi");
	iEvent.put(phi,string2);
	const std::string string3("e");
	iEvent.put(e,string3);
	const std::string string4("isHEEP");
	iEvent.put(isHEEP,string4);
	//	const std::string string5("isHEEPv50");
	//	iEvent.put(isHEEPv50,string5);
	const std::string string6("type");
	iEvent.put(type,string6);
	const std::string string7("charge");
	iEvent.put(charge,string7);
	const std::string string8("mass");
	iEvent.put(mass,string8);
	const std::string string9("pfDeltaCorrRelIso");
	iEvent.put(pfDeltaCorrRelIso,string9);
	const std::string string10("pfRhoCorrRelIso04");
	iEvent.put(pfRhoCorrRelIso04,string10);
	const std::string string11("pfRhoCorrRelIso03");
	iEvent.put(pfRhoCorrRelIso03,string11);
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
	const std::string string17("isLoose");
	iEvent.put(isLoose,string17);
        const std::string string23("isMedium");
        iEvent.put(isMedium,string23);
        const std::string string24("isTight");
        iEvent.put(isTight,string24);

	const std::string string18("SCEnergy");
	iEvent.put(SCEnergy,string18);
	const std::string string19("deltaEtaSCTracker");
	iEvent.put(deltaEtaSCTracker,string19);
	const std::string string20("deltaPhiSCTracker");
	iEvent.put(deltaPhiSCTracker,string20);
	const std::string string21("sigmaIetaIeta");
	iEvent.put(sigmaIetaIeta,string21);
	const std::string string22("sigmaIphiIphi");
	iEvent.put(sigmaIphiIphi,string22);
}


// ------------ method called once each job just before starting event loop  ------------
void 
Electron::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Electron::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
Electron::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Electron::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Electron::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Electron::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Electron::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Electron);
