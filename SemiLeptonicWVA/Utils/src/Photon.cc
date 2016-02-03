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
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

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

        void printCutFlowResult(vid::CutFlowResult &cutflow);

  // ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> > phoLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > phoTightIdFullInfoMapToken_;

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
  phoTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoTightIdMap"))),
  phoTightIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >
                               (iConfig.getParameter<edm::InputTag>("phoTightIdFullInfoMap")))
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

        const std::string string7("minPt");
        produces<std::vector<double> > (string7).setBranchAlias(string7);
        const std::string string8("phoSCEtaMultiRange");
        produces<std::vector<double> > (string8).setBranchAlias(string8);
        const std::string string9("phoSingleTowerHadOverEm");
        produces<std::vector<double> > (string9).setBranchAlias(string9);
        const std::string string10("phoFull5x5SigmaIEtaIEta");
        produces<std::vector<double> > (string10).setBranchAlias(string10);
        const std::string string11("phoAnyPFIsoWithEA");
        produces<std::vector<double> > (string11).setBranchAlias(string11);
        const std::string string12("phoAnyPFIsoWithEAAndExpoScaling");
        produces<std::vector<double> > (string12).setBranchAlias(string12);
        const std::string string13("phoAnyPFIsoWithEA1");
        produces<std::vector<double> > (string13).setBranchAlias(string13);
        const std::string string21("PassMinPt");
        produces<std::vector<bool> > (string21).setBranchAlias(string21);
        const std::string string22("PassPhoSCEtaMultiRange");
        produces<std::vector<bool> > (string22).setBranchAlias(string22);
        const std::string string23("PassPhoSingleTowerHadOverEm");
        produces<std::vector<bool> > (string23).setBranchAlias(string23);
        const std::string string24("PassPhoFull5x5SigmaIEtaIEta");
        produces<std::vector<bool> > (string24).setBranchAlias(string24);
        const std::string string25("PassPhoAnyPFIsoWithEA");
        produces<std::vector<bool> > (string25).setBranchAlias(string25);
        const std::string string26("PassPhoAnyPFIsoWithEAAndExpoScaling");
        produces<std::vector<bool> > (string26).setBranchAlias(string26);
        const std::string string27("PassPhoAnyPFIsoWithEA1");
        produces<std::vector<bool> > (string27).setBranchAlias(string27);

        const std::string string14("hasPixelSeed");
        produces<std::vector<bool> > (string14).setBranchAlias(string14);
        const std::string string15("passElectronVeto");
        produces<std::vector<bool> > (string15).setBranchAlias(string15);
        const std::string string16("photonIso");
        produces<std::vector<double> > (string16).setBranchAlias(string16);
        const std::string string17("neutralHadIso");
        produces<std::vector<double> > (string17).setBranchAlias(string17);
        const std::string string18("chargedHadIso");
        produces<std::vector<double> > (string18).setBranchAlias(string18);
        const std::string string19("puChargedHadIso");
        produces<std::vector<double> > (string19).setBranchAlias(string19);
        const std::string string20("sigmaIetaIeta");
        produces<std::vector<double> > (string20).setBranchAlias(string20);

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

        std::auto_ptr< std::vector<double> > minPt(new std::vector<double>);
        std::auto_ptr< std::vector<double> > phoSCEtaMultiRange(new std::vector<double>);
        std::auto_ptr< std::vector<double> > phoSingleTowerHadOverEm(new std::vector<double>);
        std::auto_ptr< std::vector<double> > phoFull5x5SigmaIEtaIEta(new std::vector<double>);
        std::auto_ptr< std::vector<double> > phoAnyPFIsoWithEA(new std::vector<double>);
        std::auto_ptr< std::vector<double> > phoAnyPFIsoWithEAAndExpoScaling(new std::vector<double>);
        std::auto_ptr< std::vector<double> > phoAnyPFIsoWithEA1(new std::vector<double>);
        std::auto_ptr< std::vector<bool> > PassMinPt(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > PassPhoSCEtaMultiRange(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > PassPhoSingleTowerHadOverEm(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > PassPhoFull5x5SigmaIEtaIEta(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > PassPhoAnyPFIsoWithEA(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > PassPhoAnyPFIsoWithEAAndExpoScaling(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > PassPhoAnyPFIsoWithEA1(new std::vector<bool>);

        std::auto_ptr< std::vector<bool> > passElectronVeto(new std::vector<bool>);
        std::auto_ptr< std::vector<bool> > hasPixelSeed(new std::vector<bool>);
        std::auto_ptr< std::vector<double> > photonIso(new std::vector<double>);
        std::auto_ptr< std::vector<double> > neutralHadIso(new std::vector<double>);
        std::auto_ptr< std::vector<double> > chargedHadIso(new std::vector<double>);
        std::auto_ptr< std::vector<double> > puChargedHadIso(new std::vector<double>);
        std::auto_ptr< std::vector<double> > sigmaIetaIeta(new std::vector<double>);

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
        edm::Handle<edm::ValueMap<vid::CutFlowResult> > tight_id_cutflow_data;
	iEvent.getByToken(phoLooseIdMapToken_ ,loose_id_decisions);
	iEvent.getByToken(phoMediumIdMapToken_,medium_id_decisions);
	iEvent.getByToken(phoTightIdMapToken_,tight_id_decisions);
        iEvent.getByToken(phoTightIdFullInfoMapToken_, tight_id_cutflow_data);

	if( Photons.isValid() ) {
		for(unsigned int i=0; i<Photons->size();i++)
		{
		  const auto ph = Photons->ptrAt(i);
		  bool isPassLoose = (*loose_id_decisions)[ph];
		  bool isPassMedium = (*medium_id_decisions)[ph];
		  bool isPassTight = (*tight_id_decisions)[ph];

		  bool PassElectronVeto = (pat::Photon(Photons->at(i))).passElectronVeto();
		  bool HasPixelSeed = (pat::Photon(Photons->at(i))).hasPixelSeed();

                  // Index of cut to mask (SigmaIetaIeta):
                  const int cutIndexToMask = 3;//PhoFull5x5SigmaIEtaIEtaCut_0
                  vid::CutFlowResult fullCutFlowData = (*tight_id_cutflow_data)[ph];
                  vid::CutFlowResult maskedCutFlowData = fullCutFlowData.getCutFlowResultMasking(cutIndexToMask);
                  isPassTight = maskedCutFlowData.cutFlowPassed();
                  //printf("Cut masked: %s\n", maskedCutFlowData.getNameAtIndex(cutIndexToMask).c_str());
                  //printCutFlowResult(fullCutFlowData);

                  minPt->push_back(fullCutFlowData.getValueCutUpon(0));
                  phoSCEtaMultiRange->push_back(fullCutFlowData.getValueCutUpon(1));
                  phoSingleTowerHadOverEm->push_back(fullCutFlowData.getValueCutUpon(2));
                  phoFull5x5SigmaIEtaIEta->push_back(fullCutFlowData.getValueCutUpon(3));
                  phoAnyPFIsoWithEA->push_back(fullCutFlowData.getValueCutUpon(4));
                  phoAnyPFIsoWithEAAndExpoScaling->push_back(fullCutFlowData.getValueCutUpon(5));
                  phoAnyPFIsoWithEA1->push_back(fullCutFlowData.getValueCutUpon(6));
                  PassMinPt->push_back(fullCutFlowData.getCutResultByIndex(0));
                  PassPhoSCEtaMultiRange->push_back(fullCutFlowData.getCutResultByIndex(1));
                  PassPhoSingleTowerHadOverEm->push_back(fullCutFlowData.getCutResultByIndex(2));
                  PassPhoFull5x5SigmaIEtaIEta->push_back(fullCutFlowData.getCutResultByIndex(3));
                  PassPhoAnyPFIsoWithEA->push_back(fullCutFlowData.getCutResultByIndex(4));
                  PassPhoAnyPFIsoWithEAAndExpoScaling->push_back(fullCutFlowData.getCutResultByIndex(5));
                  PassPhoAnyPFIsoWithEA1->push_back(fullCutFlowData.getCutResultByIndex(6));

                  passElectronVeto->push_back(PassElectronVeto);
                  hasPixelSeed->push_back(HasPixelSeed);
                  photonIso->push_back(Photons->at(i).photonIso());
                  neutralHadIso->push_back(Photons->at(i).neutralHadronIso());
                  chargedHadIso->push_back(Photons->at(i).chargedHadronIso());
                  puChargedHadIso->push_back(Photons->at(i).puChargedHadronIso());
                  sigmaIetaIeta->push_back(Photons->at(i).sigmaIetaIeta());

                  //std::cout << "Resulting table ... " << std::endl;
                  //printCutFlowResult(maskedCutFlowData);
                  //std::cout << "Result: " << isPassTight << std::endl;

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

        const std::string string7("minPt");
        iEvent.put(minPt,string7);
        const std::string string8("phoSCEtaMultiRange");
        iEvent.put(phoSCEtaMultiRange,string8);
        const std::string string9("phoSingleTowerHadOverEm");
        iEvent.put(phoSingleTowerHadOverEm,string9);
        const std::string string10("phoFull5x5SigmaIEtaIEta");
        iEvent.put(phoFull5x5SigmaIEtaIEta,string10);
        const std::string string11("phoAnyPFIsoWithEA");
        iEvent.put(phoAnyPFIsoWithEA,string11);
        const std::string string12("phoAnyPFIsoWithEAAndExpoScaling");
        iEvent.put(phoAnyPFIsoWithEAAndExpoScaling,string12);
        const std::string string13("phoAnyPFIsoWithEA1");
        iEvent.put(phoAnyPFIsoWithEA1,string13);
        const std::string string21("PassMinPt");
        iEvent.put(PassMinPt,string21);
        const std::string string22("PassPhoSCEtaMultiRange");
        iEvent.put(PassPhoSCEtaMultiRange,string22);
        const std::string string23("PassPhoSingleTowerHadOverEm");
        iEvent.put(PassPhoSingleTowerHadOverEm,string23);
        const std::string string24("PassPhoFull5x5SigmaIEtaIEta");
        iEvent.put(PassPhoFull5x5SigmaIEtaIEta,string24);
        const std::string string25("PassPhoAnyPFIsoWithEA");
        iEvent.put(PassPhoAnyPFIsoWithEA,string25);
        const std::string string26("PassPhoAnyPFIsoWithEAAndExpoScaling");
        iEvent.put(PassPhoAnyPFIsoWithEAAndExpoScaling,string26);
        const std::string string27("PassPhoAnyPFIsoWithEA1");
        iEvent.put(PassPhoAnyPFIsoWithEA1,string27);

        const std::string string14("hasPixelSeed");
        iEvent.put(hasPixelSeed,string14);
        const std::string string15("passElectronVeto");
        iEvent.put(passElectronVeto,string15);
        const std::string string16("photonIso");
        iEvent.put(photonIso,string16);
        const std::string string17("neutralHadIso");
        iEvent.put(neutralHadIso,string17);
        const std::string string18("chargedHadIso");
        iEvent.put(chargedHadIso,string18);
        const std::string string19("puChargedHadIso");
        iEvent.put(puChargedHadIso,string19);
        const std::string string20("sigmaIetaIeta");
        iEvent.put(sigmaIetaIeta,string20);

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

// ------------ method to print Photon ID cut flow results --------------
void Photon::printCutFlowResult(vid::CutFlowResult &cutflow){

  printf("    CutFlow name= %s    decision is %d\n", 
	 cutflow.cutFlowName().c_str(),
	 (int) cutflow.cutFlowPassed());
  int ncuts = cutflow.cutFlowSize();
  printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
  for(int icut = 0; icut<ncuts; icut++){
    printf("  %d       %50s    %d        %f          %d\n", icut,
	   cutflow.getNameAtIndex(icut).c_str(),
	   (int)cutflow.isCutMasked(icut),
	   cutflow.getValueCutUpon(icut),
	   (int)cutflow.getCutResultByIndex(icut));
  }
  printf("    WARNING: the value-cut-upon is bugged in 7.4.7, it is always 1.0\n");

}

//define this as a plug-in
DEFINE_FWK_MODULE(Photon);
