//-*- C++ -*-
//
// Package: LeptonFilter
// Class: LeptonFilter
//
//
// Original Author: Luca Brianza
// Created: Fri Jun 6 11:07:38 CDT 2015
// 
//
//
// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;


class LeptonFilter : public edm::EDFilter {
public:
  explicit LeptonFilter (const edm::ParameterSet&);
  ~LeptonFilter();
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  // ----------member data ---------------------------
  edm::InputTag electronsInputTag_;
  edm::InputTag muonsInputTag_;
  double eleFilterPtCut_;
  double muFilterPtCut_;
};

//
// constructors and destructor
//
LeptonFilter::LeptonFilter(const edm::ParameterSet& iConfig) {
  electronsInputTag_ = iConfig.getParameter<edm::InputTag> ("electronsInputTag" );
  muonsInputTag_ = iConfig.getParameter<edm::InputTag> ("muonsInputTag" );
  eleFilterPtCut_ = iConfig.getParameter<double> ("eleFilterPtCut" );
  muFilterPtCut_ = iConfig.getParameter<double> ("muFilterPtCut" );
}

LeptonFilter::~LeptonFilter() {}

void LeptonFilter::beginJob() {}

void LeptonFilter::endJob() {}

// ------------ method called to produce the data ------------

bool LeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  edm::Handle< edm::View<pat::Electron> > els_h;
  iEvent.getByLabel(electronsInputTag_,els_h);
  edm::Handle< edm::View<pat::Muon> > muon_h;
  iEvent.getByLabel(muonsInputTag_, muon_h);

  //edm::View<reco::Muon>::const_iterator muons_end = muon_h->end();
  if( els_h.isValid() ) {
    //      std::cout<<"Electrons: "<<Electrons->size()<<std::endl;                                                                                         
    for(unsigned int i=0; i<els_h->size();i++) {
      if (els_h->at(i).pt() > eleFilterPtCut_) {
	return true;
      }
    }
  }
  if( muon_h.isValid() ) {
    //      std::cout<<"Electrons: "<<Electrons->size()<<std::endl;                                                                                         
    for(unsigned int i=0; i<muon_h->size();i++) {
      if (muon_h->at(i).pt() > muFilterPtCut_) {
	return true;
      }
    }
  }

  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonFilter);
