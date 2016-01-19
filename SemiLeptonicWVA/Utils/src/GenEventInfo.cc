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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

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
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > > PUInfoToken_;
  edm::LumiReWeighting  LumiWeights_;
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
  geneventToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("geneventToken"))),
  PUInfoToken_(consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("PUInfoToken")))
{

	produces<double>("genEventWeight");
	produces<double>("PUWeight");
	produces<int>("npT");

	float PU_data_f[52] = {
	  0 ,
	  2.78298e-05 ,
	  0.000205213 ,
	  0.000278167 ,
	  0.000252931 ,
	  0.000398665 ,
	  0.000527857 ,
	  0.000885768 ,
	  0.00251079 ,
	  0.0127222 ,
	  0.0448683 ,
	  0.0894268 ,
	  0.132441 ,
	  0.163719 ,
	  0.171203 ,
	  0.150137 ,
	  0.109292 ,
	  0.0659582 ,
	  0.0332417 ,
	  0.0141861 ,
	  0.00523167 ,
	  0.00171762 ,
	  0.00052488 ,
	  0.00015846 ,
	  5.05184e-05 ,
	  1.81421e-05 ,
	  7.67416e-06 ,
	  3.78755e-06 ,
	  2.05074e-06 ,
	  1.14134e-06 ,
	  6.27311e-07 ,
	  3.34314e-07 ,
	  1.71548e-07 ,
	  8.4558e-08 ,
	  4.00085e-08 ,
	  1.81674e-08 ,
	  7.9169e-09 ,
	  3.31079e-09 ,
	  1.32868e-09 ,
	  5.11706e-10 ,
	  1.89117e-10 ,
	  6.70731e-11 ,
	  2.28284e-11 ,
	  7.45611e-12 ,
	  2.337e-12 ,
	  7.02939e-13 ,
	  2.02903e-13 ,
	  5.62046e-14 ,
	  1.49408e-14 ,
	  3.81146e-15 ,
	  9.33102e-16 ,
	  2.19227e-16 
	};

	float PU_MC_f[52] {
	  4.8551E-07,
            1.74806E-06,
            3.30868E-06,
            1.62972E-05,
            4.95667E-05,
            0.000606966,
            0.003307249,
            0.010340741,
            0.022852296,
            0.041948781,
            0.058609363,
            0.067475755,
            0.072817826,
            0.075931405,
            0.076782504,
            0.076202319,
            0.074502547,
            0.072355135,
            0.069642102,
            0.064920999,
            0.05725576,
            0.047289348,
            0.036528446,
            0.026376131,
            0.017806872,
            0.011249422,
            0.006643385,
            0.003662904,
            0.001899681,
            0.00095614,
            0.00050028,
            0.000297353,
            0.000208717,
            0.000165856,
            0.000139974,
            0.000120481,
            0.000103826,
            8.88868E-05,
            7.53323E-05,
            6.30863E-05,
            5.21356E-05,
            4.24754E-05,
            3.40876E-05,
            2.69282E-05,
            2.09267E-05,
            1.5989E-05,
            4.8551E-06,
            2.42755E-06,
            4.8551E-07,
            2.42755E-07,
            1.21378E-07,
            4.8551E-08
	    };

	std::vector<float> MC_dist;
	
	for (unsigned int iPU = 0; iPU < 52; iPU++){
	  MC_dist.push_back(PU_MC_f[iPU]);
	}

	std::vector<float> data_dist; 

	for (unsigned int iPU = 0; iPU < 52; iPU++){
	  data_dist.push_back(PU_data_f[iPU]);
	}

	LumiWeights_ = edm::LumiReWeighting(MC_dist,data_dist);
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

	edm::Handle<std::vector< PileupSummaryInfo > >  PUInfo;
	iEvent.getByToken(PUInfoToken_, PUInfo);

	double PUWeight=-1.;
	double eventWeight =1.;
	float Tnpv = -1;

	if (geneventInfo_.isValid()) {
	  eventWeight = geneventInfo_->weight();

	  std::vector<PileupSummaryInfo>::const_iterator PVI;
	  for(PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI) {
	    int BX = PVI->getBunchCrossing();
	    if(BX == 0) { 
	      Tnpv = PVI->getTrueNumInteractions();
	      continue;
	    }
	  }
	  PUWeight = LumiWeights_.weight( Tnpv );
	}

	std::auto_ptr<double> htpw(new double(PUWeight));
	iEvent.put(htpw,"PUWeight");

	std::auto_ptr<int> htpu(new int(Tnpv));
	iEvent.put(htpu,"npT");

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
