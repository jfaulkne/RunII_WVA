#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"

#include "../interface/setInputTree.h"
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"

#include "../interface/SF.h"

using namespace std;

//*****PU WEIGHT***************

vector<double> generate_weights(TH1* data_npu_estimated, int isForSynch){
  // see SimGeneral/MixingModule/python/mix_2015_25ns_Startup_PoissonOOTPU_cfi.pyy; copy and paste from there:
  const double npu_probs[52] = {
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
  
  if (isForSynch==0) { //OFFICIAL RECIPE
    vector<double> result(52);
    double s = 0.0;
    for(int npu=0; npu<52; ++npu){
      double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
      result[npu] = npu_estimated / npu_probs[npu];
      s += npu_estimated;
    }
    // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    for(int npu=0; npu<52; ++npu){
      result[npu] /= s;
    }
    return result;
  }

  else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
    vector<double> result(60);
    for(int npu=0; npu<60; ++npu){
      if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==NULL)
	result[npu] = 0.;
      else {
	double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
	result[npu] = npu_estimated;
      }
    }
    return result;
  }

}


//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  int isMC = std::atof(argv[3]);
  std::string leptonName = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string numberOfEntries = argv[8];
  int isAQGC = std::atof(argv[9]);
  std::string firstEvent = argv[10];
  int applyPUJetID = std::atof(argv[11]);

  float weight = std::atof(xSecWeight.c_str())/std::atof(numberOfEntries.c_str());
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }

  TLorentzVector W,LEP;
  TLorentzVector NU0,NU1,NU2;
  TLorentzVector JET, HADW, AK4;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector ELE,MU,PH;

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;
  std::vector<TLorentzVector> tightPh;

  std::cout<<"file: "<<(inputFolder+inputFile).c_str()<<std::endl;
  TFile *MyFile = TFile::Open((inputFolder+inputFile).c_str());
  setInputTree *ReducedTree = new setInputTree (MyFile, inputTreeName.c_str());
  if (ReducedTree->fChain == 0) return (-1);
  ReducedTree->Init();

  int n_cuts = 34;
  bool cut_flag[3][n_cuts];
  int cutEff[3][n_cuts]; for (int i = 0; i != 3; i++) for (int j = 0; j != n_cuts; j++) cutEff[i][j] = 0;

  //--------pile up file -----------------
  std::vector<double> weights_pu1; //these are made with our recipe

  TFile* pileupFile1 = TFile::Open("MyDataPileupHistogram_69mb.root");  
  TH1F *pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
  weights_pu1 = generate_weights(pileupHisto1,0);
  pileupFile1->Close();

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_")+leptonName+std::string("/")+outputFile).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);

  setOutputTree *WWTree = new setOutputTree(outTree);

  //---------start loop on events------------
  Long64_t max_evt = std::atof(numberOfEntries.c_str());
  if (max_evt == -1) max_evt = ReducedTree->fChain->GetEntries();
  if (ReducedTree->fChain->GetEntries() < max_evt) max_evt = ReducedTree->fChain->GetEntries();
  for (Long64_t jentry = 0; jentry != max_evt; jentry++) {

    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) {cout << "Error loading event " << jentry+1 << endl; break;}

    int nb = ReducedTree->fChain->GetEntry(jentry);   
    if (jentry < std::atof(firstEvent.c_str())-1) continue;

    tightPh.clear();
    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();

    if(iEntry % 25000 == 0)    
      cout << "read entry: " << iEntry << endl;

    for (int i = 0; i != 3; i++) for (int j = 0; j != n_cuts; j++) cut_flag[i][j] = false;
    WWTree->initializeVariables(); //initialize all variables
    
    if (ReducedTree->passFilterHBHELooseRerun == 0)	continue;
    for (int i = 0; i != 3; i++) {cut_flag[i][0] = true; cutEff[i][0]++;}

    if (ReducedTree->passFilterHBHEIsoRerun == 0)	continue;
    for (int i = 0; i != 3; i++) {cut_flag[i][1] = true; cutEff[i][1]++;}

    if (ReducedTree->passFilterCSCHalo == 0)		continue;
    for (int i = 0; i != 3; i++) {cut_flag[i][2] = true; cutEff[i][2]++;}

    if (ReducedTree->passFilterGoodVtx == 0)		continue;
    for (int i = 0; i != 3; i++) {cut_flag[i][3] = true; cutEff[i][3]++;}

    if (ReducedTree->passFilterEEBadSC == 0)		continue;
    for (int i = 0; i != 3; i++) {cut_flag[i][4] = true; cutEff[i][4]++;}

    WWTree->wSampleWeight = weight; //xsec/numberOfEntries
    WWTree->totalEventWeight = 1.; //temporary value
    WWTree->eff_and_pu_Weight = 1.; //temporary value
    WWTree->totalEventWeight_2 = 1.; //temporary value
    WWTree->eff_and_pu_Weight_2 = 1.; //temporary value
    
    if (ReducedTree->genEventWeight>0) WWTree->genWeight=1.;
    else if (ReducedTree->genEventWeight<0) WWTree->genWeight=-1.;

    //PILE-UP WEIGHT
    if (isMC) {
      if(ReducedTree->NVtx < weights_pu1.size()){
	WWTree->eff_and_pu_Weight = weights_pu1[ReducedTree->npT]; //official pu recipe
	WWTree->totalEventWeight *= weights_pu1[ReducedTree->npT];
      }
      else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	std::cout<<"Warning! n_pu too big"<<std::endl;
	//	throw logic_error("n_pu too big");
        WWTree->eff_and_pu_Weight = 0.;
	WWTree->totalEventWeight *= 0.;
      }    

      WWTree->eff_and_pu_Weight_2 = ReducedTree->PUWeight; //our pu recipe
      WWTree->totalEventWeight_2 *= ReducedTree->PUWeight;
    }    

    //save event variables
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi = ReducedTree->LumiBlockNum;
    WWTree->nPV  = ReducedTree->NVtx;
    WWTree->nPho = ReducedTree->PhotonsNum; 

    /////////////////THE SELECTED LEPTON
    int nTightLepton=0;
    if (strcmp(leptonName.c_str(),"el")==0) {

      if (ReducedTree->TriggerProducerTriggerPass->at(0)==0) continue; //trigger
      for (int cut = 0; cut != 3; cut++) {cut_flag[cut][5] = true; cutEff[cut][5]++;}

      float tempPt=0.;
      for (int i=0; i<ReducedTree->ElectronsNum; i++) {

	if (ReducedTree->ElectronsPt[i]<tempPt) continue;
        if (!ReducedTree->Electrons_isTight[i]) continue;
        if (!cut_flag[0][7]) for (int cut = 0; cut != 3; cut++) {cutEff[cut][7]++; cut_flag[cut][7] = true;}

        if (ReducedTree->ElectronsPt[i]<=25) continue;
        if (!cut_flag[0][8]) for (int cut = 0; cut != 3; cut++) {cutEff[cut][8]++; cut_flag[cut][8] = true;}

	ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
	if ((fabs(ELE.Eta())>=1.44 && fabs(ELE.Eta())<=1.57) || fabs(ELE.Eta())>=2.5) continue;
        if (!cut_flag[0][9]) for (int cut = 0; cut != 3; cut++) {cutEff[cut][9]++; cut_flag[cut][9] = true;}

	tightEle.push_back(ELE);
	WWTree->l_pt  = ReducedTree->ElectronsPt[i];
	WWTree->l_eta = ReducedTree->ElectronsEta[i];
	WWTree->l_phi = ReducedTree->ElectronsPhi[i];	
	WWTree->l_e= ReducedTree->ElectronsE[i];	
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    }
    else if (strcmp(leptonName.c_str(),"mu")==0) {

      if (ReducedTree->TriggerProducerTriggerPass->at(0)==0 && ReducedTree->TriggerProducerTriggerPass->at(1)==0) continue; //trigger
      for (int cut = 0; cut != 3; cut++) {cut_flag[cut][5] = true; cutEff[cut][5]++;}

      float tempPt=0.;
      for (int i=0; i<ReducedTree->MuonsNum; i++) {

	if (ReducedTree->MuonsPt[i]<tempPt) continue;
	if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
        if (!cut_flag[0][6]) for (int cut = 0; cut != 3; cut++) {cutEff[cut][6]++; cut_flag[cut][6] = true;}

        if (!ReducedTree->Muons_isTight[i]) continue;
        if (!cut_flag[0][7]) for (int cut = 0; cut != 3; cut++) {cutEff[cut][7]++; cut_flag[cut][7] = true;}

        if (ReducedTree->MuonsPt[i]<=23) continue;
        if (!cut_flag[0][8]) for (int cut = 0; cut != 3; cut++) {cutEff[cut][8]++; cut_flag[cut][8] = true;}

        if (fabs(ReducedTree->MuonsEta[i])>=2.1) continue;
        if (!cut_flag[0][9]) for (int cut = 0; cut != 3; cut++) {cutEff[cut][9]++; cut_flag[cut][9] = true;}

	MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);

	tightMuon.push_back(MU);
	WWTree->l_pt  = ReducedTree->MuonsPt[i];
	WWTree->l_eta = ReducedTree->MuonsEta[i];
	WWTree->l_phi = ReducedTree->MuonsPhi[i];
	WWTree->l_e = ReducedTree->MuonsE[i];
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    } else {cout << "Invalid lepton channel " << leptonName.c_str() << " requested" << endl; break;}

    if (nTightLepton==0) continue; //no leptons with required ID
    for (int cut = 0; cut != 3; cut++) {cutEff[cut][10]++; cut_flag[cut][10] = true;}

    //VETO ADDITIONAL LEPTONS
    int nLooseLepton=0;
    for (int i=0; i<ReducedTree->ElectronsNum; i++) {
      if (!ReducedTree->Electrons_isLoose[i]) continue;
      if (ReducedTree->ElectronsPt[i]<20) continue;       
      ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
      if ((fabs(ELE.Eta())>=1.44 && fabs(ELE.Eta())<=1.57) || fabs(ELE.Eta())>=2.5) continue;
      looseEle.push_back(ELE);      
      nLooseLepton++;
    }
    for (int i=0; i<ReducedTree->MuonsNum; i++) {
      if (!ReducedTree->Muons_isLoose[i]) continue;
      if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
      if (fabs(ReducedTree->MuonsEta[i])>=2.1) continue;
      if (ReducedTree->MuonsPt[i]<10) continue;
      MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
      looseMuon.push_back(MU);
      nLooseLepton++;
    }
    if (nLooseLepton!=1) continue; //no additional leptons
    for (int cut = 0; cut != 3; cut++) {cutEff[cut][11]++; cut_flag[cut][11] = true;}

    //preselection on jet pt and met
    if (ReducedTree->METPt <= 35) continue; 
    for (int cut = 0; cut != 3; cut++) {cutEff[cut][12]++; cut_flag[cut][12] = true;}

    //////////////THE MET

    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*

    float Wmass = 80.385;

    TLorentzVector W_mu, W_Met;

    W_mu.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    W_Met.SetPxPyPzE(ReducedTree->METPt * TMath::Cos(ReducedTree->METPhi), ReducedTree->METPt * TMath::Sin(ReducedTree->METPhi), 0., sqrt(ReducedTree->METPt*ReducedTree->METPt));

    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    METzCalculator_Run2 NeutrinoPz_run2;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_mu);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());

    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(W_mu);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());

    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    double pz1_run2 = NeutrinoPz_run2.Calculate(); 

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }

    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(W_mu);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());

    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2; 
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    if (NeutrinoPz_type2.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }

    WWTree->pfMET   = sqrt(ReducedTree->METPt*ReducedTree->METPt);
    WWTree->pfMET_Phi = ReducedTree->METPhi;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();


    /////////////////THE LEPTONIC W
    
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    NU0.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type0*WWTree->nu_pz_type0));
    NU2.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_run2*WWTree->nu_pz_run2));
    W = LEP + NU2;
    
    WWTree->v_pt = W.Pt();
    WWTree->v_eta = W.Eta();
    WWTree->v_phi = W.Phi();
    WWTree->v_mt = TMath::Sqrt(2*LEP.Et()*NU2.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2))));

    //FOR THE SYNCHORNIZATION!!! REMOVE IT FOR THE REAL ANALYSIS!!!!
    //    NU2.SetPtEtaPhiE(ReducedTree->METPt,0.,ReducedTree->METPhi,0.);
    //    W = NU2+LEP; 
    ////

    if (WWTree->v_mt<=30) continue;
    for (int cut = 0; cut != 3; cut++) {cutEff[cut][13]++; cut_flag[cut][13] = true;}

    /////////VBF and b-tag section
    bool fillVBF = true;

    std::vector<int> tmpIndexGoodJets;

    WWTree->njets=0;
    WWTree->nBTagJet_loose=0;
    WWTree->nBTagJet_medium=0;
    WWTree->nBTagJet_tight=0;

    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
      {

	if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
        for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][14]) {cutEff[cut][14]++; cut_flag[cut][14] = true;}

	if (ReducedTree->Jets_PtCorr[i]<=30 || ReducedTree->JetsPt[i]<=20) continue;
        for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][15]) {cutEff[cut][15]++; cut_flag[cut][15] = true;}

        if (fabs(ReducedTree->JetsEta[i])>=2.4)  continue;
        for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][16]) {cutEff[cut][16]++; cut_flag[cut][16] = true;}

	//fill B-Tag info
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.423)   WWTree->nBTagJet_loose++;
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.814)   WWTree->nBTagJet_medium++;
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.941)   WWTree->nBTagJet_tight++;

        if (ReducedTree->Jets_bDiscriminatorCSV[i]>=0.679)   continue;
        for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][17]) {cutEff[cut][17]++; cut_flag[cut][17] = true;}

        //Clean from PU (medium working point, 95%):
        if (applyPUJetID) {
        if (fabs(ReducedTree->JetsEta[i]) < 2. && ReducedTree->Jets_PUJetID[i] >= -0.483) continue;
        else if (fabs(ReducedTree->JetsEta[i]) >= 2. && fabs(ReducedTree->JetsEta[i]) < 2.5 && ReducedTree->Jets_PUJetID[i] >= -0.663) continue;
        else if (fabs(ReducedTree->JetsEta[i]) >= 2.5 && fabs(ReducedTree->JetsEta[i]) < 3. && ReducedTree->Jets_PUJetID[i] >= -0.445) continue;
        else if (fabs(ReducedTree->JetsEta[i]) >= 3. && fabs(ReducedTree->JetsEta[i]) < 5. && ReducedTree->Jets_PUJetID[i] >= -0.296) continue;
        }

	//CLEANING FROM LEPTON
	if (deltaR(WWTree->l_eta, WWTree->l_phi, ReducedTree->JetsEta[i], ReducedTree->JetsPhi[i]) <= 0.3) continue;
        for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][18]) {cutEff[cut][18]++; cut_flag[cut][18] = true;}

	//CLEANING FROM MET 
        if (fabs(WWTree->pfMET_Phi-ReducedTree->JetsPhi[i]) <= 0.4) continue;
        for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][19]) {cutEff[cut][19]++; cut_flag[cut][19] = true;}

	WWTree->njets++;
	tmpIndexGoodJets.push_back(i); //save index of the "good" vbf jets candidate

      }

    int JetCNT = 0;
    for (std::vector<int>::iterator i = tmpIndexGoodJets.begin(); i != tmpIndexGoodJets.end(); ++i)
       if (ReducedTree->Jets_PtCorr[*i] > 30) JetCNT++;
    if (JetCNT == 2) {cut_flag[0][20] = cut_flag[1][20] = true; cutEff[0][20]++; cutEff[1][20]++;}
    else tmpIndexGoodJets.clear();

    if (tmpIndexGoodJets.size()<2) fillVBF = false; //check if at least 2 jets are inside the collection

    // Since 'on-the-fly' JEC updates were applied to the jet collection, they are no longer pT-ordered
    // Therefore, loop over 'good' jets and reorder by pT
    std::vector<int> indexGoodJets;
    while (indexGoodJets.size() != tmpIndexGoodJets.size()) {
       double jet_pt_tmp = 0;
       int jet_i_tmp = 0;
       for (unsigned int i = 0; i != tmpIndexGoodJets.size(); i++) {
         if (tmpIndexGoodJets.at(i) >= 0)
           if (ReducedTree->Jets_PtCorr[tmpIndexGoodJets.at(i)] > jet_pt_tmp)
             {jet_i_tmp = i; jet_pt_tmp = ReducedTree->Jets_PtCorr[tmpIndexGoodJets.at(i)];}
       }
       indexGoodJets.push_back(tmpIndexGoodJets.at(jet_i_tmp));
       tmpIndexGoodJets.at(jet_i_tmp) = -1;
    }
    std::vector<int> indexGoodJetsAK8, tmpIndexGoodJetsAK8;
    if (ReducedTree->AK8JetsNum > 0)
    while (tmpIndexGoodJetsAK8.size() != ReducedTree->AK8JetsNum) {
       double jet_pt_tmp = 0;
       int jet_i_tmp = 0;
       for (int i = 0; i != ReducedTree->AK8JetsNum; i++) {
         bool skip = false;
         for (unsigned int j = 0; j != tmpIndexGoodJetsAK8.size(); j++) if (i == tmpIndexGoodJetsAK8.at(j)) skip = true;
         if (!skip && ReducedTree->AK8Jets_PtCorr[i] > jet_pt_tmp)
             {jet_i_tmp = i; jet_pt_tmp = ReducedTree->AK8Jets_PtCorr[jet_i_tmp];}
       }
       tmpIndexGoodJetsAK8.push_back(jet_i_tmp);
    }

    std::vector<int> Jet1, Jet2;
    if (fillVBF) {
       for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][21]) {cutEff[cut][21]++; cut_flag[cut][21] = true;}
       for (unsigned int i=0; i<indexGoodJets.size()-1; i++) {
          for (unsigned int ii=i+1; ii<indexGoodJets.size(); ii++) {

	       VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->Jets_ECorr[indexGoodJets.at(i)]);
	       VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(ii)],ReducedTree->JetsEta[indexGoodJets.at(ii)],ReducedTree->JetsPhi[indexGoodJets.at(ii)],ReducedTree->Jets_ECorr[indexGoodJets.at(ii)]);

               if (fabs(VBF1.Eta()-VBF2.Eta()) >= 1.4) continue;
               for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][22])
                  {cutEff[cut][22]++; cut_flag[cut][22] = true;}

               if ((VBF1+VBF2).Pt() < 70) continue;
               for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][23])
                  {cutEff[cut][23]++; cut_flag[cut][23] = true;}

	       TOT = VBF1 + VBF2;
	       if (TOT.M() > 70 && TOT.M() < 100) continue; // PLEASE CHECK THIS!!!
               for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][24])
                  {cutEff[cut][24]++; cut_flag[cut][24] = true;}

	       Jet1.push_back(indexGoodJets.at(i));
	       Jet2.push_back(indexGoodJets.at(ii));

          }
       }

    }

    // Check for boosted jets:
    if (tmpIndexGoodJetsAK8.size() > 0) {
      JetCNT = 0;
      for (unsigned int i = 0; i != tmpIndexGoodJetsAK8.size(); i++) {

         if (!ReducedTree->AK8Jets_AK8isLooseJetId[tmpIndexGoodJetsAK8.at(i)]) continue;
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][14]) {cutEff[cut][14]++; cut_flag[cut][14] = true;}

         if (ReducedTree->AK8Jets_PtCorr[tmpIndexGoodJetsAK8.at(i)] <= 80) continue;
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][15]) {cutEff[cut][15]++; cut_flag[cut][15] = true;}

         if (ReducedTree->AK8JetsEta[tmpIndexGoodJetsAK8.at(i)] >= 2.4) continue;
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][16]) {cutEff[cut][16]++; cut_flag[cut][16] = true;}

         if (ReducedTree->AK8Jets_bDiscriminatorCSV[tmpIndexGoodJetsAK8.at(i)]>=0.679)   continue;
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][17]) {cutEff[cut][17]++; cut_flag[cut][17] = true;}

         if (deltaR(WWTree->l_eta, WWTree->l_phi, ReducedTree->AK8JetsEta[tmpIndexGoodJetsAK8.at(i)], ReducedTree->AK8JetsPhi[tmpIndexGoodJetsAK8.at(i)]) <= 0.3) continue;
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][18]) {cutEff[cut][18]++; cut_flag[cut][18] = true;}

         if (fabs(WWTree->pfMET_Phi-ReducedTree->AK8JetsPhi[tmpIndexGoodJetsAK8.at(i)]) <= 0.4) continue;
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][19]) {cutEff[cut][19]++; cut_flag[cut][19] = true;}

         JetCNT++;
      }

      if (JetCNT == 1) {
      for (unsigned int i = 0; i != tmpIndexGoodJetsAK8.size(); i++) {

         if (ReducedTree->AK8Jets_PtCorr[tmpIndexGoodJetsAK8.at(i)] < 200) continue;
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][20]) {cutEff[cut][20]++; cut_flag[cut][20] = true;}

         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][21]) {cutEff[cut][21]++; cut_flag[cut][21] = true;}

         if (ReducedTree->AK8Jets_tau2[tmpIndexGoodJetsAK8.at(i)]/ReducedTree->AK8Jets_tau1[tmpIndexGoodJetsAK8.at(i)] >= 0.6) continue;
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][22]) {cutEff[cut][22]++; cut_flag[cut][22] = true;}

         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][23]) {cutEff[cut][23]++; cut_flag[cut][23] = true;}

	 if (ReducedTree->AK8Jets_prunedMass[tmpIndexGoodJetsAK8.at(i)] > 70 && ReducedTree->AK8Jets_prunedMass[tmpIndexGoodJetsAK8.at(i)] < 100) continue; // PLEASE CHECK THIS!!!
         for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][24]) {cutEff[cut][24]++; cut_flag[cut][24] = true;}

         indexGoodJetsAK8.push_back(tmpIndexGoodJetsAK8.at(i)); //save index of the "good" vbf jets candidate

       }
       }
    }

    if (Jet1.size() == 0 && indexGoodJetsAK8.size() == 0) continue;
    if (Jet1.size() != 0) for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][25])
       {cut_flag[cut][25] = true; cutEff[cut][25]++;}
    if (indexGoodJetsAK8.size() != 0) for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][25])
       {cut_flag[cut][25] = true; cutEff[cut][25]++;}


    int nVBF1=-1, nVBF2=-1; //index of the two vbf jets
    int nBJ=-1;
    /////////////////////////The Selected Photon:
    int nTightPhoton=0;
    float tempPt_ph = 0.;
    for (unsigned i=0; i != ReducedTree->PhotonsNum; i++) {

        if (ReducedTree->PhotonsPt[i]<tempPt_ph) continue;

        if (ReducedTree->Photons_hasPixelSeed[i]) continue;

        if (ReducedTree->Photons_isTight[i]==false) continue;
        if (!cut_flag[0][26]) for (int cut = 0; cut != 3; cut++) if (cut_flag[cut][25])
           {cutEff[cut][26]++; cut_flag[cut][26] = true;}

        if (ReducedTree->PhotonsPt[i]<=30) continue;
        if (!cut_flag[0][27]) for (int cut = 0; cut != 3; cut++) if (cut_flag[cut][26])
           {cutEff[cut][27]++; cut_flag[cut][27] = true;}

	if (fabs(ReducedTree->PhotonsEta[i])>=1.4421) continue;
        if (!cut_flag[0][28]) for (int cut = 0; cut != 3; cut++) if (cut_flag[cut][27])
           {cutEff[cut][28]++; cut_flag[cut][28] = true;}

        PH.SetPtEtaPhiE(ReducedTree->PhotonsPt[i],ReducedTree->PhotonsEta[i],ReducedTree->PhotonsPhi[i],
		ReducedTree->PhotonsE[i]);
        if (deltaR(PH.Eta(),PH.Phi(),LEP.Eta(),LEP.Phi())<=0.5) continue;
        if (!cut_flag[0][29]) for (int cut = 0; cut != 3; cut++) if (cut_flag[cut][28])
           {cutEff[cut][29]++; cut_flag[cut][29] = true;}

	if (strcmp(leptonName.c_str(),"el")==0) {
            if (fabs((LEP+PH).M()-91.1876)<=10.) continue;
	}
        if (!cut_flag[0][30]) for (int cut = 0; cut != 3; cut++) if (cut_flag[cut][29])
           {cutEff[cut][30]++; cut_flag[cut][30] = true;}

	int nVBF1_=-1, nVBF2_=-1, nBJ_=-1;
	if (Jet2.size() > 0) {

           unsigned int j = -1;
           while (nVBF1_ == -1 && j != Jet1.size()-1) {
               j++;
               VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[Jet1.at(j)],ReducedTree->JetsEta[Jet1.at(j)],ReducedTree->JetsPhi[Jet1.at(j)],ReducedTree->Jets_ECorr[Jet1.at(j)]);
               VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[Jet2.at(j)],ReducedTree->JetsEta[Jet2.at(j)],ReducedTree->JetsPhi[Jet2.at(j)],ReducedTree->Jets_ECorr[Jet2.at(j)]);

               if (deltaR(PH.Eta(),PH.Phi(),VBF1.Eta(),VBF1.Phi())<=0.5) continue;
	       if (deltaR(PH.Eta(),PH.Phi(),VBF2.Eta(),VBF2.Phi())<=0.5) continue;
	       nVBF1_ = Jet1.at(j);
	       nVBF2_ = Jet2.at(j);
	   }// End loop over AK4 jets

        }

	if (indexGoodJetsAK8.size() > 0) {
           unsigned int j = -1;
           while (nBJ_ == -1 && j != indexGoodJetsAK8.size()-1) {
               j++;
               VBF1.SetPtEtaPhiE(ReducedTree->AK8Jets_PtCorr[indexGoodJetsAK8.at(j)],ReducedTree->AK8JetsEta[indexGoodJetsAK8.at(j)],ReducedTree->AK8JetsPhi[indexGoodJetsAK8.at(j)],ReducedTree->AK8Jets_ECorr[indexGoodJetsAK8.at(j)]);
               if (deltaR(PH.Eta(),PH.Phi(),VBF1.Eta(),VBF1.Phi())<=0.5) continue;
               nBJ_ = indexGoodJetsAK8.at(j);

           }// End loop over AK8 jets

	}

        if (nVBF2_ == -1 && nBJ_ == -1) continue;

	nVBF1 = nVBF1_; nVBF2 = nVBF2_; nBJ = nBJ_;

        tightPh.push_back(PH);
        WWTree->photon_pt  = PH.Pt();
        WWTree->photon_eta = PH.Eta();
        WWTree->photon_phi = PH.Phi();
        WWTree->photon_e= PH.E();
        tempPt_ph = PH.Pt();
        nTightPhoton++;
    }
    if (nTightPhoton==0) continue;
    for (int cut = 0; cut != 3; cut++) if (!cut_flag[cut][31] && cut_flag[cut][30])
       {cutEff[cut][31]++; cut_flag[cut][31] = true;}

    /////////////////MC Infos
    if (isMC)
      {
	TLorentzVector temp, temp2;
	//	std::cout<<"entry: "<<iEntry<<" "<<GenNuNum<<std::endl;
	double deltaPhiOld=100.;
	for (int i=0; i<ReducedTree->GenBosonNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenBosonPhi[i],WWTree->v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;
	  //	  std::cout<<"bosone: "<<i<<" "<<ReducedTree->GenBosonPhi[i]<<" "<<v_phi<<std::endl;
	  temp.SetPtEtaPhiE(ReducedTree->GenBosonPt[i],ReducedTree->GenBosonEta[i],ReducedTree->GenBosonPhi[i],ReducedTree->GenBosonE[i]);
	  WWTree->W_pt_gen = ReducedTree->GenBosonPt[i];
	  WWTree->W_pz_gen = temp.Pz();
	  WWTree->W_rap_gen = temp.Rapidity();
	  deltaPhiOld = deltaPhi;
	}	
	if (ReducedTree->GenBosonNum==2) {
	  temp.SetPtEtaPhiE(ReducedTree->GenBosonPt[0],ReducedTree->GenBosonEta[0],ReducedTree->GenBosonPhi[0],ReducedTree->GenBosonE[0]);
	  temp2.SetPtEtaPhiE(ReducedTree->GenBosonPt[1],ReducedTree->GenBosonEta[1],ReducedTree->GenBosonPhi[1],ReducedTree->GenBosonE[1]);
	  WWTree->genGravMass=(temp+temp2).M();	
	}

	deltaPhiOld=100.;
       	for (int i=0; i<ReducedTree->GenNuNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenNuPhi[i],WWTree->v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;	  
	  temp.SetPtEtaPhiE(ReducedTree->GenNuPt[i],ReducedTree->GenNuEta[i],ReducedTree->GenNuPhi[i],ReducedTree->GenNuE[i]);
	  WWTree->nu_pz_gen=temp.Pz();	  
	  WWTree->nu_pt_gen=temp.Pt();	  
	  WWTree->nu_phi_gen=temp.Phi();	  
	  WWTree->nu_eta_gen=temp.Eta();
	  deltaPhiOld = deltaPhi;
	}		
      }

    if (WWTree->v_pt <= 200 && nVBF2 != -1) {

       VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF1],ReducedTree->JetsEta[nVBF1],ReducedTree->JetsPhi[nVBF1],ReducedTree->Jets_ECorr[nVBF1]);
       VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF2],ReducedTree->JetsEta[nVBF2],ReducedTree->JetsPhi[nVBF2],ReducedTree->Jets_ECorr[nVBF2]);

       WWTree->vbf_maxpt_j2_pt = VBF2.Pt();
       WWTree->vbf_maxpt_j2_eta = VBF2.Eta();
       WWTree->vbf_maxpt_j2_phi = VBF2.Phi();
       WWTree->vbf_maxpt_j2_e = VBF2.E();
       WWTree->vbf_maxpt_j1_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorCSV[nVBF1];
       WWTree->vbf_maxpt_j2_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorCSV[nVBF2];

       TOT = VBF1+VBF2;
       WWTree->vbf_maxpt_jj_m = TOT.M();	
       WWTree->isBoostedSignal = 0;

    } else if (WWTree->v_pt > 200 && nBJ != -1) {

       VBF1.SetPtEtaPhiE(ReducedTree->AK8Jets_PtCorr[nBJ],ReducedTree->AK8JetsEta[nBJ],ReducedTree->AK8JetsPhi[nBJ],ReducedTree->AK8Jets_ECorr[nBJ]);
       WWTree->vbf_maxpt_j1_bDiscriminatorCSV = ReducedTree->AK8Jets_bDiscriminatorCSV[nBJ];
       WWTree->jet_tau2tau1 = ReducedTree->AK8Jets_tau2[nBJ]/ReducedTree->AK8Jets_tau1[nBJ];

       TOT = VBF1;
       WWTree->vbf_maxpt_jj_m = ReducedTree->AK8Jets_prunedMass[nBJ];	
       WWTree->isBoostedSignal = 1;

    } else continue;

    if (WWTree->isBoostedSignal == 0) for (int cut = 0; cut != 2; cut++) if (!cut_flag[cut][32] && cut_flag[cut][31])
       {cutEff[cut][32]++; cut_flag[cut][32] = true;}
    if (WWTree->isBoostedSignal == 1) for (int cut = 0; cut < 3; cut+=2) if (!cut_flag[cut][32] && cut_flag[cut][31])
       {cutEff[cut][32]++; cut_flag[cut][32] = true;}

    WWTree->vbf_maxpt_j1_pt = VBF1.Pt();
    WWTree->vbf_maxpt_j1_eta = VBF1.Eta();
    WWTree->vbf_maxpt_j1_phi = VBF1.Phi();
    WWTree->vbf_maxpt_j1_e = VBF1.E();

    WWTree->vbf_maxpt_jj_pt = TOT.Pt();
    WWTree->vbf_maxpt_jj_eta = TOT.Eta();
    WWTree->vbf_maxpt_jj_phi = TOT.Phi();

    //Pass along the AQGC weights:
    if (isAQGC) {
    WWTree->originalWeight = ReducedTree->originalWeight;
    WWTree->AQGCweights = *(ReducedTree->AQGCweights);
    }

    // Pass along SF's:
    if (isMC) {
       WWTree->scaleFactor *= SF_Photon(WWTree->photon_eta, WWTree->photon_pt);
       if (strcmp(leptonName.c_str(),"el")==0) WWTree->scaleFactor *= SF_Elec(WWTree->l_eta, WWTree->l_pt);
       else WWTree->scaleFactor *= SF_Muon(WWTree->l_eta, WWTree->l_pt);
    }

    //fill the tree
    outTree->Fill();
    cutEff[0][33]++; cut_flag[0][33] = true;
  }

  std::cout<<"\nProcessed/Total Events:	"<<max_evt<<"/"<<ReducedTree->fChain->GetEntries()<<std::endl
           <<"\nEvent Selection"<<"			Tot.	Res.	Boosted"<<std::endl
           <<"HBHE Loose Rerun Filter:	"<<cutEff[0][0]<<"	"<<cutEff[1][0]<<"	"<<cutEff[2][0]<<std::endl
           <<"HBHE Iso Rerun Filter:		"<<cutEff[0][1]<<"	"<<cutEff[1][1]<<"	"<<cutEff[2][1]<<std::endl
           <<"CSC Halo Filter:		"<<cutEff[0][2]<<"	"<<cutEff[1][2]<<"	"<<cutEff[2][2]<<std::endl
           <<"Good Vertex Fiter:		"<<cutEff[0][3]<<"	"<<cutEff[1][3]<<"	"<<cutEff[2][3]<<std::endl
           <<"EE Bad SC Filter:		"<<cutEff[0][4]<<"	"<<cutEff[1][4]<<"	"<<cutEff[2][4]<<std::endl
           <<"Lepton Triggers:		"<<cutEff[0][5]<<"	"<<cutEff[1][5]<<"	"<<cutEff[2][5]<<std::endl
           <<"Muon Tracks:			"<<cutEff[0][6]<<"	"<<cutEff[1][6]<<"	"<<cutEff[2][6]<<std::endl
           <<"Lepton ID:			"<<cutEff[0][7]<<"	"<<cutEff[1][7]<<"	"<<cutEff[2][7]<<std::endl
	   <<"Lepton pT:			"<<cutEff[0][8]<<"	"<<cutEff[1][8]<<"	"<<cutEff[2][8]<<std::endl
	   <<"Lepton Eta:			"<<cutEff[0][9]<<"	"<<cutEff[1][9]<<"	"<<cutEff[2][9]<<std::endl
	   <<"Single Lep. Veto:		"<<cutEff[0][10]<<"	"<<cutEff[1][10]<<"	"<<cutEff[2][10]<<std::endl
	   <<"Second Lep. Veto:		"<<cutEff[0][11]<<"	"<<cutEff[1][11]<<"	"<<cutEff[2][11]<<std::endl
	   <<"MET pT:				"<<cutEff[0][12]<<"	"<<cutEff[1][12]<<"	"<<cutEff[2][12]<<std::endl
	   <<"W mT:				"<<cutEff[0][13]<<"	"<<cutEff[1][13]<<"	"<<cutEff[2][13]<<std::endl
	   <<"Jet ID:				"<<cutEff[0][14]<<"	"<<cutEff[1][14]<<"	"<<cutEff[2][14]<<std::endl
           <<"Jet pT:				"<<cutEff[0][15]<<"	"<<cutEff[1][15]<<"	"<<cutEff[2][15]<<std::endl
           <<"Jet Eta:			"<<cutEff[0][16]<<"	"<<cutEff[1][16]<<"	"<<cutEff[2][16]<<std::endl
           <<"b Jet Veto:			"<<cutEff[0][17]<<"	"<<cutEff[1][17]<<"	"<<cutEff[2][17]<<std::endl
           <<"Jet-Lep. iso:			"<<cutEff[0][18]<<"	"<<cutEff[1][18]<<"	"<<cutEff[2][18]<<std::endl
           <<"Jet-MET iso:			"<<cutEff[0][19]<<"	"<<cutEff[1][19]<<"	"<<cutEff[2][19]<<std::endl
	   <<"Jet pT Veto:			"<<cutEff[0][20]<<"	"<<cutEff[1][20]<<"	"<<cutEff[2][20]<<std::endl
           <<"Jet Count Veto:			"<<cutEff[0][21]<<"	"<<cutEff[1][21]<<"	"<<cutEff[2][21]<<std::endl
           <<"Jet substructure:		"<<cutEff[0][22]<<"	"<<cutEff[1][22]<<"	"<<cutEff[2][22]<<std::endl
           <<"Jet(s) pT:			"<<cutEff[0][23]<<"	"<<cutEff[1][23]<<"	"<<cutEff[2][23]<<std::endl
           <<"Jet(s) Mass:			"<<cutEff[0][24]<<"	"<<cutEff[1][24]<<"	"<<cutEff[2][24]<<std::endl
           <<"Final Jet Veto:			"<<cutEff[0][25]<<"	"<<cutEff[1][25]<<"	"<<cutEff[2][25]<<std::endl
           <<"Tight Photon:			"<<cutEff[0][26]<<"	"<<cutEff[1][26]<<"	"<<cutEff[2][26]<<std::endl
           <<"Photon pT:			"<<cutEff[0][27]<<"	"<<cutEff[1][27]<<"	"<<cutEff[2][27]<<std::endl
           <<"Photon Eta:			"<<cutEff[0][28]<<"	"<<cutEff[1][28]<<"	"<<cutEff[2][28]<<std::endl
           <<"Photon-Lep. iso:		"<<cutEff[0][29]<<"	"<<cutEff[1][29]<<"	"<<cutEff[2][29]<<std::endl
           <<"Fake e- cut:			"<<cutEff[0][30]<<"	"<<cutEff[1][30]<<"	"<<cutEff[2][30]<<std::endl
           <<"Jet-Photon iso:			"<<cutEff[0][31]<<"	"<<cutEff[1][31]<<"	"<<cutEff[2][31]<<std::endl
           <<"Lep. W pT:			"<<cutEff[0][32]<<"	"<<cutEff[1][32]<<"	"<<cutEff[2][32]<<std::endl
           <<"All cuts:			"<<cutEff[0][33]<<std::endl<<std::endl;

  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}

