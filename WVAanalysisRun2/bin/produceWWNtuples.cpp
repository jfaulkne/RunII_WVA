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

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  bool isMC = argv[3];
  std::string leptonName = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string numberOfEntries = argv[8];
  float weight = std::atof(xSecWeight.c_str())/std::atof(numberOfEntries.c_str());
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }

  TLorentzVector W,MET,LEP;
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

  int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_")+leptonName+std::string("/")+outputFile).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);

  setOutputTree *WWTree = new setOutputTree(outTree);

  //---------start loop on events------------
  for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++) {

    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;

    int nb = ReducedTree->fChain->GetEntry(jentry);   
    // if (Cut(ientry) < 0) continue;                    

    tightPh.clear();
    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();

    if(iEntry % 1000 == 0)    
      cout << "read entry: " << iEntry << endl;

    WWTree->initializeVariables(); //initialize all variables
    
    WWTree->issignal = 0;
    WWTree->wSampleWeight = weight; //xsec/numberOfEntries
    WWTree->totalEventWeight = 1.; //temporary value
    WWTree->eff_and_pu_Weight = 1.; //temporary value
    
    //save event variables
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi = ReducedTree->LumiBlockNum;
    WWTree->nPV  = ReducedTree->NVtx;
    WWTree->nPho = ReducedTree->PhotonsNum; 
    
    /////////////////THE SELECTED LEPTON
    int nTightLepton=0;
    if (strcmp(leptonName.c_str(),"el")==0) {
      float tempPt=0.;
      for (int i=0; i<ReducedTree->ElectronsNum; i++) {
	//if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
        if (!ReducedTree->Electrons_isTight[i]) continue;
        if (ReducedTree->ElectronsPt[i]<=30) continue;
	if (ReducedTree->ElectronsPt[i]<tempPt) continue;
	ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
	if ((fabs(ELE.Eta())>=1.44 && fabs(ELE.Eta())<=1.57) || fabs(ELE.Eta())>=2.5) continue;
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
      float tempPt=0.;
      for (int i=0; i<ReducedTree->MuonsNum; i++) {
	//if (ReducedTree->Muons_isHighPt[i]==false) continue;
        if (!ReducedTree->Muons_isTight[i]) continue;
	if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
        if (ReducedTree->MuonsPt[i]<=25) continue;
        if (fabs(ReducedTree->MuonsEta[i])>=2.1) continue;
	MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
	tightMuon.push_back(MU);
	if (ReducedTree->MuonsPt[i]<tempPt) continue;
	WWTree->l_pt  = ReducedTree->MuonsPt[i];
	WWTree->l_eta = ReducedTree->MuonsEta[i];
	WWTree->l_phi = ReducedTree->MuonsPhi[i];
	WWTree->l_e = ReducedTree->MuonsE[i];
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    }
    if (nTightLepton==0) continue; //no leptons with required ID
    cutEff[0]++;

    //VETO ADDITIONAL LEPTONS
    int nLooseLepton=0;
    for (int i=0; i<ReducedTree->ElectronsNum; i++) {
      //if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
      if (!ReducedTree->Electrons_isTight[i]) continue;
      if (ReducedTree->ElectronsPt[i]<20) continue;       
      ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
      if ((fabs(ELE.Eta())>=1.44 && fabs(ELE.Eta())<=1.57) || fabs(ELE.Eta())>=2.5) continue;
      looseEle.push_back(ELE);      
      nLooseLepton++;
    }
    for (int i=0; i<ReducedTree->MuonsNum; i++) {
      //if (ReducedTree->Muons_isHighPt[i]==false) continue;
      if (!ReducedTree->Muons_isTight[i]) continue;
      if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
      if (fabs(ReducedTree->MuonsEta[i])>=2.1) continue;
      if (ReducedTree->MuonsPt[i]<10) continue;
      MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
      looseMuon.push_back(MU);
      nLooseLepton++;
    }
    if (nLooseLepton!=1) continue; //no additional leptons
    cutEff[1]++;


    //preselection on jet pt and met
    if (ReducedTree->METPt <= 35) continue; 
    cutEff[2]++;

    MET.SetPtEtaPhiE(ReducedTree->METPt,0.,ReducedTree->METPhi,0.);
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);

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
    double pz2_type0 = NeutrinoPz_type0.getOther(); // Default one

    double pz1_run2 = NeutrinoPz_run2.Calculate(); 

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    // chenge the neutrino pT in case of complex solution in order to make it real
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
    double pz2_type2 = NeutrinoPz_type2.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    // chenge the neutrino pT in case of complex solution in order to make it real
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
    cutEff[3]++;

    /////////VBF and b-tag section
    bool fillVBF = true;

    std::vector<int> indexGoodJets;
    indexGoodJets.clear();
    if (indexGoodJets.size()!=0)  fillVBF=false;

    WWTree->njets=0;
    WWTree->nBTagJet_loose=0;
    WWTree->nBTagJet_medium=0;
    WWTree->nBTagJet_tight=0;

    float oldDeltaR = 1000.;
    float oldDeltaRLep = 1000.;
    int indexCloserJet = -1;
    int indexCloserJetLep = -1;

    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
      {
	bool isCleanedJet = true;
	if (ReducedTree->Jets_PtCorr[i]<=30 || ReducedTree->JetsPt[i]<=20 || fabs(ReducedTree->JetsEta[i])>=2.4)  continue;
	if (ReducedTree->Jets_isLooseJetId[i]==false) continue;

	//CLEANING FROM LEPTONS
	if (deltaR(WWTree->l_eta, WWTree->l_phi, ReducedTree->JetsEta[i], ReducedTree->JetsPhi[i]) <= 0.3){
	    isCleanedJet = false;
	}
	if (isCleanedJet==false) continue;

	//fill B-Tag info
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.423)   WWTree->nBTagJet_loose++;
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.814)   WWTree->nBTagJet_medium++;
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.941)   WWTree->nBTagJet_tight++;

        if (ReducedTree->Jets_bDiscriminatorCSV[i]>=0.679)   continue;
	WWTree->njets++;
	indexGoodJets.push_back(i); //save index of the "good" vbf jets candidate
      }
    if (indexGoodJets.size()<2){ fillVBF=false; /*continue;*/ } //check if at least 2 jets are inside the collection

    float tempPtMax=0.;
    std::vector<int> Jet1, Jet2;
    if (fillVBF) {
       cutEff[4]++;
       for (unsigned int i=0; i<indexGoodJets.size()-1; i++) {
          for (unsigned int ii=i+1; ii<indexGoodJets.size(); ii++) {
	       VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->Jets_ECorr[indexGoodJets.at(i)]);
	       VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(ii)],ReducedTree->JetsEta[indexGoodJets.at(ii)],ReducedTree->JetsPhi[indexGoodJets.at(ii)],ReducedTree->Jets_ECorr[indexGoodJets.at(ii)]);
               if (fabs(VBF1.Eta()-VBF2.Eta()) >= 1.4) continue;
	       TOT = VBF1 + VBF2;
	       if (TOT.M() <= 70 || TOT.M() >= 100) continue; // PLEASE CHECK THIS!!!
	       Jet1.push_back(indexGoodJets.at(i));
	       Jet2.push_back(indexGoodJets.at(ii));
	       if (TOT.Pt() < tempPtMax) continue;
	       tempPtMax = TOT.Pt(); //take the jet pair with largest Pt
          }
       }
    }

    // Check for boosted jets - IF no VBF jets:
    if (Jet1.size() == 0 && (ReducedTree->AK8JetsNum < 1 || ReducedTree->AK8JetsPt[0] < 30)) continue;
    cutEff[5]++;

    int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
    /////////////////////////The Selected Photon:
    int nTightPhoton=0;
    float tempPt_ph = 0.;
    for (unsigned i=0; i<ReducedTree->PhotonsNum; i++) {
        if (ReducedTree->Photons_isTight[i]==false) continue;
        if (ReducedTree->PhotonsPt[i]<=30) continue;
        if (ReducedTree->PhotonsPt[i]<tempPt_ph) continue;
	if (fabs(ReducedTree->PhotonsEta[i])>=1.4421) continue;

        PH.SetPtEtaPhiE(ReducedTree->PhotonsPt[i],ReducedTree->PhotonsEta[i],ReducedTree->PhotonsPhi[i],
		ReducedTree->PhotonsE[i]);
//        if ((fabs(PH.Eta())>=1.44 && fabs(PH.Eta())<=1.57) || fabs(PH.Eta())>=2.5) continue;
        if (deltaR(PH.Eta(),PH.Phi(),LEP.Eta(),LEP.Phi())<=0.5) continue;
	if (strcmp(leptonName.c_str(),"el")==0) {
            if (fabs((LEP+PH).M()-91.1876)<=10.) continue;
	}

	float tempPt_jets = 0.;
	int nVBF1_=-1, nVBF2_=-1;
	if (Jet2.size() > 0) {

	   for (unsigned int j = 0; j<Jet1.size(); j++) {
               VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[Jet1.at(j)],ReducedTree->JetsEta[Jet1.at(j)],ReducedTree->JetsPhi[Jet1.at(j)],ReducedTree->Jets_ECorr[Jet1.at(j)]);
               VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[Jet2.at(j)],ReducedTree->JetsEta[Jet2.at(j)],ReducedTree->JetsPhi[Jet2.at(j)],ReducedTree->Jets_ECorr[Jet2.at(j)]);

               if (deltaR(PH.Eta(),PH.Phi(),VBF1.Eta(),VBF1.Phi())<=0.5) continue;
	       if (deltaR(PH.Eta(),PH.Phi(),VBF2.Eta(),VBF2.Phi())<=0.5) continue;
	       TOT = VBF1 + VBF2;
	       if (TOT.M() < tempPt_jets) continue;
	       tempPt_jets = TOT.M();
	       nVBF1_ = Jet1.at(j);
	       nVBF2_ = Jet2.at(j);
	   }// End loop over AK4 jets

        }

	if (nVBF1_ == -1 && ReducedTree->AK8JetsNum > 0 && ReducedTree->AK8JetsPt[0] > 30) {

           tempPt_jets = 0.;
           for (unsigned int j = 0; j<ReducedTree->AK8JetsNum; j++) {

	       if (ReducedTree->AK8JetsPt[j] < 30 || ReducedTree->AK8JetsEta[j] > 4.7) continue;
	       if (deltaR(WWTree->l_eta, WWTree->l_phi, ReducedTree->AK8JetsEta[j], ReducedTree->AK8JetsPhi[j]) <= 0.3) continue;
               if (ReducedTree->AK8JetsPt[j] <= tempPt_jets) continue;
               VBF1.SetPtEtaPhiE(ReducedTree->AK8JetsPt[j],ReducedTree->AK8JetsEta[j],ReducedTree->AK8JetsPhi[j],ReducedTree->AK8JetsE[j]);
               if (deltaR(PH.Eta(),PH.Phi(),VBF1.Eta(),VBF1.Phi())<=0.5) continue;
               if (VBF1.M() <= 70. || VBF1.M() >= 100.) continue;
               if (VBF1.M() < tempPt_jets) continue;
               tempPt_jets = VBF1.M();
               nVBF1_ = (int)j;

           }// End loop over AK8 jets

	}

        if (nVBF1_ == -1) continue;
	nVBF1 = nVBF1_; nVBF2 = nVBF2_;

        tightPh.push_back(PH);
        WWTree->photon_pt  = PH.Pt();
        WWTree->photon_eta = PH.Eta();
        WWTree->photon_phi = PH.Phi();
        WWTree->photon_e= PH.E();
        tempPt_ph = PH.Pt();
        nTightPhoton++;
    }
    if (nTightPhoton==0) continue;
    cutEff[6]++;

    if (nVBF2 > 0) {
       VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF1],ReducedTree->JetsEta[nVBF1],ReducedTree->JetsPhi[nVBF1],ReducedTree->Jets_ECorr[nVBF1]);
       VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF2],ReducedTree->JetsEta[nVBF2],ReducedTree->JetsPhi[nVBF2],ReducedTree->Jets_ECorr[nVBF2]);
       WWTree->vbf_maxpt_j2_pt = ReducedTree->Jets_PtCorr[nVBF2];
       WWTree->vbf_maxpt_j2_eta = ReducedTree->JetsEta[nVBF2];
       WWTree->vbf_maxpt_j2_phi = ReducedTree->JetsPhi[nVBF2];
       WWTree->vbf_maxpt_j2_e = ReducedTree->Jets_ECorr[nVBF2];
       WWTree->vbf_maxpt_j2_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nVBF2];
    } 
    else VBF1.SetPtEtaPhiE(ReducedTree->AK8JetsPt[nVBF1],ReducedTree->AK8JetsEta[nVBF1],ReducedTree->AK8JetsPhi[nVBF1],ReducedTree->AK8JetsE[nVBF1]);
    WWTree->vbf_maxpt_j1_pt = VBF1.Pt();
    WWTree->vbf_maxpt_j1_eta = VBF1.Eta();
    WWTree->vbf_maxpt_j1_phi = VBF1.Phi();
    WWTree->vbf_maxpt_j1_e = VBF1.E();
    WWTree->vbf_maxpt_j1_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nVBF1];

    TOT = VBF1;
    if (nVBF2 > 0) TOT += VBF2;

    WWTree->vbf_maxpt_jj_pt = TOT.Pt();
    WWTree->vbf_maxpt_jj_eta = TOT.Eta();
    WWTree->vbf_maxpt_jj_phi = TOT.Phi();
    WWTree->vbf_maxpt_jj_m = TOT.M();	

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
    cutEff[7]++;
    
    //fill the tree
    outTree->Fill();
  }

  std::cout<<std::endl<<"lepton eff: "<<cutEff[0]<<std::endl
	   <<"single lep eff: "<<cutEff[1]<<std::endl
	   <<"met eff:    "<<cutEff[2]<<std::endl
	   <<"W eff:      "<<cutEff[3]<<std::endl
	   <<"2 AK4 jets found:  "<<cutEff[4]<<std::endl
	   <<"Jet sec:	  "<<cutEff[5]<<std::endl
	   <<"Photon eff:    "<<cutEff[6]<<std::endl
	   <<"MC eff:    "<<cutEff[7]<<std::endl;

  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}

