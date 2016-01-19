/*//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

Channel: WW+gamma semileptonic Muon channel (W -> uv, W -> jj, gamma)

Usage: Reads in multiple ROOT files (RD Trees) representing data, MC
       signal, and background.  The background normalization is done
       using theoretical cross sections, observed luminosity, and
       sample sizes (user can also add K-factors, etc.).  User
       specifies kinematical distribution to observe, as well as
       histogram range and binning.  User can also specify the
       Systematic/Statistical uncertainty.

Output: ROOT file containing 1D histograms for:
	(1) Observed data (data_obs)
	(2) SM backgrounds + SM signal (background)
	(3) Background + uncertainty (background_mudijet_backshapeUp)
	(4) Background - uncertainty (background_mudijet_backshapeDown)
	(5-?) Anomalous signal - SM signal (signal_??)

Syntax: root -l -b -q mkROOTsm.C

///////////////////////////////////////////////////////////////////////
*//////////////////////////////////////////////////////////////////////


/////////////////////
// Load header files:
#include <iostream>
#include "TString.h"
#include "TLatex.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TF1.h"


///////////////////////////////////////////
// Define Type to store histogram settings:
struct plotVar_t {
  char* plotvar;
  double MINRange;
  double MAXRange;
  int    NBINS;
  int    slog;
  char* xlabel;
};

void cmspre(double COM, double intlumifbinv)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);

  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.9,0.925,Form("#sqrt{s} = %.0f TeV", (float)COM));
  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.7,0.925,Form("#int #font[12]{L} dt = %.1f pb^{-1}", (float)intlumifbinv));

  latex.SetTextAlign(11); // align left
  latex.DrawLatex(0.1,0.925,"CMS preliminary");

}


///////////////////////
// Begin Main function:
void mkROOTsm_ControlPlots(){

  //////////////////////////////////////////////////
  // Provide channel-specific int. Lumi in inv. pb.:
  const double intLUMI[2] = {2168.45, 2151.91};// [pb^-1]

  ////////////////////////////////
  // Specify event selection cuts:
  // Added high jet pT cut for boosted regime
  TCut the_cut("((totalEventWeight==0)? 1:totalEventWeight)*genWeight*(photon_pt>0&&isBoostedSignal==0&&(vbf_maxpt_jj_m<70||vbf_maxpt_jj_m>100))");
  TCut the_cut_LPT("((totalEventWeight==0)? 1:totalEventWeight)*genWeight*(photon_pt<500&&isBoostedSignal==0&&(vbf_maxpt_jj_m<70||vbf_maxpt_jj_m>100))");
  TCut the_cut_HPT("((totalEventWeight==0)? 1:totalEventWeight)*genWeight*(photon_pt>=500&&isBoostedSignal==0&&(vbf_maxpt_jj_m<70||vbf_maxpt_jj_m>100))");

  TString FakeRate = TString("((photon_pt<30)? .33:(photon_pt<40)? .246:(photon_pt<50)? .214:(photon_pt<75)? .151:(photon_pt<100)? .087:(photon_pt<135)? 0.06:0.063)");
  //TString FakeRate = TString("((photon_pt<50)? .12:(photon_pt<75)? .14:(photon_pt<90)? .23:(photon_pt<135)? .22:.39)");

  TCut the_cut_FP(FakeRate+TString("*((totalEventWeight==0)? 1:totalEventWeight)*genWeight*(photon_pt>0&&isBoostedSignal==0&&(vbf_maxpt_jj_m<70||vbf_maxpt_jj_m>100))"));
  TCut the_cut_FP_LPT(FakeRate+TString("*((totalEventWeight==0)? 1:totalEventWeight)*genWeight*(photon_pt<500&&isBoostedSignal==0&&(vbf_maxpt_jj_m<70||vbf_maxpt_jj_m>100))"));
  TCut the_cut_FP_HPT(FakeRate+TString("*((totalEventWeight==0)? 1:totalEventWeight)*genWeight*(photon_pt>=500&&isBoostedSignal==0&&(vbf_maxpt_jj_m<70||vbf_maxpt_jj_m>100))"));
//  TCut the_cut("vbf_maxpt_j1_pt > 250");
//  TCut the_cut_FP(TString("((photon_pt<50)? .12:(photon_pt<75)? .14:(photon_pt<90)? .23:(photon_pt<135)? .22:.39)*(vbf_maxpt_j1_pt>250)"));


  ///////////////////////////
  // Create output ROOT file:
  TFile f("RunII_WVA_ControlPlots.root", "RECREATE");


  //////////////////////////////////////////////////
  // Create file pointers for each sample ROOT file:
  TFile * wwaShape_file[2], * wzaShape_file[2], * wajetsShape_file[2], * wajetsHPTShape_file[2], * zajets_file[2];
  TFile * ttajets_file[2], * tajets_file[2];
  //TFile * sTops_file[2], * sTopt_file[2], * sToptw_file[2], * sTopat_file[2], * sTopatw_file[2];
  TFile * data_file[2];


  //////////////////////////////
  // Open each sample ROOT file:

  ////////////////////
  // Muon channel = 0:

  data_file[0] = new TFile("../WVAanalysisRun2/output/output_mu/RD_Run2015CD_SingleMuon.root");
  wwaShape_file[0] = new TFile("../WVAanalysisRun2/output/output_mu/RD_WWG_mu.root");
  wzaShape_file[0] = new TFile("../WVAanalysisRun2/output/output_mu/RD_WZG_mu.root");
  wajetsShape_file[0] = new TFile("../WVAanalysisRun2/output/output_mu/RD_WGjets_mu.root");// 25ns off.
  wajetsHPTShape_file[0] = new TFile("../WVAanalysisRun2/output/output_mu/RD_WGjets_PtG500_mu.root");
  zajets_file[0] = new TFile("../WVAanalysisRun2/output/output_mu/RD_ZGjets_mu.root");// 25ns off.
  ttajets_file[0] = new TFile("../WVAanalysisRun2/output/output_mu/RD_TTGJets_mu.root");// 25ns off.
  tajets_file[0] = new TFile("../WVAanalysisRun2/output/output_mu/RD_TGJets_mu.root");// 25ns off.


  ////////////////////////
  // Electron channel = 1:

  data_file[1] = new TFile("../WVAanalysisRun2/output/output_el/RD_Run2015CD_SingleElectron.root");
  wwaShape_file[1] = new TFile("../WVAanalysisRun2/output/output_el/RD_WWG_el.root");
  wzaShape_file[1] = new TFile("../WVAanalysisRun2/output/output_el/RD_WZG_el.root");
  wajetsShape_file[1] = new TFile("../WVAanalysisRun2/output/output_el/RD_WGjets_el.root");// 25ns off.
  wajetsHPTShape_file[1] = new TFile("../WVAanalysisRun2/output/output_el/RD_WGjets_PtG500_el.root");
  zajets_file[1] = new TFile("../WVAanalysisRun2/output/output_el/RD_ZGjets_el.root");// 25ns off.
  ttajets_file[1] = new TFile("../WVAanalysisRun2/output/output_el/RD_TTGJets_el.root");// 25ns off.
  tajets_file[1] = new TFile("../WVAanalysisRun2/output/output_el/RD_TGJets_el.root");// 25ns off.


  //////////////////////////////
  // Create histogram directory:
  TDirectory * chDir[2];
  TDirectory * chSubDir[1]; int plot = 0;


  ///////////////////////////////////////////////////////////////////////////////////
  // Specify what kinematical distribution to observe, as well as histogram settings:
  vector<plotVar_t> PV = {
        {"pfMET", 35, 155, 12, 3, "PF MET (GeV)"},
        {"pfMET_Phi", -3.5, 3.5, 14, 3, "PF MET #phi (Radians)"},
        {"l_pt", 25, 155, 13, 3, "Lepton p_{T} (GeV)"},
        {"l_eta", -3, 3, 12, 3, "Lepton #eta"},
        {"l_phi", -3.5, 3.5, 14, 3, "Lepton #phi"},
        {"photon_pt", 30, 450, 10, 3, "Photon E_{T} (GeV)"},
        {"photon_eta", -1.5, 1.5, 6, 3, "Photon #eta"},
        {"photon_phi", -3.5, 3.5, 14, 3, "Photon #phi (Radians)"},
        {"vbf_maxpt_j1_pt", 30, 190, 8, 3, "Leading Jet p_{T} (GeV)"},
        {"vbf_maxpt_j1_eta", -5, 5, 20, 3, "Leading Jet #eta"},
        {"vbf_maxpt_j1_phi", -3.5, 3.5, 14, 3, "Leading Jet #phi (Radians)"},
        {"vbf_maxpt_j2_pt", 30, 190, 8, 3, "Second Leading Jet p_{T} (GeV)"},
        {"vbf_maxpt_j2_eta", -5, 5, 20, 3, "Second Leading Jet #eta"},
        {"vbf_maxpt_j2_phi", -3.5, 3.5, 14, 3, "Second Leading Jet #phi (Radians)"},
        {"vbf_maxpt_jj_pt", 0, 550, 11, 3, "Dijet p_{T} (GeV)"},
        {"vbf_maxpt_jj_eta", -5, 5, 20, 3, "Dijet #eta"},
        {"vbf_maxpt_jj_phi", -3.5, 3.5, 14, 3, "Dijet #phi (Radians)"},
        {"vbf_maxpt_jj_m", 0, 400, 40, 3, "Dijet Mass (GeV)"},
        {"jet_tau2tau1", 0, 1., 5, 1, "Merged Jet #tau_{2}:#tau_{1}"}
  };


  /////////////////////////////////////
  // Loop over muon + electron channel:
  for (int ch = 0; ch != 2; ch++) {
  chDir[ch] = f.mkdir((ch==0)? "Muon":"Electron");


  //////////////////////////////////////
  // Set MC normalization scale factors:

// WA+jets norm - data driven : 1.09896 +- 0.091
//  const double WAJets_scale = 1.09896 * 129.6 * intLUMI[ch]/3423715;
//  const double ZAJets_scale  = 18.2 * intLUMI[ch]/2836163;

  const double WAJets_norm[2] = {2.31, 2.18};//2.35, 2.17};// {mu, el}
  const double FakePhotons_norm[2] = {1.0, 1.0};//1.05733, 1.05899};// {mu, el}

  const double WAJets_scale      = 405.271     * intLUMI[ch]/6099599;//25ns official prod
  const double WAJets_HPT_scale  =   0.0117887 * intLUMI[ch]/1392843;//25ns official prod
  const double ZAJets_scale      = 117.864     * intLUMI[ch]/4455517;// 25ns official prod
  const double TTAJets_scale     =   3.697     * intLUMI[ch]/4843746;//25ns official prod
  const double TAJets_scale      =   2.967     * intLUMI[ch]/400000;//25ns official prod

  // PT dependent k-factor included (same for WWA,WZA aQGC)
  const double WWA_scale     = 0.2147  * intLUMI[ch]/1000000;
  const double WZA_scale     = 0.04123 * intLUMI[ch]/1000000;


  ////////////////////////////
  // Loop over plot variables:
  for (vector<plotVar_t>::iterator pv = PV.begin(); pv != PV.end(); ++pv) {
  std::cout << TString(pv->plotvar) << "\t"<<pv->MINRange<<"\t" << pv->MAXRange<<"\t" << pv->NBINS<<"\tTHE CUT " << endl;


  ///////////////////////////////////////////////////
  // Retrieve ROOT tree with kinematic distributions:
  TTree* treedata    = (TTree*)  data_file[ch]->Get("otree");
  TTree* treewwa     = (TTree*)  wwaShape_file[ch]->Get("otree");
  TTree* treewza     = (TTree*)  wzaShape_file[ch]->Get("otree");
  TTree* treewaj     = (TTree*)  wajetsShape_file[ch]->Get("otree");
  TTree* treewaj_HPT = (TTree*)  wajetsHPTShape_file[ch]->Get("otree");
  TTree* treezaj     = (TTree*)  zajets_file[ch]->Get("otree");
  TTree* treettaj    = (TTree*)  ttajets_file[ch]->Get("otree");
  TTree* treetaj     = (TTree*)  tajets_file[ch]->Get("otree");


  ////////////////////////////////////////////////////////////
  // Create kinematic-distribution histograms for each sample:
  TH1* th1data = new TH1D("th1data", "Run2015CD 13 TeV 25 ns 3.8 T Data", pv->NBINS, pv->MINRange, pv->MAXRange);
  TH1* th1wwa = new TH1D("th1wwa", "SM WWA (Signal)", pv->NBINS,pv->MINRange,pv->MAXRange);
  TH1* th1wza = new TH1D("th1wza", "SM WZA (Signal)", pv->NBINS,pv->MINRange,pv->MAXRange);
  TH1* th1wajets  = new TH1D("th1wajets",  "SM WA+Jets (BKG)",  pv->NBINS ,pv->MINRange,pv->MAXRange);
  TH1* th1wajets_HPT  = new TH1D("th1wajets_HPT",  "SM WA+Jets (BKG) p_{T} > 500 GeV",  pv->NBINS ,pv->MINRange,pv->MAXRange);
  TH1* th1zajets = new TH1D("th1zajets", "SM ZA+Jets (BKG)", pv->NBINS,pv->MINRange,pv->MAXRange);
  TH1* th1ttajets = new TH1D("th1ttajets", "SM TTA+Jets (BKG)", pv->NBINS, pv->MINRange, pv->MAXRange);
  TH1* th1tajets = new TH1D("th1tajets", "SM TA+Jets (BKG)", pv->NBINS, pv->MINRange, pv->MAXRange);

  TH1* th1FakePhotons_ = new TH1D("th1FakePhotons_", "Jet #rightarrow Photon (Fake Photon)", pv->NBINS, pv->MINRange, pv->MAXRange);
  TH1* th1FakePhotons_1 = new TH1D("th1FakePhotons_1", "Jet #rightarrow Photon (Fake Photon)", pv->NBINS, pv->MINRange, pv->MAXRange);
  TH1* th1FakePhotons_2 = new TH1D("th1FakePhotons_2", "Jet #rightarrow Photon (Fake Photon)", pv->NBINS, pv->MINRange, pv->MAXRange);
  TH1* th1FakePhotons_3 = new TH1D("th1FakePhotons_3", "Jet #rightarrow Photon (Fake Photon)", pv->NBINS, pv->MINRange, pv->MAXRange);
  TH1* th1FakePhotons_4 = new TH1D("th1FakePhotons_4", "Jet #rightarrow Photon (Fake Photon)", pv->NBINS, pv->MINRange, pv->MAXRange);
  TH1* th1FakePhotons_5 = new TH1D("th1FakePhotons_5", "Jet #rightarrow Photon (Fake Photon)", pv->NBINS, pv->MINRange, pv->MAXRange);
  TH1* th1FakePhotons_6 = new TH1D("th1FakePhotons_6", "Jet #rightarrow Photon (Fake Photon)", pv->NBINS, pv->MINRange, pv->MAXRange);


  /////////////////////////////////////////////////////////////////////////
  // Specify histograms to store Sum of Squares of Weights for each sample:
  th1data->Sumw2();
  th1wwa->Sumw2();
  th1wza->Sumw2();
  th1wajets->Sumw2();
  th1wajets_HPT->Sumw2();
  th1zajets->Sumw2();
  th1ttajets->Sumw2();
  th1tajets->Sumw2();


  ///////////////////////////////////////////////////////////////////////////////////
  // Fill kinematical distribution for each sample according to event selection cuts:
  std::cout<<"Fill Data Histogram..."<<std::endl;
  treedata->Draw(TString(pv->plotvar)+TString(">>th1data"), the_cut, "goff");
  th1data->AddBinContent(pv->NBINS,th1data->GetBinContent(pv->NBINS+1));th1data->SetBinContent(pv->NBINS+1,0.);

  std::cout<<"Fill SM WWA Histogram..."<<std::endl;
  treewwa->Draw(TString(pv->plotvar)+TString(">>th1wwa"), the_cut, "goff");
  th1wwa->AddBinContent(pv->NBINS,th1wwa->GetBinContent(pv->NBINS+1));th1wwa->SetBinContent(pv->NBINS+1,0.);

  std::cout<<"Fill SM WZA Histogram..."<<std::endl;
  treewza->Draw(TString(pv->plotvar)+TString(">>th1wza"), the_cut, "goff");
  th1wza->AddBinContent(pv->NBINS,th1wza->GetBinContent(pv->NBINS+1));th1wza->SetBinContent(pv->NBINS+1,0.);

  std::cout<<"Fill WA+Jets Histogram..."<<std::endl;
  treewaj->Draw(TString(pv->plotvar)+TString(">>th1wajets"), the_cut_LPT, "goff");
  th1wajets->AddBinContent(pv->NBINS,th1wajets->GetBinContent(pv->NBINS+1));th1wajets->SetBinContent(pv->NBINS+1,0.);

  std::cout<<"Fill WA+Jets Hight PT Histogram..."<<std::endl;
  treewaj_HPT->Draw(TString(pv->plotvar)+TString(">>th1wajets_HPT"), the_cut_HPT, "goff");
  th1wajets_HPT->AddBinContent(pv->NBINS,th1wajets_HPT->GetBinContent(pv->NBINS+1));th1wajets_HPT->SetBinContent(pv->NBINS+1,0.);

  std::cout<<"Fill ZA+Jets Histogram..."<<std::endl;
  treezaj->Draw(TString(pv->plotvar)+TString(">>th1zajets"), the_cut, "goff");
  th1zajets->AddBinContent(pv->NBINS,th1zajets->GetBinContent(pv->NBINS+1));th1zajets->SetBinContent(pv->NBINS+1,0.);

  std::cout<<"Fill TTA+Jets Histogram..."<<std::endl;
  treettaj->Draw(TString(pv->plotvar)+TString(">>th1ttajets"), the_cut, "goff");
  th1ttajets->AddBinContent(pv->NBINS,th1ttajets->GetBinContent(pv->NBINS+1));th1ttajets->SetBinContent(pv->NBINS+1,0.);

  std::cout<<"Fill TA+Jets Histogram..."<<std::endl;
  treetaj->Draw(TString(pv->plotvar)+TString(">>th1tajets"), the_cut, "goff");
  th1tajets->AddBinContent(pv->NBINS,th1tajets->GetBinContent(pv->NBINS+1));th1tajets->SetBinContent(pv->NBINS+1,0.);

  std::cout<<"Fill Fake Photons Histogram..."<<std::endl;
  treedata->Draw(TString(pv->plotvar)+TString(">>th1FakePhotons_"), the_cut_FP, "goff");
/*
  treewwa->Draw(TString(pv->plotvar)+TString(">>th1FakePhotons_"), the_cut_FP*Form("%f",WWA_scale), "goff");
  treewza->Draw(TString(pv->plotvar)+TString(">>th1FakePhotons_1"), the_cut_FP*Form("%f",WZA_scale), "goff");
  treewaj->Draw(TString(pv->plotvar)+TString(">>th1FakePhotons_2"), the_cut_FP_LPT*Form("%f",WAJets_scale), "goff");
  treewaj_HPT->Draw(TString(pv->plotvar)+TString(">>th1FakePhotons_3"), the_cut_FP_HPT*Form("%f",WAJets_HPT_scale), "goff");
  treezaj->Draw(TString(pv->plotvar)+TString(">>th1FakePhotons_4"), the_cut_FP*Form("%f",ZAJets_scale), "goff");
  treettaj->Draw(TString(pv->plotvar)+TString(">>th1FakePhotons_5"), the_cut_FP*Form("%f",TTAJets_scale), "goff");
  treetaj->Draw(TString(pv->plotvar)+TString(">>th1FakePhotons_6"), the_cut_FP*Form("%f",TAJets_scale), "goff");
  th1FakePhotons_->Add(th1FakePhotons_1,1); th1FakePhotons_->Add(th1FakePhotons_2,1);
  th1FakePhotons_->Add(th1FakePhotons_3,1); th1FakePhotons_->Add(th1FakePhotons_4,1);
  th1FakePhotons_->Add(th1FakePhotons_5,1); th1FakePhotons_->Add(th1FakePhotons_6,1);
*/
  th1FakePhotons_->AddBinContent(pv->NBINS,th1FakePhotons_->GetBinContent(pv->NBINS+1));
  th1FakePhotons_->SetBinContent(pv->NBINS+1,0.);
  TH1D* th1FakePhotons = (TH1D*)th1wwa->Clone("th1FakePhotons");
  th1FakePhotons->Reset();
  for (int bin = 0; bin < pv->NBINS; bin++) th1FakePhotons->SetBinContent(bin+1, th1FakePhotons_->GetBinContent(bin+1));
  th1FakePhotons->SetTitle("Jet #rightarrow Photon (Fake Photon)");


  /////////////////////////
  // Normalize each sample:
  std::cout<<"\nScale Histograms..."<<std::endl;
  th1wwa->Scale(WWA_scale);
  th1wza->Scale(WZA_scale);
  th1FakePhotons->Scale(FakePhotons_norm[ch]);
  th1wajets->Scale(WAJets_scale);
  th1wajets_HPT->Scale(WAJets_HPT_scale);
  th1zajets->Scale(ZAJets_scale);
  th1ttajets->Scale(TTAJets_scale);
  th1tajets->Scale(TAJets_scale);

  ///////////////////
  // Combine Samples:
  TH1D * th1WVA = (TH1D*)th1wwa->Clone("th1WVA");
  th1WVA->Add(th1wza,1);
  th1wajets->Add(th1wajets_HPT,1);
  th1wajets->Scale(WAJets_norm[ch]);


  /////////////////////
  // Display integrals:
  std::cout << "\nSample Contribution ("+TString((ch==0)? "Mu":"El")+"):" << std::endl;
  std::cout << "WVA: " << th1WVA->Integral() << std::endl;
  std::cout << "WWA: " << th1wwa->Integral() << std::endl;
  std::cout << "WZA: " << th1wza->Integral() << std::endl;
  std::cout << "WA+Jets: " << th1wajets->Integral() << std::endl;
  std::cout << "ZA+Jets: " << th1zajets->Integral() << std::endl;
  std::cout << "TTA+Jets: " << th1ttajets->Integral() << std::endl;
  std::cout << "Single Top: " << th1tajets->Integral() << std::endl;
  std::cout << "Fake Photons: " << th1FakePhotons->Integral() << std::endl;
  std::cout << "Data: " << th1data->Integral() << std::endl;

  ///////////////////////////
  // Combine Signal and BKGs:
  TH1D *signal_SM = (TH1D*)th1WVA->Clone("signal_SM");
  TH1D *data_obs = (TH1D*)th1data->Clone("data_obs");


  ///////////////////////////
  // Create 50x Signal Hist:
  th1WVA->Scale(50.);


  ////////////////////////
  // Set histogram labels:
  const double BINWIDTH = ((pv->MAXRange-pv->MINRange)/pv->NBINS);
  char tmpc[100];    sprintf(tmpc,"Events / %.1f GeV",BINWIDTH);
  if (pv->slog==1)    sprintf(tmpc,"Events / %.1f",BINWIDTH);
  if (pv->slog==2)    sprintf(tmpc,"Events / %.2f",BINWIDTH);
  if (pv->slog==3)    sprintf(tmpc,"Events / %.0f GeV",BINWIDTH);
  if (pv->slog==6)    sprintf(tmpc,"Events / %.2f rad",BINWIDTH);

  th1wwa->SetYTitle(tmpc);
  th1wwa->GetXaxis()->SetTitle(pv->xlabel);
  th1wwa->GetYaxis()->CenterTitle(true);
  th1wza->SetYTitle(tmpc);
  th1wza->GetXaxis()->SetTitle(pv->xlabel);
  th1wza->GetYaxis()->CenterTitle(true);
  th1wajets->SetYTitle(tmpc);
  th1wajets->GetXaxis()->SetTitle(pv->xlabel);
  th1wajets->GetYaxis()->CenterTitle(true);
  th1zajets->SetYTitle(tmpc);
  th1zajets->GetXaxis()->SetTitle(pv->xlabel);
  th1zajets->GetYaxis()->CenterTitle(true);
  th1ttajets->SetYTitle(tmpc);
  th1ttajets->GetXaxis()->SetTitle(pv->xlabel);
  th1ttajets->GetYaxis()->CenterTitle(true);
  th1tajets->SetYTitle(tmpc);
  th1tajets->GetXaxis()->SetTitle(pv->xlabel);
  th1tajets->GetYaxis()->CenterTitle(true);

  th1FakePhotons->SetYTitle(tmpc);
  th1FakePhotons->GetXaxis()->SetTitle(pv->xlabel);
  th1FakePhotons->GetYaxis()->CenterTitle(true);

  signal_SM->SetYTitle(tmpc);
  signal_SM->GetXaxis()->SetTitle(pv->xlabel);
  signal_SM->GetYaxis()->CenterTitle(true);
  data_obs->SetYTitle(tmpc);
  data_obs->GetXaxis()->SetTitle(pv->xlabel);
  data_obs->GetYaxis()->CenterTitle(true);


  // Create Stacked Histogram:
  th1wwa->SetFillColor(kGreen);
  th1wwa->SetLineColor(kGreen);
  th1wwa->SetLineWidth(0);
  th1wza->SetFillColor(kGreen);
  th1wza->SetLineColor(kGreen);
  th1wza->SetLineWidth(0);

  th1wajets->SetFillColor(kRed);
  th1wajets->SetLineColor(kRed);
  th1wajets->SetLineWidth(0);

  th1zajets->SetFillColor(kYellow);
  th1zajets->SetLineColor(kYellow);
  th1zajets->SetLineWidth(0);

  th1ttajets->SetFillColor(kCyan);
  th1ttajets->SetLineColor(kCyan);
  th1ttajets->SetLineWidth(0);

  th1tajets->SetFillColor(kMagenta);
  th1tajets->SetLineColor(kMagenta);
  th1tajets->SetLineWidth(0);

  th1FakePhotons->SetFillColor(kBlue);
  th1FakePhotons->SetLineColor(kBlue);
  th1FakePhotons->SetLineWidth(0);

  TCanvas* c1 = new TCanvas(pv->plotvar,pv->plotvar,10,10, 800, 800);
  TPad * TopPad = new TPad("TopPad", "", 0., 0.2, 1., 1., 21);
  TPad * BottomPad = new TPad("BottomPad", "", 0., 0., 1., 0.2, 22);
  TopPad->SetFillColor(0);
  TopPad->SetFrameBorderMode(0);
  TopPad->SetBorderMode(0);
  TopPad->SetBorderSize(0);
  TopPad->SetBottomMargin(0);
  BottomPad->SetFillColor(0);
  BottomPad->SetFrameBorderMode(0);
  BottomPad->SetBorderMode(0);
  BottomPad->SetBorderSize(0);
  BottomPad->SetTopMargin(0);
  BottomPad->SetBottomMargin(.3);
  TopPad->Draw();
  BottomPad->Draw();
  TopPad->cd();
  gPad->SetLogy();
  gPad->SetTickx(1);
  gPad->SetTicky(1);

  THStack * hs = new THStack("hs",";"+TString(pv->xlabel)+";"+tmpc);
  hs->Add(th1wwa);
  hs->Add(th1wza);
  hs->Add(th1tajets);
  hs->Add(th1ttajets);
  hs->Add(th1zajets);
  hs->Add(th1wajets);
  hs->Add(th1FakePhotons);

  double maxY = th1data->GetBinContent(th1data->GetMaximumBin());
  double minY = 0.001;
  if (hs->GetMaximum() > maxY) maxY = hs->GetMaximum();
  hs->SetMinimum(minY);
  hs->SetMaximum(maxY*1.1);
  hs->Draw("HIST");
  hs->GetHistogram()->GetYaxis()->CenterTitle();
  hs->GetHistogram()->GetYaxis()->SetTitleSize(0.03);
  hs->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
  hs->GetHistogram()->GetXaxis()->SetLabelOffset(999);
  hs->GetHistogram()->GetXaxis()->SetLabelSize(0);
  hs->GetHistogram()->GetXaxis()->SetTitleOffset(999);
  hs->GetHistogram()->GetXaxis()->SetTitleSize(0);
  cmspre(13, intLUMI[ch]);

  // Add in 50x Signal Hist:
  th1WVA->SetLineStyle(2);
  th1WVA->SetLineColor(kBlack);
  th1WVA->SetLineWidth(2);
  th1WVA->Draw("SAME HIST");

  th1data->SetMarkerStyle(20);
  th1data->Draw("esame");

  TLegend * leg = new TLegend(0.65, 0.65, 0.85, 0.85);
  leg->AddEntry(th1WVA, "WV#gamma x 50", "l");
  leg->AddEntry(th1wwa, "WV#gamma", "f");
  leg->AddEntry(th1tajets, "Single Top+#gamma+Jets", "f");
  leg->AddEntry(th1ttajets, "t#bar{t}#gamma+Jets", "f");
  leg->AddEntry(th1zajets, "Z#gamma+Jets", "f");
  leg->AddEntry(th1wajets, "W#gamma+Jets", "f");
  leg->AddEntry(th1FakePhotons, "Jet #rightarrow #gamma", "f");
  leg->AddEntry(th1data, TString((ch==0)? "Mu":"El")+" Data", "lep");
  leg->Draw();

  gPad->RedrawAxis();

  BottomPad->cd();
  gPad->SetTickx(1);
  gPad->SetTicky(1);

  TH1D * Ratio = (TH1D*)th1data->Clone("Ratio");
  Ratio->Reset();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  Ratio->SetYTitle("Data/MC");
  Ratio->GetYaxis()->SetNdivisions(5);
  Ratio->GetYaxis()->CenterTitle();
  Ratio->GetYaxis()->SetTitleOffset(0.35);
  Ratio->GetYaxis()->SetTitleSize(0.12);
  Ratio->GetYaxis()->SetLabelSize(0.12);
  Ratio->GetYaxis()->SetLabelColor(kRed);
  Ratio->GetXaxis()->SetTitle(pv->xlabel);
  Ratio->GetXaxis()->CenterTitle();
  Ratio->GetXaxis()->SetTitleSize(0.12);
  Ratio->GetXaxis()->SetLabelSize(0.12);
  Ratio->Add(th1data,1);
  Ratio->Divide((TH1D*)hs->GetStack()->Last());
  minY = Ratio->GetMinimum();
  maxY = Ratio->GetMaximum();
  Ratio->SetMinimum(0.); Ratio->SetMaximum(2.0);
  //if (fabs(1-minY) > fabs(1-maxY)){Ratio->SetMinimum(1-1.1*fabs(1-minY)); Ratio->SetMaximum(1+1.1*fabs(1-minY));}
  //else {Ratio->SetMinimum(1-1.1*fabs(1-maxY)); Ratio->SetMaximum(1+1.1*fabs(1-maxY));}
  Ratio->Draw();

  TLine * unity = new TLine(pv->MINRange,1,pv->MAXRange,1);
  unity->SetLineWidth(3);
  unity->SetLineColor(kYellow);
  unity->Draw();
  Ratio->Draw("SAME");
  gPad->RedrawAxis();

  TString outfile = TString((ch==0)? "mu":"el")+"_"+TString(pv->plotvar);
  c1->Print("Plots/"+outfile+".png");
  c1->SaveAs("Plots/"+outfile+".pdf");
  //c1->Modified();
  //c1->Update();

  //////////////////////////////////////////////////////////
  // Save Observed Data, Background (+ SM WWA), Background +
  // Uncertainty, Background - Uncertainty, Anomalous Signal
  // histograms to output ROOT file:
  std::cout<<"Save Histograms..."<<std::endl;
  chDir[ch]->cd();
  chSubDir[plot] = chDir[ch]->mkdir(pv->plotvar);
  chSubDir[plot]->cd();
  th1wwa->Write();
  th1wza->Write();
  th1wajets->Write();
  th1zajets->Write();
  th1ttajets->Write();
  th1tajets->Write();
  th1FakePhotons->Write();
  signal_SM->Write();
  data_obs->Write();
  f.cd();


  /////////////////////////////////////
  // Free memory for next plot variable
  delete th1data;
  delete th1wwa;
  delete th1wza;
  delete th1WVA;
  delete th1wajets;
  delete th1wajets_HPT;
  delete th1zajets;
  delete th1ttajets;
  delete th1tajets;
  delete th1FakePhotons;
  delete th1FakePhotons_;
  delete th1FakePhotons_1;
  delete th1FakePhotons_2;
  delete th1FakePhotons_3;
  delete th1FakePhotons_4;
  delete th1FakePhotons_5;
  delete th1FakePhotons_6;
  delete signal_SM;
  delete data_obs;
  delete c1;
  delete Ratio;


  }// End loop over plot variables


  }// End loop over channels


  ////////////////////
  // Close ROOT files:
  for (int ch = 0; ch != 2; ch++) {
    wwaShape_file[ch]->Close();
    wzaShape_file[ch]->Close();
    wajetsShape_file[ch]->Close();
    wajetsHPTShape_file[ch]->Close();
    zajets_file[ch]->Close();
    ttajets_file[ch]->Close();
    tajets_file[ch]->Close();
    data_file[ch]->Close();
  }
  f.Close();


}// End Main function
