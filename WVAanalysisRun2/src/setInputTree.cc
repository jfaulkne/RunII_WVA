#include "../interface/setInputTree.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

setInputTree::setInputTree(TFile* inputFile, std::string inputTreeName)
{
  if(inputFile == 0) {
    TFile* f = new TFile("/gwteray/users/brianza/WWNtupleRun2/ReducedTree/ReducedSelection_TTbar.root");
    fChain = (TTree*) f -> Get("WJet");
  }
  else fChain = (TTree*) inputFile -> Get(inputTreeName.c_str());
  Init();
}

setInputTree::~setInputTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t setInputTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain)  return 0; 
   return fChain->GetEntry(entry);
}
Long64_t setInputTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void setInputTree::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
  //   if (!tree) return;
  //   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
   fChain->SetBranchAddress("LumiBlockNum", &LumiBlockNum, &b_LumiBlockNum);
   fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("BTags", &BTags, &b_BTags);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
   fChain->SetBranchAddress("npT", &npT, &b_npT);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("genEventWeight", &genEventWeight, &b_genEventWeight);
   fChain->SetBranchAddress("originalWeight", &originalWeight, &b_originalWeight);
   fChain->SetBranchAddress("GenEventInfoAQGCweights", &AQGCweights, &b_AQGCweights);
   fChain->SetBranchAddress("passFilterHBHE", &passFilterHBHE, &b_passFilterHBHE);
   fChain->SetBranchAddress("passFilterHBHEIso", &passFilterHBHEIso, &b_passFilterHBHEIso);
   fChain->SetBranchAddress("passFilterCSCHalo", &passFilterCSCHalo, &b_passFilterCSCHalo);
   fChain->SetBranchAddress("passFilterGoodVtx", &passFilterGoodVtx, &b_passFilterGoodVtx);
   fChain->SetBranchAddress("passFilterEEBadSC", &passFilterEEBadSC, &b_passFilterEEBadSC);
   fChain->SetBranchAddress("passFilterHBHELooseRerun", &passFilterHBHELooseRerun, &b_passFilterHBHELooseRerun);
   fChain->SetBranchAddress("passFilterHBHETightRerun", &passFilterHBHETightRerun, &b_passFilterHBHETightRerun);
   fChain->SetBranchAddress("passFilterHBHEIsoRerun", &passFilterHBHEIsoRerun, &b_passFilterHBHEIsoRerun);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("METPt", &METPt, &b_METPt);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("METPtRaw", &METPtRaw, &b_METPtRaw);
   fChain->SetBranchAddress("METPhiRaw", &METPhiRaw, &b_METPhiRaw);
   fChain->SetBranchAddress("METPtDefault", &METPtDefault, &b_METPtDefault);
   fChain->SetBranchAddress("METPhiDefault", &METPhiDefault, &b_METPhiDefault);
   fChain->SetBranchAddress("METPtType1", &METPtType1, &b_METPtType1);
   fChain->SetBranchAddress("METPhiType1", &METPhiType1, &b_METPhiType1);
   fChain->SetBranchAddress("METPtType1XYSmear", &METPtType1XYSmear, &b_METPtType1XYSmear);
   fChain->SetBranchAddress("METPhiType1XYSmear", &METPhiType1XYSmear, &b_METPhiType1XYSmear);
   fChain->SetBranchAddress("METPtType1Smear", &METPtType1Smear, &b_METPtType1Smear);
   fChain->SetBranchAddress("METPhiType1Smear", &METPhiType1Smear, &b_METPhiType1Smear);
   fChain->SetBranchAddress("METPtType1XY", &METPtType1XY, &b_METPtType1XY);
   fChain->SetBranchAddress("METPhiType1XY", &METPhiType1XY, &b_METPhiType1XY);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("DeltaPhi1", &DeltaPhi1, &b_DeltaPhi1);
   fChain->SetBranchAddress("DeltaPhi2", &DeltaPhi2, &b_DeltaPhi2);
   fChain->SetBranchAddress("DeltaPhi3", &DeltaPhi3, &b_DeltaPhi3);
   fChain->SetBranchAddress("IsolatedTracksNum", &IsolatedTracksNum, &b_IsolatedTracksNum);
   fChain->SetBranchAddress("IsolatedTracksPt", IsolatedTracksPt, &b_IsolatedTracksPt);
   fChain->SetBranchAddress("IsolatedTracksEta", IsolatedTracksEta, &b_IsolatedTracksEta);
   fChain->SetBranchAddress("IsolatedTracksPhi", IsolatedTracksPhi, &b_IsolatedTracksPhi);
   fChain->SetBranchAddress("IsolatedTracksE", IsolatedTracksE, &b_IsolatedTracksE);
   fChain->SetBranchAddress("IsolatedTracksTLorentzVector", IsolatedTracksTLorentzVector, &b_IsolatedTracksTLorentzVector);
   fChain->SetBranchAddress("GenPhotonNum", &GenPhotonNum, &b_GenPhotonNum);
   fChain->SetBranchAddress("GenPhotonPt", GenPhotonPt, &b_GenPhotonPt);
   fChain->SetBranchAddress("GenPhotonEta", GenPhotonEta, &b_GenPhotonEta);
   fChain->SetBranchAddress("GenPhotonPhi", GenPhotonPhi, &b_GenPhotonPhi);
   fChain->SetBranchAddress("GenPhotonE", GenPhotonE, &b_GenPhotonE);
   fChain->SetBranchAddress("GenPhotonTLorentzVector", GenPhotonTLorentzVector, &b_GenPhotonTLorentzVector);
   fChain->SetBranchAddress("GenPhoton_GenPhotonIsPromt", GenPhoton_GenPhotonIsPrompt, &b_GenPhoton_GenPhotonIsPrompt);
   fChain->SetBranchAddress("GenPhoton_GenPhotonStatus", GenPhoton_GenPhotonStatus, &b_GenPhoton_GenPhotonStatus);
   fChain->SetBranchAddress("GenBosonNum", &GenBosonNum, &b_GenBosonNum);
   fChain->SetBranchAddress("GenBosonPt", GenBosonPt, &b_GenBosonPt);
   fChain->SetBranchAddress("GenBosonEta", GenBosonEta, &b_GenBosonEta);
   fChain->SetBranchAddress("GenBosonPhi", GenBosonPhi, &b_GenBosonPhi);
   fChain->SetBranchAddress("GenBosonE", GenBosonE, &b_GenBosonE);
   fChain->SetBranchAddress("GenBosonTLorentzVector", GenBosonTLorentzVector, &b_GenBosonTLorentzVector);
   fChain->SetBranchAddress("GenBoson_GenBosonPDGId", GenBoson_GenBosonPDGId, &b_GenBoson_GenBosonPDGId);
   fChain->SetBranchAddress("GenMuNum", &GenMuNum, &b_GenMuNum);
   fChain->SetBranchAddress("GenMuPt", GenMuPt, &b_GenMuPt);
   fChain->SetBranchAddress("GenMuEta", GenMuEta, &b_GenMuEta);
   fChain->SetBranchAddress("GenMuPhi", GenMuPhi, &b_GenMuPhi);
   fChain->SetBranchAddress("GenMuE", GenMuE, &b_GenMuE);
   fChain->SetBranchAddress("GenMuTLorentzVector", GenMuTLorentzVector, &b_GenMuTLorentzVector);
   fChain->SetBranchAddress("GenMu_GenMuFromTau", GenMu_GenMuFromTau, &b_GenMu_GenMuFromTau);
   fChain->SetBranchAddress("GenElecNum", &GenElecNum, &b_GenElecNum);
   fChain->SetBranchAddress("GenElecPt", GenElecPt, &b_GenElecPt);
   fChain->SetBranchAddress("GenElecEta", GenElecEta, &b_GenElecEta);
   fChain->SetBranchAddress("GenElecPhi", GenElecPhi, &b_GenElecPhi);
   fChain->SetBranchAddress("GenElecE", GenElecE, &b_GenElecE);
   fChain->SetBranchAddress("GenElecTLorentzVector", GenElecTLorentzVector, &b_GenElecTLorentzVector);
   fChain->SetBranchAddress("GenElec_GenElecFromTau", GenElec_GenElecFromTau, &b_GenElec_GenElecFromTau);
   fChain->SetBranchAddress("GenTauNum", &GenTauNum, &b_GenTauNum);
   fChain->SetBranchAddress("GenTauPt", GenTauPt, &b_GenTauPt);
   fChain->SetBranchAddress("GenTauEta", GenTauEta, &b_GenTauEta);
   fChain->SetBranchAddress("GenTauPhi", GenTauPhi, &b_GenTauPhi);
   fChain->SetBranchAddress("GenTauE", GenTauE, &b_GenTauE);
   fChain->SetBranchAddress("GenTauTLorentzVector", GenTauTLorentzVector, &b_GenTauTLorentzVector);
   fChain->SetBranchAddress("GenTau_GenTauHad", GenTau_GenTauHad, &b_GenTau_GenTauHad);
   fChain->SetBranchAddress("GenNuNum", &GenNuNum, &b_GenNuNum);
   fChain->SetBranchAddress("GenNuPt", GenNuPt, &b_GenNuPt);
   fChain->SetBranchAddress("GenNuEta", GenNuEta, &b_GenNuEta);
   fChain->SetBranchAddress("GenNuPhi", GenNuPhi, &b_GenNuPhi);
   fChain->SetBranchAddress("GenNuE", GenNuE, &b_GenNuE);
   fChain->SetBranchAddress("GenNuTLorentzVector", GenNuTLorentzVector, &b_GenNuTLorentzVector);
   fChain->SetBranchAddress("JetsNum", &JetsNum, &b_JetsNum);
   fChain->SetBranchAddress("JetsPt", JetsPt, &b_JetsPt);
   fChain->SetBranchAddress("JetsEta", JetsEta, &b_JetsEta);
   fChain->SetBranchAddress("JetsPhi", JetsPhi, &b_JetsPhi);
   fChain->SetBranchAddress("JetsE", JetsE, &b_JetsE);
   fChain->SetBranchAddress("JetsTLorentzVector", JetsTLorentzVector, &b_JetsTLorentzVector);
   fChain->SetBranchAddress("Jets_bDiscriminatorICSV", Jets_bDiscriminatorICSV, &b_Jets_bDiscriminatorICSV);
   fChain->SetBranchAddress("Jets_bDiscriminatorCSV", Jets_bDiscriminatorCSV, &b_Jets_bDiscriminatorCSV);
   fChain->SetBranchAddress("Jets_chargedEmEnergyFraction", Jets_chargedEmEnergyFraction, &b_Jets_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jets_chargedHadronEnergyFraction", Jets_chargedHadronEnergyFraction, &b_Jets_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jets_chargedHadronMultiplicity", Jets_chargedHadronMultiplicity, &b_Jets_chargedHadronMultiplicity);
   fChain->SetBranchAddress("Jets_electronMultiplicity", Jets_electronMultiplicity, &b_Jets_electronMultiplicity);
   fChain->SetBranchAddress("Jets_jetArea", Jets_jetArea, &b_Jets_jetArea);
   fChain->SetBranchAddress("Jets_muonEnergyFraction", Jets_muonEnergyFraction, &b_Jets_muonEnergyFraction);
   fChain->SetBranchAddress("Jets_muonMultiplicity", Jets_muonMultiplicity, &b_Jets_muonMultiplicity);
   fChain->SetBranchAddress("Jets_neutralEmEnergyFraction", Jets_neutralEmEnergyFraction, &b_Jets_neutralEmEnergyFraction);
   fChain->SetBranchAddress("Jets_neutralHadronMultiplicity", Jets_neutralHadronMultiplicity, &b_Jets_neutralHadronMultiplicity);
   fChain->SetBranchAddress("Jets_photonEnergyFraction", Jets_photonEnergyFraction, &b_Jets_photonEnergyFraction);
   fChain->SetBranchAddress("Jets_photonMultiplicity", Jets_photonMultiplicity, &b_Jets_photonMultiplicity);
   fChain->SetBranchAddress("Jets_isLooseJetId", Jets_isLooseJetId, &b_Jets_isLooseJetId);
   fChain->SetBranchAddress("Jets_PtCorr", Jets_PtCorr, &b_Jets_PtCorr);
   fChain->SetBranchAddress("Jets_EtaCorr", Jets_EtaCorr, &b_Jets_EtaCorr);
   fChain->SetBranchAddress("Jets_PhiCorr", Jets_PhiCorr, &b_Jets_PhiCorr);
   fChain->SetBranchAddress("Jets_ECorr", Jets_ECorr, &b_Jets_ECorr);
   fChain->SetBranchAddress("Jets_PUJetID", Jets_PUJetID, &b_Jets_PUJetID);
   fChain->SetBranchAddress("AK8JetsNum", &AK8JetsNum, &b_AK8JetsNum);
   fChain->SetBranchAddress("AK8JetsPt", AK8JetsPt, &b_AK8JetsPt);
   fChain->SetBranchAddress("AK8JetsEta", AK8JetsEta, &b_AK8JetsEta);
   fChain->SetBranchAddress("AK8JetsPhi", AK8JetsPhi, &b_AK8JetsPhi);
   fChain->SetBranchAddress("AK8JetsE", AK8JetsE, &b_AK8JetsE);
   fChain->SetBranchAddress("AK8JetsTLorentzVector", AK8JetsTLorentzVector, &b_AK8JetsTLorentzVector);
   fChain->SetBranchAddress("AK8Jets_bDiscriminatorICSV", AK8Jets_bDiscriminatorICSV, &b_AK8Jets_bDiscriminatorICSV);
   fChain->SetBranchAddress("AK8Jets_bDiscriminatorCSV", AK8Jets_bDiscriminatorCSV, &b_AK8Jets_bDiscriminatorCSV);
   fChain->SetBranchAddress("AK8Jets_chargedEmEnergyFraction", AK8Jets_chargedEmEnergyFraction, &b_AK8Jets_chargedEmEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_chargedHadronEnergyFraction", AK8Jets_chargedHadronEnergyFraction, &b_AK8Jets_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_chargedHadronMultiplicity", AK8Jets_chargedHadronMultiplicity, &b_AK8Jets_chargedHadronMultiplicity);
   fChain->SetBranchAddress("AK8Jets_electronMultiplicity", AK8Jets_electronMultiplicity, &b_AK8Jets_electronMultiplicity);
   fChain->SetBranchAddress("AK8Jets_jetArea", AK8Jets_jetArea, &b_AK8Jets_jetArea);
   fChain->SetBranchAddress("AK8Jets_muonEnergyFraction", AK8Jets_muonEnergyFraction, &b_AK8Jets_muonEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_muonMultiplicity", AK8Jets_muonMultiplicity, &b_AK8Jets_muonMultiplicity);
   fChain->SetBranchAddress("AK8Jets_neutralEmEnergyFraction", AK8Jets_neutralEmEnergyFraction, &b_AK8Jets_neutralEmEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_neutralHadronMultiplicity", AK8Jets_neutralHadronMultiplicity, &b_AK8Jets_neutralHadronMultiplicity);
   fChain->SetBranchAddress("AK8Jets_photonEnergyFraction", AK8Jets_photonEnergyFraction, &b_AK8Jets_photonEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_photonMultiplicity", AK8Jets_photonMultiplicity, &b_AK8Jets_photonMultiplicity);
   fChain->SetBranchAddress("AK8Jets_prunedMass", AK8Jets_prunedMass, &b_AK8Jets_prunedMass);
   fChain->SetBranchAddress("AK8Jets_softDropMass", AK8Jets_softDropMass, &b_AK8Jets_softDropMass);
   fChain->SetBranchAddress("AK8Jets_trimmedMass", AK8Jets_trimmedMass, &b_AK8Jets_trimmedMass);
   fChain->SetBranchAddress("AK8Jets_filteredMass", AK8Jets_filteredMass, &b_AK8Jets_filteredMass);
   fChain->SetBranchAddress("AK8Jets_tau1", AK8Jets_tau1, &b_AK8Jets_tau1);
   fChain->SetBranchAddress("AK8Jets_tau2", AK8Jets_tau2, &b_AK8Jets_tau2);
   fChain->SetBranchAddress("AK8Jets_tau3", AK8Jets_tau3, &b_AK8Jets_tau3);
   fChain->SetBranchAddress("AK8Jets_AK8isLooseJetId", AK8Jets_AK8isLooseJetId, &b_AK8Jets_AK8isLooseJetId);
   fChain->SetBranchAddress("AK8Jets_PtCorr", AK8Jets_PtCorr, &b_AK8Jets_PtCorr);
   fChain->SetBranchAddress("AK8Jets_EtaCorr", AK8Jets_EtaCorr, &b_AK8Jets_EtaCorr);
   fChain->SetBranchAddress("AK8Jets_PhiCorr", AK8Jets_PhiCorr, &b_AK8Jets_PhiCorr);
   fChain->SetBranchAddress("AK8Jets_ECorr", AK8Jets_ECorr, &b_AK8Jets_ECorr);
   fChain->SetBranchAddress("ElectronsNum", &ElectronsNum, &b_ElectronsNum);
   fChain->SetBranchAddress("ElectronsPt", ElectronsPt, &b_ElectronsPt);
   fChain->SetBranchAddress("ElectronsEta", ElectronsEta, &b_ElectronsEta);
   fChain->SetBranchAddress("ElectronsPhi", ElectronsPhi, &b_ElectronsPhi);
   fChain->SetBranchAddress("ElectronsE", ElectronsE, &b_ElectronsE);
   fChain->SetBranchAddress("ElectronsTLorentzVector", ElectronsTLorentzVector, &b_ElectronsTLorentzVector);
   fChain->SetBranchAddress("Electrons_charge", Electrons_charge, &b_Electrons_charge);
   fChain->SetBranchAddress("Electrons_isHEEP", Electrons_isHEEP, &b_Electrons_isHEEP);
   fChain->SetBranchAddress("Electrons_type", Electrons_type, &b_Electrons_type);
   fChain->SetBranchAddress("Electrons_mass", Electrons_mass, &b_Electrons_mass);
   fChain->SetBranchAddress("Electrons_pfDeltaCorrRelIso", Electrons_pfDeltaCorrRelIso, &b_Electrons_pfDeltaCorrRelIso);
   fChain->SetBranchAddress("Electrons_pfRhoCorrRelIso04", Electrons_pfRhoCorrRelIso04, &b_Electrons_pfRhoCorrRelIso04);
   fChain->SetBranchAddress("Electrons_pfRhoCorrRelIso03", Electrons_pfRhoCorrRelIso03, &b_Electrons_pfRhoCorrRelIso03);
   fChain->SetBranchAddress("Electrons_pfRelIso", Electrons_pfRelIso, &b_Electrons_pfRelIso);
   fChain->SetBranchAddress("Electrons_photonIso", Electrons_photonIso, &b_Electrons_photonIso);
   fChain->SetBranchAddress("Electrons_neutralHadIso", Electrons_neutralHadIso, &b_Electrons_neutralHadIso);
   fChain->SetBranchAddress("Electrons_chargedHadIso", Electrons_chargedHadIso, &b_Electrons_chargedHadIso);
   fChain->SetBranchAddress("Electrons_trackIso", Electrons_trackIso, &b_Electrons_trackIso);
   fChain->SetBranchAddress("Electrons_isLoose", Electrons_isLoose, &b_Electrons_isLoose);
   fChain->SetBranchAddress("Electrons_isMedium", Electrons_isMedium, &b_Electrons_isMedium);
   fChain->SetBranchAddress("Electrons_isTight", Electrons_isTight, &b_Electrons_isTight);
   fChain->SetBranchAddress("MuonsNum", &MuonsNum, &b_MuonsNum);
   fChain->SetBranchAddress("MuonsPt", MuonsPt, &b_MuonsPt);
   fChain->SetBranchAddress("MuonsEta", MuonsEta, &b_MuonsEta);
   fChain->SetBranchAddress("MuonsPhi", MuonsPhi, &b_MuonsPhi);
   fChain->SetBranchAddress("MuonsE", MuonsE, &b_MuonsE);
   fChain->SetBranchAddress("MuonsTLorentzVector", MuonsTLorentzVector, &b_MuonsTLorentzVector);
   fChain->SetBranchAddress("Muons_charge", Muons_charge, &b_Muons_charge);
   fChain->SetBranchAddress("Muons_isHighPt", Muons_isHighPt, &b_Muons_isHighPt);
   fChain->SetBranchAddress("Muons_type", Muons_type, &b_Muons_type);
   fChain->SetBranchAddress("Muons_mass", Muons_mass, &b_Muons_mass);
   fChain->SetBranchAddress("Muons_pfDeltaCorrRelIso", Muons_pfDeltaCorrRelIso, &b_Muons_pfDeltaCorrRelIso);
   fChain->SetBranchAddress("Muons_pfRelIso", Muons_pfRelIso, &b_Muons_pfRelIso);
   fChain->SetBranchAddress("Muons_photonIso", Muons_photonIso, &b_Muons_photonIso);
   fChain->SetBranchAddress("Muons_neutralHadIso", Muons_neutralHadIso, &b_Muons_neutralHadIso);
   fChain->SetBranchAddress("Muons_chargedHadIso", Muons_chargedHadIso, &b_Muons_chargedHadIso);
   fChain->SetBranchAddress("Muons_trackIso", Muons_trackIso, &b_Muons_trackIso);
   fChain->SetBranchAddress("Muons_isLoose", Muons_isLoose, &b_Muons_isLoose);
   fChain->SetBranchAddress("Muons_isMedium", Muons_isMedium, &b_Muons_isMedium);
   fChain->SetBranchAddress("Muons_isTight", Muons_isTight, &b_Muons_isTight);
   fChain->SetBranchAddress("Muons_isPFMuon", Muons_isPFMuon, &b_Muons_isPFMuon);

   fChain->SetBranchAddress("PhotonsNum", &PhotonsNum, &b_PhotonsNum);
   fChain->SetBranchAddress("PhotonsPt", PhotonsPt, &b_PhotonsPt);
   fChain->SetBranchAddress("PhotonsEta", PhotonsEta, &b_PhotonsEta);
   fChain->SetBranchAddress("PhotonsPhi", PhotonsPhi, &b_PhotonsPhi);
   fChain->SetBranchAddress("PhotonsE", PhotonsE, &b_PhotonsE);
   fChain->SetBranchAddress("PhotonsTLorentzVector", PhotonsTLorentzVector, &b_PhotonsTLorentzVector);
   fChain->SetBranchAddress("Photons_isLoose", Photons_isLoose, &b_Photons_isLoose);
   fChain->SetBranchAddress("Photons_isMedium", Photons_isMedium, &b_Photons_isMedium);
   fChain->SetBranchAddress("Photons_isTight", Photons_isTight, &b_Photons_isTight);
   fChain->SetBranchAddress("Photons_minPt", Photons_minPt, &b_Photons_minPt);
   fChain->SetBranchAddress("Photons_phoSCEtaMultiRange", Photons_phoSCEtaMultiRange, &b_Photons_phoSCEtaMultiRange);
   fChain->SetBranchAddress("Photons_phoSingleTowerHadOverEm", Photons_phoSingleTowerHadOverEm, &b_Photons_phoSingleTowerHadOverEm);
   fChain->SetBranchAddress("Photons_phoFull5x5SigmaIEtaIEta", Photons_phoFull5x5SigmaIEtaIEta, &b_Photons_phoFull5x5SigmaIEtaIEta);
   fChain->SetBranchAddress("Photons_phoAnyPFIsoWithEA", Photons_phoAnyPFIsoWithEA, &b_Photons_phoAnyPFIsoWithEA);
   fChain->SetBranchAddress("Photons_phoAnyPFIsoWithEAAndExpoScaling", Photons_phoAnyPFIsoWithEAAndExpoScaling, &b_Photons_phoAnyPFIsoWithEAAndExpoScaling);
   fChain->SetBranchAddress("Photons_phoAnyPFIsoWithEA1", Photons_phoAnyPFIsoWithEA1, &b_Photons_phoAnyPFIsoWithEA1);
   fChain->SetBranchAddress("Photons_PassMinPt", Photons_PassMinPt, &b_Photons_PassMinPt);
   fChain->SetBranchAddress("Photons_PassPhoSCEtaMultiRange", Photons_PassPhoSCEtaMultiRange, &b_Photons_PassPhoSCEtaMultiRange);
   fChain->SetBranchAddress("Photons_PassPhoSingleTowerHadOverEm", Photons_PassPhoSingleTowerHadOverEm, &b_Photons_PassPhoSingleTowerHadOverEm);
   fChain->SetBranchAddress("Photons_PassPhoFull5x5SigmaIEtaIEta", Photons_PassPhoFull5x5SigmaIEtaIEta, &b_Photons_PassPhoFull5x5SigmaIEtaIEta);
   fChain->SetBranchAddress("Photons_PassPhoAnyPFIsoWithEA", Photons_PassPhoAnyPFIsoWithEA, &b_Photons_PassPhoAnyPFIsoWithEA);
   fChain->SetBranchAddress("Photons_PassPhoAnyPFIsoWithEAAndExpoScaling", Photons_PassPhoAnyPFIsoWithEAAndExpoScaling, &b_Photons_PassPhoAnyPFIsoWithEAAndExpoScaling);
   fChain->SetBranchAddress("Photons_PassPhoAnyPFIsoWithEA1", Photons_PassPhoAnyPFIsoWithEA1, &b_Photons_PassPhoAnyPFIsoWithEA1);
   fChain->SetBranchAddress("Photons_hasPixelSeed", Photons_hasPixelSeed, &b_Photons_hasPixelSeed);
   fChain->SetBranchAddress("Photons_passElectronVeto", Photons_passElectronVeto, &b_Photons_passElectronVeto);
   fChain->SetBranchAddress("Photons_photonIso", Photons_photonIso, &b_Photons_photonIso);
   fChain->SetBranchAddress("Photons_neutralHadIso", Photons_neutralHadIso, &b_Photons_neutralHadIso);
   fChain->SetBranchAddress("Photons_chargedHadIso", Photons_chargedHadIso, &b_Photons_chargedHadIso);
   fChain->SetBranchAddress("Photons_puChargedHadIso", Photons_puChargedHadIso, &b_Photons_puChargedHadIso);
   fChain->SetBranchAddress("Photons_sigmaIetaIeta", Photons_sigmaIetaIeta, &b_Photons_sigmaIetaIeta);

   fChain->SetBranchAddress("TriggerProducerTriggerPrescales", &TriggerProducerTriggerPrescales, &b_TriggerProducerTriggerPrescales);
   fChain->SetBranchAddress("TriggerProducerTriggerPass", &TriggerProducerTriggerPass, &b_TriggerProducerTriggerPass);
   fChain->SetBranchAddress("TriggerProducerTriggerNames", &TriggerProducerTriggerNames, &b_TriggerProducerTriggerNames);
   //   Notify();
}

Bool_t setInputTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void setInputTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t setInputTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
