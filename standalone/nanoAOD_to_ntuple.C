#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

// Structure to hold GenParticle information
struct GenParticle {
    float pt, eta, phi, mass;
    int pdgId;
    int statusFlags;
    std::vector<int> dauIdx;
    GenParticle* genW = nullptr;
    GenParticle* genB = nullptr;
    std::vector<GenParticle*> daus;
    
    TLorentzVector p4() const {
        TLorentzVector vec;
        vec.SetPtEtaPhiM(pt, eta, phi, mass);
        return vec;
    }
};

// Structure to hold FatJet with gen-matching info
struct FatJetInfo {
    float pt, eta, phi, mass;
    GenParticle* genH = nullptr;
    GenParticle* genZ = nullptr;
    GenParticle* genW = nullptr;
    GenParticle* genT = nullptr;
    GenParticle* genLepT = nullptr;
    float dr_H = 99;
    float dr_Z = 99;
    float dr_W = 99;
    float dr_T = 99;
    float dr_LepT = 99;
    
    TLorentzVector p4() const {
        TLorentzVector vec;
        vec.SetPtEtaPhiM(pt, eta, phi, mass);
        return vec;
    }
};

// Calculate deltaR
float deltaR(const TLorentzVector& v1, const TLorentzVector& v2) {
    float deta = v1.Eta() - v2.Eta();
    float dphi = TVector2::Phi_mpi_pi(v1.Phi() - v2.Phi());
    return sqrt(deta*deta + dphi*dphi);
}

float deltaR(const FatJetInfo& fj, const GenParticle* gp) {
    if (!gp) return 99;
    return deltaR(fj.p4(), gp->p4());
}

float deltaR(const FatJetInfo& fj, const GenParticle& gp) {
    return deltaR(fj.p4(), gp.p4());
}

// Check if particle decays hadronically
bool isHadronic(const GenParticle* gp, const std::vector<GenParticle>& genparts) {
    if (!gp || gp->dauIdx.size() < 2) return false;
    for (int idx : gp->dauIdx) {
        int pdgId = abs(genparts[idx].pdgId);
        // Check if any daughter is a quark (1-5)
        if (pdgId >= 1 && pdgId <= 5) return true;
    }
    return false;
}

// Get final state particle (follow decay chain)
GenParticle* getFinal(GenParticle* gp, std::vector<GenParticle>& genparts) {
    if (!gp) return nullptr;
    if (gp->dauIdx.empty()) return gp;
    
    // Look for same-flavor daughter
    for (int idx : gp->dauIdx) {
        if (abs(genparts[idx].pdgId) == abs(gp->pdgId)) {
            return getFinal(&genparts[idx], genparts);
        }
    }
    return gp;
}

// Find closest particle
std::pair<GenParticle*, float> closest(const FatJetInfo& fj, const std::vector<GenParticle*>& particles) {
    float minDR = 999;
    GenParticle* closest_gp = nullptr;
    
    for (auto* gp : particles) {
        if (!gp) continue;
        float dr = deltaR(fj, gp);
        if (dr < minDR) {
            minDR = dr;
            closest_gp = gp;
        }
    }
    
    return {closest_gp, minDR};
}

void nanoAOD_to_ntuple(const char* inputFile, const char* outputFile) {
    
    // Prepend XRootD redirector if input file starts with "/store"
    std::string inputPath(inputFile);
    if (inputPath.substr(0, 7) == "/store/") {
        inputPath = "root://xrootd-cms.infn.it/" + inputPath;
        std::cout << "Using XRootD path: " << inputPath << std::endl;
    }
    
    // Open input file
    TFile *inFile = TFile::Open(inputPath.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return;
    }
    
    TTree *inTree = (TTree*)inFile->Get("Events");
    if (!inTree) {
        std::cerr << "Cannot find Events tree in input file" << std::endl;
        inFile->Close();
        return;
    }
    
    // Set up input branches for FatJet
    Int_t nFatJet;
    Float_t FatJet_pt[100], FatJet_eta[100], FatJet_phi[100], FatJet_mass[100];
    Float_t FatJet_msoftdrop[100];
    Float_t FatJet_particleNetLegacy_QCD[100], FatJet_particleNetLegacy_Xbb[100];
    Float_t FatJet_particleNetLegacy_Xcc[100], FatJet_particleNetLegacy_Xqq[100];
    Float_t FatJet_particleNetLegacy_mass[100];
    Float_t FatJet_particleNet_QCD[100], FatJet_particleNet_QCD0HF[100];
    Float_t FatJet_particleNet_QCD1HF[100], FatJet_particleNet_QCD2HF[100];
    Float_t FatJet_particleNet_WVsQCD[100], FatJet_particleNet_XbbVsQCD[100];
    Float_t FatJet_particleNet_XccVsQCD[100], FatJet_particleNet_XggVsQCD[100];
    Float_t FatJet_particleNet_XqqVsQCD[100], FatJet_particleNet_XteVsQCD[100];
    Float_t FatJet_particleNet_XtmVsQCD[100], FatJet_particleNet_XttVsQCD[100];
    Float_t FatJet_particleNet_massCorr[100];
    Float_t FatJet_particleNetWithMass_H4qvsQCD[100], FatJet_particleNetWithMass_HbbvsQCD[100];
    Float_t FatJet_particleNetWithMass_HccvsQCD[100], FatJet_particleNetWithMass_QCD[100];
    Float_t FatJet_particleNetWithMass_TvsQCD[100], FatJet_particleNetWithMass_WvsQCD[100];
    Float_t FatJet_particleNetWithMass_ZvsQCD[100];
    Float_t FatJet_globalParT3_QCD[100], FatJet_globalParT3_TopbWev[100];
    Float_t FatJet_globalParT3_TopbWmv[100], FatJet_globalParT3_TopbWq[100];
    Float_t FatJet_globalParT3_TopbWqq[100], FatJet_globalParT3_TopbWtauhv[100];
    Float_t FatJet_globalParT3_WvsQCD[100], FatJet_globalParT3_XWW3q[100];
    Float_t FatJet_globalParT3_XWW4q[100], FatJet_globalParT3_XWWqqev[100];
    Float_t FatJet_globalParT3_XWWqqmv[100], FatJet_globalParT3_Xbb[100];
    Float_t FatJet_globalParT3_Xcc[100], FatJet_globalParT3_Xcs[100];
    Float_t FatJet_globalParT3_Xqq[100], FatJet_globalParT3_Xtauhtaue[100];
    Float_t FatJet_globalParT3_Xtauhtauh[100], FatJet_globalParT3_Xtauhtaum[100];
    Float_t FatJet_globalParT3_massCorrGeneric[100], FatJet_globalParT3_massCorrX2p[100];
    Float_t FatJet_globalParT3_withMassTopvsQCD[100], FatJet_globalParT3_withMassWvsQCD[100];
    Float_t FatJet_globalParT3_withMassZvsQCD[100];
    
    inTree->SetBranchAddress("nFatJet", &nFatJet);
    inTree->SetBranchAddress("FatJet_pt", FatJet_pt);
    inTree->SetBranchAddress("FatJet_eta", FatJet_eta);
    inTree->SetBranchAddress("FatJet_phi", FatJet_phi);
    inTree->SetBranchAddress("FatJet_mass", FatJet_mass);
    inTree->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop);
    inTree->SetBranchAddress("FatJet_particleNetLegacy_QCD", FatJet_particleNetLegacy_QCD);
    inTree->SetBranchAddress("FatJet_particleNetLegacy_Xbb", FatJet_particleNetLegacy_Xbb);
    inTree->SetBranchAddress("FatJet_particleNetLegacy_Xcc", FatJet_particleNetLegacy_Xcc);
    inTree->SetBranchAddress("FatJet_particleNetLegacy_Xqq", FatJet_particleNetLegacy_Xqq);
    inTree->SetBranchAddress("FatJet_particleNetLegacy_mass", FatJet_particleNetLegacy_mass);
    inTree->SetBranchAddress("FatJet_particleNet_QCD", FatJet_particleNet_QCD);
    inTree->SetBranchAddress("FatJet_particleNet_QCD0HF", FatJet_particleNet_QCD0HF);
    inTree->SetBranchAddress("FatJet_particleNet_QCD1HF", FatJet_particleNet_QCD1HF);
    inTree->SetBranchAddress("FatJet_particleNet_QCD2HF", FatJet_particleNet_QCD2HF);
    inTree->SetBranchAddress("FatJet_particleNet_WVsQCD", FatJet_particleNet_WVsQCD);
    inTree->SetBranchAddress("FatJet_particleNet_XbbVsQCD", FatJet_particleNet_XbbVsQCD);
    inTree->SetBranchAddress("FatJet_particleNet_XccVsQCD", FatJet_particleNet_XccVsQCD);
    inTree->SetBranchAddress("FatJet_particleNet_XggVsQCD", FatJet_particleNet_XggVsQCD);
    inTree->SetBranchAddress("FatJet_particleNet_XqqVsQCD", FatJet_particleNet_XqqVsQCD);
    inTree->SetBranchAddress("FatJet_particleNet_XteVsQCD", FatJet_particleNet_XteVsQCD);
    inTree->SetBranchAddress("FatJet_particleNet_XtmVsQCD", FatJet_particleNet_XtmVsQCD);
    inTree->SetBranchAddress("FatJet_particleNet_XttVsQCD", FatJet_particleNet_XttVsQCD);
    inTree->SetBranchAddress("FatJet_particleNet_massCorr", FatJet_particleNet_massCorr);
    inTree->SetBranchAddress("FatJet_particleNetWithMass_H4qvsQCD", FatJet_particleNetWithMass_H4qvsQCD);
    inTree->SetBranchAddress("FatJet_particleNetWithMass_HbbvsQCD", FatJet_particleNetWithMass_HbbvsQCD);
    inTree->SetBranchAddress("FatJet_particleNetWithMass_HccvsQCD", FatJet_particleNetWithMass_HccvsQCD);
    inTree->SetBranchAddress("FatJet_particleNetWithMass_QCD", FatJet_particleNetWithMass_QCD);
    inTree->SetBranchAddress("FatJet_particleNetWithMass_TvsQCD", FatJet_particleNetWithMass_TvsQCD);
    inTree->SetBranchAddress("FatJet_particleNetWithMass_WvsQCD", FatJet_particleNetWithMass_WvsQCD);
    inTree->SetBranchAddress("FatJet_particleNetWithMass_ZvsQCD", FatJet_particleNetWithMass_ZvsQCD);
    inTree->SetBranchAddress("FatJet_globalParT3_QCD", FatJet_globalParT3_QCD);
    inTree->SetBranchAddress("FatJet_globalParT3_TopbWev", FatJet_globalParT3_TopbWev);
    inTree->SetBranchAddress("FatJet_globalParT3_TopbWmv", FatJet_globalParT3_TopbWmv);
    inTree->SetBranchAddress("FatJet_globalParT3_TopbWq", FatJet_globalParT3_TopbWq);
    inTree->SetBranchAddress("FatJet_globalParT3_TopbWqq", FatJet_globalParT3_TopbWqq);
    inTree->SetBranchAddress("FatJet_globalParT3_TopbWtauhv", FatJet_globalParT3_TopbWtauhv);
    inTree->SetBranchAddress("FatJet_globalParT3_WvsQCD", FatJet_globalParT3_WvsQCD);
    inTree->SetBranchAddress("FatJet_globalParT3_XWW3q", FatJet_globalParT3_XWW3q);
    inTree->SetBranchAddress("FatJet_globalParT3_XWW4q", FatJet_globalParT3_XWW4q);
    inTree->SetBranchAddress("FatJet_globalParT3_XWWqqev", FatJet_globalParT3_XWWqqev);
    inTree->SetBranchAddress("FatJet_globalParT3_XWWqqmv", FatJet_globalParT3_XWWqqmv);
    inTree->SetBranchAddress("FatJet_globalParT3_Xbb", FatJet_globalParT3_Xbb);
    inTree->SetBranchAddress("FatJet_globalParT3_Xcc", FatJet_globalParT3_Xcc);
    inTree->SetBranchAddress("FatJet_globalParT3_Xcs", FatJet_globalParT3_Xcs);
    inTree->SetBranchAddress("FatJet_globalParT3_Xqq", FatJet_globalParT3_Xqq);
    inTree->SetBranchAddress("FatJet_globalParT3_Xtauhtaue", FatJet_globalParT3_Xtauhtaue);
    inTree->SetBranchAddress("FatJet_globalParT3_Xtauhtauh", FatJet_globalParT3_Xtauhtauh);
    inTree->SetBranchAddress("FatJet_globalParT3_Xtauhtaum", FatJet_globalParT3_Xtauhtaum);
    inTree->SetBranchAddress("FatJet_globalParT3_massCorrGeneric", FatJet_globalParT3_massCorrGeneric);
    inTree->SetBranchAddress("FatJet_globalParT3_massCorrX2p", FatJet_globalParT3_massCorrX2p);
    inTree->SetBranchAddress("FatJet_globalParT3_withMassTopvsQCD", FatJet_globalParT3_withMassTopvsQCD);
    inTree->SetBranchAddress("FatJet_globalParT3_withMassWvsQCD", FatJet_globalParT3_withMassWvsQCD);
    inTree->SetBranchAddress("FatJet_globalParT3_withMassZvsQCD", FatJet_globalParT3_withMassZvsQCD);
    
    // Set up input branches for GenPart
    Int_t nGenPart;
    Float_t GenPart_pt[500], GenPart_eta[500], GenPart_phi[500], GenPart_mass[500];
    Int_t GenPart_pdgId[500];
    UShort_t GenPart_statusFlags[500];
    Short_t GenPart_genPartIdxMother[500];
    
    inTree->SetBranchAddress("nGenPart", &nGenPart);
    inTree->SetBranchAddress("GenPart_pt", GenPart_pt);
    inTree->SetBranchAddress("GenPart_eta", GenPart_eta);
    inTree->SetBranchAddress("GenPart_phi", GenPart_phi);
    inTree->SetBranchAddress("GenPart_mass", GenPart_mass);
    inTree->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    inTree->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags);
    inTree->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    
    // Create output file and tree
    TFile *outFile = new TFile(outputFile, "RECREATE");
    TTree *outTree = new TTree("FatJets", "FatJet n-tuple with gen-matching");
    
    // Output branches - FatJet basic variables
    Float_t out_pt, out_eta, out_phi, out_mass;
    Float_t out_msoftdrop;
    Float_t out_particleNetLegacy_QCD, out_particleNetLegacy_Xbb;
    Float_t out_particleNetLegacy_Xcc, out_particleNetLegacy_Xqq;
    Float_t out_particleNetLegacy_mass;
    Float_t out_particleNet_QCD, out_particleNet_QCD0HF;
    Float_t out_particleNet_QCD1HF, out_particleNet_QCD2HF;
    Float_t out_particleNet_WVsQCD, out_particleNet_XbbVsQCD;
    Float_t out_particleNet_XccVsQCD, out_particleNet_XggVsQCD;
    Float_t out_particleNet_XqqVsQCD, out_particleNet_XteVsQCD;
    Float_t out_particleNet_XtmVsQCD, out_particleNet_XttVsQCD;
    Float_t out_particleNet_massCorr;
    Float_t out_particleNetWithMass_H4qvsQCD, out_particleNetWithMass_HbbvsQCD;
    Float_t out_particleNetWithMass_HccvsQCD, out_particleNetWithMass_QCD;
    Float_t out_particleNetWithMass_TvsQCD, out_particleNetWithMass_WvsQCD;
    Float_t out_particleNetWithMass_ZvsQCD;
    Float_t out_globalParT3_QCD, out_globalParT3_TopbWev;
    Float_t out_globalParT3_TopbWmv, out_globalParT3_TopbWq;
    Float_t out_globalParT3_TopbWqq, out_globalParT3_TopbWtauhv;
    Float_t out_globalParT3_WvsQCD, out_globalParT3_XWW3q;
    Float_t out_globalParT3_XWW4q, out_globalParT3_XWWqqev;
    Float_t out_globalParT3_XWWqqmv, out_globalParT3_Xbb;
    Float_t out_globalParT3_Xcc, out_globalParT3_Xcs;
    Float_t out_globalParT3_Xqq, out_globalParT3_Xtauhtaue;
    Float_t out_globalParT3_Xtauhtauh, out_globalParT3_Xtauhtaum;
    Float_t out_globalParT3_massCorrGeneric, out_globalParT3_massCorrX2p;
    Float_t out_globalParT3_withMassTopvsQCD, out_globalParT3_withMassWvsQCD;
    Float_t out_globalParT3_withMassZvsQCD;
    
    // Output branches - gen-matching info
    Float_t out_dr_H, out_dr_H_daus, out_H_pt;
    Int_t out_H_decay;
    
    Float_t out_dr_Z, out_dr_Z_daus, out_Z_pt;
    Int_t out_Z_decay;
    
    Float_t out_dr_W, out_dr_W_daus, out_W_pt;
    Int_t out_W_decay;
    
    Float_t out_dr_T, out_dr_T_b, out_dr_T_Wq_max, out_dr_T_Wq_min, out_T_pt;
    Int_t out_T_Wq_max_pdgId, out_T_Wq_min_pdgId;
    
    Float_t out_dr_LepT, out_dr_LepT_b, out_dr_LepT_l, out_LepT_pt;
    Int_t out_dr_LepT_l_pdgId;
    
    // Set output branches
    outTree->Branch("pt", &out_pt);
    outTree->Branch("eta", &out_eta);
    outTree->Branch("phi", &out_phi);
    outTree->Branch("mass", &out_mass);
    outTree->Branch("msoftdrop", &out_msoftdrop);
    outTree->Branch("particleNetLegacy_QCD", &out_particleNetLegacy_QCD);
    outTree->Branch("particleNetLegacy_Xbb", &out_particleNetLegacy_Xbb);
    outTree->Branch("particleNetLegacy_Xcc", &out_particleNetLegacy_Xcc);
    outTree->Branch("particleNetLegacy_Xqq", &out_particleNetLegacy_Xqq);
    outTree->Branch("particleNetLegacy_mass", &out_particleNetLegacy_mass);
    outTree->Branch("particleNet_QCD", &out_particleNet_QCD);
    outTree->Branch("particleNet_QCD0HF", &out_particleNet_QCD0HF);
    outTree->Branch("particleNet_QCD1HF", &out_particleNet_QCD1HF);
    outTree->Branch("particleNet_QCD2HF", &out_particleNet_QCD2HF);
    outTree->Branch("particleNet_WVsQCD", &out_particleNet_WVsQCD);
    outTree->Branch("particleNet_XbbVsQCD", &out_particleNet_XbbVsQCD);
    outTree->Branch("particleNet_XccVsQCD", &out_particleNet_XccVsQCD);
    outTree->Branch("particleNet_XggVsQCD", &out_particleNet_XggVsQCD);
    outTree->Branch("particleNet_XqqVsQCD", &out_particleNet_XqqVsQCD);
    outTree->Branch("particleNet_XteVsQCD", &out_particleNet_XteVsQCD);
    outTree->Branch("particleNet_XtmVsQCD", &out_particleNet_XtmVsQCD);
    outTree->Branch("particleNet_XttVsQCD", &out_particleNet_XttVsQCD);
    outTree->Branch("particleNet_massCorr", &out_particleNet_massCorr);
    outTree->Branch("particleNetWithMass_H4qvsQCD", &out_particleNetWithMass_H4qvsQCD);
    outTree->Branch("particleNetWithMass_HbbvsQCD", &out_particleNetWithMass_HbbvsQCD);
    outTree->Branch("particleNetWithMass_HccvsQCD", &out_particleNetWithMass_HccvsQCD);
    outTree->Branch("particleNetWithMass_QCD", &out_particleNetWithMass_QCD);
    outTree->Branch("particleNetWithMass_TvsQCD", &out_particleNetWithMass_TvsQCD);
    outTree->Branch("particleNetWithMass_WvsQCD", &out_particleNetWithMass_WvsQCD);
    outTree->Branch("particleNetWithMass_ZvsQCD", &out_particleNetWithMass_ZvsQCD);
    outTree->Branch("globalParT3_QCD", &out_globalParT3_QCD);
    outTree->Branch("globalParT3_TopbWev", &out_globalParT3_TopbWev);
    outTree->Branch("globalParT3_TopbWmv", &out_globalParT3_TopbWmv);
    outTree->Branch("globalParT3_TopbWq", &out_globalParT3_TopbWq);
    outTree->Branch("globalParT3_TopbWqq", &out_globalParT3_TopbWqq);
    outTree->Branch("globalParT3_TopbWtauhv", &out_globalParT3_TopbWtauhv);
    outTree->Branch("globalParT3_WvsQCD", &out_globalParT3_WvsQCD);
    outTree->Branch("globalParT3_XWW3q", &out_globalParT3_XWW3q);
    outTree->Branch("globalParT3_XWW4q", &out_globalParT3_XWW4q);
    outTree->Branch("globalParT3_XWWqqev", &out_globalParT3_XWWqqev);
    outTree->Branch("globalParT3_XWWqqmv", &out_globalParT3_XWWqqmv);
    outTree->Branch("globalParT3_Xbb", &out_globalParT3_Xbb);
    outTree->Branch("globalParT3_Xcc", &out_globalParT3_Xcc);
    outTree->Branch("globalParT3_Xcs", &out_globalParT3_Xcs);
    outTree->Branch("globalParT3_Xqq", &out_globalParT3_Xqq);
    outTree->Branch("globalParT3_Xtauhtaue", &out_globalParT3_Xtauhtaue);
    outTree->Branch("globalParT3_Xtauhtauh", &out_globalParT3_Xtauhtauh);
    outTree->Branch("globalParT3_Xtauhtaum", &out_globalParT3_Xtauhtaum);
    outTree->Branch("globalParT3_massCorrGeneric", &out_globalParT3_massCorrGeneric);
    outTree->Branch("globalParT3_massCorrX2p", &out_globalParT3_massCorrX2p);
    outTree->Branch("globalParT3_withMassTopvsQCD", &out_globalParT3_withMassTopvsQCD);
    outTree->Branch("globalParT3_withMassWvsQCD", &out_globalParT3_withMassWvsQCD);
    outTree->Branch("globalParT3_withMassZvsQCD", &out_globalParT3_withMassZvsQCD);
    
    outTree->Branch("dr_H", &out_dr_H);
    outTree->Branch("dr_H_daus", &out_dr_H_daus);
    outTree->Branch("H_pt", &out_H_pt);
    outTree->Branch("H_decay", &out_H_decay);
    
    outTree->Branch("dr_Z", &out_dr_Z);
    outTree->Branch("dr_Z_daus", &out_dr_Z_daus);
    outTree->Branch("Z_pt", &out_Z_pt);
    outTree->Branch("Z_decay", &out_Z_decay);
    
    outTree->Branch("dr_W", &out_dr_W);
    outTree->Branch("dr_W_daus", &out_dr_W_daus);
    outTree->Branch("W_pt", &out_W_pt);
    outTree->Branch("W_decay", &out_W_decay);
    
    outTree->Branch("dr_T", &out_dr_T);
    outTree->Branch("dr_T_b", &out_dr_T_b);
    outTree->Branch("dr_T_Wq_max", &out_dr_T_Wq_max);
    outTree->Branch("dr_T_Wq_min", &out_dr_T_Wq_min);
    outTree->Branch("T_Wq_max_pdgId", &out_T_Wq_max_pdgId);
    outTree->Branch("T_Wq_min_pdgId", &out_T_Wq_min_pdgId);
    outTree->Branch("T_pt", &out_T_pt);
    
    outTree->Branch("dr_LepT", &out_dr_LepT);
    outTree->Branch("dr_LepT_b", &out_dr_LepT_b);
    outTree->Branch("dr_LepT_l", &out_dr_LepT_l);
    outTree->Branch("dr_LepT_l_pdgId", &out_dr_LepT_l_pdgId);
    outTree->Branch("LepT_pt", &out_LepT_pt);
    
    // Loop over events
    Long64_t nEntries = inTree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        if (iEntry % 1000 == 0) {
            std::cout << "Processing event " << iEntry << " / " << nEntries << std::endl;
        }
        
        inTree->GetEntry(iEntry);
        
        // Build GenParticle collection
        std::vector<GenParticle> genparts(nGenPart);
        for (Int_t i = 0; i < nGenPart; i++) {
            genparts[i].pt = GenPart_pt[i];
            genparts[i].eta = GenPart_eta[i];
            genparts[i].phi = GenPart_phi[i];
            genparts[i].mass = GenPart_mass[i];
            genparts[i].pdgId = GenPart_pdgId[i];
            genparts[i].statusFlags = GenPart_statusFlags[i];
        }
        
        // Build mother-daughter relationships
        for (Int_t i = 0; i < nGenPart; i++) {
            int motherIdx = GenPart_genPartIdxMother[i];
            if (motherIdx >= 0 && motherIdx < (int)nGenPart) {
                genparts[motherIdx].dauIdx.push_back(i);
            }
        }
        
        // Find hadronic/leptonic gen particles
        std::vector<GenParticle*> lepGenTops, hadGenTops, hadGenWs, hadGenZs, hadGenHs;
        
        for (auto& gp : genparts) {
            // Check if isLastCopy (bit 13)
            if ((gp.statusFlags & (1 << 13)) == 0) continue;
            
            if (abs(gp.pdgId) == 6) {  // top quark
                for (int idx : gp.dauIdx) {
                    auto& dau = genparts[idx];
                    if (abs(dau.pdgId) == 24) {  // W boson
                        GenParticle* genW = getFinal(&dau, genparts);
                        gp.genW = genW;
                        if (isHadronic(genW, genparts)) {
                            hadGenTops.push_back(&gp);
                        } else {
                            lepGenTops.push_back(&gp);
                        }
                    } else if (abs(dau.pdgId) >= 1 && abs(dau.pdgId) <= 5) {  // b quark
                        gp.genB = &dau;
                    }
                }
            } else if (abs(gp.pdgId) == 24) {  // W boson
                if (isHadronic(&gp, genparts)) {
                    hadGenWs.push_back(&gp);
                }
            } else if (abs(gp.pdgId) == 23) {  // Z boson
                if (isHadronic(&gp, genparts)) {
                    hadGenZs.push_back(&gp);
                }
            } else if (abs(gp.pdgId) == 25) {  // Higgs boson
                if (isHadronic(&gp, genparts)) {
                    hadGenHs.push_back(&gp);
                }
            }
        }
        
        // Print gen particle counts for this event
        // std::cout << "Event " << iEntry << ": "
        //           << "lepGenTops=" << lepGenTops.size() << ", "
        //           << "hadGenTops=" << hadGenTops.size() << ", "
        //           << "hadGenWs=" << hadGenWs.size() << ", "
        //           << "hadGenZs=" << hadGenZs.size() << ", "
        //           << "hadGenHs=" << hadGenHs.size() << std::endl;
        
        // Set up daughter collections for tops
        for (auto* parton : lepGenTops) {
            if (parton->genB && parton->genW && parton->genW->dauIdx.size() >= 2) {
                parton->daus = {parton->genB, 
                               &genparts[parton->genW->dauIdx[0]], 
                               &genparts[parton->genW->dauIdx[1]]};
                parton->genW->daus = {parton->daus[1], parton->daus[2]};
            }
        }
        
        for (auto* parton : hadGenTops) {
            if (parton->genB && parton->genW && parton->genW->dauIdx.size() >= 2) {
                parton->daus = {parton->genB, 
                               &genparts[parton->genW->dauIdx[0]], 
                               &genparts[parton->genW->dauIdx[1]]};
                parton->genW->daus = {parton->daus[1], parton->daus[2]};
            }
        }
        
        // Set up daughter collections for W/Z/H
        for (auto* parton : hadGenWs) {
            if (parton->dauIdx.size() >= 2) {
                parton->daus = {&genparts[parton->dauIdx[0]], 
                               &genparts[parton->dauIdx[1]]};
            }
        }
        
        for (auto* parton : hadGenZs) {
            if (parton->dauIdx.size() >= 2) {
                parton->daus = {&genparts[parton->dauIdx[0]], 
                               &genparts[parton->dauIdx[1]]};
            }
        }
        
        for (auto* parton : hadGenHs) {
            if (parton->dauIdx.size() >= 2) {
                parton->daus = {&genparts[parton->dauIdx[0]], 
                               &genparts[parton->dauIdx[1]]};
            }
        }
        
        // Loop over FatJets
        for (Int_t iFJ = 0; iFJ < nFatJet; iFJ++) {
            FatJetInfo fj;
            fj.pt = FatJet_pt[iFJ];
            fj.eta = FatJet_eta[iFJ];
            fj.phi = FatJet_phi[iFJ];
            fj.mass = FatJet_mass[iFJ];
            
            // Find closest gen particles
            auto [genH, dr_H] = closest(fj, hadGenHs);
            fj.genH = genH;
            fj.dr_H = dr_H;
            
            auto [genZ, dr_Z] = closest(fj, hadGenZs);
            fj.genZ = genZ;
            fj.dr_Z = dr_Z;
            
            auto [genW, dr_W] = closest(fj, hadGenWs);
            fj.genW = genW;
            fj.dr_W = dr_W;
            
            auto [genT, dr_T] = closest(fj, hadGenTops);
            fj.genT = genT;
            fj.dr_T = dr_T;
            
            auto [genLepT, dr_LepT] = closest(fj, lepGenTops);
            fj.genLepT = genLepT;
            fj.dr_LepT = dr_LepT;
            
            // Fill basic FatJet variables
            out_pt = fj.pt;
            out_eta = fj.eta;
            out_phi = fj.phi;
            out_mass = fj.mass;
            out_msoftdrop = FatJet_msoftdrop[iFJ];
            out_particleNetLegacy_QCD = FatJet_particleNetLegacy_QCD[iFJ];
            out_particleNetLegacy_Xbb = FatJet_particleNetLegacy_Xbb[iFJ];
            out_particleNetLegacy_Xcc = FatJet_particleNetLegacy_Xcc[iFJ];
            out_particleNetLegacy_Xqq = FatJet_particleNetLegacy_Xqq[iFJ];
            out_particleNetLegacy_mass = FatJet_particleNetLegacy_mass[iFJ];
            out_particleNet_QCD = FatJet_particleNet_QCD[iFJ];
            out_particleNet_QCD0HF = FatJet_particleNet_QCD0HF[iFJ];
            out_particleNet_QCD1HF = FatJet_particleNet_QCD1HF[iFJ];
            out_particleNet_QCD2HF = FatJet_particleNet_QCD2HF[iFJ];
            out_particleNet_WVsQCD = FatJet_particleNet_WVsQCD[iFJ];
            out_particleNet_XbbVsQCD = FatJet_particleNet_XbbVsQCD[iFJ];
            out_particleNet_XccVsQCD = FatJet_particleNet_XccVsQCD[iFJ];
            out_particleNet_XggVsQCD = FatJet_particleNet_XggVsQCD[iFJ];
            out_particleNet_XqqVsQCD = FatJet_particleNet_XqqVsQCD[iFJ];
            out_particleNet_XteVsQCD = FatJet_particleNet_XteVsQCD[iFJ];
            out_particleNet_XtmVsQCD = FatJet_particleNet_XtmVsQCD[iFJ];
            out_particleNet_XttVsQCD = FatJet_particleNet_XttVsQCD[iFJ];
            out_particleNet_massCorr = FatJet_particleNet_massCorr[iFJ];
            out_particleNetWithMass_H4qvsQCD = FatJet_particleNetWithMass_H4qvsQCD[iFJ];
            out_particleNetWithMass_HbbvsQCD = FatJet_particleNetWithMass_HbbvsQCD[iFJ];
            out_particleNetWithMass_HccvsQCD = FatJet_particleNetWithMass_HccvsQCD[iFJ];
            out_particleNetWithMass_QCD = FatJet_particleNetWithMass_QCD[iFJ];
            out_particleNetWithMass_TvsQCD = FatJet_particleNetWithMass_TvsQCD[iFJ];
            out_particleNetWithMass_WvsQCD = FatJet_particleNetWithMass_WvsQCD[iFJ];
            out_particleNetWithMass_ZvsQCD = FatJet_particleNetWithMass_ZvsQCD[iFJ];
            out_globalParT3_QCD = FatJet_globalParT3_QCD[iFJ];
            out_globalParT3_TopbWev = FatJet_globalParT3_TopbWev[iFJ];
            out_globalParT3_TopbWmv = FatJet_globalParT3_TopbWmv[iFJ];
            out_globalParT3_TopbWq = FatJet_globalParT3_TopbWq[iFJ];
            out_globalParT3_TopbWqq = FatJet_globalParT3_TopbWqq[iFJ];
            out_globalParT3_TopbWtauhv = FatJet_globalParT3_TopbWtauhv[iFJ];
            out_globalParT3_WvsQCD = FatJet_globalParT3_WvsQCD[iFJ];
            out_globalParT3_XWW3q = FatJet_globalParT3_XWW3q[iFJ];
            out_globalParT3_XWW4q = FatJet_globalParT3_XWW4q[iFJ];
            out_globalParT3_XWWqqev = FatJet_globalParT3_XWWqqev[iFJ];
            out_globalParT3_XWWqqmv = FatJet_globalParT3_XWWqqmv[iFJ];
            out_globalParT3_Xbb = FatJet_globalParT3_Xbb[iFJ];
            out_globalParT3_Xcc = FatJet_globalParT3_Xcc[iFJ];
            out_globalParT3_Xcs = FatJet_globalParT3_Xcs[iFJ];
            out_globalParT3_Xqq = FatJet_globalParT3_Xqq[iFJ];
            out_globalParT3_Xtauhtaue = FatJet_globalParT3_Xtauhtaue[iFJ];
            out_globalParT3_Xtauhtauh = FatJet_globalParT3_Xtauhtauh[iFJ];
            out_globalParT3_Xtauhtaum = FatJet_globalParT3_Xtauhtaum[iFJ];
            out_globalParT3_massCorrGeneric = FatJet_globalParT3_massCorrGeneric[iFJ];
            out_globalParT3_massCorrX2p = FatJet_globalParT3_massCorrX2p[iFJ];
            out_globalParT3_withMassTopvsQCD = FatJet_globalParT3_withMassTopvsQCD[iFJ];
            out_globalParT3_withMassWvsQCD = FatJet_globalParT3_withMassWvsQCD[iFJ];
            out_globalParT3_withMassZvsQCD = FatJet_globalParT3_withMassZvsQCD[iFJ];
            
            // Fill H matching info
            out_dr_H = fj.dr_H;
            if (fj.genH && fj.genH->daus.size() >= 2) {
                float maxDR = std::max(deltaR(fj, fj.genH->daus[0]), 
                                      deltaR(fj, fj.genH->daus[1]));
                out_dr_H_daus = maxDR;
                out_H_pt = fj.genH->pt;
                out_H_decay = abs(fj.genH->daus[0]->pdgId);
            } else {
                out_dr_H_daus = 99;
                out_H_pt = -1;
                out_H_decay = 0;
            }
            
            // Fill Z matching info
            out_dr_Z = fj.dr_Z;
            if (fj.genZ && fj.genZ->daus.size() >= 2) {
                float maxDR = std::max(deltaR(fj, fj.genZ->daus[0]), 
                                      deltaR(fj, fj.genZ->daus[1]));
                out_dr_Z_daus = maxDR;
                out_Z_pt = fj.genZ->pt;
                out_Z_decay = abs(fj.genZ->daus[0]->pdgId);
            } else {
                out_dr_Z_daus = 99;
                out_Z_pt = -1;
                out_Z_decay = 0;
            }
            
            // Fill W matching info
            out_dr_W = fj.dr_W;
            if (fj.genW && fj.genW->daus.size() >= 2) {
                float maxDR = std::max(deltaR(fj, fj.genW->daus[0]), 
                                      deltaR(fj, fj.genW->daus[1]));
                out_dr_W_daus = maxDR;
                out_W_pt = fj.genW->pt;
                int decay1 = abs(fj.genW->daus[0]->pdgId);
                int decay2 = abs(fj.genW->daus[1]->pdgId);
                out_W_decay = std::max(decay1, decay2);
            } else {
                out_dr_W_daus = 99;
                out_W_pt = -1;
                out_W_decay = 0;
            }
            
            // Fill hadronic top matching info
            out_dr_T = fj.dr_T;
            if (fj.genT && fj.genT->genB && fj.genT->genW && fj.genT->genW->daus.size() >= 2) {
                out_dr_T_b = deltaR(fj, fj.genT->genB);
                
                float drwq1 = deltaR(fj, fj.genT->genW->daus[0]);
                float drwq2 = deltaR(fj, fj.genT->genW->daus[1]);
                int wq1_pdgId = fj.genT->genW->daus[0]->pdgId;
                int wq2_pdgId = fj.genT->genW->daus[1]->pdgId;
                
                if (drwq1 < drwq2) {
                    std::swap(drwq1, drwq2);
                    std::swap(wq1_pdgId, wq2_pdgId);
                }
                
                out_dr_T_Wq_max = drwq1;
                out_dr_T_Wq_min = drwq2;
                out_T_Wq_max_pdgId = wq1_pdgId;
                out_T_Wq_min_pdgId = wq2_pdgId;
                out_T_pt = fj.genT->pt;
            } else {
                out_dr_T_b = 99;
                out_dr_T_Wq_max = 99;
                out_dr_T_Wq_min = 99;
                out_T_Wq_max_pdgId = 0;
                out_T_Wq_min_pdgId = 0;
                out_T_pt = -1;
            }
            
            // Fill leptonic top matching info
            out_dr_LepT = fj.dr_LepT;
            if (fj.genLepT && fj.genLepT->genB && fj.genLepT->genW && fj.genLepT->genW->daus.size() >= 2) {
                out_dr_LepT_b = deltaR(fj, fj.genLepT->genB);
                
                // Find the lepton (e, mu, tau)
                GenParticle* lepton = nullptr;
                if (abs(fj.genLepT->genW->daus[0]->pdgId) == 11 || 
                    abs(fj.genLepT->genW->daus[0]->pdgId) == 13 || 
                    abs(fj.genLepT->genW->daus[0]->pdgId) == 15) {
                    lepton = fj.genLepT->genW->daus[0];
                } else if (abs(fj.genLepT->genW->daus[1]->pdgId) == 11 || 
                           abs(fj.genLepT->genW->daus[1]->pdgId) == 13 || 
                           abs(fj.genLepT->genW->daus[1]->pdgId) == 15) {
                    lepton = fj.genLepT->genW->daus[1];
                }
                
                if (lepton) {
                    out_dr_LepT_l = deltaR(fj, lepton);
                    out_dr_LepT_l_pdgId = lepton->pdgId;
                } else {
                    out_dr_LepT_l = 99;
                    out_dr_LepT_l_pdgId = 0;
                }
                
                out_LepT_pt = fj.genLepT->pt;
            } else {
                out_dr_LepT_b = 99;
                out_dr_LepT_l = 99;
                out_dr_LepT_l_pdgId = 0;
                out_LepT_pt = -1;
            }
            
            // Fill the output tree
            outTree->Fill();
        }
    }
    Int_t nOutputEntries = outTree->GetEntries();
    
    // Write and close output file
    outFile->cd();
    outTree->Write();
    outFile->Close();
    
    std::cout << "Output written to " << outputFile << std::endl;
    std::cout << "Total FatJet entries: " << nOutputEntries << std::endl;

    inFile->Close();
}

