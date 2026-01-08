// makePileupDataHistWithVar.C
// Usage:
//   root -l -b -q 'makePileupDataHistWithVar.C()'
//
// Produces: pileupHistogram-..._WithVar.root
// with histograms: pileup (nom), pileup_plus (up), pileup_minus (down)

#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TClass.h"
#include "TString.h"
#include <vector>
#include <iostream>

TH1* GetFirstTH1(TFile* f) {
  if (!f || f->IsZombie()) return nullptr;

  TIter nextKey(f->GetListOfKeys());
  while (TKey* key = (TKey*)nextKey()) {
    TClass* cl = gROOT->GetClass(key->GetClassName());
    if (!cl) continue;
    if (cl->InheritsFrom("TH1")) {
      TH1* h = (TH1*)key->ReadObj();
      return h;
    }
  }
  return nullptr;
}

TH1* GetPUHist(TFile* f) {
  if (!f || f->IsZombie()) return nullptr;

  // Common histogram names in pileupCalc outputs
  const std::vector<TString> candidates = {
    "pileup",
    "pileup_hist",
    "pileupHistogram",
    "hPileup",
    "pu",
    "PU"
  };

  for (const auto& name : candidates) {
    TH1* h = dynamic_cast<TH1*>(f->Get(name));
    if (h) return h;
  }

  // Fallback: first TH1 in file
  return GetFirstTH1(f);
}

void WriteClone(TFile* out, TH1* hin, const char* outName) {
  if (!out || !hin) return;
  out->cd();

  TH1* h = (TH1*)hin->Clone(outName);
  h->SetDirectory(out);   // ownership by output file
  h->Sumw2(false);
  h->Write();             // write ONCE, no overwrite
}

void makePileupDataHistWithVar() {
  // ---- EDIT THESE PATHS ----
//   const char* fDownPath = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/PileUp/BCD/pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-66000ub-99bins.root";
//   const char* fNomPath  = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/PileUp/BCD/pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-69200ub-99bins.root";
//   const char* fUpPath   = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/PileUp/BCD/pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-72400ub-99bins.root";

//   const char* outPath   = "pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV_WithVar.root";

  // special path for 2024
  const char* fDownPath = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/PileUp/dataPileupHistogram-2024BCDEFGHI-66000ub.root";
  const char* fNomPath  = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/PileUp/dataPileupHistogram-2024BCDEFGHI-69200ub.root";
  const char* fUpPath   = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/PileUp/dataPileupHistogram-2024BCDEFGHI-72400ub.root";

  const char* outPath   = "pileupHistogram-Cert_Collisions2024BCDEFGHI_GoldenJson-13p6TeV_WithVar.root";
  // --------------------------

  TFile* fDown = TFile::Open(fDownPath, "READ");
  TFile* fNom  = TFile::Open(fNomPath,  "READ");
  TFile* fUp   = TFile::Open(fUpPath,   "READ");

  if (!fDown || fDown->IsZombie() || !fNom || fNom->IsZombie() || !fUp || fUp->IsZombie()) {
    std::cerr << "ERROR: Could not open one or more input files.\n";
    return;
  }

  TH1* hDown = GetPUHist(fDown);
  TH1* hNom  = GetPUHist(fNom);
  TH1* hUp   = GetPUHist(fUp);

  if (!hDown || !hNom || !hUp) {
    std::cerr << "ERROR: Could not find a TH1 in one or more input files.\n";
    std::cerr << "Try inspecting with: root -l <file>.root and then <file>.ls()\n";
    return;
  }

  // Basic sanity checks
  if (hDown->GetNbinsX() != hNom->GetNbinsX() || hUp->GetNbinsX() != hNom->GetNbinsX()) {
    std::cerr << "WARNING: Different bin counts among histograms.\n";
  }
  if (hDown->GetXaxis()->GetXmin() != hNom->GetXaxis()->GetXmin() ||
      hDown->GetXaxis()->GetXmax() != hNom->GetXaxis()->GetXmax() ||
      hUp->GetXaxis()->GetXmin()   != hNom->GetXaxis()->GetXmin()   ||
      hUp->GetXaxis()->GetXmax()   != hNom->GetXaxis()->GetXmax()) {
    std::cerr << "WARNING: Different x-axis ranges among histograms.\n";
  }

  TFile* fout = TFile::Open(outPath, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: Could not create output file: " << outPath << "\n";
    return;
  }

  // Convention: nominal, plus(up), minus(down)
  WriteClone(fout, hNom,  "pileup");
  WriteClone(fout, hUp,   "pileup_plus");
  WriteClone(fout, hDown, "pileup_minus");

  fout->Close();

  fDown->Close();
  fNom->Close();
  fUp->Close();

  std::cout << "Wrote: " << outPath << "\n";
  std::cout << "Contains: pileup, pileup_plus, pileup_minus\n";
}
