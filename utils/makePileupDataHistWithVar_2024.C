// makePileupDataHistWithVar_2024.C
// Build a 2024 "WithVar" PU file by summing eras B..I for each minBiasXsec.
// Output contains exactly one key per histogram: pileup, pileup_plus, pileup_minus.

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
    if (cl->InheritsFrom("TH1")) return (TH1*)key->ReadObj();
  }
  return nullptr;
}

TH1* GetPUHist(TFile* f) {
  if (!f || f->IsZombie()) return nullptr;

  // Try common names first
  const std::vector<TString> candidates = {
    "pileup", "pileup_plus", "pileup_minus",
    "pileup_hist", "pileupHistogram", "hPileup", "pu", "PU"
  };

  for (const auto& name : candidates) {
    if (TH1* h = dynamic_cast<TH1*>(f->Get(name))) return h;
  }

  // Fallback: first TH1 in file
  return GetFirstTH1(f);
}

// Sum a list of ROOT files into a single TH1 (clone of first + Add the rest).
// Returns a heap object owned by caller (SetDirectory(nullptr)).
TH1* SumEraFiles(const std::vector<TString>& paths, const char* sumName) {
  TH1* hsum = nullptr;

  for (const auto& path : paths) {
    TFile* fin = TFile::Open(path, "READ");
    if (!fin || fin->IsZombie()) {
      std::cerr << "ERROR: Cannot open " << path << "\n";
      if (fin) fin->Close();
      continue;
    }

    TH1* h = GetPUHist(fin);
    if (!h) {
      std::cerr << "ERROR: No TH1 found in " << path << "\n";
      fin->Close();
      continue;
    }

    // Detach histogram from input file (clone) so it's safe after closing fin
    if (!hsum) {
      hsum = (TH1*)h->Clone(sumName);
      hsum->SetDirectory(nullptr);
    } else {
      // Optional safety checks (binning/range)
      if (h->GetNbinsX() != hsum->GetNbinsX() ||
          h->GetXaxis()->GetXmin() != hsum->GetXaxis()->GetXmin() ||
          h->GetXaxis()->GetXmax() != hsum->GetXaxis()->GetXmax()) {
        std::cerr << "WARNING: Binning/range mismatch in " << path
                  << " (still adding as-is)\n";
      }
      hsum->Add(h);
    }

    fin->Close();
  }

  if (!hsum) {
    std::cerr << "ERROR: SumEraFiles produced null histogram for " << sumName << "\n";
  }
  return hsum;
}

void WriteOnce(TFile* out, TH1* hin, const char* outName) {
  if (!out || !hin) return;
  out->cd();

  // ensure only one cycle even if you re-run in interactive session on same open file
  out->Delete(Form("%s;*", outName));

  TH1* h = (TH1*)hin->Clone(outName);
  h->SetDirectory(out);
  h->Sumw2(false);
  h->Write(); // write exactly once (no fout->Write() later)
}

void makePileupDataHistWithVar_2024() {
  const std::vector<TString> eras = {"B","C","D","E","F","G","H","I"};

  auto buildPaths = [&](int xsec_ub) {
    std::vector<TString> paths;
    paths.reserve(eras.size());
    for (const auto& era : eras) {
      TString p = Form("/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/PileUp/dataPileupHistogram-2024%s_Golden-%dub.root", era.Data(), xsec_ub);
      paths.push_back(p);
    }
    return paths;
  };

  // Down/Nom/Up sets
  auto pathsDown = buildPaths(66000);
  auto pathsNom  = buildPaths(69200);
  auto pathsUp   = buildPaths(72400);

  // Sum each set across eras B..I
  TH1* hDown = SumEraFiles(pathsDown, "sum_66000ub");
  TH1* hNom  = SumEraFiles(pathsNom,  "sum_69200ub");
  TH1* hUp   = SumEraFiles(pathsUp,   "sum_72400ub");

  if (!hDown || !hNom || !hUp) {
    std::cerr << "ERROR: Failed to build one or more summed histograms.\n";
    return;
  }

  // Output file
  const char* outPath = "dataPileupHistogram-2024BCDEFGHI_Golden-13p6TeV_WithVar.root";
  TFile* fout = TFile::Open(outPath, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: Could not create output file: " << outPath << "\n";
    return;
  }

  // Convention: nominal / up / down
  WriteOnce(fout, hNom,  "pileup");
  WriteOnce(fout, hUp,   "pileup_plus");
  WriteOnce(fout, hDown, "pileup_minus");

  fout->Close();

  // clean up heap histograms we created
  delete hDown;
  delete hNom;
  delete hUp;

  std::cout << "Wrote: " << outPath << "\n"
            << "Keys: pileup, pileup_plus, pileup_minus\n";
}
