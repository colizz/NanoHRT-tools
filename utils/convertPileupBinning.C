// File: convertPileupBinning.C
// Usage from ROOT:
//   root -l -q 'convertPileupBinning.C()'

#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"
#include "TString.h"
#include <iostream>
#include <cmath>

static bool AlmostEqual(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) < eps;
}

static void AssertInputBinning(const TH1D* h, const char* name) {
  if (!h) {
    std::cerr << "ERROR: histogram '" << name << "' is null.\n";
    throw std::runtime_error("Missing histogram");
  }

  const int nb = h->GetNbinsX();
  const double xmin = h->GetXaxis()->GetXmin();
  const double xmax = h->GetXaxis()->GetXmax();
  const double bw   = h->GetXaxis()->GetBinWidth(1);

  // Expect exactly 99 bins covering [0,99] with width 1
  if (nb != 99 ||
      !AlmostEqual(xmin, 0.0) ||
      !AlmostEqual(xmax, 99.0) ||
      !AlmostEqual(bw, 1.0)) {
    std::cerr << "ERROR: histogram '" << name << "' has unexpected binning.\n"
              << "  Got: nbins=" << nb << " range=[" << xmin << "," << xmax << "]"
              << " binWidth=" << bw << "\n"
              << "  Expected: nbins=99 range=[0,99] binWidth=1\n";
    throw std::runtime_error("Unexpected binning");
  }

  // Also sanity check: every bin edge should align to integers
  for (int i = 1; i <= nb + 1; ++i) {
    const double edge = h->GetXaxis()->GetBinLowEdge(i);
    if (!AlmostEqual(edge, std::round(edge))) {
      std::cerr << "ERROR: histogram '" << name << "' has non-integer bin edge at i="
                << i << " edge=" << edge << "\n";
      throw std::runtime_error("Non-integer bin edge");
    }
  }
}

static TH1D* MakeExtended(const TH1D* hin, const char* outName) {
  // New binning: [0,100] with width 1 => 100 bins
  TH1D* hout = new TH1D(outName, hin->GetTitle(), 100, 0.0, 100.0);

  // Copy bin contents/errors for bins 1..99 (0-1 ... 98-99)
  for (int i = 1; i <= 99; ++i) {
    hout->SetBinContent(i, hin->GetBinContent(i));
  }

  // Last bin (bin 100 corresponds to [99,100)) should be 0
  hout->SetBinContent(100, 0.0);

  return hout;
}

void convertPileupBinning() {

  const char* inputFilePath = "/afs/cern.ch/user/c/coli/work/hcc/nano/nanov15/CMSSW_15_0_10/src/PhysicsTools/NanoHRTTools/python/postprocessing/data/pileup/pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV_WithVar.root";
  const char* outputFilePath = "/afs/cern.ch/user/c/coli/work/hcc/nano/nanov15/CMSSW_15_0_10/src/PhysicsTools/NanoHRTTools/python/postprocessing/data/pileup/pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV_100bins_WithVar.root";

  // Open input
  TFile* fin = TFile::Open(inputFilePath, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << inputFilePath << "\n";
    throw std::runtime_error("Failed to open input file");
  }

  const char* names[3] = {"pileup", "pileup_plus", "pileup_minus"};

  TH1D* hIn[3] = {nullptr, nullptr, nullptr};
  for (int i = 0; i < 3; ++i) {
    fin->GetObject(names[i], hIn[i]);
    if (!hIn[i]) {
      std::cerr << "ERROR: cannot find TH1D named '" << names[i]
                << "' in file " << inputFilePath << "\n";
      fin->Close();
      throw std::runtime_error("Histogram not found");
    }
    // Ensure it is really TH1D
    if (std::string(hIn[i]->ClassName()) != "TH1D") {
      std::cerr << "ERROR: object '" << names[i] << "' is not TH1D (it is "
                << hIn[i]->ClassName() << ")\n";
      fin->Close();
      throw std::runtime_error("Wrong histogram type");
    }

    AssertInputBinning(hIn[i], names[i]);
  }

  // Create output file
  TFile* fout = TFile::Open(outputFilePath, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: cannot create output file: " << outputFilePath << "\n";
    fin->Close();
    throw std::runtime_error("Failed to create output file");
  }

  // Create and write extended histograms
  for (int i = 0; i < 3; ++i) {
    fout->cd();
    TH1D* hOut = MakeExtended(hIn[i], names[i]);
    hOut->Write(names[i], TObject::kOverwrite);
    delete hOut; // file now owns a written copy; safe to delete local object
  }

  fout->Close();
  fin->Close();

  std::cout << "Wrote converted histograms to: " << outputFilePath << "\n";
}
