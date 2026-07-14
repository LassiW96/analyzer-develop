#ifndef FADCSATANDALONEANALYZER_H
#define FADCSATANDALONEANALYZER_H

#include <cmath>
#include <vector>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TMath.h"

//_________________Data structures which change only per-module
//basis________________ Configurations
struct FADCConfig {
  // configuration parameters - taken from db
  Int_t fNPedestalSamples;
  Double_t fSampThreshold;
  Int_t fNSB;
  Int_t fNSA;
  Int_t fNSAT;
  Int_t fMaxNPulses;
};

// Calibrations
struct FADCCalib {
  // Calibration parameters - these can be also taken from db at run time
  Double_t ChanTomV;
  Double_t pC_Conv;

  // These change per channel basis
  // Double_t tcal;
  // Double_t gain;
};

// Output Pulse Data
struct FADCPulse {
  Double_t fSampPed;
  Bool_t fHasMulti;
  Int_t fNSampPulses;
  Int_t fNPeakSamples;
  Double_t fPeakPedestalRatio;

  // Output Pulse Data
  std::vector<Double_t> fSampPulseInt;
  std::vector<Double_t> fSampPulseAmp;
  std::vector<Double_t> fSampPulseTime;
  std::vector<Double_t> fSampPulseIntPedSub;
  std::vector<Double_t> fSampPulseIntMOLLERRaw;
  std::vector<Double_t> fSampPulseIntMOLLERVal;
  std::vector<Double_t> fSampPulseTOT;
  std::vector<Int_t> fSampPulsePileup;
};

class FadcStandaloneAnalyzer {
public:
  FadcStandaloneAnalyzer();
  virtual ~FadcStandaloneAnalyzer() = default;

  // Configuration setters
  void setConfig(const FADCConfig &config) { fConfig = config; }
  void setCalib(const FADCCalib &calib) { fCalib = calib; }

  // Getters
  FADCConfig getConfig() const { return fConfig; }
  FADCCalib getCalib() const { return fCalib; }

  // Main analysis function
  FADCPulse Analyze(const std::vector<Double_t> &samples, Double_t gain,
                    Double_t tcal) const;

private:
  FADCConfig fConfig;
  FADCCalib fCalib;

  // Any helper functions needed, to help the main analysis function ...
  // Pedestal calculation function
  // Pulse pillup finding function
  // NSAT passing findng function

  // ClassDef(FADCStandaloneAnalyzer, 0);
}; // end FADCStandaloneAnalyzer

#endif // #ifndef FADCSATANDALONEANALYZER_H