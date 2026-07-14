#ifndef WAVEFORMGENERATOR_H
#define WAVEFORMGENERATOR_H

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TMath.h"
#include <random>
#include <vector>

#include "../Podd/FadcStandaloneAnalyzer.h"

struct WaveformData {
  Int_t TotNSamples;
  Int_t InitPedSamples;

  Int_t Pedestal;
  Int_t PedRange;

  Int_t NPulses;
  Int_t PulseAmp;
  Int_t PulseWidth;
  Int_t PulseSeparation;
};

struct WaveformInfo {
  std::vector<Double_t>
      PedSubIntegral; // Double to match analyzer's floating-point formula
  std::vector<Int_t> PulseAmp;
  std::vector<Int_t> PulseTime;
};

class waveformGenerator {
public:
  waveformGenerator(); // Seed the random number generator. New seed everytime
                       // run.
  virtual ~waveformGenerator() = default;

  // Waveform Data Setters
  void setWaveformData(const WaveformData &waveformData) {
    fWaveformData = waveformData;
  }
  void setAnalyzer(const FadcStandaloneAnalyzer &analyzer) {
    saAnalyzer = analyzer;
  }

  // Waveform generators
  std::vector<Int_t> noiseWaveform() const;
  std::vector<Int_t> singlePulseWaveform() const;
  std::vector<Int_t> multiplePulseWaveform() const;

  // Calc pulse parameters
  // std::vector<Int_t> calcPulseIntegral(const std::vector<Int_t> wave) const;
  void calcPulseParameters(const std::vector<Int_t> wave);

  // Accessors for pulse parameters computed by calcPulseParameters()
  const WaveformInfo &getWaveformInfo() const { return fWaveformInfo; }
  void clearWaveformInfo() { fWaveformInfo = {}; }

private:
  WaveformData fWaveformData;
  WaveformInfo fWaveformInfo;

  FadcStandaloneAnalyzer saAnalyzer;

  // Per-instance RNG — thread-safe (no shared global state)
  mutable std::mt19937 fRng;
};

#endif // WAVEFORMGENERATOR_H