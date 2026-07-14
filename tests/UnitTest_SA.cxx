#include "../Podd/FadcStandaloneAnalyzer.h"
#include "waveformGenerator.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>

// ─── Helpers ─────────────────────────────────────────────────────────────────

// waveformGenerator returns Int_t samples; FADCStandaloneAnalyzer::Analyze()
// takes Double_t. This helper does the element-wise cast.
std::vector<Double_t> ToDoubleWave(const std::vector<Int_t> &wave) {
  return std::vector<Double_t>(wave.begin(), wave.end());
}

// Shared analyzer configuration.
FADCConfig GetTestConfig() {
  FADCConfig conf;
  conf.fNPedestalSamples = 4;
  conf.fSampThreshold = 50.0;
  conf.fNSB = 1;
  conf.fNSA = 2;
  conf.fNSAT = 1;
  conf.fMaxNPulses = 3;
  return conf;
}

// ─── Tests ───────────────────────────────────────────────────────────────────

TEST_CASE("WaveformAnalyzer correctly processes fADC data", "[analyzer]") {

  FadcStandaloneAnalyzer analyzer;
  analyzer.setConfig(GetTestConfig());
  FADCCalib calib = {0.2441, 4.0 / 50.0};
  analyzer.setCalib(calib);

  WaveformData wfData;
  wfData.TotNSamples = 100;
  wfData.InitPedSamples = 4;
  wfData.Pedestal = 100;
  wfData.PedRange = 4;
  wfData.NPulses = 4;
  wfData.PulseAmp = 300;
  wfData.PulseWidth = 5;
  wfData.PulseSeparation = 10;

  waveformGenerator gen;
  gen.setWaveformData(wfData);
  gen.setAnalyzer(analyzer); // inject same config so singleWaveIntegral uses
                             // identical FADCConfig

  SECTION("Empty waveform returns no pulses and doesn't crash") {
    std::vector<Double_t> empty_wave;
    auto results = analyzer.Analyze(empty_wave, 1.0, 4.0);
    REQUIRE(results.fNSampPulses == 0);
  }

  SECTION("Noise waveform (below threshold) returns no pulses") {
    // Generate a pure noise waveform and confirm the analyzer finds nothing.
    std::vector<Double_t> wave = ToDoubleWave(gen.noiseWaveform());
    auto results = analyzer.Analyze(wave, 1.0, 4.0);
    REQUIRE(results.fNSampPulses == 0);
  }

  SECTION("Single pulse waveform is detected with a positive "
          "pedestal-subtracted integral") {
    std::vector<Int_t> int_wave = gen.singlePulseWaveform();
    std::vector<Double_t> wave = ToDoubleWave(int_wave);

    // --- FADCStandaloneAnalyzer result ---
    auto results = analyzer.Analyze(wave, 1.0, 4.0);

    // At least one pulse must be detected.
    REQUIRE(results.fNSampPulses != 0);

    // Pulse parameters using waveformGenerator class (for comparison)
    gen.clearWaveformInfo();
    gen.calcPulseParameters(int_wave);
    const WaveformInfo &wfInfo = gen.getWaveformInfo();
    REQUIRE(results.fSampPulseIntPedSub[0] ==
            Catch::Approx(static_cast<Double_t>(wfInfo.PedSubIntegral[0]))
                .margin(wfData.PedRange));
    REQUIRE(results.fSampPulseAmp[0] ==
            Catch::Approx(static_cast<Double_t>(wfInfo.PulseAmp[0]))
                .margin(wfData.PedRange));

    // std::cout << "\nSingle Pulse:\nPulse Time: From Waveform generator " <<
    // wfInfo.PulseTime[0] << " from FADCStandaloneAnalyzer " <<
    // results.fSampPulseTime[0] << "\n";
    REQUIRE(results.fSampPulseTime[0] ==
            Catch::Approx(static_cast<Double_t>(wfInfo.PulseTime[0]))
                .margin(wfData.PedRange));
  }

  SECTION("Multiple pulse waveform is detected with a positive "
          "pedestal-subtracted integral") {
    std::vector<Int_t> int_wave = gen.multiplePulseWaveform();
    std::vector<Double_t> wave = ToDoubleWave(int_wave);

    // --- FADCStandaloneAnalyzer result ---
    auto results = analyzer.Analyze(wave, 1.0, 4.0);

    // Output shouldn't be empty
    REQUIRE(results.fNSampPulses >= 2);

    // Pulse parameters using waveformGenerator class (for comparison)
    gen.clearWaveformInfo();
    gen.calcPulseParameters(int_wave);
    const WaveformInfo &wfInfo = gen.getWaveformInfo();

    // std::cout << "\nCalculated Integrals and times for Multiple Pulses:\n";
    for (int i = 0; i < (int)wfInfo.PedSubIntegral.size(); i++) {
      // std::cout << "\tPulse " << i + 1 << ": " << wfInfo.PedSubIntegral[i] <<
      // std::endl;
      REQUIRE(results.fSampPulseIntPedSub[i] ==
              Catch::Approx(static_cast<Double_t>(wfInfo.PedSubIntegral[i]))
                  .margin(wfData.PedRange));
      REQUIRE(results.fSampPulseAmp[i] ==
              Catch::Approx(static_cast<Double_t>(wfInfo.PulseAmp[i]))
                  .margin(wfData.PedRange));

      // std::cout << "\tPulse Time: From Waveform generator " <<
      // wfInfo.PulseTime[i] << " from FADCStandaloneAnalyzer " <<
      // results.fSampPulseTime[i] << "\n\n";
      REQUIRE(results.fSampPulseTime[i] ==
              Catch::Approx(static_cast<Double_t>(wfInfo.PulseTime[i]))
                  .margin(wfData.PedRange));
    }
  }
}