#include "waveformGenerator.h"
#include "../Podd/FadcStandaloneAnalyzer.h"
#include "RtypesCore.h"
#include <vector>

waveformGenerator::waveformGenerator() : fRng(std::random_device{}()) {}

std::vector<Int_t> waveformGenerator::noiseWaveform() const {
  int min = fWaveformData.Pedestal - fWaveformData.PedRange / 2;
  int max = fWaveformData.Pedestal + fWaveformData.PedRange / 2 + 1;
  std::uniform_real_distribution<double> rndm(0.0, 1.0);
  std::vector<Int_t> wave;
  for (int i = 0; i < fWaveformData.TotNSamples; i++) {
    wave.push_back(min + (max - min) * rndm(fRng));
  }
  return wave;
}

std::vector<Int_t> waveformGenerator::singlePulseWaveform() const {
  int min_ped = fWaveformData.Pedestal - fWaveformData.PedRange / 2;
  int max_ped = fWaveformData.Pedestal + fWaveformData.PedRange / 2 + 1;
  std::uniform_real_distribution<double> rndm(0.0, 1.0);
  std::vector<Int_t> wave;

  int pulse_min = fWaveformData.InitPedSamples;
  int pulse_max = fWaveformData.TotNSamples - 1;
  int random_idx = pulse_min + int(pulse_max - pulse_min) * rndm(fRng);
  for (int i = 0; i < fWaveformData.TotNSamples; i++) {
    if (i >= random_idx && i < random_idx + fWaveformData.PulseWidth) {
      int j = i - random_idx;
      wave.push_back(min_ped + fWaveformData.PulseAmp / TMath::Power(2, j));
    } else {
      wave.push_back(min_ped + (max_ped - min_ped) * rndm(fRng));
    }
  }
  return wave;
}

std::vector<Int_t> waveformGenerator::multiplePulseWaveform() const {
  int min_ped = fWaveformData.Pedestal - fWaveformData.PedRange / 2;
  int max_ped = fWaveformData.Pedestal + fWaveformData.PedRange / 2 + 1;
  std::uniform_real_distribution<double> rndm(0.0, 1.0);
  std::vector<Int_t> wave;
  std::vector<std::pair<Int_t, Int_t>> pulse_info;

  int num_pulses = 2 + int((fWaveformData.NPulses - 2) * rndm(fRng));

  // Maximum pulse space required
  int max_space = (num_pulses * fWaveformData.PulseWidth) +
                  ((num_pulses - 1) * fWaveformData.PulseSeparation);

  // Maximum available space for pulse start
  int max_first_start = fWaveformData.TotNSamples - max_space;

  if (max_first_start < fWaveformData.InitPedSamples) {
    max_first_start = fWaveformData.InitPedSamples;
  }

  // Current starting point
  int current_start =
      fWaveformData.InitPedSamples +
      int((max_first_start - fWaveformData.InitPedSamples) * rndm(fRng));

  // Random amps and separations
  for (int p = 0; p < num_pulses; p++) {
    int min_amp = fWaveformData.PulseAmp / 2;
    int current_amp =
        min_amp + int((fWaveformData.PulseAmp - min_amp) * rndm(fRng));

    pulse_info.push_back({current_start, current_amp});

    int current_separation = fWaveformData.PulseSeparation / 2 +
                             int((fWaveformData.PulseSeparation -
                                  fWaveformData.PulseSeparation / 2.) *
                                 rndm(fRng));
    current_start += fWaveformData.PulseWidth + current_separation;
  }

  // Build the waveform
  for (int i = 0; i < fWaveformData.TotNSamples; i++) {
    Bool_t in_pulse = false;
    Int_t pulse_relative_idx = 0;
    Int_t active_amp = 0;

    for (auto &p : pulse_info) {
      if (i >= p.first && i < p.first + fWaveformData.PulseWidth) {
        in_pulse = true;
        pulse_relative_idx = i - p.first;
        active_amp = p.second;
        break;
      }
    }
    if (in_pulse) {
      wave.push_back(min_ped +
                     active_amp / TMath::Power(2, pulse_relative_idx));
    } else {
      wave.push_back(min_ped + (max_ped - min_ped) * rndm(fRng));
    }
  }
  return wave;
}

// std::vector<Int_t> waveformGenerator::calcPulseIntegral(const
// std::vector<Int_t> wave) const {
//     FADCConfig fConfig = saAnalyzer.getConfig();
//     Int_t NPedSamps = fConfig.fNPedestalSamples;
//     Int_t NSB       = fConfig.fNSB;
//     Int_t NSA       = fConfig.fNSA;
//     Int_t Threshold = fConfig.fSampThreshold;
//     Int_t NSamples  = (Int_t)wave.size();

//     std::vector<Int_t> AvgIntegrals = {};

//     // Pedestal: raw sum of first NPedSamps samples (matches analyzer's
//     fSampPed) Int_t PedSum = 0; for (int i = 0; i < NPedSamps; i++)
//         PedSum += wave[i];
//     Int_t AvgPed = PedSum / NPedSamps;

//     Bool_t PulseFound = false;

//     for(int i = 0; i < wave.size(); i++) {
//         if(wave[i] - AvgPed > Threshold) {
//             if(!PulseFound) {
//                 // Integration window: [FirstSamp - NSB, FirstSamp + NSA -
//                 1], clamped to valid range Int_t lo = TMath::Max(i - NSB, 0);
//                 Int_t hi = TMath::Min(i + NSA - 1, NSamples - 1);
//                 // Pedestal-subtracted integral, using the same formula as
//                 the analyzer:
//                 //   PedSubInt = RawInt - PedSum * (window_width / NPedSamps)
//                 //             = RawInt - AvgPed * window_width
//                 Int_t WindowWidth = hi - lo + 1;
//                 Int_t AvgIntegral = 0;
//                 Int_t RawIntegral = 0;
//                 for (int j = lo; j <= hi; j++) {
//                     RawIntegral += wave[j];
//                 }
//                 AvgIntegral = RawIntegral - AvgPed * WindowWidth;
//                 PulseFound = true;
//                 AvgIntegrals.push_back(AvgIntegral);
//             }
//         } else {
//             PulseFound = false;
//         }
//     }
//     return AvgIntegrals;
// }

void waveformGenerator::calcPulseParameters(const std::vector<Int_t> wave) {
  FADCConfig fConfig = saAnalyzer.getConfig();
  Int_t NPedSamps = fConfig.fNPedestalSamples;
  Int_t NSB = fConfig.fNSB;
  Int_t NSA = fConfig.fNSA;
  Int_t Threshold = fConfig.fSampThreshold;
  Int_t NSamples = (Int_t)wave.size();

  std::vector<Int_t> AvgIntegrals = {};

  // Pedestal: raw sum of first NPedSamps samples (matches analyzer's fSampPed)
  Double_t PedSum = 0.0;
  for (int i = 0; i < NPedSamps; i++)
    PedSum += wave[i];
  // Use the same formula as the analyzer:
  //   AvgPed (float) = PedSum / NPedSamps  — NOT integer division
  Double_t AvgPed = PedSum / NPedSamps;

  Bool_t PulseFound = false;

  for (int i = 0; i < wave.size(); i++) {
    if (wave[i] - AvgPed > Threshold) {
      if (!PulseFound) {
        // Integration window: [FirstSamp - NSB, FirstSamp + NSA - 1], clamped
        // to valid range
        Int_t lo = TMath::Max(i - NSB, 0);
        Int_t hi = TMath::Min(i + NSA - 1, NSamples - 1);
        // Pedestal-subtracted integral, using the same formula as the analyzer:
        //   PedSubInt = RawInt - PedSum * (window_width / NPedSamps)
        //             = RawInt - AvgPed * window_width
        Int_t WindowWidth = hi - lo + 1;
        Double_t AvgIntegral = 0.0;
        Int_t RawIntegral = 0;
        for (int j = lo; j <= hi; j++) {
          RawIntegral += wave[j];
        }
        // Match analyzer: PedSubInt = RawInt - PedSum * (WindowWidth /
        // NPedSamps)
        AvgIntegral = RawIntegral -
                      PedSum * (static_cast<Double_t>(WindowWidth) / NPedSamps);
        PulseFound = true;
        fWaveformInfo.PedSubIntegral.push_back(AvgIntegral);

        Int_t pulse_amp = 0;
        Int_t pulse_Time = 0;

        Int_t PeakBin = 0;
        Double_t PeakVal = wave[i] - AvgPed;
        for (Int_t nt = i + 1; nt < TMath::Min((i + NSA), int(NSamples));
             nt++) {
          if ((wave[nt] - AvgPed) < PeakVal && PeakBin == 0) {
            PeakBin = nt - 1;
          } else if (PeakBin == 0) {
            PeakVal = wave[nt] - AvgPed;
          }
        }
        if (PeakBin > 0) {
          pulse_amp = wave[PeakBin];
          Int_t Time = i * 64;
          Double_t VMid = (wave[PeakBin] - AvgPed) / 2.;
          for (Int_t nt = TMath::Max(i - fConfig.fNSB, 0);
               nt < TMath::Min(PeakBin, int(NSamples - 1)); nt++) {
            if (VMid >= (wave[nt] - AvgPed) && VMid < (wave[nt + 1] - AvgPed)) {
              Time = 64 * nt +
                     int(64 * (VMid - (wave[nt] - AvgPed)) /
                         ((wave[nt + 1] - AvgPed) - (wave[nt] - AvgPed)));
            }
          }
          pulse_Time = Time;
        } else {
          pulse_amp = wave[i];
          pulse_Time = 64 * i;
        }

        fWaveformInfo.PulseAmp.push_back(pulse_amp);
        fWaveformInfo.PulseTime.push_back(pulse_Time);
      }
    } else {
      PulseFound = false;
    }
  }
}
