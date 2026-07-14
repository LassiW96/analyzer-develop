#include "FadcStandaloneAnalyzer.h"

FadcStandaloneAnalyzer::FadcStandaloneAnalyzer() {
  fConfig = {/*safety config parameters*/
             4, 50, 2, 5, 2, 3};

  fCalib = {/*safety calibration parameters*/
            0.2441, 4.0 / 50.0};
}

FADCPulse FadcStandaloneAnalyzer::Analyze(const std::vector<Double_t> &samples,
                                          Double_t gain, Double_t tcal) const {
  Int_t fNSamples = samples.size();
  // std::vector<FADCPulse> pulses;
  FADCPulse data;

  data.fSampPulseInt.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseAmp.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseTime.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseTOT.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulsePileup.assign(fConfig.fMaxNPulses, 0);
  data.fSampPulseIntPedSub.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseIntMOLLERRaw.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseIntMOLLERVal.assign(fConfig.fMaxNPulses, 0.0);

  if (fNSamples <= 0 || fConfig.fNPedestalSamples <= 0)
    return data;
  auto GetIntegral = [&](Int_t lo, Int_t hi) {
    Double_t sum = 0.0;
    lo = TMath::Max(lo, 0);
    hi = TMath::Min(hi, fNSamples - 1);
    for (Int_t i = lo; i <= hi; i++)
      sum += samples[i];
    return sum;
  };

  auto GetSampleRaw = [&](Int_t i) { return samples[i]; };

  auto GetSample = [&](Int_t i) {
    return samples[i] - data.fSampPed / fConfig.fNPedestalSamples;
  };

  data.fHasMulti = kTRUE;
  data.fSampPed = GetIntegral(0, fConfig.fNPedestalSamples - 1);

  Int_t NS = fConfig.fNPedestalSamples - 1;
  data.fNSampPulses = 0;

  if (fConfig.fSampThreshold == 0) {
    data.fSampPulseInt[data.fNSampPulses] =
        GetIntegral(TMath::Max(NS - fConfig.fNSB, 0),
                    TMath::Min((NS + fConfig.fNSA - 1), int(fNSamples - 1)));
    data.fNPeakSamples =
        TMath::Min((NS + fConfig.fNSA - 1), int(fNSamples - 1)) -
        TMath::Max(NS - fConfig.fNSB, 0) + 1;
    data.fPeakPedestalRatio =
        1.0 * data.fNPeakSamples / fConfig.fNPedestalSamples;
    data.fSampPulseIntPedSub[data.fNSampPulses] =
        data.fSampPulseInt[data.fNSampPulses] -
        data.fSampPed * data.fPeakPedestalRatio;
    data.fSampPulseIntMOLLERRaw[data.fNSampPulses] =
        data.fSampPulseInt[data.fNSampPulses] * fCalib.ChanTomV *
        fCalib.pC_Conv;

    data.fSampPulseIntMOLLERVal[data.fNSampPulses] =
        data.fSampPulseIntPedSub[data.fNSampPulses] * fCalib.ChanTomV *
        fCalib.pC_Conv * gain;

    data.fSampPulseAmp[data.fNSampPulses] = GetSampleRaw(NS);
    data.fSampPulseTime[data.fNSampPulses] = 64 * NS;
    data.fNSampPulses = 1;
  } else {
    Bool_t CheckSampBelowThres = kFALSE;
    Int_t LastFourSampPed =
        GetIntegral(fNSamples - fConfig.fNPedestalSamples, fNSamples - 1);
    if (LastFourSampPed < data.fSampPed) {
      data.fSampPed = LastFourSampPed;
      CheckSampBelowThres = kTRUE;
    }

    while (NS < int(fNSamples) && data.fNSampPulses < fConfig.fMaxNPulses) {
      if (CheckSampBelowThres) {
        if (GetSample(NS) < fConfig.fSampThreshold)
          CheckSampBelowThres = kFALSE;
      } else {
        Int_t ns_found = 0;
        for (Int_t nt = NS; nt < TMath::Min(NS + fConfig.fNSAT, int(fNSamples));
             nt++) {
          if (GetSample(nt) > fConfig.fSampThreshold)
            ns_found++;
        }
        if (ns_found == fConfig.fNSAT) {
          data.fSampPulseInt[data.fNSampPulses] = GetIntegral(
              TMath::Max(NS - fConfig.fNSB, 0),
              TMath::Min(NS + fConfig.fNSA - 1, int(fNSamples - 1)));
          data.fNPeakSamples =
              TMath::Min((NS + fConfig.fNSA - 1), int(fNSamples - 1)) -
              TMath::Max(NS - fConfig.fNSB, 0) + 1;
          data.fPeakPedestalRatio =
              1.0 * data.fNPeakSamples / fConfig.fNPedestalSamples;

          data.fSampPulseIntPedSub[data.fNSampPulses] =
              data.fSampPulseInt[data.fNSampPulses] -
              data.fSampPed * data.fPeakPedestalRatio;

          data.fSampPulseIntMOLLERRaw[data.fNSampPulses] =
              data.fSampPulseInt[data.fNSampPulses] * fCalib.ChanTomV *
              fCalib.pC_Conv;

          data.fSampPulseIntMOLLERVal[data.fNSampPulses] =
              data.fSampPulseIntPedSub[data.fNSampPulses] * fCalib.ChanTomV *
              fCalib.pC_Conv * gain;

          // Peak and time calculation
          data.fSampPulseAmp[data.fNSampPulses] = 0;
          data.fSampPulseTime[data.fNSampPulses] = 0;

          Int_t PeakBin = 0;
          Double_t PeakVal = GetSample(NS);
          for (Int_t nt = NS + 1;
               nt < TMath::Min((NS + fConfig.fNSA), int(fNSamples)); nt++) {
            if (GetSample(nt) < PeakVal && PeakBin == 0) {
              PeakBin = nt - 1;
            } else if (PeakBin == 0) {
              PeakVal = GetSample(nt);
            }
          }
          if (PeakBin > 0) {
            data.fSampPulseAmp[data.fNSampPulses] = GetSampleRaw(PeakBin);
            Int_t Time = NS * 64;
            Double_t VMid = (GetSample(PeakBin)) / 2.;
            for (Int_t nt = TMath::Max(NS - fConfig.fNSB, 0);
                 nt < TMath::Min(PeakBin, int(fNSamples - 1)); nt++) {
              if (VMid >= GetSample(nt) && VMid < GetSample(nt + 1)) {
                Time = 64 * nt + int(64 * (VMid - GetSample(nt)) /
                                     (GetSample(nt + 1) - GetSample(nt)));
              }
            }
            data.fSampPulseTime[data.fNSampPulses] = Time;
          } else {
            data.fSampPulseAmp[data.fNSampPulses] = GetSampleRaw(NS);
            data.fSampPulseTime[data.fNSampPulses] = 64 * NS;
          }

          // Pileup calculations
          Int_t hi = TMath::Min((NS + fConfig.fNSA - 1), int(fNSamples - 1));
          Int_t tot = 0;
          for (Int_t nt = NS; nt <= hi; nt++) {
            if (GetSample(nt) > fConfig.fSampThreshold)
              tot++;
          }
          data.fSampPulseTOT[data.fNSampPulses] = 64.0 * tot;

          Int_t pileup = 0;
          for (Int_t nt = NS + 1; nt < hi; nt++) {
            if (GetSample(nt) > fConfig.fSampThreshold &&
                GetSample(nt) > GetSample(nt - 1) &&
                GetSample(nt) > GetSample(nt + 1)) {
              pileup = 1;
              break;
            }
          }
          data.fSampPulsePileup[data.fNSampPulses] = pileup;

          data.fNSampPulses++;
          NS = NS + fConfig.fNSA;
          CheckSampBelowThres = kTRUE;
        }
      }
      NS++;
    }
  }

  // if(data.fNSampPulses > 0) pulses.push_back(data);
  return data;
}