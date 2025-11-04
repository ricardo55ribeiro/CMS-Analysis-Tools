# RooFit_Fits

This directory contains the RooFit-based fitting methods developed during a CMS Internship at LIP, analyzing data from the RUN 3 of the CMS experiment at the LHC (CERN).

The scripts perform data and Monte Carlo fits for B mesons and charmonium states and include tools for evaluating systematic uncertainties.

---

## Folder Structure

```
RooFit_Fits/
├── Data_Fits/                  # Fits to real CMS data (B⁺, B⁰, ψ(2S), X(3872))
│   ├── data_fit_Bplus_erfc.C
│   ├── data_fit_Bplus_without_erfc.C
│   ├── data_fit_Bzero.C
│   └── data_fit_PSI_and_X.C
│
├── MonteCarlo_Fits/            # Fits to Monte Carlo simulated samples
│   └── mc_signal_fits.C
│
└── Systematic_Variations/      # Scripts for uncertainty estimation and model variations
    ├── Bplus_variations.C
    └── uncertainty_variations.C
```

---

## Description

### Data_Fits
Performs unbinned maximum-likelihood fits to CMS data using various signal and background models to extract particle yields and parameters.

### MonteCarlo_Fits
Fits Monte Carlo signal samples to determine resolution and shape parameters applied in data fits.

### Systematic_Variations
Evaluates systematic uncertainties by varying model assumptions such as background shape, signal model, or fit range.

---

## Notes

- Scripts follow CMS RUN 3 conventions for B-meson reconstruction.  
- ROOT files referenced in the macros are not included due to data-access restrictions.  
- The code is provided for educational and methodological purposes only.
