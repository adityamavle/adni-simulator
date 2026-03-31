# adni-simulator

Synthetic Alzheimer's biomarker simulation project with ADNI-informed assumptions, not raw ADNI blood biomarker calibration.

## Overview

Generates a mixed longitudinal cohort of `100` patients across years `0-8`:
- `75` typical-progressing patients
- `25` fast/tau-early patients

The simulator generates smooth biomarker trajectories and then tests whether the mixed subtypes can still be separated from the resulting blood-marker patterns.

## Theory

This project is a rule-based synthetic simulation, not a model trained on raw ADNI blood biomarker tables.

The current idea is:
- `Abeta42` starts higher and tends to decrease over time
- `p_tau` starts lower and tends to increase over time
- fast or tau-early patients show an earlier and steeper `p_tau` rise

The ranges are chosen to stay in a biologically plausible synthetic regime inspired by ADNI-related literature:
- `Abeta42` is constrained to roughly `200-800 pg/mL`
- `p_tau` is constrained to roughly `20-110 pg/mL`

Each patient follows a smooth logistic-style progression curve with added patient-level variation and visit-level noise, so the groups overlap rather than separating perfectly.

## Implementation

The simulator is implemented in plain Python in `scripts/generate_synthetic_biomarkers.py`.

It uses:
- logistic-style progression curves for smooth longitudinal change
- random patient-level heterogeneity
- visit-level measurement noise
- a simple `p_tau` slope threshold as a first subtype-recovery check

## Outputs

- `output/synthetic_biomarkers.csv`
- `output/synthetic_patient_metrics.csv`
- `output/synthetic_summary.txt`

## Run

```powershell
python scripts/generate_synthetic_biomarkers.py
```

## Current Results

In the current synthetic experiment:
- typical group mean `p_tau`: `27.45 -> 55.79`
- fast group mean `p_tau`: `27.34 -> 70.75`
- simple `p_tau` slope classifier: accuracy `0.81`, sensitivity `0.84`, specificity `0.80`

This is a proof-of-concept synthetic result, not validation against raw ADNI blood data.

## ADNI Note

The local `DATA_ADNI` folder in this repo contains imaging metadata exports, not the blood biomarker tables needed to fit `Abeta42` or `p_tau` directly from ADNI. Because of that, the current pipeline is ADNI-inspired rather than ADNI-learned.
