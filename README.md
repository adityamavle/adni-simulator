# adni-simulator

Synthetic Alzheimer's biomarker simulation project with ADNI-informed assumptions, not raw ADNI blood biomarker calibration.

## What It Does

Generates a mixed longitudinal cohort of `100` patients across years `0-8`:
- `75` typical-progressing patients
- `25` fast/tau-early patients

Outputs:
- `output/synthetic_biomarkers.csv`
- `output/synthetic_patient_metrics.csv`
- `output/synthetic_summary.txt`

## Run

```powershell
python scripts/generate_synthetic_biomarkers.py
```

## Current Results

In the current synthetic experiment:
- typical group mean `p_tau`: `25.48 -> 55.21`
- fast group mean `p_tau`: `27.39 -> 91.50`
- simple `p_tau` slope classifier: accuracy `1.0`, sensitivity `1.0`, specificity `1.0`

This is a proof-of-concept synthetic result, not validation against raw ADNI blood data.
