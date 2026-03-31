import argparse
import csv
import math
import os
import random
import statistics


def logistic(t, midpoint, slope):
    return 1.0 / (1.0 + math.exp(-slope * (t - midpoint)))


def clamp(value, lower, upper):
    return max(lower, min(upper, value))


def simulate_patient(patient_id, subtype, years):
    rows = []

    abeta_start = random.uniform(610.0, 790.0)
    abeta_drop = random.uniform(170.0, 330.0)
    abeta_midpoint = random.uniform(3.0, 5.8)
    abeta_slope = random.uniform(0.4, 0.85)
    abeta_patient_shift = random.gauss(0.0, 15.0)

    if subtype == "fast":
        ptau_start = random.uniform(20.0, 30.0)
        ptau_rise = random.uniform(34.0, 62.0)
        ptau_midpoint = random.uniform(2.8, 4.8)
        ptau_slope = random.uniform(0.55, 1.0)
    else:
        ptau_start = random.uniform(20.0, 30.0)
        ptau_rise = random.uniform(22.0, 50.0)
        ptau_midpoint = random.uniform(3.6, 6.2)
        ptau_slope = random.uniform(0.4, 0.9)

    ptau_patient_shift = random.gauss(0.0, 4.5)
    visit_noise_scale = random.uniform(2.5, 6.5)

    for year in years:
        abeta_signal = (
            abeta_start
            - abeta_drop * logistic(year, abeta_midpoint, abeta_slope)
            + abeta_patient_shift
        )
        ptau_signal = (
            ptau_start
            + ptau_rise * logistic(year, ptau_midpoint, ptau_slope)
            + ptau_patient_shift
        )

        year_specific_ptau_noise = random.gauss(0.0, visit_noise_scale)
        year_specific_abeta_noise = random.gauss(0.0, 22.0)

        abeta_value = clamp(abeta_signal + year_specific_abeta_noise, 200.0, 800.0)
        ptau_value = clamp(ptau_signal + year_specific_ptau_noise, 20.0, 110.0)

        rows.append(
            {
                "patient_id": patient_id,
                "subtype": subtype,
                "year": year,
                "abeta42_pg_ml": round(abeta_value, 2),
                "p_tau_pg_ml": round(ptau_value, 2),
            }
        )

    return rows


def build_dataset(seed):
    random.seed(seed)
    years = list(range(9))
    rows = []

    for patient_index in range(1, 76):
        patient_id = f"T{patient_index:03d}"
        rows.extend(simulate_patient(patient_id, "typical", years))

    for patient_index in range(1, 26):
        patient_id = f"F{patient_index:03d}"
        rows.extend(simulate_patient(patient_id, "fast", years))

    return rows


def summarize(rows):
    by_subtype = {"typical": [], "fast": []}

    for row in rows:
        by_subtype[row["subtype"]].append(row)

    summary = {}
    for subtype, values in by_subtype.items():
        baseline_ptau = [row["p_tau_pg_ml"] for row in values if row["year"] == 0]
        final_ptau = [row["p_tau_pg_ml"] for row in values if row["year"] == 8]
        baseline_abeta = [row["abeta42_pg_ml"] for row in values if row["year"] == 0]
        final_abeta = [row["abeta42_pg_ml"] for row in values if row["year"] == 8]

        summary[subtype] = {
            "patients": len({row["patient_id"] for row in values}),
            "baseline_ptau_mean": round(statistics.mean(baseline_ptau), 2),
            "final_ptau_mean": round(statistics.mean(final_ptau), 2),
            "baseline_abeta_mean": round(statistics.mean(baseline_abeta), 2),
            "final_abeta_mean": round(statistics.mean(final_abeta), 2),
        }

    return summary


def patient_level_metrics(rows):
    patients = {}

    for row in rows:
        patient_id = row["patient_id"]
        patients.setdefault(
            patient_id,
            {
                "patient_id": patient_id,
                "subtype": row["subtype"],
                "year_to_ptau": {},
                "year_to_abeta": {},
            },
        )
        patients[patient_id]["year_to_ptau"][row["year"]] = row["p_tau_pg_ml"]
        patients[patient_id]["year_to_abeta"][row["year"]] = row["abeta42_pg_ml"]

    metrics = []
    for patient in patients.values():
        baseline_ptau = patient["year_to_ptau"][0]
        final_ptau = patient["year_to_ptau"][8]
        baseline_abeta = patient["year_to_abeta"][0]
        final_abeta = patient["year_to_abeta"][8]

        ptau_delta = round(final_ptau - baseline_ptau, 2)
        abeta_delta = round(final_abeta - baseline_abeta, 2)
        ptau_slope = round(ptau_delta / 8.0, 2)
        abeta_slope = round(abeta_delta / 8.0, 2)

        metrics.append(
            {
                "patient_id": patient["patient_id"],
                "subtype": patient["subtype"],
                "baseline_p_tau_pg_ml": baseline_ptau,
                "final_p_tau_pg_ml": final_ptau,
                "p_tau_delta_pg_ml": ptau_delta,
                "p_tau_slope_pg_ml_per_year": ptau_slope,
                "baseline_abeta42_pg_ml": baseline_abeta,
                "final_abeta42_pg_ml": final_abeta,
                "abeta42_delta_pg_ml": abeta_delta,
                "abeta42_slope_pg_ml_per_year": abeta_slope,
            }
        )

    return sorted(metrics, key=lambda item: item["patient_id"])


def classify_by_ptau_slope(metrics):
    typical_slopes = [
        item["p_tau_slope_pg_ml_per_year"] for item in metrics if item["subtype"] == "typical"
    ]
    fast_slopes = [
        item["p_tau_slope_pg_ml_per_year"] for item in metrics if item["subtype"] == "fast"
    ]
    threshold = (statistics.mean(typical_slopes) + statistics.mean(fast_slopes)) / 2.0

    counts = {"tp": 0, "tn": 0, "fp": 0, "fn": 0}
    classified = []
    for item in metrics:
        predicted = "fast" if item["p_tau_slope_pg_ml_per_year"] >= threshold else "typical"
        actual = item["subtype"]
        classified.append({**item, "predicted_subtype": predicted})

        if actual == "fast" and predicted == "fast":
            counts["tp"] += 1
        elif actual == "typical" and predicted == "typical":
            counts["tn"] += 1
        elif actual == "typical" and predicted == "fast":
            counts["fp"] += 1
        else:
            counts["fn"] += 1

    total = len(metrics)
    accuracy = (counts["tp"] + counts["tn"]) / total
    sensitivity = counts["tp"] / (counts["tp"] + counts["fn"])
    specificity = counts["tn"] / (counts["tn"] + counts["fp"])

    results = {
        "threshold_p_tau_slope_pg_ml_per_year": round(threshold, 2),
        "accuracy": round(accuracy, 3),
        "sensitivity": round(sensitivity, 3),
        "specificity": round(specificity, 3),
        "true_positive": counts["tp"],
        "true_negative": counts["tn"],
        "false_positive": counts["fp"],
        "false_negative": counts["fn"],
    }

    return classified, results


def write_csv(rows, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fieldnames = ["patient_id", "subtype", "year", "abeta42_pg_ml", "p_tau_pg_ml"]

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_patient_metrics(metrics, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fieldnames = [
        "patient_id",
        "subtype",
        "baseline_p_tau_pg_ml",
        "final_p_tau_pg_ml",
        "p_tau_delta_pg_ml",
        "p_tau_slope_pg_ml_per_year",
        "baseline_abeta42_pg_ml",
        "final_abeta42_pg_ml",
        "abeta42_delta_pg_ml",
        "abeta42_slope_pg_ml_per_year",
        "predicted_subtype",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metrics)


def write_summary(summary, classification_results, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("Synthetic biomarker simulation summary\n")
        handle.write("\n")
        for subtype, values in summary.items():
            handle.write(
                f"{subtype}: patients={values['patients']}, "
                f"baseline_p_tau_mean={values['baseline_ptau_mean']}, "
                f"final_p_tau_mean={values['final_ptau_mean']}, "
                f"baseline_abeta_mean={values['baseline_abeta_mean']}, "
                f"final_abeta_mean={values['final_abeta_mean']}\n"
            )

        handle.write("\n")
        handle.write("Simple subtype recovery using p_tau slope threshold\n")
        handle.write(
            f"threshold_p_tau_slope_pg_ml_per_year="
            f"{classification_results['threshold_p_tau_slope_pg_ml_per_year']}\n"
        )
        handle.write(f"accuracy={classification_results['accuracy']}\n")
        handle.write(f"sensitivity={classification_results['sensitivity']}\n")
        handle.write(f"specificity={classification_results['specificity']}\n")
        handle.write(f"true_positive={classification_results['true_positive']}\n")
        handle.write(f"true_negative={classification_results['true_negative']}\n")
        handle.write(f"false_positive={classification_results['false_positive']}\n")
        handle.write(f"false_negative={classification_results['false_negative']}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a synthetic longitudinal Alzheimer's biomarker dataset."
    )
    parser.add_argument(
        "--output",
        default=os.path.join("output", "synthetic_biomarkers.csv"),
        help="Path to the output CSV file.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible output.",
    )
    parser.add_argument(
        "--patient-metrics-output",
        default=os.path.join("output", "synthetic_patient_metrics.csv"),
        help="Path to the patient-level metrics CSV file.",
    )
    parser.add_argument(
        "--summary-output",
        default=os.path.join("output", "synthetic_summary.txt"),
        help="Path to the text summary output.",
    )
    args = parser.parse_args()

    rows = build_dataset(args.seed)
    write_csv(rows, args.output)
    summary = summarize(rows)
    metrics = patient_level_metrics(rows)
    classified_metrics, classification_results = classify_by_ptau_slope(metrics)
    write_patient_metrics(classified_metrics, args.patient_metrics_output)
    write_summary(summary, classification_results, args.summary_output)

    print(f"Wrote {len(rows)} rows to {args.output}")
    for subtype, values in summary.items():
        print(
            f"{subtype}: "
            f"patients={values['patients']}, "
            f"baseline_p_tau_mean={values['baseline_ptau_mean']}, "
            f"final_p_tau_mean={values['final_ptau_mean']}, "
            f"baseline_abeta_mean={values['baseline_abeta_mean']}, "
            f"final_abeta_mean={values['final_abeta_mean']}"
        )
    print(
        "classification: "
        f"threshold_p_tau_slope_pg_ml_per_year="
        f"{classification_results['threshold_p_tau_slope_pg_ml_per_year']}, "
        f"accuracy={classification_results['accuracy']}, "
        f"sensitivity={classification_results['sensitivity']}, "
        f"specificity={classification_results['specificity']}"
    )


if __name__ == "__main__":
    main()
