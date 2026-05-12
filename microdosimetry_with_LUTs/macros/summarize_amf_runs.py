#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path


FIELDNAMES = [
    "run_directory",
    "quantity",
    "phase_space",
    "phase_space_header",
    "world_half_length_cm",
    "scoring_half_length_x_mm",
    "scoring_half_length_y_mm",
    "scoring_half_length_z_mm",
    "scoring_trans_x_mm",
    "scoring_trans_y_mm",
    "scoring_trans_z_mm",
    "domain_radius_um",
    "nucleus_radius_um",
    "beta_ref_per_gy2",
    "electron_cut_m",
    "stopping_power",
    "step_calculator",
    "scorer_output",
    "scorer_output_exists",
    "spectra_csv",
    "spectra_csv_exists",
    "all_expected_outputs_found",
    "stdout_log",
    "stdout_log_exists",
    "stderr_log",
    "stderr_log_exists",
    "manifest_file",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize expected AMF output files across staged runs."
    )
    parser.add_argument(
        "--staged-root",
        default="amf_runtime/staged_runs",
        help="Root directory containing AMF staged run folders.",
    )
    parser.add_argument(
        "--output-csv",
        default="amf_runtime/amf_run_summary.csv",
        help="CSV path where the run summary will be written.",
    )
    parser.add_argument(
        "--quantity",
        choices=["AMFSpectra", "AMF_yD", "AMF_yS"],
        default=None,
        help="Only summarize staged runs for this AMF quantity.",
    )
    parser.add_argument(
        "--print-missing",
        action="store_true",
        help="Print runs with missing expected outputs to stderr.",
    )
    return parser.parse_args()


def strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1]
    return value


def parse_manifest(manifest_path: Path):
    data = {}
    with manifest_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("AMF staged run manifest"):
                continue
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            data[key.strip()] = strip_quotes(value)
    return data


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def expected_files(run_directory: Path, manifest: dict):
    output_file = manifest.get("output_file", "")
    quantity = manifest.get("quantity", "")

    scorer_output = run_directory / f"{output_file}.csv"
    spectra_csv = ""
    if quantity == "AMFSpectra":
        spectra_csv = run_directory / f"{output_file}_MicrodosimetricSpectra.csv"

    return scorer_output, spectra_csv


def summarize_manifest(manifest_path: Path):
    manifest = parse_manifest(manifest_path)
    run_directory = Path(
        manifest.get("run_directory", str(manifest_path.parent))
    )

    scorer_output, spectra_csv = expected_files(run_directory, manifest)
    scorer_exists = scorer_output.exists()
    spectra_exists = True if spectra_csv == "" else spectra_csv.exists()

    stdout_log = run_directory / "topas_stdout.log"
    stderr_log = run_directory / "topas_stderr.log"

    return {
        "run_directory": str(run_directory),
        "quantity": manifest.get("quantity", ""),
        "phase_space": manifest.get("phase_space", ""),
        "phase_space_header": manifest.get("phase_space_header", ""),
        "world_half_length_cm": manifest.get("world_half_length_cm", ""),
        "scoring_half_length_x_mm": manifest.get(
            "scoring_half_length_x_mm",
            manifest.get("scoring_half_length_mm", ""),
        ),
        "scoring_half_length_y_mm": manifest.get(
            "scoring_half_length_y_mm",
            manifest.get("scoring_half_length_mm", ""),
        ),
        "scoring_half_length_z_mm": manifest.get(
            "scoring_half_length_z_mm",
            manifest.get("scoring_half_length_mm", ""),
        ),
        "scoring_trans_x_mm": manifest.get("scoring_trans_x_mm", ""),
        "scoring_trans_y_mm": manifest.get("scoring_trans_y_mm", ""),
        "scoring_trans_z_mm": manifest.get("scoring_trans_z_mm", ""),
        "domain_radius_um": manifest.get("domain_radius_um", ""),
        "nucleus_radius_um": manifest.get("nucleus_radius_um", ""),
        "beta_ref_per_gy2": manifest.get("beta_ref_per_gy2", ""),
        "electron_cut_m": manifest.get("electron_cut_m", ""),
        "stopping_power": manifest.get("stopping_power", ""),
        "step_calculator": manifest.get("step_calculator", ""),
        "scorer_output": str(scorer_output),
        "scorer_output_exists": bool_text(scorer_exists),
        "spectra_csv": "" if spectra_csv == "" else str(spectra_csv),
        "spectra_csv_exists": "" if spectra_csv == "" else bool_text(spectra_exists),
        "all_expected_outputs_found": bool_text(scorer_exists and spectra_exists),
        "stdout_log": str(stdout_log),
        "stdout_log_exists": bool_text(stdout_log.exists()),
        "stderr_log": str(stderr_log),
        "stderr_log_exists": bool_text(stderr_log.exists()),
        "manifest_file": str(manifest_path),
    }


def main():
    args = parse_args()
    staged_root = Path(args.staged_root)
    output_csv = Path(args.output_csv)

    if not staged_root.is_dir():
        raise NotADirectoryError(f"Staged root not found: {staged_root}")

    manifests = sorted(staged_root.glob("*/amf_run_manifest.txt"))
    rows = []
    for manifest_path in manifests:
        row = summarize_manifest(manifest_path)
        if args.quantity is not None and row["quantity"] != args.quantity:
            continue
        rows.append(row)

    if not rows:
        raise ValueError(f"No AMF run manifests found under {staged_root}")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    missing_rows = [
        row for row in rows if row["all_expected_outputs_found"] != "true"
    ]

    print(f"Wrote AMF run summary: {output_csv}")
    print(f"Summarized runs: {len(rows)}")
    print(f"Runs with missing expected outputs: {len(missing_rows)}")

    if args.print_missing and missing_rows:
        print("Missing expected AMF outputs:", file=sys.stderr)
        for row in missing_rows:
            print(
                f"  {row['run_directory']} ({row['quantity']}): "
                f"{row['scorer_output']}",
                file=sys.stderr,
            )
            if row["spectra_csv"]:
                print(f"    {row['spectra_csv']}", file=sys.stderr)

    return 1 if missing_rows else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
