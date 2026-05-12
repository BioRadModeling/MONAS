#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def infer_amf_y_midpoints(bin_count: int):
    """AMFSpectra uses 400 log bins from 10^-3 to 10^5 keV/um."""
    y_edges = []
    y_power = -3.0
    y_step = 0.02

    for _ in range(bin_count + 1):
        y_edges.append(10.0 ** y_power)
        y_power += y_step

    return [
        0.5 * (y_edges[i] + y_edges[i + 1])
        for i in range(bin_count)
    ]


def read_amf_spectrum(csv_path: Path, voxel_row_index: int):
    rows = list(csv.reader(csv_path.open("r", newline="")))
    if len(rows) < 4:
        raise ValueError("AMF spectra CSV must contain at least four rows.")

    midpoint_row = rows[2]
    data_rows = rows[3:]

    if voxel_row_index < 0 or voxel_row_index >= len(data_rows):
        raise ValueError(
            f"voxel_row_index={voxel_row_index} is out of range. "
            f"Available data rows: 0..{len(data_rows) - 1}"
        )

    data_row = data_rows[voxel_row_index]
    voxel_index = tuple(data_row[:3])
    yd_y = [float(value) for value in data_row[3:]]

    try:
        y_midpoints = [float(value) for value in midpoint_row[3:]]
    except ValueError:
        y_midpoints = []

    if (
        len(y_midpoints) != len(yd_y)
        or not y_midpoints
        or all(value == 0.0 for value in y_midpoints)
    ):
        y_midpoints = infer_amf_y_midpoints(len(yd_y))

    return y_midpoints, yd_y, voxel_index


def main():
    parser = argparse.ArgumentParser(
        description="Plot yd(y) vs y from an AMFSpectra _MicrodosimetricSpectra.csv file."
    )
    parser.add_argument(
        "csv_file",
        help="Path to *_MicrodosimetricSpectra.csv",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output image path. Default: <csv stem>_ydy_vs_y.png",
    )
    parser.add_argument(
        "--voxel-row-index",
        type=int,
        default=0,
        help="Zero-based data row to plot when the file contains multiple voxels.",
    )
    parser.add_argument(
        "--title",
        default="AMF yd(y) vs y",
        help="Plot title.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image resolution.",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv_file).resolve()
    if not csv_path.is_file():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")

    output_path = (
        Path(args.output).resolve()
        if args.output
        else csv_path.with_name(f"{csv_path.stem}_ydy_vs_y.png")
    )

    y_midpoints, yd_y, voxel_index = read_amf_spectrum(
        csv_path,
        args.voxel_row_index,
    )

    plt.figure(figsize=(8, 5.5))
    plt.plot(y_midpoints, yd_y, linewidth=1.8)
    plt.xscale("log")
    plt.xlabel("y (keV/um)")
    plt.ylabel("yd(y)")
    plt.title(f"{args.title} - voxel {voxel_index}")
    plt.grid(True, which="both", alpha=0.25)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=args.dpi)
    plt.close()

    print(f"Saved plot to: {output_path}")


if __name__ == "__main__":
    main()
