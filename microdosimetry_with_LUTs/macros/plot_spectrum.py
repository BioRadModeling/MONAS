#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


ALLOWED_COLUMNS = {
    "y_keV_per_um",
    "f_y",
    "yf_y",
    "d_y",
    "yd_y",
}


def read_csv(csv_path: Path, xcol: str, ycol: str):
    x = []
    y = []

    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)

        missing = [c for c in (xcol, ycol) if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"Missing required column(s): {missing}. "
                f"Available columns: {reader.fieldnames}"
            )

        for row in reader:
            try:
                xv = float(row[xcol])
                yv = float(row[ycol])
            except Exception:
                continue

            x.append(xv)
            y.append(yv)

    if not x:
        raise ValueError("No valid numeric rows were read from the CSV.")

    return x, y


def main():
    parser = argparse.ArgumentParser(
        description="Plot a selected spectrum quantity from poly_spectrum.csv"
    )
    parser.add_argument("csv_file", help="Path to poly_spectrum.csv")
    parser.add_argument("--x", default="y_keV_per_um",
                        help="Column to use on x-axis")
    parser.add_argument("--y", required=True,
                        help="Column to use on y-axis")
    parser.add_argument("--output-dir", required=True,
                        help="Directory where JPEG will be saved")
    parser.add_argument("--title", default=None,
                        help="Optional custom plot title")
    parser.add_argument("--logx", action="store_true",
                        help="Use logarithmic x-axis")
    parser.add_argument("--logy", action="store_true",
                        help="Use logarithmic y-axis")
    parser.add_argument("--dpi", type=int, default=300,
                        help="JPEG resolution")
    args = parser.parse_args()

    xcol = args.x
    ycol = args.y

    if xcol not in ALLOWED_COLUMNS:
        raise ValueError(f"Unsupported x column: {xcol}")
    if ycol not in ALLOWED_COLUMNS:
        raise ValueError(f"Unsupported y column: {ycol}")

    csv_path = Path(args.csv_file).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    x, y = read_csv(csv_path, xcol, ycol)

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, linewidth=1.2)

    plt.xlabel(xcol)
    plt.ylabel(ycol)

    if args.title:
        plt.title(args.title)
    else:
        plt.title(f"{ycol} vs {xcol}")

    if args.logx:
        plt.xscale("log")
    if args.logy:
        plt.yscale("log")

    plt.tight_layout()

    out_name = f"{ycol}_vs_{xcol}.jpg"
    out_path = output_dir / out_name
    plt.savefig(out_path, format="jpeg", dpi=args.dpi)
    plt.close()

    print(f"Saved plot to: {out_path}")


if __name__ == "__main__":
    main()