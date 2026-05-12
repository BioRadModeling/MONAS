#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path


VALID_QUANTITIES = {"AMFSpectra", "AMF_yD", "AMF_yS"}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run or stage AMF TOPAS replay jobs for every .phsp/.header pair "
            "matching a pattern."
        )
    )
    parser.add_argument(
        "--app",
        default="./build/microdosimetry_with_LUTs",
        help="Path to the microdosimetry_with_LUTs executable.",
    )
    parser.add_argument(
        "--topas",
        default=None,
        help="Path to a TOPAS executable compiled with the AMF extension.",
    )
    parser.add_argument(
        "--lookup-root",
        default="lookup_tables",
        help="Lookup root containing AMF/tsed.dat.",
    )
    parser.add_argument(
        "--input-dir",
        default="input",
        help="Directory containing .phsp/.header phase-space pairs.",
    )
    parser.add_argument(
        "--staged-root",
        default="amf_runtime/staged_runs",
        help="Directory where per-depth AMF run folders will be created.",
    )
    parser.add_argument(
        "--quantity",
        required=True,
        choices=sorted(VALID_QUANTITIES),
        help="AMF quantity to score.",
    )
    parser.add_argument(
        "--pattern",
        default="*.phsp",
        help="Glob pattern used inside --input-dir.",
    )
    parser.add_argument(
        "--stage-only",
        action="store_true",
        help="Only create staged TOPAS run folders; do not execute TOPAS.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue to later phase-space pairs if one run fails.",
    )
    parser.add_argument(
        "--auto-slab-geometry",
        action="store_true",
        help=(
            "Derive a macroscopic slab geometry from each phase-space file: "
            "center Z on the phase-space Z midpoint and use a 1 mm slab "
            "thickness by default."
        ),
    )
    parser.add_argument(
        "--slab-half-x-mm",
        type=float,
        default=50.0,
        help="Half-width of the auto slab in X. Default is 50 mm for a 10 cm slab.",
    )
    parser.add_argument(
        "--slab-half-y-mm",
        type=float,
        default=50.0,
        help="Half-width of the auto slab in Y. Default is 50 mm for a 10 cm slab.",
    )
    parser.add_argument(
        "--slab-half-z-mm",
        type=float,
        default=0.5,
        help="Half-thickness of the auto slab in Z. Default is 0.5 mm for a 1 mm slab.",
    )
    parser.add_argument(
        "--auto-world-margin-mm",
        type=float,
        default=100.0,
        help="Extra world half-length margin around the auto slab. Default is 100 mm.",
    )

    args, amf_options = parser.parse_known_args()
    if amf_options and amf_options[0] == "--":
        amf_options = amf_options[1:]
    args.amf_options = amf_options
    return args


def phase_space_base(phsp_file: Path) -> Path:
    return phsp_file.with_suffix("")


def require_existing_file(path: Path, label: str):
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def discover_phase_spaces(input_dir: Path, pattern: str):
    phsp_files = sorted(input_dir.glob(pattern))
    pairs = []
    missing_headers = []

    for phsp_file in phsp_files:
        base = phase_space_base(phsp_file)
        header = base.with_suffix(".header")
        if header.is_file():
            pairs.append(base)
        else:
            missing_headers.append(header)

    return pairs, missing_headers


def phase_space_bounds_cm(phase_space_file: Path):
    bounds = {
        "min_x": None,
        "max_x": None,
        "min_y": None,
        "max_y": None,
        "min_z": None,
        "max_z": None,
    }
    row_count = 0

    with phase_space_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                x_cm, y_cm, z_cm = (float(parts[0]), float(parts[1]), float(parts[2]))
            except ValueError:
                continue

            row_count += 1
            if bounds["min_x"] is None:
                bounds.update(
                    {
                        "min_x": x_cm,
                        "max_x": x_cm,
                        "min_y": y_cm,
                        "max_y": y_cm,
                        "min_z": z_cm,
                        "max_z": z_cm,
                    }
                )
                continue

            bounds["min_x"] = min(bounds["min_x"], x_cm)
            bounds["max_x"] = max(bounds["max_x"], x_cm)
            bounds["min_y"] = min(bounds["min_y"], y_cm)
            bounds["max_y"] = max(bounds["max_y"], y_cm)
            bounds["min_z"] = min(bounds["min_z"], z_cm)
            bounds["max_z"] = max(bounds["max_z"], z_cm)

    if row_count == 0:
        raise ValueError(f"No numeric phase-space rows found: {phase_space_file}")

    return bounds


def auto_slab_options(base: Path, args):
    bounds = phase_space_bounds_cm(base.with_suffix(".phsp"))
    z_mid_mm = 10.0 * (bounds["min_z"] + bounds["max_z"]) / 2.0
    world_half_length_mm = max(
        args.slab_half_x_mm,
        args.slab_half_y_mm,
        abs(z_mid_mm) + args.slab_half_z_mm,
    ) + args.auto_world_margin_mm
    world_half_length_cm = world_half_length_mm / 10.0

    return [
        "--world-half-length-cm",
        f"{world_half_length_cm:.6g}",
        "--scoring-half-length-x-mm",
        f"{args.slab_half_x_mm:.6g}",
        "--scoring-half-length-y-mm",
        f"{args.slab_half_y_mm:.6g}",
        "--scoring-half-length-z-mm",
        f"{args.slab_half_z_mm:.6g}",
        "--scoring-x-mm",
        "0",
        "--scoring-y-mm",
        "0",
        "--scoring-z-mm",
        f"{z_mid_mm:.6g}",
    ]


def run_command(command, dry_run: bool):
    print(" ".join(str(part) for part in command), flush=True)
    if dry_run:
        return 0
    return subprocess.run(command, check=False).returncode


def main():
    args = parse_args()

    app = Path(args.app)
    lookup_root = Path(args.lookup_root)
    input_dir = Path(args.input_dir)
    staged_root = Path(args.staged_root)

    require_existing_file(app, "microdosimetry executable")
    require_existing_file(lookup_root / "AMF" / "tsed.dat", "AMF tsed.dat")

    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory not found: {input_dir}")

    if not args.stage_only:
        if args.topas is None:
            raise ValueError("--topas is required unless --stage-only is used")
        require_existing_file(Path(args.topas), "TOPAS executable")

    phase_space_bases, missing_headers = discover_phase_spaces(
        input_dir, args.pattern
    )

    if missing_headers:
        print("Skipping .phsp files with missing .header files:", file=sys.stderr)
        for header in missing_headers:
            print(f"  {header}", file=sys.stderr)

    if not phase_space_bases:
        raise ValueError(
            f"No usable .phsp/.header pairs found in {input_dir} "
            f"with pattern '{args.pattern}'"
        )

    failures = []
    for base in phase_space_bases:
        staged_run_dir = staged_root / f"{base.name}_{args.quantity}"
        if args.stage_only:
            command = [
                str(app),
                "AMF-stage",
                str(lookup_root),
                str(base),
                str(staged_run_dir),
                args.quantity,
            ]
        else:
            command = [
                str(app),
                "AMF-run",
                str(args.topas),
                str(lookup_root),
                str(base),
                str(staged_run_dir),
                args.quantity,
            ]

        if args.auto_slab_geometry:
            command.extend(auto_slab_options(base, args))

        command.extend(args.amf_options)
        exit_code = run_command(command, args.dry_run)

        if exit_code != 0:
            failures.append((base, exit_code))
            if not args.continue_on_error:
                break

    if failures:
        print("AMF batch completed with failures:", file=sys.stderr)
        for base, exit_code in failures:
            print(f"  {base}: exit code {exit_code}", file=sys.stderr)
        return 1

    print(f"AMF batch completed successfully for {len(phase_space_bases)} pair(s).")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
