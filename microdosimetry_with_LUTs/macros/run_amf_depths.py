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
        "--source-topas",
        required=True,
        help="Resolved TOPAS input file that generated the phase-space pairs.",
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
        "--scoring-radius-mm",
        required=True,
        type=float,
        help="Radius of the generated water AMF scoring sphere.",
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
    source_topas = Path(args.source_topas)
    staged_root = Path(args.staged_root)

    require_existing_file(app, "microdosimetry executable")
    require_existing_file(lookup_root / "AMF" / "tsed.dat", "AMF tsed.dat")
    require_existing_file(source_topas, "source TOPAS input file")

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
                str(source_topas),
                str(staged_run_dir),
                args.quantity,
                "--scoring-radius-mm",
                f"{args.scoring_radius_mm:g}",
            ]
        else:
            command = [
                str(app),
                "AMF-run",
                str(args.topas),
                str(lookup_root),
                str(base),
                str(source_topas),
                str(staged_run_dir),
                args.quantity,
                "--scoring-radius-mm",
                f"{args.scoring_radius_mm:g}",
            ]

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
