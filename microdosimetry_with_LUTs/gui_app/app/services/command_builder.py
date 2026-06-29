from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from shlex import quote

from app.state import AppState


@dataclass(frozen=True)
class CommandSpec:
    label: str
    argv: list[str]
    workdir: Path | None = None

    def render(self) -> str:
        return " ".join(quote(part) for part in self.argv)


def build_command_specs(state: AppState) -> list[CommandSpec]:
    executable = str(state.executable_path)
    lookup_root = str(state.lookup_root)
    phase_space_file = str(state.phase_space_file)
    output_dir = str(state.output_dir)

    if state.approach == "amf":
        return _build_amf_command_specs(state)

    if state.approach == "spectrum":
        if state.spectrum_family == "Cartechini":
            return [
                CommandSpec(
                    label="Build spectrum (Cartechini)",
                    argv=[
                        executable,
                        "build-spectrum",
                        lookup_root,
                        "Cartechini",
                        state.cartechini_particle,
                        state.cartechini_radius,
                        phase_space_file,
                        output_dir,
                        "--rebin-samples",
                        str(state.rebin_samples),
                        "--rebin-seed",
                        str(state.rebin_seed),
                    ],
                )
            ]

        return [
            CommandSpec(
                label="Build spectrum (DeCunha)",
                argv=[
                    executable,
                    "build-spectrum",
                    lookup_root,
                    "DeCunha",
                    state.decunha_voxel_size,
                    state.decunha_energy_grid,
                    phase_space_file,
                    output_dir,
                    "--rebin-samples",
                    str(state.rebin_samples),
                    "--rebin-seed",
                    str(state.rebin_seed),
                ],
            )
        ]

    if state.approach == "means":
        commands: list[CommandSpec] = []

        if state.enable_let:
            commands.append(
                CommandSpec(
                    label="LET summary",
                    argv=[executable, "LET", lookup_root, phase_space_file, output_dir],
                )
            )

        if state.enable_magini:
            commands.append(
                CommandSpec(
                    label="Magini summary",
                    argv=[executable, "Magini", lookup_root, phase_space_file, output_dir],
                )
            )

        if state.enable_inaniwa:
            commands.append(
                CommandSpec(
                    label="Inaniwa summary",
                    argv=[executable, "Inaniwa", lookup_root, phase_space_file, output_dir],
                )
            )

        return commands

    return []


def _build_amf_command_specs(state: AppState) -> list[CommandSpec]:
    commands: list[CommandSpec] = []

    commands.append(
        CommandSpec(
            label="Create AMF output directory",
            argv=["/bin/mkdir", "-p", str(state.amf_output_dir)],
            workdir=state.build_workdir,
        )
    )

    if state.amf_stopping_power == "ExternalTable":
        _append_copy_command(
            commands,
            "Stage AMF external stopping-power table",
            state.amf_external_stopping_power_file,
            state.amf_output_dir / "StoppingPower.txt",
            state.build_workdir,
        )

    amf_args = [
        str(state.executable_path),
        "AMF-run",
        str(state.lookup_root),
        str(state.amf_phase_space_base),
        str(state.amf_source_topas_file),
        str(state.amf_output_dir),
        state.amf_quantity,
        "--scoring-radius-mm",
        f"{state.amf_scoring_radius_mm:g}",
        "--domain-radius",
        f"{state.amf_domain_radius_um:g}",
        "--stopping-power",
        state.amf_stopping_power,
        "--output-file",
        amf_output_stem(state),
    ]
    if state.amf_disable_phase_space_precheck:
        amf_args.append("--no-phase-space-precheck")
    else:
        amf_args.append("--phase-space-precheck")

    commands.append(
        CommandSpec(
            label="Run AMF phase-space replay",
            argv=amf_args,
            workdir=state.build_workdir,
        )
    )
    return commands


def _append_copy_command(
    commands: list[CommandSpec],
    label: str,
    source: Path,
    destination: Path,
    workdir: Path,
) -> None:
    if _same_path(source, destination):
        return
    commands.append(
        CommandSpec(
            label=label,
            argv=["/bin/cp", str(source), str(destination)],
            workdir=workdir,
        )
    )


def _same_path(left: Path, right: Path) -> bool:
    return left.resolve(strict=False) == right.resolve(strict=False)


def render_command_preview(state: AppState) -> str:
    commands = build_command_specs(state)
    if not commands:
        return "Choose an analysis approach to preview the exact commands."
    return "\n".join(command.render() for command in commands)


def amf_output_stem(state: AppState) -> str:
    return f"{state.amf_phase_space_base.name}_{state.amf_quantity}"


def planned_output_files(state: AppState) -> list[Path]:
    if state.approach == "amf":
        return _planned_amf_output_files(state)

    output_dir = state.output_dir

    if state.approach is None:
        return []

    if state.approach == "spectrum":
        return [
            output_dir / "particle_matches.csv",
            output_dir / "poly_spectrum.csv",
            output_dir / "poly_spectrum_moments.csv",
        ]

    files: list[Path] = []
    if state.enable_let:
        files.append(output_dir / "let_summary.csv")
    if state.enable_magini:
        files.append(output_dir / "magini_summary.csv")
    if state.enable_inaniwa:
        files.append(output_dir / "inaniwa_summary.csv")
    return files


def _planned_amf_output_files(state: AppState) -> list[Path]:
    run_dir = state.amf_output_dir
    output_stem = amf_output_stem(state)
    files = [
        run_dir / "amf_run_manifest.txt",
        run_dir / "replay_amf.txt",
        run_dir / "topas_stdout.log",
        run_dir / "topas_stderr.log",
    ]
    if state.amf_quantity == "AMFSpectra":
        files.append(run_dir / f"{output_stem}_MicrodosimetricSpectra.csv")
        files.append(run_dir / f"{output_stem}_MicrodosimetricMoments.csv")
    else:
        files.append(run_dir / f"{output_stem}.csv")
    if state.amf_stopping_power == "ExternalTable":
        files.append(run_dir / "StoppingPower.txt")
    return files
