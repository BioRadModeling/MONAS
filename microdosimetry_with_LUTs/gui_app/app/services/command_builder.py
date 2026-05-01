from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from shlex import quote

from app.state import AppState


@dataclass(frozen=True)
class CommandSpec:
    label: str
    argv: list[str]

    def render(self) -> str:
        return " ".join(quote(part) for part in self.argv)


def build_command_specs(state: AppState) -> list[CommandSpec]:
    executable = str(state.executable_path)
    lookup_root = str(state.lookup_root)
    phase_space_file = str(state.phase_space_file)
    output_dir = str(state.output_dir)

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


def render_command_preview(state: AppState) -> str:
    commands = build_command_specs(state)
    if not commands:
        return "Choose an analysis approach to preview the exact commands."
    return "\n".join(command.render() for command in commands)


def planned_output_files(state: AppState) -> list[Path]:
    output_dir = state.output_dir

    if state.approach is None:
        return []

    if state.approach == "spectrum":
        return [
            output_dir / "proton_matches.csv",
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
