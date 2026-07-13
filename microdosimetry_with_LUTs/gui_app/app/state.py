from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path


@dataclass
class AppState:
    phase_space_file: Path
    lookup_root: Path
    output_dir: Path
    executable_path: Path
    build_workdir: Path
    approach: str | None = None
    spectrum_family: str = "DeCunha"
    decunha_voxel_size: str = "1mm"
    decunha_energy_grid: str = "log"
    cartechini_particle: str = "proton"
    cartechini_radius: str = "0.5um"
    rebin_samples: int = 1_000_000
    rebin_seed: int = 1_511_504_868
    enable_let: bool = True
    enable_magini: bool = True
    enable_inaniwa: bool = True
    enable_at: bool = True
    at_particle: str = "proton"
    amf_quantity: str = "AMFSpectra"
    amf_domain_radius_um: float = 0.28
    amf_stopping_power: str = "Topas"
    amf_phase_space_base: Path = Path()
    amf_source_topas_file: Path = Path()
    amf_output_dir: Path = Path()
    amf_external_stopping_power_file: Path = Path()
    amf_disable_phase_space_precheck: bool = True
    amf_scoring_radius_mm: float = 6.35

    def clone(self) -> "AppState":
        return replace(self)


def default_app_state(project_root: Path) -> AppState:
    return AppState(
        phase_space_file=project_root / "input" / "PhaseSpace_33mm.phsp",
        lookup_root=project_root / "lookup_tables",
        output_dir=project_root / "output" / "gui_run_33mm",
        executable_path=project_root / "build" / "microdosimetry_with_LUTs",
        build_workdir=project_root / "build",
        amf_phase_space_base=project_root / "input" / "PhaseSpace_curved_33mm",
        amf_source_topas_file=(
            project_root
            / "amf_runtime"
            / "comparison_seed1_100k"
            / "01_generate_phase_space.txt"
        ),
        amf_output_dir=project_root / "amf_runtime" / "staged_runs" / "gui_amf_run",
        amf_external_stopping_power_file=project_root / "amf_runtime" / "staged_runs" / "gui_amf_run" / "StoppingPower.txt",
    )
