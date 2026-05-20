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
    cartechini_radius: str = "0.5um"
    rebin_samples: int = 1_000_000
    rebin_seed: int = 1_511_504_868
    enable_let: bool = True
    enable_magini: bool = True
    enable_inaniwa: bool = True
    topas_executable_path: Path = Path("topas")
    amf_run_mode: str = "replay"
    amf_quantity: str = "AMFSpectra"
    amf_detector: str = "water"
    amf_domain_radius_um: float = 0.28
    amf_stopping_power: str = "Topas"
    amf_step_calculator: str = "MidStep"
    amf_phase_space_base: Path = Path()
    amf_staged_run_dir: Path = Path()
    amf_full_simulation_file: Path = Path()
    amf_external_stopping_power_file: Path = Path()
    amf_disable_phase_space_precheck: bool = True
    amf_world_half_length_cm: float = 148.79
    amf_scoring_x_mm: float = 0.0
    amf_scoring_y_mm: float = 0.0
    amf_scoring_z_mm: float = -1387.4

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
        amf_staged_run_dir=project_root / "amf_runtime" / "staged_runs" / "gui_amf_run",
        amf_full_simulation_file=project_root / "amf_runtime" / "full_simulation_amf.txt",
        amf_external_stopping_power_file=project_root / "amf_runtime" / "staged_runs" / "gui_amf_run" / "StoppingPower.txt",
    )
