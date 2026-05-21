from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
import shutil

from app.state import AppState


LET_ELEMENTS = ("H", "He", "Li", "Be", "B", "C", "N", "O")
INANIWA_ATOMIC_NUMBERS = tuple(range(1, 11))
PROTON_FILENAME_PATTERN = re.compile(r"^Proton_([0-9]+(?:\.[0-9]+)?)_MeV\.csv$")
CARTECHINI_FILENAME_PATTERN = re.compile(
    r"^H_E([0-9]+(?:\.[0-9]+)?)_R(?:0\.5|8(?:\.0)?)(?:_[^.]+)?\.txt$"
)


@dataclass
class PhaseSpaceScan:
    charged_rows: int = 0
    proton_rows: int = 0
    supported_ion_rows: int = 0
    unsupported_charged_rows: int = 0
    min_energy_mev: float | None = None
    max_energy_mev: float | None = None
    species_counts: dict[str, int] = field(default_factory=dict)
    proton_energies: list[float] = field(default_factory=list)


@dataclass
class InputCheckResult:
    status: str = "warning"
    status_text: str = "Not checked"
    summary_text: str = "No input check results yet."
    species_text: str = ""
    warnings_text: str = ""


@dataclass(frozen=True)
class ParticleIdentity:
    is_charged: bool
    atomic_number: int
    mass_number: int


class InputChecker:
    def __init__(self) -> None:
        self._phase_space_cache: dict[Path, PhaseSpaceScan] = {}
        self._spectrum_range_cache: dict[Path, tuple[float, float] | None] = {}

    def run(self, state: AppState) -> InputCheckResult:
        if state.approach == "amf":
            return self._run_amf_check(state)

        errors: list[str] = []
        warnings: list[str] = []

        if not state.phase_space_file.exists():
            errors.append(f"Phase-space file not found: {state.phase_space_file}")
        if not state.lookup_root.exists():
            errors.append(f"Lookup root not found: {state.lookup_root}")
        if not state.output_dir.parent.exists():
            errors.append(
                f"Output parent directory does not exist: {state.output_dir.parent}"
            )

        if errors:
            return InputCheckResult(
                status="error",
                status_text="Errors found",
                summary_text="Unable to run input check until the required paths exist.",
                warnings_text="\n\n".join(errors),
            )

        scan = self._scan_phase_space(state.phase_space_file)
        species_text = self._format_species_text(scan)

        if scan.charged_rows == 0:
            warnings.append("No charged rows were found in the selected phase-space file.")
        if scan.proton_rows == 0:
            warnings.append("No proton rows were found in the selected phase-space file.")

        if state.approach == "spectrum":
            self._add_spectrum_warnings(state, scan, warnings, errors)
        elif state.approach == "means":
            self._add_means_warnings(state, scan, warnings)
        else:
            warnings.append("Choose an analysis approach to enable approach-specific checks.")

        if errors:
            return InputCheckResult(
                status="error",
                status_text="Errors found",
                summary_text=self._format_summary_text(scan),
                species_text=species_text,
                warnings_text="\n\n".join(errors + warnings),
            )

        return InputCheckResult(
            status="warning" if warnings else "ready",
            status_text="Warnings present" if warnings else "Ready",
            summary_text=self._format_summary_text(scan),
            species_text=species_text,
            warnings_text="\n\n".join(warnings) if warnings else "No blocking issues detected.",
        )

    def _run_amf_check(self, state: AppState) -> InputCheckResult:
        errors: list[str] = []
        warnings: list[str] = []

        if not state.lookup_root.exists():
            errors.append(f"Lookup root not found: {state.lookup_root}")
        if not (state.lookup_root / "AMF" / "tsed.dat").exists():
            errors.append(f"AMF tsed.dat not found: {state.lookup_root / 'AMF' / 'tsed.dat'}")
        if not state.executable_path.exists():
            errors.append(f"MONAS executable not found: {state.executable_path}")
        if not self._is_executable_available(state.topas_executable_path):
            errors.append(f"TOPAS executable not found: {state.topas_executable_path}")
        if not 0.0015 <= state.amf_domain_radius_um <= 0.5:
            errors.append("AMF domain radius must be between 0.0015 um and 0.5 um.")

        if state.amf_stopping_power == "ExternalTable":
            if not state.amf_external_stopping_power_file.exists():
                errors.append(
                    f"External stopping-power table not found: {state.amf_external_stopping_power_file}"
                )
            elif state.amf_external_stopping_power_file.name != "StoppingPower.txt":
                warnings.append(
                    "The selected external stopping-power table will be copied as StoppingPower.txt in the TOPAS run directory."
                )

        scan: PhaseSpaceScan | None = None
        species_text = ""

        if state.amf_run_mode == "replay":
            phase_space_file = state.amf_phase_space_base.with_suffix(".phsp")
            header_file = state.amf_phase_space_base.with_suffix(".header")
            if not phase_space_file.exists():
                errors.append(f"AMF phase-space file not found: {phase_space_file}")
            if not header_file.exists():
                errors.append(f"AMF phase-space header not found: {header_file}")
            if not state.amf_staged_run_dir.parent.exists():
                errors.append(
                    f"AMF staged run parent directory does not exist: {state.amf_staged_run_dir.parent}"
                )

            if phase_space_file.exists():
                scan = self._scan_phase_space(phase_space_file)
                species_text = self._format_species_text(scan)
                if scan.charged_rows == 0:
                    warnings.append("No charged rows were found in the selected AMF phase-space file.")
            detector_warning = self._detector_phase_space_warning(state)
            if detector_warning:
                warnings.append(detector_warning)
            if state.amf_disable_phase_space_precheck:
                warnings.append("TOPAS phase-space precheck will be disabled for replay.")
        else:
            if not state.amf_full_simulation_file.exists():
                errors.append(
                    f"AMF full simulation TOPAS file not found: {state.amf_full_simulation_file}"
                )
            if not state.amf_full_simulation_file.parent.exists():
                errors.append(
                    f"AMF full simulation directory does not exist: {state.amf_full_simulation_file.parent}"
                )
            warnings.append(
                "Full simulation runs the selected TOPAS file directly; the file should already contain the AMF detector and scorer block."
            )

        if errors:
            return InputCheckResult(
                status="error",
                status_text="Errors found",
                summary_text=self._format_amf_summary_text(state, scan),
                species_text=species_text,
                warnings_text="\n\n".join(errors + warnings),
            )

        return InputCheckResult(
            status="warning" if warnings else "ready",
            status_text="Warnings present" if warnings else "Ready",
            summary_text=self._format_amf_summary_text(state, scan),
            species_text=species_text,
            warnings_text="\n\n".join(warnings) if warnings else "No blocking issues detected.",
        )

    def _scan_phase_space(self, phase_space_file: Path) -> PhaseSpaceScan:
        cached = self._phase_space_cache.get(phase_space_file)
        if cached is not None:
            return cached

        scan = PhaseSpaceScan(species_counts={element: 0 for element in LET_ELEMENTS})

        with phase_space_file.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.strip():
                    continue

                fields = line.split()
                if len(fields) < 8:
                    continue

                try:
                    energy_mev = float(fields[5])
                    pdg_code = int(float(fields[7]))
                except ValueError:
                    continue

                identity = decode_particle_identity(pdg_code)
                if not identity.is_charged:
                    continue

                scan.charged_rows += 1
                scan.min_energy_mev = energy_mev if scan.min_energy_mev is None else min(scan.min_energy_mev, energy_mev)
                scan.max_energy_mev = energy_mev if scan.max_energy_mev is None else max(scan.max_energy_mev, energy_mev)

                if pdg_code == 2212:
                    scan.proton_rows += 1
                    scan.proton_energies.append(energy_mev)

                let_family = let_family_for_pdg(pdg_code)
                if let_family in LET_ELEMENTS:
                    scan.species_counts[let_family] += 1
                else:
                    scan.unsupported_charged_rows += 1

                if 1 <= identity.atomic_number <= 10 and identity.mass_number > 0:
                    scan.supported_ion_rows += 1

        self._phase_space_cache[phase_space_file] = scan
        return scan

    def _add_spectrum_warnings(
        self,
        state: AppState,
        scan: PhaseSpaceScan,
        warnings: list[str],
        errors: list[str],
    ) -> None:
        spectrum_dir = selected_spectrum_directory(state)
        if not spectrum_dir.exists():
            errors.append(f"Spectrum LUT directory not found: {spectrum_dir}")
            return

        lut_range = self._read_spectrum_range(spectrum_dir)
        if lut_range is None:
            errors.append(
                f"Could not determine proton energy range from LUT files in: {spectrum_dir}"
            )
            return

        low, high = lut_range
        out_of_range = sum(
            1
            for energy in scan.proton_energies
            if energy < low or energy > high
        )
        warnings.append(
            f"Spectrum LUT range: {low:.6g} to {high:.6g} MeV in {spectrum_dir.name}."
        )
        if out_of_range > 0:
            warnings.append(
                f"{out_of_range} proton rows sit outside the selected spectrum LUT range. If Cartechini selected, fallback to DeCunha will occur. If DeCunha selected, clamping to DeCunha high endpoint will occur."
            )

    def _add_means_warnings(
        self,
        state: AppState,
        scan: PhaseSpaceScan,
        warnings: list[str],
    ) -> None:
        if not any((state.enable_let, state.enable_magini, state.enable_inaniwa)):
            warnings.append("No mean-value calculation is enabled.")

        if state.enable_let:
            missing = [
                element
                for element in LET_ELEMENTS
                if not (state.lookup_root / "LET" / f"{element}_water.txt").exists()
            ]
            if missing:
                warnings.append("Missing LET lookup tables: " + ", ".join(missing))

        if state.enable_inaniwa:
            missing = [
                str(number)
                for number in INANIWA_ATOMIC_NUMBERS
                if not (state.lookup_root / "Inaniwa" / f"Zp_{number}.csv").exists()
            ]
            if missing:
                warnings.append("Missing Inaniwa lookup tables: Zp_" + ", Zp_".join(missing))

        if (state.enable_let or state.enable_inaniwa) and scan.unsupported_charged_rows > 0:
            warnings.append(
                f"{scan.unsupported_charged_rows} charged rows do not map to the currently supported LET or ion mean-value lookup families and would be skipped."
            )

    def _read_spectrum_range(self, spectrum_dir: Path) -> tuple[float, float] | None:
        if spectrum_dir in self._spectrum_range_cache:
            return self._spectrum_range_cache[spectrum_dir]

        energies: list[float] = []
        for path in spectrum_dir.iterdir():
            if not path.is_file():
                continue

            match = PROTON_FILENAME_PATTERN.match(path.name)
            if match is None:
                match = CARTECHINI_FILENAME_PATTERN.match(path.name)
            if match is None:
                continue

            try:
                energies.append(float(match.group(1)))
            except ValueError:
                continue

        if not energies:
            self._spectrum_range_cache[spectrum_dir] = None
            return None

        result = (min(energies), max(energies))
        self._spectrum_range_cache[spectrum_dir] = result
        return result

    @staticmethod
    def _format_summary_text(scan: PhaseSpaceScan) -> str:
        if scan.min_energy_mev is None or scan.max_energy_mev is None:
            energy_text = "n/a"
        else:
            energy_text = f"{scan.min_energy_mev:.6g} - {scan.max_energy_mev:.6g} MeV"

        return (
            f"Charged rows: {scan.charged_rows:,}\n"
            f"Proton rows: {scan.proton_rows:,}\n"
            f"Supported ion rows: {scan.supported_ion_rows:,}\n"
            f"Energy range: {energy_text}"
        )

    @staticmethod
    def _format_species_text(scan: PhaseSpaceScan) -> str:
        parts = [f"{element} {scan.species_counts.get(element, 0):,}" for element in LET_ELEMENTS]
        parts.append(f"Unsupported charged species {scan.unsupported_charged_rows:,}")
        return "   ".join(parts)

    @staticmethod
    def _format_amf_summary_text(
        state: AppState,
        scan: PhaseSpaceScan | None,
    ) -> str:
        lines = [
            f"AMF mode: {'phase-space replay' if state.amf_run_mode == 'replay' else 'full simulation'}",
            f"Quantity: {state.amf_quantity}",
            f"Detector: {state.amf_detector}",
            f"Domain radius: {state.amf_domain_radius_um:g} um",
            f"Stopping power: {state.amf_stopping_power}",
            f"Scoring position: ({state.amf_scoring_x_mm:g}, {state.amf_scoring_y_mm:g}, {state.amf_scoring_z_mm:g}) mm",
        ]
        if scan is not None:
            lines.extend(
                [
                    f"Charged rows: {scan.charged_rows:,}",
                    f"Proton rows: {scan.proton_rows:,}",
                ]
            )
        return "\n".join(lines)

    @staticmethod
    def _is_executable_available(path: Path) -> bool:
        if path.is_absolute() or path.parent != Path("."):
            return path.exists()
        return shutil.which(str(path)) is not None

    @staticmethod
    def _detector_phase_space_warning(state: AppState) -> str:
        token = state.amf_phase_space_base.name.lower()
        if state.amf_detector == "silicon" and not any(
            marker in token for marker in ("soi", "silicon", "si")
        ):
            return (
                "Silicon replay is using a phase-space base whose name does not look like an SOI/silicon phase space. "
                "The manual SOI comparison uses PhaseSpace_seed1_100k_soi; replaying a different phase space through the silicon detector can change the high-y tail."
            )
        if state.amf_detector == "TEgas" and not any(
            marker in token for marker in ("tegas", "te_gas", "tepc", "gas")
        ):
            return (
                "TEgas replay is using a phase-space base whose name does not look like a TE gas or TEPC phase space."
            )
        return ""


def decode_particle_identity(pdg_code: int) -> ParticleIdentity:
    abs_pdg = abs(pdg_code)

    if abs_pdg in (11, 13, 15):
        return ParticleIdentity(True, 0, 0)
    if abs_pdg == 2212:
        return ParticleIdentity(True, 1, 1)
    if abs_pdg in (211, 321):
        return ParticleIdentity(True, 0, 0)
    if abs_pdg >= 1_000_000_000:
        atomic_number = (abs_pdg // 10_000) % 1_000
        mass_number = (abs_pdg // 10) % 1_000
        if atomic_number > 0 and mass_number > 0:
            return ParticleIdentity(True, atomic_number, mass_number)
    return ParticleIdentity(False, 0, 0)


def let_family_for_pdg(pdg_code: int) -> str | None:
    match pdg_code:
        case 2212 | 1000010010 | 1000010020 | 1000010030:
            return "H"
        case 1000020030 | 1000020040:
            return "He"
        case 1000030060 | 1000030070:
            return "Li"
        case 1000040070 | 1000040090:
            return "Be"
        case 1000050100 | 1000050110:
            return "B"
        case 1000060110 | 1000060120 | 1000060130:
            return "C"
        case 1000070130 | 1000070140 | 1000070150:
            return "N"
        case 1000080150 | 1000080160:
            return "O"
        case 1000090190:
            return "F"
        case 1000100200 | 1000100220:
            return "Ne"
        case _:
            return None


def selected_spectrum_directory(state: AppState) -> Path:
    if state.spectrum_family == "Cartechini":
        folder = "R0.5" if state.cartechini_radius == "0.5um" else "R8.0"
        return state.lookup_root / "Cartechini" / folder

    suffix = {
        ("1mm", "linear"): "1mm_linear",
        ("1mm", "log"): "1mm_logarithmic",
        ("5um", "linear"): "5um_linear",
        ("5um", "log"): "5um_logarithmic",
    }[(state.decunha_voxel_size, state.decunha_energy_grid)]
    return state.lookup_root / "DeCunha" / suffix
