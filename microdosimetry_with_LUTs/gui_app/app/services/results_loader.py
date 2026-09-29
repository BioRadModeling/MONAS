from __future__ import annotations

from dataclasses import dataclass, field
import csv
from pathlib import Path
import re

from app.state import AppState
from app.services.command_builder import amf_output_stem


@dataclass
class SpectrumPoint:
    y_keV_per_um: float
    f_y: float
    yf_y: float
    d_y: float
    yd_y: float


@dataclass
class LoadedResults:
    spectrum_summary_text: str = "Spectrum results have not been loaded yet."
    spectrum_points: list[SpectrumPoint] = field(default_factory=list)
    spectrum_file: Path | None = None
    preferred_file: Path | None = None
    preferred_metric_files: dict[str, Path] = field(default_factory=dict)
    let_headers: list[str] = field(default_factory=list)
    let_rows: list[list[str]] = field(default_factory=list)
    magini_pairs: list[tuple[str, str]] = field(default_factory=list)
    inaniwa_pairs: list[tuple[str, str]] = field(default_factory=list)
    at_pairs: list[tuple[str, str]] = field(default_factory=list)
    amf_pairs: list[tuple[str, str]] = field(default_factory=list)
    file_paths: list[Path] = field(default_factory=list)
    file_previews: dict[str, str] = field(default_factory=dict)


class ResultsLoader:
    def load(self, state: AppState) -> LoadedResults:
        output_dir = self._effective_output_dir(state)
        loaded = LoadedResults()

        if not output_dir.exists():
            return loaded

        phase_space_token = self._phase_space_token(state.phase_space_file)
        loaded.file_paths = sorted(
            [path for path in output_dir.iterdir() if path.is_file()],
            key=lambda path: self._file_sort_key(path, state, phase_space_token),
        )
        loaded.file_previews = {
            str(path): self._preview_file(path) for path in loaded.file_paths
        }
        loaded.preferred_metric_files = self._preferred_metric_files(loaded.file_paths, state, phase_space_token)

        if state.approach == "amf":
            self._load_amf_results(loaded, output_dir, state)
            return loaded

        spectrum_path = self._resolve_output_file(output_dir, "poly_spectrum", state, phase_space_token)
        let_path = self._resolve_output_file(output_dir, "let_summary", state, phase_space_token)
        magini_path = self._resolve_output_file(output_dir, "magini_summary", state, phase_space_token)
        inaniwa_path = self._resolve_output_file(output_dir, "inaniwa_summary", state, phase_space_token)
        at_path = self._resolve_output_file(output_dir, "at_summary", state, phase_space_token)

        loaded.spectrum_file = spectrum_path if spectrum_path.exists() else None
        loaded.spectrum_summary_text, loaded.spectrum_points = self._load_spectrum(spectrum_path)
        loaded.let_headers, loaded.let_rows = self._load_csv_table(let_path)
        loaded.magini_pairs = self._load_key_value_summary(magini_path)
        loaded.inaniwa_pairs = self._load_key_value_summary(inaniwa_path)
        loaded.at_pairs = self._load_key_value_summary(at_path)
        loaded.preferred_file = self._preferred_default_file(
            loaded.file_paths,
            loaded.preferred_metric_files,
            spectrum_path if spectrum_path.exists() else None,
            let_path if let_path.exists() else None,
            magini_path if magini_path.exists() else None,
            inaniwa_path if inaniwa_path.exists() else None,
            at_path if at_path.exists() else None,
            state,
        )
        return loaded

    def _load_amf_results(
        self,
        loaded: LoadedResults,
        output_dir: Path,
        state: AppState,
    ) -> None:
        output_stem = amf_output_stem(state)
        allow_name_fallback = False

        if state.amf_quantity == "AMFSpectra":
            spectrum_path = self._resolve_amf_spectrum_file(
                output_dir,
                output_stem,
                allow_name_fallback,
            )
            moments_path = self._resolve_amf_spectrum_moments_file(
                output_dir,
                output_stem,
                allow_name_fallback,
            )
            loaded.spectrum_file = spectrum_path if spectrum_path.exists() else None
            loaded.spectrum_summary_text, loaded.spectrum_points = self._load_amf_spectrum(spectrum_path)
            loaded.preferred_file = spectrum_path if spectrum_path.exists() else None
            loaded.amf_pairs = [
                ("Quantity", state.amf_quantity),
                ("Domain radius [um]", f"{state.amf_domain_radius_um:g}"),
                ("Spectrum file", str(spectrum_path) if spectrum_path.exists() else "Not found"),
                ("Moments file", str(moments_path) if moments_path.exists() else "Not found"),
            ]
            return

        scalar_path = self._resolve_amf_scalar_file(
            output_dir,
            output_stem,
            state.amf_quantity,
            allow_name_fallback,
        )
        value = self._load_amf_scalar(scalar_path)
        loaded.preferred_file = scalar_path if scalar_path.exists() else None
        loaded.amf_pairs = [
            ("Quantity", state.amf_quantity),
            ("Domain radius [um]", f"{state.amf_domain_radius_um:g}"),
            ("Result file", str(scalar_path) if scalar_path.exists() else "Not found"),
            ("Value", value if value is not None else "No numeric result found"),
        ]

    def _effective_output_dir(self, state: AppState) -> Path:
        if state.approach != "amf":
            return state.output_dir
        return state.amf_output_dir

    def _resolve_amf_spectrum_file(
        self,
        output_dir: Path,
        output_stem: str,
        allow_name_fallback: bool,
    ) -> Path:
        exact = output_dir / f"{output_stem}_MicrodosimetricSpectra.csv"
        if exact.exists() or not allow_name_fallback:
            return exact

        named = output_dir / "AMF_Spectra_MicrodosimetricSpectra.csv"
        if named.exists():
            return named

        matches = sorted(output_dir.glob("*MicrodosimetricSpectra.csv"))
        if matches:
            return matches[0]
        return exact

    def _resolve_amf_spectrum_moments_file(
        self,
        output_dir: Path,
        output_stem: str,
        allow_name_fallback: bool,
    ) -> Path:
        exact = output_dir / f"{output_stem}_MicrodosimetricMoments.csv"
        if exact.exists() or not allow_name_fallback:
            return exact

        named = output_dir / "AMF_Spectra_MicrodosimetricMoments.csv"
        if named.exists():
            return named

        matches = sorted(output_dir.glob("*MicrodosimetricMoments.csv"))
        if matches:
            return matches[0]
        return exact

    def _resolve_amf_scalar_file(
        self,
        output_dir: Path,
        output_stem: str,
        quantity: str,
        allow_name_fallback: bool,
    ) -> Path:
        exact = output_dir / f"{output_stem}.csv"
        if exact.exists() or not allow_name_fallback:
            return exact

        named = output_dir / f"{quantity}.csv"
        if named.exists():
            return named

        matches = sorted(output_dir.glob(f"*{quantity}*.csv"))
        if matches:
            return matches[0]
        return exact

    def _resolve_output_file(
        self,
        output_dir: Path,
        stem: str,
        state: AppState,
        phase_space_token: str,
    ) -> Path:
        exact = output_dir / f"{stem}.csv"
        if exact.exists():
            return exact

        matches = sorted(
            output_dir.glob(f"{stem}*.csv"),
            key=lambda path: self._candidate_score(path, state, phase_space_token, stem),
        )
        if matches:
            return matches[0]
        return exact

    def _load_spectrum(self, path: Path) -> tuple[str, list[SpectrumPoint]]:
        if not path.exists():
            return ("No `poly_spectrum.csv` found in the selected output directory.", [])

        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            rows = list(csv.DictReader(handle))

        if not rows:
            return ("`poly_spectrum.csv` exists but contains no spectrum rows.", [])

        points: list[SpectrumPoint] = []
        for row in rows:
            try:
                points.append(
                    SpectrumPoint(
                        y_keV_per_um=float(row["y_keV_per_um"]),
                        f_y=float(row["f_y"]),
                        yf_y=float(row["yf_y"]),
                        d_y=float(row["d_y"]),
                        yd_y=float(row["yd_y"]),
                    )
                )
            except (KeyError, ValueError):
                continue

        first = rows[0]
        last = rows[-1]
        summary = (
            "Spectrum summary\n\n"
            f"Bins: {len(rows):,}\n"
            f"y range: {first.get('y_keV_per_um', 'n/a')} to {last.get('y_keV_per_um', 'n/a')} keV/um\n"
            f"First yf(y): {first.get('yf_y', 'n/a')}\n"
            f"First yd(y): {first.get('yd_y', 'n/a')}\n"
        )
        return summary, points

    def _load_amf_spectrum(self, path: Path) -> tuple[str, list[SpectrumPoint]]:
        if not path.exists():
            return ("No AMF microdosimetric spectra CSV found in the selected output directory.", [])

        rows = self._read_csv_rows(path)
        if len(rows) < 2:
            scorer_total = self._load_amf_scalar(path.with_name(path.name.replace("_MicrodosimetricSpectra", "")))
            total_text = f" Scorer total: {scorer_total}." if scorer_total is not None else ""
            return (
                "AMF spectra CSV exists but contains no spectra data row."
                f"{total_text} Check that the replay detector overlaps the phase-space plane.",
                [],
            )

        y_headers = rows[0][3:]
        value_row = rows[1][3:]
        points: list[SpectrumPoint] = []
        for y_text, value_text in zip(y_headers, value_row):
            try:
                y_value = float(y_text)
                yd_value = float(value_text)
            except ValueError:
                continue
            points.append(
                SpectrumPoint(
                    y_keV_per_um=y_value,
                    f_y=0.0,
                    yf_y=0.0,
                    d_y=0.0,
                    yd_y=yd_value,
                )
            )

        if not points:
            return ("AMF spectra CSV was found, but no numeric y bins could be loaded.", [])

        summary = (
            "AMF spectrum summary\n\n"
            f"Bins: {len(points):,}\n"
            f"y range: {points[0].y_keV_per_um:.6g} to {points[-1].y_keV_per_um:.6g} keV/um\n"
            f"First yd(y): {points[0].yd_y:.6g}\n"
        )
        return summary, points

    def _load_amf_scalar(self, path: Path) -> str | None:
        if not path.exists():
            return None

        try:
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    parts = [part.strip() for part in re.split(r"[,\s]+", stripped) if part.strip()]
                    return parts[-1] if parts else stripped
        except OSError:
            return None
        return None

    def _load_csv_table(self, path: Path) -> tuple[list[str], list[list[str]]]:
        if not path.exists():
            return [], []

        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle)
            rows = list(reader)

        if not rows:
            return [], []

        headers = rows[0]
        body = rows[1:]
        return headers, body

    def _read_csv_rows(self, path: Path) -> list[list[str]]:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            return [row for row in csv.reader(handle) if row and not row[0].startswith("#")]

    def _load_key_value_summary(self, path: Path) -> list[tuple[str, str]]:
        if not path.exists():
            return []

        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            rows = list(csv.DictReader(handle))

        if not rows:
            return []

        row = rows[0]
        return [(key, value) for key, value in row.items()]

    def _preview_file(self, path: Path, max_lines: int = 40) -> str:
        if path.suffix.lower() in {".jpg", ".jpeg", ".png"}:
            return f"Binary image file:\n{path.name}"

        lines: list[str] = []
        try:
            with path.open("r", encoding="utf-8", errors="replace") as handle:
                for index, line in enumerate(handle):
                    if index >= max_lines:
                        lines.append("...")
                        break
                    lines.append(line.rstrip("\n"))
        except OSError as error:
            return f"Unable to read {path.name}: {error}"

        return "\n".join(lines) if lines else f"{path.name} is empty."

    def _phase_space_token(self, phase_space_file: Path) -> str:
        stem = phase_space_file.stem
        match = re.search(r"(\d+(?:\.\d+)?mm)", stem, flags=re.IGNORECASE)
        if match:
            return match.group(1).lower()

        normalized = stem.lower().replace("phasespace_", "").replace("phase_space_", "")
        return normalized

    def _candidate_score(
        self,
        path: Path,
        state: AppState,
        phase_space_token: str,
        stem: str,
    ) -> tuple[int, int, int, str]:
        name = path.name.lower()
        score = 0

        if path.name == f"{stem}.csv":
            score -= 100

        if phase_space_token and phase_space_token in name:
            score -= 40

        if state.approach == "spectrum":
            family = state.spectrum_family.lower()
            if family in name:
                score -= 25

            if family == "cartechini":
                particle_token = state.cartechini_particle.lower()
                if particle_token in name:
                    score -= 10
                radius_token = state.cartechini_radius.replace("um", "").lower()
                if radius_token in name:
                    score -= 10

        return (
            score,
            len(name),
            0 if name.endswith(".csv") else 1,
            name,
        )

    def _file_sort_key(
        self,
        path: Path,
        state: AppState,
        phase_space_token: str,
    ) -> tuple[int, int, int, str]:
        name = path.name.lower()
        group_score = 3

        if state.approach == "spectrum" and (
            name.startswith("poly_spectrum")
            or name.startswith("poly_spectrum_moments")
            or name.startswith("particle_matches")
            or name.startswith("proton_matches")
            or "_vs_y_" in name
        ):
            group_score = 0
        elif state.approach == "amf" and (
            name.startswith("amf_")
            or name.startswith("replay_amf")
            or name.startswith("topas_")
            or "amfspectra" in name
            or "amf_yd" in name
            or "amf_ys" in name
            or name == "tsed.dat"
        ):
            group_score = 0
        elif state.approach == "means" and (
            name.startswith("let_summary")
            or name.startswith("magini_summary")
            or name.startswith("inaniwa_summary")
            or name.startswith("at_summary")
            or name.startswith("at_diagnostics")
        ):
            group_score = 0
        elif name.endswith(".csv"):
            group_score = 1
        elif name.endswith((".jpg", ".jpeg", ".png")):
            group_score = 2

        candidate_score = self._candidate_score(
            path,
            state,
            phase_space_token,
            self._logical_stem(path.name),
        )
        return (group_score, *candidate_score)

    def _logical_stem(self, filename: str) -> str:
        lower = filename.lower()
        for stem in (
            "poly_spectrum_moments",
            "poly_spectrum",
            "particle_matches",
            "proton_matches",
            "let_summary",
            "magini_summary",
            "inaniwa_summary",
            "at_summary",
            "at_diagnostics",
        ):
            if lower.startswith(stem):
                return stem
        return Path(filename).stem

    def _preferred_metric_files(
        self,
        file_paths: list[Path],
        state: AppState,
        phase_space_token: str,
    ) -> dict[str, Path]:
        if state.approach != "spectrum":
            return {}

        metric_prefixes = {
            "yd_y": "yd_y_vs_y",
            "yf_y": "yf_y_vs_y",
        }
        preferred: dict[str, Path] = {}
        for metric_key, prefix in metric_prefixes.items():
            matches = [
                path
                for path in file_paths
                if path.suffix.lower() in {".jpg", ".jpeg", ".png"} and path.name.lower().startswith(prefix)
            ]
            if not matches:
                continue
            preferred[metric_key] = sorted(
                matches,
                key=lambda path: self._candidate_score(path, state, phase_space_token, prefix),
            )[0]
        return preferred

    def _preferred_default_file(
        self,
        file_paths: list[Path],
        preferred_metric_files: dict[str, Path],
        spectrum_path: Path | None,
        let_path: Path | None,
        magini_path: Path | None,
        inaniwa_path: Path | None,
        at_path: Path | None,
        state: AppState,
    ) -> Path | None:
        if not file_paths:
            return None

        if state.approach == "spectrum":
            return preferred_metric_files.get("yd_y") or spectrum_path or file_paths[0]

        if state.approach == "means":
            for path in (let_path, magini_path, inaniwa_path, at_path):
                if path is not None:
                    return path

        return file_paths[0]
