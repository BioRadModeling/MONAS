import csv
import math
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "build" / "microdosimetry_with_LUTs"
CSV_WRITER_CPP = ROOT / "src" / "CsvWriter.cpp"
AMF_RESULT_PARSER_CPP = ROOT / "src" / "AmfResultParser.cpp"
AMF_SPECTRA_SCORER_CC = ROOT / "amf_extension" / "ScoreAMFSpectra.cc"

AMF_SPECTRA = (
    ROOT
    / "amf_runtime"
    / "staged_runs"
    / "PhaseSpace_curved_33mm_AMFSpectra"
    / "PhaseSpace_curved_33mm_AMFSpectra_MicrodosimetricSpectra.csv"
)
AMF_YD = (
    ROOT
    / "amf_runtime"
    / "staged_runs"
    / "PhaseSpace_curved_33mm_AMF_yD"
    / "PhaseSpace_curved_33mm_AMF_yD.csv"
)
LUT_SPECTRUM = ROOT / "output" / "curved_surface_33mm" / "poly_spectrum.csv"
LUT_MOMENTS = ROOT / "output" / "curved_surface_33mm" / "poly_spectrum_moments.csv"


def read_csv_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return [row for row in csv.reader(handle) if row and not row[0].startswith("#")]


def read_amf_direct_yd(path=AMF_YD):
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            return float(line)
    raise AssertionError(f"No numeric AMF_yD value found in {path}")


def amf_y_grid():
    edges = []
    ypower = -3.0
    for _ in range(401):
        edges.append(10.0**ypower)
        ypower += 0.02
    return [0.5 * (edges[i] + edges[i + 1]) for i in range(400)]


def logarithmic_probability_factor(y_values):
    return math.log(10.0) * (math.log10(y_values[1]) - math.log10(y_values[0]))


def read_amf_spectrum_values(path=AMF_SPECTRA):
    rows = read_csv_rows(path)
    if len(rows) < 2:
        raise AssertionError(f"No AMF spectra data row found in {path}")
    return [float(value) for value in rows[1][3:]]


def integrated_yd_mean(y_values, yd_values, low=0.0, high=math.inf):
    factor = logarithmic_probability_factor(y_values)
    numerator = 0.0
    probability = 0.0
    for y, yd in zip(y_values, yd_values):
        if low <= y <= high:
            probability += factor * yd
            numerator += factor * y * yd
    return numerator, probability


def ratio_standard_error(samples):
    numerator_sum = sum(n for n, _ in samples)
    denominator_sum = sum(d for _, d in samples)
    ratio = numerator_sum / denominator_sum
    residual_sum = sum((n - ratio * d) ** 2 for n, d in samples)
    return math.sqrt(
        len(samples)
        / (len(samples) - 1)
        * residual_sum
        / denominator_sum**2
    )


class AmfDiagnosticsTest(unittest.TestCase):
    def test_ratio_standard_error_formula_matches_manual_calculation(self):
        samples = [(10.0, 2.0), (12.0, 3.0), (28.0, 4.0)]
        numerator_sum = sum(n for n, _ in samples)
        denominator_sum = sum(d for _, d in samples)
        ratio = numerator_sum / denominator_sum

        self.assertAlmostEqual(ratio, 50.0 / 9.0)
        self.assertAlmostEqual(
            ratio_standard_error(samples),
            math.sqrt(
                3.0
                / 2.0
                * sum((n - ratio * d) ** 2 for n, d in samples)
                / denominator_sum**2
            ),
        )

    def test_poly_spectrum_moments_schema_includes_standard_error(self):
        writer_source = CSV_WRITER_CPP.read_text(encoding="utf-8")
        self.assertIn("mean_standard_error_keV_per_um", writer_source)

    def test_amf_spectra_moments_output_is_wired(self):
        parser_source = AMF_RESULT_PARSER_CPP.read_text(encoding="utf-8")
        scorer_source = AMF_SPECTRA_SCORER_CC.read_text(encoding="utf-8")

        self.assertIn("_MicrodosimetricMoments.csv", parser_source)
        self.assertIn("_MicrodosimetricMoments.csv", scorer_source)
        self.assertIn("mean_standard_error_keV_per_um", scorer_source)

    def test_amf_status_printer_reports_moments_csv(self):
        main_source = (ROOT / "src" / "main.cpp").read_text(encoding="utf-8")

        self.assertIn("Expected spectra moments CSV:", main_source)

    def test_amf_stage_accepts_rotated_curved_surface_without_component(self):
        if not BINARY.exists():
            self.skipTest(f"Missing built executable: {BINARY}")
        if not (ROOT / "lookup_tables" / "AMF" / "tsed.dat").exists():
            self.skipTest("Missing AMF tsed.dat")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            source = temp / "source.txt"
            phase_space_base = temp / "TEgasSV_EntranceHemisphere_PhaseSpace_15cm"
            run_dir = temp / "staged"

            source.write_text(
                """
s:Ge/World/Material = "Air"
d:Ge/World/HLX = 5.0 m
d:Ge/World/HLY = 5.0 m
d:Ge/World/HLZ = 5.0 m
s:Ge/Phantom/Type = "TsBox"
s:Ge/Phantom/Parent = "World"
s:Ge/Phantom/Material = "G4_WATER"
d:Ge/Phantom/HLX = 7.0 cm
d:Ge/Phantom/HLY = 7.0 cm
d:Ge/Phantom/HLZ = 15.0 cm
d:Ge/Phantom/TransX = 0. cm
d:Ge/Phantom/TransY = 0. cm
d:Ge/Phantom/TransZ = -16. cm
d:Ge/Phantom/RotX = 180. deg
d:Ge/Phantom/RotY = 0. deg
d:Ge/Phantom/RotZ = 0. deg
s:Ge/AlShellOut/Parent = "Phantom"
s:Ge/AlShellOut/Type = "TsSphere"
s:Ge/AlShellOut/Material = "Aluminum"
d:Ge/AlShellOut/RMin = 0. mm
d:Ge/AlShellOut/RMax = 10.0 mm
d:Ge/AlShellOut/TransX = 0. mm
d:Ge/AlShellOut/TransY = 0. mm
d:Ge/AlShellOut/TransZ = -115.0 mm
s:Ge/AlShellIn/Parent = "AlShellOut"
s:Ge/AlShellIn/Type = "TsSphere"
s:Ge/AlShellIn/Material = "TEgas"
d:Ge/AlShellIn/RMin = 0. mm
d:Ge/AlShellIn/RMax = 9.82 mm
s:Ge/A150/Parent = "AlShellIn"
s:Ge/A150/Type = "TsSphere"
s:Ge/A150/Material = "G4_A-150_TISSUE"
d:Ge/A150/RMin = 0. mm
d:Ge/A150/RMax = 7.62 mm
s:Ge/TEgasSV/Parent = "A150"
s:Ge/TEgasSV/Type = "TsSphere"
s:Ge/TEgasSV/Material = "TEgas"
d:Ge/TEgasSV/RMin = 0. mm
d:Ge/TEgasSV/RMax = 6.35 mm
s:Sc/TEgasSVEntrancePS/Quantity = "PhaseSpace"
s:Sc/TEgasSVEntrancePS/Surface = "TEgasSV/OuterCurvedSurface"
s:Sc/TEgasSVEntrancePS/OutputType = "ASCII"
s:Sc/TEgasSVEntrancePS/OutputFile = "TEgasSV_EntranceHemisphere_PhaseSpace"
""".lstrip(),
                encoding="utf-8",
            )
            (phase_space_base.with_suffix(".header")).write_text(
                "TOPAS ASCII Phase Space\n", encoding="utf-8"
            )
            (phase_space_base.with_suffix(".phsp")).write_text(
                "\n".join(
                    [
                        "0.0 0.0 -4.49999 0 0 1 1 2212 0 1",
                        "0.0 0.0 -3.86500 0 0 1 1 2212 0 1",
                        "0.635 0.0 -4.2 0 0 1 1 2212 0 1",
                        "-0.635 0.0 -4.2 0 0 1 1 2212 0 1",
                        "0.0 0.635 -4.2 0 0 1 1 2212 0 1",
                        "0.0 -0.635 -4.2 0 0 1 1 2212 0 1",
                        "0.1 0.1 -4.2 0 0 0 1 11 0 0",
                    ]
                ),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    str(BINARY),
                    "AMF-stage",
                    str(ROOT / "lookup_tables"),
                    str(phase_space_base),
                    str(source),
                    str(run_dir),
                    "AMFSpectra",
                    "--scoring-radius-mm",
                    "6.35",
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(
                result.returncode,
                0,
                msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
            )
            manifest = (run_dir / "amf_run_manifest.txt").read_text(
                encoding="utf-8"
            )
            replay = (run_dir / "replay_amf.txt").read_text(encoding="utf-8")
            self.assertIn("phase_space_component = TEgasSV", manifest)
            self.assertIn("phase_space_energy_floored_rows = 1", manifest)
            self.assertIn("scoring_trans_z_mm = -45", manifest)
            self.assertIn("d:Ge/Phantom/TransZ = -160", replay)
            self.assertIn("d:Ge/AMFScoringVolume/TransZ = 115", replay)
            staged_phase_space = (
                run_dir / "TEgasSV_EntranceHemisphere_PhaseSpace_15cm.phsp"
            ).read_text(encoding="utf-8")
            self.assertIn("0.1 0.1 -4.2 0 0 1e-09 1 11 0 0", staged_phase_space)

    def test_lut_poly_spectrum_recomputes_reported_dose_mean(self):
        if not LUT_SPECTRUM.exists() or not LUT_MOMENTS.exists():
            self.skipTest("Missing LUT spectrum fixtures")

        rows = read_csv_rows(LUT_SPECTRUM)
        header = rows[0]
        y_index = header.index("y_keV_per_um")
        d_index = header.index("d_y")
        y_values = [float(row[y_index]) for row in rows[1:]]
        d_values = [float(row[d_index]) for row in rows[1:]]

        factor = logarithmic_probability_factor(y_values)
        probability = sum(factor * y * d for y, d in zip(y_values, d_values))
        mean = (
            sum(factor * y * d * y for y, d in zip(y_values, d_values))
            / probability
        )

        moments = {
            row[0]: float(row[1])
            for row in read_csv_rows(LUT_MOMENTS)[1:]
        }
        self.assertLess(abs(probability - 1.0), 5.0e-5)
        self.assertAlmostEqual(mean, moments["dose"], places=4)

    def test_amf_direct_yd_matches_amf_spectrum_integral(self):
        if not AMF_SPECTRA.exists() or not AMF_YD.exists():
            self.skipTest("Missing AMF spectra/direct yD fixtures")

        y_values = amf_y_grid()
        yd_values = read_amf_spectrum_values()

        mean, probability = integrated_yd_mean(y_values, yd_values)
        direct_yd = read_amf_direct_yd()

        self.assertAlmostEqual(probability, 1.0, places=4)
        self.assertLess(abs(mean - direct_yd) / direct_yd, 0.01)

    def test_amf_lut_difference_is_dominated_by_high_y_tail(self):
        if not AMF_SPECTRA.exists() or not LUT_MOMENTS.exists():
            self.skipTest("Missing AMF or LUT moments fixtures")

        y_values = amf_y_grid()
        yd_values = read_amf_spectrum_values()
        full_mean, _ = integrated_yd_mean(y_values, yd_values)
        below_100_mean, below_100_probability = integrated_yd_mean(
            y_values, yd_values, high=100.0
        )
        above_100_mean, above_100_probability = integrated_yd_mean(
            y_values, yd_values, low=100.0
        )
        lut_dose_mean = {
            row[0]: float(row[1])
            for row in read_csv_rows(LUT_MOMENTS)[1:]
        }["dose"]

        self.assertGreater(full_mean / lut_dose_mean, 4.0)
        self.assertLess((below_100_mean / below_100_probability) / lut_dose_mean, 1.5)
        self.assertLess(above_100_probability, 0.02)
        self.assertGreater(above_100_mean / full_mean, 0.70)

    def test_saved_amf_spectra_header_contains_nonzero_y_centers(self):
        if not AMF_SPECTRA.exists():
            self.skipTest(f"Missing fixture: {AMF_SPECTRA}")

        rows = read_csv_rows(AMF_SPECTRA)
        y_headers = [float(value) for value in rows[0][3:]]
        self.assertEqual(len(y_headers), 400)
        self.assertGreater(sum(1 for y in y_headers if y > 0.0), 390)


if __name__ == "__main__":
    unittest.main()
