from __future__ import annotations

from datetime import datetime
from pathlib import Path
import sys

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QHBoxLayout,
    QMainWindow,
    QMessageBox,
    QScrollArea,
    QSplitter,
    QWidget,
)

from app.services.command_builder import build_command_specs, planned_output_files
from app.services.input_checker import InputChecker
from app.services.process_runner import ProcessRunner
from app.services.results_loader import LoadedResults, ResultsLoader
from app.state import AppState, default_app_state
from app.widgets.results_panel import OverviewData, ResultsPanel, RunHistoryEntry
from app.widgets.setup_panel import SetupPanel


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.project_root = Path(__file__).resolve().parents[2]
        self.state = default_app_state(self.project_root)
        self.runner = ProcessRunner(self)
        self.input_checker = InputChecker()
        self.results_loader = ResultsLoader()
        self.run_history: list[RunHistoryEntry] = []
        self._loaded_results = LoadedResults()
        self._run_status = "Ready to configure a run."

        self.setWindowTitle("MONAS Fast Microdosimetry GUI")
        self.resize(1400, 920)

        central = QWidget()
        layout = QHBoxLayout(central)
        layout.setContentsMargins(12, 12, 12, 12)

        splitter = QSplitter()
        self.setup_panel = SetupPanel(self.state)
        self.setup_scroll = QScrollArea()
        self.setup_scroll.setWidgetResizable(True)
        self.setup_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        self.setup_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setup_scroll.setWidget(self.setup_panel)
        self.results_panel = ResultsPanel()
        splitter.addWidget(self.setup_scroll)
        splitter.addWidget(self.results_panel)
        splitter.setSizes([640, 760])
        layout.addWidget(splitter)
        self.setCentralWidget(central)

        self._apply_stylesheet()
        self._wire_signals()
        self._refresh_input_check()
        self._refresh_loaded_results()

    def _wire_signals(self) -> None:
        self.setup_panel.stateChanged.connect(self._handle_state_changed)
        self.setup_panel.runRequested.connect(self._run_analysis)
        self.setup_panel.resetRequested.connect(self._handle_reset)
        self.setup_panel.inputCheckRequested.connect(self._refresh_input_check)

        self.runner.commandStarted.connect(self._handle_command_started)
        self.runner.outputReady.connect(self.results_panel.append_log)
        self.runner.runFinished.connect(self._handle_run_finished)

    def _handle_state_changed(self, state: AppState) -> None:
        self.state = state
        if not self.runner.is_running:
            self._run_status = "Selections updated. Review the command preview or start a run."
        self._refresh_input_check()
        self._refresh_loaded_results()

    def _handle_reset(self) -> None:
        self.state = self.setup_panel.current_state()
        self.results_panel.reset_to_overview()
        self.results_panel.clear_log()
        self._run_status = "Configuration reset to the default sample paths."
        self._refresh_input_check()
        self._refresh_loaded_results()

    def _sync_overview(self, state: AppState) -> None:
        outputs = planned_output_files(state)
        self.results_panel.set_planned_outputs(outputs)

        detail_lines: list[str] = []
        approach_line = "No analysis approach selected yet."
        phase_space_display = str(state.phase_space_file)
        output_display = str(state.output_dir)
        if state.approach == "spectrum":
            approach_line = f"Full microdosimetric spectrum via {state.spectrum_family}."
            if state.spectrum_family == "DeCunha":
                detail_lines.append(
                    f"DeCunha settings: voxel={state.decunha_voxel_size}, grid={state.decunha_energy_grid}"
                )
            else:
                detail_lines.append(
                    f"Cartechini settings: radius={state.cartechini_radius}"
                )
        elif state.approach == "means":
            approach_line = "Mean values only."
            enabled = []
            if state.enable_let:
                enabled.append("LET")
            if state.enable_magini:
                enabled.append("Magini")
            if state.enable_inaniwa:
                enabled.append("Inaniwa")
            detail_lines.append(
                "Enabled mean-value modes: " + (", ".join(enabled) if enabled else "none")
            )
        elif state.approach == "amf":
            if state.amf_run_mode == "replay":
                approach_line = "AMF phase-space replay."
                phase_space_display = str(state.amf_phase_space_base)
                output_display = str(state.amf_staged_run_dir)
            else:
                approach_line = "AMF full simulation."
                phase_space_display = str(state.amf_full_simulation_file)
                output_display = str(state.amf_full_simulation_file.parent)
            detail_lines.extend(
                [
                    f"Quantity={state.amf_quantity}, detector={state.amf_detector}, domain radius={state.amf_domain_radius_um:g} um",
                    f"Stopping power={state.amf_stopping_power}, step calculator={state.amf_step_calculator}",
                    f"Scoring position=({state.amf_scoring_x_mm:g}, {state.amf_scoring_y_mm:g}, {state.amf_scoring_z_mm:g}) mm",
                ]
            )

        self.results_panel.set_overview_data(
            OverviewData(
                run_status=self._run_status,
                approach=approach_line,
                details=detail_lines,
                phase_space_file=phase_space_display,
                lookup_root=str(state.lookup_root),
                output_dir=output_display,
                planned_outputs=[path.name for path in outputs],
                available_outputs=[path.name for path in self._loaded_results.file_paths],
                history=list(reversed(self.run_history)),
            )
        )

    def _run_analysis(self) -> None:
        input_result = self.input_checker.run(self.state)
        self.setup_panel.set_input_check_result(input_result)
        if input_result.status == "error":
            self._run_status = "Run blocked by input-check errors."
            self._sync_overview(self.state)
            QMessageBox.critical(
                self,
                "Input check failed",
                "Fix the input-check errors before starting the run.\n\n"
                + input_result.warnings_text,
            )
            return

        commands = build_command_specs(self.state)
        if not commands:
            QMessageBox.information(
                self,
                "Nothing to run",
                "Choose an analysis approach before starting the run.",
            )
            return

        self.results_panel.reset_to_overview()
        self.results_panel.clear_log()
        self._run_status = f"Run in progress. {len(commands)} command(s) queued."
        self._sync_overview(self.state)
        self.results_panel.append_log("Starting run...\n\n")
        self.runner.run_commands(commands, self.state.build_workdir)

    def _handle_command_started(self, rendered_command: str) -> None:
        self.results_panel.append_log(f"$ {rendered_command}\n")

    def _handle_run_finished(self, success: bool) -> None:
        if success:
            self._run_status = "Run completed successfully."
            self.results_panel.append_log("Run completed successfully.\n")
            self._refresh_loaded_results()
            self.results_panel.reset_to_overview()
        else:
            self._run_status = "Run stopped before completion."
            self.results_panel.append_log("Run stopped before completion.\n")
            self._sync_overview(self.state)

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.run_history.append(
            RunHistoryEntry(
                timestamp=timestamp,
                status="success" if success else "stopped",
                description=self._history_description(self.state),
            )
        )
        self.run_history = self.run_history[-12:]
        self._sync_overview(self.state)

    def _refresh_input_check(self) -> None:
        result = self.input_checker.run(self.state)
        self.setup_panel.set_input_check_result(result)

    def _refresh_loaded_results(self) -> None:
        self._loaded_results = self.results_loader.load(self.state)
        self.results_panel.set_loaded_results(self._loaded_results)
        self._sync_overview(self.state)

    def _history_description(self, state: AppState) -> str:
        if state.approach == "spectrum":
            if state.spectrum_family == "DeCunha":
                return (
                    "Spectrum build with DeCunha "
                    f"({state.decunha_voxel_size}, {state.decunha_energy_grid})"
                )
            return f"Spectrum build with Cartechini ({state.cartechini_radius})"

        if state.approach == "means":
            enabled: list[str] = []
            if state.enable_let:
                enabled.append("LET")
            if state.enable_magini:
                enabled.append("Magini")
            if state.enable_inaniwa:
                enabled.append("Inaniwa")
            return "Mean values only: " + (", ".join(enabled) if enabled else "no summary mode selected")

        if state.approach == "amf":
            mode = "replay" if state.amf_run_mode == "replay" else "full simulation"
            return (
                f"AMF {mode}: {state.amf_quantity}, "
                f"{state.amf_detector}, {state.amf_domain_radius_um:g} um"
            )

        return "No analysis approach selected"

    def _apply_stylesheet(self) -> None:
        self.setStyleSheet(
            """
            QMainWindow {
                background: #f4efe7;
            }
            QGroupBox {
                font-weight: 700;
                border: 1px solid #d8d0c4;
                border-radius: 10px;
                margin-top: 10px;
                padding-top: 10px;
                background: #fffaf2;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 4px 0 4px;
            }
            QFrame#OverviewCard {
                background: #fffaf2;
                border: 1px solid #d8d0c4;
                border-radius: 10px;
            }
            QLabel#OverviewCardTitle {
                color: #5d6c71;
                font-size: 12px;
                font-weight: 700;
                text-transform: uppercase;
            }
            QPlainTextEdit, QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QListWidget {
                background: #fffdf8;
                border: 1px solid #d8d0c4;
                border-radius: 8px;
                padding: 6px;
            }
            QComboBox, QSpinBox, QDoubleSpinBox {
                padding-top: 4px;
                padding-bottom: 4px;
            }
            QPushButton, QToolButton {
                background: #ffffff;
                border: 1px solid #d8d0c4;
                border-radius: 8px;
                padding: 8px 12px;
            }
            QPushButton:hover, QToolButton:hover {
                background: #f0ece4;
            }
            QLabel#EyebrowLabel {
                color: #5d6c71;
                text-transform: uppercase;
            }
            QLabel#PanelTitle {
                font-size: 24px;
                font-weight: 700;
                color: #19343b;
            }
            QLabel[status="ready"] {
                color: #2b7a55;
                font-weight: 700;
            }
            QLabel[status="warning"] {
                color: #b66a2f;
                font-weight: 700;
            }
            QLabel[status="error"] {
                color: #a13b2d;
                font-weight: 700;
            }
            """
        )


def launch() -> int:
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()
