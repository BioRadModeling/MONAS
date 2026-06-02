from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFontDatabase
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QRadioButton,
    QSlider,
    QSpinBox,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from app.services.command_builder import render_command_preview
from app.services.input_checker import InputCheckResult
from app.state import AppState, default_app_state


AMF_DOMAIN_RADIUS_MIN_UM = 0.0015
AMF_DOMAIN_RADIUS_MAX_UM = 0.5
AMF_DOMAIN_RADIUS_SLIDER_STEPS = 100_000


class SetupPanel(QWidget):
    stateChanged = Signal(object)
    runRequested = Signal()
    resetRequested = Signal()
    inputCheckRequested = Signal()

    def __init__(self, initial_state: AppState, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._state = initial_state.clone()
        self._project_root = initial_state.phase_space_file.parents[1]

        root_layout = QVBoxLayout(self)
        root_layout.setContentsMargins(0, 0, 0, 0)
        root_layout.setSpacing(12)

        root_layout.addWidget(self._build_header())
        root_layout.addWidget(self._build_settings_group())
        root_layout.addWidget(self._build_input_check_group())
        root_layout.addWidget(self._build_approach_group())
        root_layout.addWidget(self._build_amf_group())
        root_layout.addWidget(self._build_spectrum_group())
        root_layout.addWidget(self._build_means_group())
        root_layout.addWidget(self._build_command_group())
        root_layout.addStretch(1)

        self._load_state(self._state)
        self._refresh_visibility()
        self._update_command_preview()

    def current_state(self) -> AppState:
        return self._state.clone()

    def _build_header(self) -> QWidget:
        frame = QFrame()
        layout = QHBoxLayout(frame)
        layout.setContentsMargins(0, 0, 0, 0)

        title_box = QVBoxLayout()
        eyebrow = QLabel("Version 1.0")
        eyebrow.setObjectName("EyebrowLabel")
        title = QLabel("FAST MICRODOSIMETRY")
        title.setObjectName("PanelTitle")
        title_box.addWidget(eyebrow)
        title_box.addWidget(title)

        self.settings_button = QToolButton()
        self.settings_button.setText("Path settings")
        self.settings_button.setCheckable(True)
        self.settings_button.toggled.connect(self._toggle_settings_group)

        layout.addLayout(title_box, 1)
        layout.addWidget(self.settings_button, 0, Qt.AlignTop)
        return frame

    def _build_settings_group(self) -> QGroupBox:
        group = QGroupBox("Paths")
        group.setVisible(False)
        layout = QVBoxLayout(group)

        form = QFormLayout()
        self.phase_space_edit = QLineEdit()
        self.phase_space_edit.setReadOnly(True)
        self.lookup_root_edit = QLineEdit()
        self.lookup_root_edit.setReadOnly(True)
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setReadOnly(True)

        form.addRow(
            "Phase-space file",
            self._row_with_browse(self.phase_space_edit, self._choose_phase_space),
        )
        form.addRow(
            "Lookup root",
            self._row_with_browse(self.lookup_root_edit, self._choose_lookup_root),
        )
        form.addRow(
            "Output directory",
            self._row_with_browse(self.output_dir_edit, self._choose_output_dir),
        )
        layout.addLayout(form)

        buttons = QHBoxLayout()
        preview_button = QPushButton("Preview commands")
        preview_button.clicked.connect(self._show_preview_hint)
        refresh_button = QPushButton("Refresh input check")
        refresh_button.clicked.connect(self.inputCheckRequested.emit)
        buttons.addWidget(preview_button)
        buttons.addWidget(refresh_button)
        buttons.addStretch(1)
        layout.addLayout(buttons)

        self.settings_group = group
        return group

    def _build_input_check_group(self) -> QGroupBox:
        group = QGroupBox("Input Check")
        layout = QVBoxLayout(group)

        self.input_check_status = QLabel("Not checked")
        layout.addWidget(self.input_check_status)

        self.input_check_summary = QPlainTextEdit()
        self.input_check_summary.setReadOnly(True)
        self.input_check_summary.setMaximumHeight(110)
        layout.addWidget(self.input_check_summary)

        self.input_check_species = QPlainTextEdit()
        self.input_check_species.setReadOnly(True)
        self.input_check_species.setMaximumHeight(72)
        layout.addWidget(self.input_check_species)

        self.input_check_warnings = QPlainTextEdit()
        self.input_check_warnings.setReadOnly(True)
        self.input_check_warnings.setMaximumHeight(170)
        layout.addWidget(self.input_check_warnings)
        return group

    def _build_approach_group(self) -> QGroupBox:
        group = QGroupBox("Analysis Approach")
        layout = QVBoxLayout(group)
        layout.addWidget(QLabel("Choose one path first. Spectrum build is no longer the default."))

        self.approach_group = QButtonGroup(self)
        self.no_approach_radio = QRadioButton("No approach selected yet")
        self.amf_radio = QRadioButton("AMF")
        self.spectrum_radio = QRadioButton("Full microdosimetric spectrum")
        self.means_radio = QRadioButton("Mean values only")

        self.approach_group.addButton(self.no_approach_radio)
        self.approach_group.addButton(self.amf_radio)
        self.approach_group.addButton(self.spectrum_radio)
        self.approach_group.addButton(self.means_radio)

        for widget in (
            self.no_approach_radio,
            self.amf_radio,
            self.spectrum_radio,
            self.means_radio,
        ):
            widget.toggled.connect(self._sync_state_from_widgets)
            layout.addWidget(widget)
        return group

    def _build_amf_group(self) -> QGroupBox:
        group = QGroupBox("AMF Configuration")
        layout = QVBoxLayout(group)

        grid = QGridLayout()
        self.amf_quantity_combo = QComboBox()
        self.amf_quantity_combo.addItems(["AMFSpectra", "AMF_yD", "AMF_yS"])
        self.amf_domain_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.amf_domain_radius_slider.setRange(0, AMF_DOMAIN_RADIUS_SLIDER_STEPS)
        self.amf_domain_radius_slider.setSingleStep(100)
        self.amf_domain_radius_slider.setPageStep(2_500)
        self.amf_domain_radius_input = QDoubleSpinBox()
        self.amf_domain_radius_input.setRange(
            AMF_DOMAIN_RADIUS_MIN_UM,
            AMF_DOMAIN_RADIUS_MAX_UM,
        )
        self.amf_domain_radius_input.setDecimals(4)
        self.amf_domain_radius_input.setSingleStep(0.01)
        self.amf_domain_radius_input.setSuffix(" um")
        self.amf_domain_radius_input.setMinimumWidth(120)
        self.amf_stopping_power_combo = QComboBox()
        self.amf_stopping_power_combo.addItems(["Topas", "ExternalTable"])
        self.amf_precheck_checkbox = QCheckBox("Disable phase-space precheck")

        grid.addWidget(QLabel("Quantity"), 0, 0)
        grid.addWidget(self.amf_quantity_combo, 0, 1)
        grid.addWidget(QLabel("Domain radius"), 1, 0)
        grid.addWidget(self._build_domain_radius_control(), 1, 1)
        grid.addWidget(QLabel("Stopping power"), 2, 0)
        grid.addWidget(self.amf_stopping_power_combo, 2, 1)
        grid.addWidget(self.amf_precheck_checkbox, 3, 1)
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid)

        paths_form = QFormLayout()
        self.amf_phase_space_base_edit = QLineEdit()
        self.amf_phase_space_base_edit.setReadOnly(True)
        self.amf_source_topas_file_edit = QLineEdit()
        self.amf_source_topas_file_edit.setReadOnly(True)
        self.amf_output_dir_edit = QLineEdit()
        self.amf_external_stopping_power_edit = QLineEdit()
        self.amf_external_stopping_power_edit.setReadOnly(True)

        self.amf_replay_base_row = self._row_with_browse(
            self.amf_phase_space_base_edit,
            self._choose_amf_phase_space_base,
        )
        self.amf_source_topas_file_row = self._row_with_browse(
            self.amf_source_topas_file_edit,
            self._choose_amf_source_topas_file,
        )
        self.amf_output_dir_row = self._row_with_browse(
            self.amf_output_dir_edit,
            self._choose_amf_output_dir,
        )
        self.amf_external_stopping_power_row = self._row_with_browse(
            self.amf_external_stopping_power_edit,
            self._choose_amf_external_stopping_power_file,
        )

        self.amf_paths_form = paths_form
        paths_form.addRow("Phase-space base", self.amf_replay_base_row)
        paths_form.addRow("Source TOPAS simulation txt", self.amf_source_topas_file_row)
        paths_form.addRow("Output directory", self.amf_output_dir_row)
        paths_form.addRow("StoppingPower.txt", self.amf_external_stopping_power_row)
        self.amf_external_stopping_power_label = paths_form.labelForField(
            self.amf_external_stopping_power_row
        )
        layout.addLayout(paths_form)

        for widget in (
            self.amf_quantity_combo,
            self.amf_stopping_power_combo,
        ):
            widget.currentIndexChanged.connect(self._sync_state_from_widgets)
        self.amf_domain_radius_slider.valueChanged.connect(self._handle_domain_radius_slider_changed)
        self.amf_domain_radius_input.valueChanged.connect(self._handle_domain_radius_input_changed)
        self.amf_precheck_checkbox.toggled.connect(self._sync_state_from_widgets)
        self.amf_output_dir_edit.editingFinished.connect(self._sync_state_from_widgets)

        self.amf_group_box = group
        return group

    def _build_domain_radius_control(self) -> QWidget:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.amf_domain_radius_slider, 1)
        layout.addWidget(self.amf_domain_radius_input)
        return row

    def _build_spectrum_group(self) -> QGroupBox:
        group = QGroupBox("Spectrum Build")
        layout = QVBoxLayout(group)
        layout.addWidget(QLabel("Choose the full-spectrum family and required specifications."))

        family_row = QHBoxLayout()
        self.decunha_radio = QRadioButton("DeCunha")
        self.cartechini_radio = QRadioButton("Cartechini")
        self.decunha_radio.toggled.connect(self._sync_state_from_widgets)
        self.cartechini_radio.toggled.connect(self._sync_state_from_widgets)
        family_row.addWidget(self.decunha_radio)
        family_row.addWidget(self.cartechini_radio)
        family_row.addStretch(1)
        layout.addLayout(family_row)

        grid = QGridLayout()
        self.decunha_voxel_combo = QComboBox()
        self.decunha_voxel_combo.addItems(["1mm", "5um"])
        self.decunha_grid_combo = QComboBox()
        self.decunha_grid_combo.addItems(["linear", "log"])
        self.cartechini_radius_combo = QComboBox()
        self.cartechini_radius_combo.addItems(["0.5um", "8um"])
        self.rebin_samples_spin = QSpinBox()
        self.rebin_samples_spin.setRange(1, 100_000_000)
        self.rebin_seed_spin = QSpinBox()
        self.rebin_seed_spin.setRange(0, 2_000_000_000)
        for widget in (
            self.decunha_voxel_combo,
            self.decunha_grid_combo,
            self.cartechini_radius_combo,
            self.rebin_samples_spin,
            self.rebin_seed_spin,
        ):
            widget.setMinimumHeight(32)

        grid.addWidget(QLabel("DeCunha voxel size"), 0, 0)
        grid.addWidget(self.decunha_voxel_combo, 0, 1)
        grid.addWidget(QLabel("DeCunha energy grid"), 1, 0)
        grid.addWidget(self.decunha_grid_combo, 1, 1)
        grid.addWidget(QLabel("Cartechini radius"), 2, 0)
        grid.addWidget(self.cartechini_radius_combo, 2, 1)
        grid.addWidget(QLabel("Rebin samples"), 3, 0)
        grid.addWidget(self.rebin_samples_spin, 3, 1)
        grid.addWidget(QLabel("Rebin seed"), 4, 0)
        grid.addWidget(self.rebin_seed_spin, 4, 1)
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid)

        for widget in (
            self.decunha_voxel_combo,
            self.decunha_grid_combo,
            self.cartechini_radius_combo,
        ):
            widget.currentIndexChanged.connect(self._sync_state_from_widgets)
        self.rebin_samples_spin.valueChanged.connect(self._sync_state_from_widgets)
        self.rebin_seed_spin.valueChanged.connect(self._sync_state_from_widgets)

        self.spectrum_group_box = group
        return group

    def _build_means_group(self) -> QGroupBox:
        group = QGroupBox("Mean Value Calculations")
        layout = QVBoxLayout(group)
        layout.addWidget(QLabel("Enable one or more modes to calculate mean microdosimetric values without generating the full spectrum."))

        self.let_checkbox = QCheckBox("LET")
        self.magini_checkbox = QCheckBox("Magini")
        self.inaniwa_checkbox = QCheckBox("Inaniwa")

        for widget in (self.let_checkbox, self.magini_checkbox, self.inaniwa_checkbox):
            widget.toggled.connect(self._sync_state_from_widgets)
            layout.addWidget(widget)

        #note = QLabel(
        #    "Current modes reflect the repository's available mean-value-only LUTs "
        #    ": LET, Magini, and Inaniwa."
        #)
        #note.setWordWrap(True)
        #layout.addWidget(note)

        self.means_group_box = group
        return group

    def _build_command_group(self) -> QGroupBox:
        group = QGroupBox("Command Preview")
        layout = QVBoxLayout(group)

        self.command_preview = QPlainTextEdit()
        self.command_preview.setReadOnly(True)
        self.command_preview.setMinimumHeight(120)
        self.command_preview.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.command_preview.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        layout.addWidget(self.command_preview)

        button_row = QHBoxLayout()
        self.run_button = QPushButton("Run analysis")
        self.reset_button = QPushButton("Reset")
        self.run_button.clicked.connect(self.runRequested.emit)
        self.reset_button.clicked.connect(self._reset_state)
        button_row.addWidget(self.run_button)
        button_row.addWidget(self.reset_button)
        button_row.addStretch(1)
        layout.addLayout(button_row)
        return group

    def _row_with_browse(self, edit: QLineEdit, callback) -> QWidget:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        button = QPushButton("Browse")
        button.clicked.connect(callback)
        layout.addWidget(edit, 1)
        layout.addWidget(button)
        return row

    def _toggle_settings_group(self, checked: bool) -> None:
        self.settings_group.setVisible(checked)

    def _choose_phase_space(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Choose phase-space file",
            str(self._state.phase_space_file.parent),
            "Phase-space files (*.phsp *.txt);;All files (*)",
        )
        if path:
            self.phase_space_edit.setText(path)
            self._sync_state_from_widgets()

    def _choose_lookup_root(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Choose lookup root",
            str(self._state.lookup_root),
        )
        if path:
            self.lookup_root_edit.setText(path)
            self._sync_state_from_widgets()

    def _choose_output_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Choose output directory",
            str(self._state.output_dir),
        )
        if path:
            self.output_dir_edit.setText(path)
            self._sync_state_from_widgets()

    def _choose_amf_phase_space_base(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Choose AMF phase-space file",
            str(self._state.amf_phase_space_base.parent),
            "Phase-space files (*.phsp *.header);;All files (*)",
        )
        if path:
            selected = Path(path)
            if selected.suffix in {".phsp", ".header"}:
                selected = selected.with_suffix("")
            self.amf_phase_space_base_edit.setText(str(selected))
            self._sync_state_from_widgets()

    def _choose_amf_source_topas_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Choose source TOPAS simulation txt file",
            str(self._state.amf_source_topas_file.parent),
            "TOPAS parameter files (*.txt);;All files (*)",
        )
        if path:
            self.amf_source_topas_file_edit.setText(path)
            self._sync_state_from_widgets()

    def _choose_amf_output_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Choose AMF output directory",
            str(self._state.amf_output_dir.parent),
        )
        if path:
            self.amf_output_dir_edit.setText(path)
            self._sync_state_from_widgets()

    def _choose_amf_external_stopping_power_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Choose external stopping-power table",
            str(self._state.amf_external_stopping_power_file.parent),
            "Text files (*.txt);;All files (*)",
        )
        if path:
            self.amf_external_stopping_power_edit.setText(path)
            self._sync_state_from_widgets()

    def _show_preview_hint(self) -> None:
        QMessageBox.information(
            self,
            "Preview commands",
            "The command preview below updates from the current selections.",
        )

    def _reset_state(self) -> None:
        fresh_state = default_app_state(self._project_root)
        self._load_state(fresh_state)
        self.resetRequested.emit()

    def _load_state(self, state: AppState) -> None:
        self._state = state.clone()
        self.phase_space_edit.setText(str(state.phase_space_file))
        self.lookup_root_edit.setText(str(state.lookup_root))
        self.output_dir_edit.setText(str(state.output_dir))

        self.no_approach_radio.setChecked(state.approach is None)
        self.amf_radio.setChecked(state.approach == "amf")
        self.spectrum_radio.setChecked(state.approach == "spectrum")
        self.means_radio.setChecked(state.approach == "means")

        self.decunha_radio.setChecked(state.spectrum_family == "DeCunha")
        self.cartechini_radio.setChecked(state.spectrum_family == "Cartechini")
        self.decunha_voxel_combo.setCurrentText(state.decunha_voxel_size)
        self.decunha_grid_combo.setCurrentText(state.decunha_energy_grid)
        self.cartechini_radius_combo.setCurrentText(state.cartechini_radius)
        self.rebin_samples_spin.setValue(state.rebin_samples)
        self.rebin_seed_spin.setValue(state.rebin_seed)

        self.let_checkbox.setChecked(state.enable_let)
        self.magini_checkbox.setChecked(state.enable_magini)
        self.inaniwa_checkbox.setChecked(state.enable_inaniwa)

        self.amf_quantity_combo.setCurrentText(state.amf_quantity)
        self._set_domain_radius_controls(state.amf_domain_radius_um)
        self.amf_stopping_power_combo.setCurrentText(state.amf_stopping_power)
        self.amf_phase_space_base_edit.setText(str(state.amf_phase_space_base))
        self.amf_source_topas_file_edit.setText(str(state.amf_source_topas_file))
        self.amf_output_dir_edit.setText(str(state.amf_output_dir))
        self.amf_external_stopping_power_edit.setText(str(state.amf_external_stopping_power_file))
        self.amf_precheck_checkbox.setChecked(state.amf_disable_phase_space_precheck)

        self._refresh_visibility()
        self._update_command_preview()
        self.stateChanged.emit(self.current_state())

    def _sync_state_from_widgets(self) -> None:
        if self.amf_radio.isChecked():
            approach = "amf"
        elif self.spectrum_radio.isChecked():
            approach = "spectrum"
        elif self.means_radio.isChecked():
            approach = "means"
        else:
            approach = None

        spectrum_family = "Cartechini" if self.cartechini_radio.isChecked() else "DeCunha"

        self._state.phase_space_file = Path(self.phase_space_edit.text())
        self._state.lookup_root = Path(self.lookup_root_edit.text())
        self._state.output_dir = Path(self.output_dir_edit.text())
        self._state.approach = approach
        self._state.spectrum_family = spectrum_family
        self._state.decunha_voxel_size = self.decunha_voxel_combo.currentText()
        self._state.decunha_energy_grid = self.decunha_grid_combo.currentText()
        self._state.cartechini_radius = self.cartechini_radius_combo.currentText()
        self._state.rebin_samples = self.rebin_samples_spin.value()
        self._state.rebin_seed = self.rebin_seed_spin.value()
        self._state.enable_let = self.let_checkbox.isChecked()
        self._state.enable_magini = self.magini_checkbox.isChecked()
        self._state.enable_inaniwa = self.inaniwa_checkbox.isChecked()
        self._state.amf_quantity = self.amf_quantity_combo.currentText()
        self._state.amf_domain_radius_um = self.amf_domain_radius_input.value()
        self._state.amf_stopping_power = self.amf_stopping_power_combo.currentText()
        self._state.amf_phase_space_base = Path(self.amf_phase_space_base_edit.text())
        self._state.amf_source_topas_file = Path(self.amf_source_topas_file_edit.text())
        self._state.amf_output_dir = Path(self.amf_output_dir_edit.text())
        self._state.amf_external_stopping_power_file = Path(self.amf_external_stopping_power_edit.text())
        self._state.amf_disable_phase_space_precheck = self.amf_precheck_checkbox.isChecked()

        self._refresh_visibility()
        self._update_command_preview()
        self.stateChanged.emit(self.current_state())

    def _refresh_visibility(self) -> None:
        self.amf_group_box.setVisible(self._state.approach == "amf")
        self.spectrum_group_box.setVisible(self._state.approach == "spectrum")
        self.means_group_box.setVisible(self._state.approach == "means")

        decunha_selected = self._state.spectrum_family == "DeCunha"
        self.decunha_voxel_combo.setEnabled(decunha_selected)
        self.decunha_grid_combo.setEnabled(decunha_selected)
        self.cartechini_radius_combo.setEnabled(not decunha_selected)

        external_selected = self._state.amf_stopping_power == "ExternalTable"
        self.amf_external_stopping_power_row.setVisible(external_selected)
        self.amf_external_stopping_power_label.setVisible(external_selected)

    def _update_command_preview(self) -> None:
        self.command_preview.setPlainText(render_command_preview(self._state))

    def _handle_domain_radius_slider_changed(self, value: int) -> None:
        radius_um = self._slider_to_domain_radius(value)
        self.amf_domain_radius_input.blockSignals(True)
        self.amf_domain_radius_input.setValue(radius_um)
        self.amf_domain_radius_input.blockSignals(False)
        self._sync_state_from_widgets()

    def _handle_domain_radius_input_changed(self, value: float) -> None:
        self.amf_domain_radius_slider.blockSignals(True)
        self.amf_domain_radius_slider.setValue(self._domain_radius_to_slider(value))
        self.amf_domain_radius_slider.blockSignals(False)
        self._sync_state_from_widgets()

    def _set_domain_radius_controls(self, radius_um: float) -> None:
        self.amf_domain_radius_slider.blockSignals(True)
        self.amf_domain_radius_input.blockSignals(True)
        self.amf_domain_radius_slider.setValue(self._domain_radius_to_slider(radius_um))
        self.amf_domain_radius_input.setValue(self._slider_to_domain_radius(self.amf_domain_radius_slider.value()))
        self.amf_domain_radius_input.blockSignals(False)
        self.amf_domain_radius_slider.blockSignals(False)

    def _slider_to_domain_radius(self, value: int) -> float:
        fraction = value / AMF_DOMAIN_RADIUS_SLIDER_STEPS
        radius_um = AMF_DOMAIN_RADIUS_MIN_UM + (
            AMF_DOMAIN_RADIUS_MAX_UM - AMF_DOMAIN_RADIUS_MIN_UM
        ) * fraction
        return round(radius_um, 4)

    def _domain_radius_to_slider(self, radius_um: float) -> int:
        bounded = max(AMF_DOMAIN_RADIUS_MIN_UM, min(AMF_DOMAIN_RADIUS_MAX_UM, radius_um))
        fraction = (
            (bounded - AMF_DOMAIN_RADIUS_MIN_UM)
            / (AMF_DOMAIN_RADIUS_MAX_UM - AMF_DOMAIN_RADIUS_MIN_UM)
        )
        return round(fraction * AMF_DOMAIN_RADIUS_SLIDER_STEPS)

    def set_input_check_result(self, result: InputCheckResult) -> None:
        self.input_check_status.setText(result.status_text)
        self.input_check_status.setProperty("status", result.status)
        self.style().unpolish(self.input_check_status)
        self.style().polish(self.input_check_status)
        self.input_check_summary.setPlainText(result.summary_text)
        self.input_check_species.setPlainText(result.species_text)
        self.input_check_warnings.setPlainText(result.warnings_text)
