from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from PySide6.QtCharts import QChart, QChartView, QLineSeries, QLogValueAxis, QValueAxis
from PySide6.QtCore import QPointF, Qt, QUrl
from PySide6.QtGui import QDesktopServices, QFontDatabase, QPainter, QPixmap
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QStackedWidget,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from app.services.results_loader import LoadedResults, SpectrumPoint


@dataclass
class RunHistoryEntry:
    timestamp: str
    status: str
    description: str


@dataclass
class OverviewData:
    run_status: str = "Ready to configure a run."
    approach: str = "No analysis approach selected yet."
    details: list[str] = field(default_factory=list)
    phase_space_file: str = ""
    lookup_root: str = ""
    output_dir: str = ""
    planned_outputs: list[str] = field(default_factory=list)
    available_outputs: list[str] = field(default_factory=list)
    history: list[RunHistoryEntry] = field(default_factory=list)


class ResultsPanel(QWidget):
    _SPECTRUM_METRICS = {
        "yf_y": ("yf(y)", "yf(y)"),
        "yd_y": ("yd(y)", "yd(y)"),
    }

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        self.tabs.addTab(self._build_overview_tab(), "Overview")
        self.tabs.addTab(self._build_spectrum_tab(), "Spectrum")

        self.let_table = QTableWidget()
        self.tabs.addTab(self.let_table, "LET")

        self.magini_table = self._build_key_value_table()
        self.tabs.addTab(self.magini_table, "Magini")

        self.inaniwa_table = self._build_key_value_table()
        self.tabs.addTab(self.inaniwa_table, "Inaniwa")

        self.tabs.addTab(self._build_files_tab(), "Files")

        self.run_log = QPlainTextEdit()
        self.run_log.setReadOnly(True)
        self.run_log.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.tabs.addTab(self.run_log, "Run Log")

        self._planned_outputs: list[Path] = []
        self._file_previews: dict[str, str] = {}
        self._selected_file: Path | None = None
        self._selected_pixmap: QPixmap | None = None
        self._spectrum_points: list[SpectrumPoint] = []
        self._spectrum_file: Path | None = None
        self._preferred_file: Path | None = None
        self._preferred_metric_files: dict[str, Path] = {}
        self._programmatic_file_selection = False
        self._auto_selected_file = False

    def _build_overview_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)

        cards = QGridLayout()
        cards.setHorizontalSpacing(12)
        cards.setVerticalSpacing(12)

        self.status_value = QLabel("Ready to configure a run.")
        self.approach_value = QLabel("No analysis approach selected yet.")
        self.details_value = QLabel("Choose a workflow to see calculation-specific details.")
        self.paths_value = QLabel("")
        self.outputs_value = QLabel("Planned outputs will appear here.")

        for label in (
            self.status_value,
            self.approach_value,
            self.details_value,
            self.paths_value,
            self.outputs_value,
        ):
            label.setWordWrap(True)
            label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        cards.addWidget(self._build_card("Run Status", self.status_value), 0, 0)
        cards.addWidget(self._build_card("Approach", self.approach_value), 0, 1)
        cards.addWidget(self._build_card("Current Details", self.details_value), 1, 0)
        cards.addWidget(self._build_card("Paths", self.paths_value), 1, 1)
        cards.addWidget(self._build_card("Outputs", self.outputs_value), 2, 0, 1, 2)

        self.history_list = QListWidget()
        history_card = self._build_card("Recent Runs", self.history_list)

        layout.addLayout(cards)
        layout.addWidget(history_card, 1)
        return tab

    def _build_spectrum_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Metric"))
        self.spectrum_metric_combo = QComboBox()
        self.spectrum_metric_combo.addItem("yf(y)", "yf_y")
        self.spectrum_metric_combo.addItem("yd(y)", "yd_y")
        self.spectrum_metric_combo.setMinimumContentsLength(8)
        self.spectrum_metric_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToContents
        )
        self.spectrum_metric_combo.setMinimumWidth(120)
        self.spectrum_metric_combo.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Fixed,
        )
        self.spectrum_metric_combo.view().setMinimumWidth(120)
        self.spectrum_metric_combo.setCurrentIndex(1)
        self.spectrum_metric_combo.currentIndexChanged.connect(self._refresh_spectrum_chart)
        self.spectrum_metric_combo.currentIndexChanged.connect(self._sync_preferred_file_with_metric)
        controls.addWidget(self.spectrum_metric_combo)

        self.export_plot_button = QPushButton("Export plot")
        self.export_plot_button.clicked.connect(self._export_plot)
        controls.addWidget(self.export_plot_button)
        controls.addStretch(1)
        layout.addLayout(controls)

        self.spectrum_source_label = QLabel("No spectrum file loaded.")
        self.spectrum_source_label.setWordWrap(True)
        self.spectrum_source_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.spectrum_source_label)

        self.spectrum_summary = QPlainTextEdit()
        self.spectrum_summary.setReadOnly(True)
        self.spectrum_summary.setMaximumHeight(120)
        self.spectrum_summary.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        layout.addWidget(self.spectrum_summary)

        self.spectrum_chart = QChartView()
        self.spectrum_chart.setRenderHint(QPainter.RenderHint.Antialiasing)
        layout.addWidget(self.spectrum_chart, 1)
        return tab

    def _build_files_tab(self) -> QWidget:
        self.files_list = QListWidget()
        self.files_list.currentItemChanged.connect(self._handle_file_selection)

        self.file_selection_hint = QLabel("Select an output file to preview it here.")
        self.file_selection_hint.setWordWrap(True)

        self.file_path_label = QLabel("Select an output file to preview it here.")
        self.file_path_label.setWordWrap(True)
        self.file_path_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self.file_text_preview = QPlainTextEdit()
        self.file_text_preview.setReadOnly(True)
        self.file_text_preview.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))

        self.file_image_label = QLabel("Image preview will appear here.")
        self.file_image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.file_image_label.setMinimumHeight(260)
        self.file_image_label.setStyleSheet("border: 1px solid #d8d0c4; border-radius: 8px;")

        self.file_preview_stack = QStackedWidget()
        self.file_preview_stack.addWidget(self.file_text_preview)
        self.file_preview_stack.addWidget(self.file_image_label)

        preview_panel = QWidget()
        preview_layout = QVBoxLayout(preview_panel)
        preview_layout.setContentsMargins(0, 0, 0, 0)
        preview_layout.addWidget(self.file_selection_hint)
        preview_layout.addWidget(self.file_path_label)
        preview_layout.addWidget(self.file_preview_stack, 1)

        self.open_file_button = QPushButton("Open selected file")
        self.open_file_button.clicked.connect(self._open_selected_file)
        self.open_file_button.setEnabled(False)

        preview_layout.addWidget(self.open_file_button, 0, Qt.AlignmentFlag.AlignLeft)

        file_splitter = QSplitter()
        file_splitter.addWidget(self.files_list)
        file_splitter.addWidget(preview_panel)
        file_splitter.setSizes([260, 520])

        files_tab = QWidget()
        files_layout = QVBoxLayout(files_tab)
        files_layout.setContentsMargins(0, 0, 0, 0)
        files_layout.addWidget(file_splitter, 1)
        return files_tab

    def _build_card(self, title: str, content_widget: QWidget) -> QFrame:
        card = QFrame()
        card.setObjectName("OverviewCard")

        layout = QVBoxLayout(card)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(8)

        title_label = QLabel(title)
        title_label.setObjectName("OverviewCardTitle")
        layout.addWidget(title_label)
        layout.addWidget(content_widget, 1)
        return card

    def _build_key_value_table(self) -> QTableWidget:
        table = QTableWidget(0, 2)
        table.setHorizontalHeaderLabels(["Metric", "Value"])
        table.horizontalHeader().setStretchLastSection(True)
        return table

    def reset_to_overview(self) -> None:
        self.tabs.setCurrentIndex(0)

    def clear_log(self) -> None:
        self.run_log.clear()

    def append_log(self, text: str) -> None:
        cursor = self.run_log.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        cursor.insertText(text)
        self.run_log.setTextCursor(cursor)
        self.run_log.ensureCursorVisible()

    def set_planned_outputs(self, paths: list[Path]) -> None:
        self._planned_outputs = paths[:]
        if not self._selected_file and not self.files_list.count():
            self._show_planned_outputs_placeholder()

    def set_overview_data(self, overview: OverviewData) -> None:
        self.status_value.setText(overview.run_status)
        self.approach_value.setText(overview.approach)
        self.details_value.setText("\n".join(overview.details) if overview.details else "No workflow details yet.")
        self.paths_value.setText(
            "\n".join(
                (
                    f"Phase-space: {overview.phase_space_file}",
                    f"Lookup root: {overview.lookup_root}",
                    f"Output directory: {overview.output_dir}",
                )
            )
        )

        outputs_lines: list[str] = []
        if overview.planned_outputs:
            outputs_lines.extend(self._format_output_lines("Planned", overview.planned_outputs))
        if overview.available_outputs:
            if outputs_lines:
                outputs_lines.append("")
            outputs_lines.extend(
                self._format_output_lines(
                    "Available in output directory",
                    overview.available_outputs,
                )
            )
        self.outputs_value.setText("\n".join(outputs_lines) if outputs_lines else "No output files planned yet.")

        self.history_list.clear()
        if overview.history:
            for entry in overview.history[:8]:
                item = QListWidgetItem(f"{entry.timestamp}  [{entry.status}]  {entry.description}")
                self.history_list.addItem(item)
        else:
            self.history_list.addItem("No runs recorded in this session yet.")

    def set_loaded_results(self, loaded: LoadedResults) -> None:
        self.spectrum_summary.setPlainText(loaded.spectrum_summary_text)
        self.spectrum_source_label.setText(
            f"Spectrum source: {loaded.spectrum_file}" if loaded.spectrum_file else "No spectrum file loaded."
        )
        self._spectrum_points = loaded.spectrum_points
        self._spectrum_file = loaded.spectrum_file
        self._preferred_file = loaded.preferred_file
        self._preferred_metric_files = loaded.preferred_metric_files
        self._refresh_spectrum_chart()
        self._populate_table(self.let_table, loaded.let_headers, loaded.let_rows)
        self._populate_key_value_table(self.magini_table, loaded.magini_pairs)
        self._populate_key_value_table(self.inaniwa_table, loaded.inaniwa_pairs)
        self.set_loaded_files(loaded.file_paths, loaded.file_previews)

    def set_loaded_files(self, paths: list[Path], previews: dict[str, str]) -> None:
        self.files_list.clear()
        self._file_previews = previews
        self._selected_file = None
        self._selected_pixmap = None
        self.open_file_button.setEnabled(False)
        self._auto_selected_file = False
        self.file_selection_hint.setText("Select an output file to preview it here.")

        for path in paths:
            item = QListWidgetItem(path.name)
            item.setData(Qt.ItemDataRole.UserRole, str(path))
            self.files_list.addItem(item)

        if paths:
            self._select_preferred_file()
        else:
            self._show_planned_outputs_placeholder()

    def _refresh_spectrum_chart(self) -> None:
        metric_key = self.spectrum_metric_combo.currentData()
        series_name, axis_title = self._SPECTRUM_METRICS.get(metric_key, ("yd(y)", "yd(y)"))

        chart = QChart()
        chart.setTitle(f"{series_name} vs log(y)")
        chart.legend().hide()

        if not self._spectrum_points:
            chart.addSeries(QLineSeries())
            self.spectrum_chart.setChart(chart)
            self.export_plot_button.setEnabled(False)
            return

        self.export_plot_button.setEnabled(True)

        series = QLineSeries()
        series.setName(series_name)
        y_values: list[float] = []

        for point in self._spectrum_points:
            value = getattr(point, metric_key)
            y_values.append(value)
            series.append(QPointF(point.y_keV_per_um, value))

        chart.addSeries(series)

        x_values = [point.y_keV_per_um for point in self._spectrum_points if point.y_keV_per_um > 0.0]
        x_axis = QLogValueAxis()
        x_axis.setBase(10.0)
        x_axis.setTitleText("y (keV/um, log scale)")
        x_axis.setLabelFormat("%.3g")
        x_min = min(x_values) if x_values else 1e-3
        x_max = max(x_values) if x_values else 1.0
        if x_min == x_max:
            x_max = x_min * 10.0
        x_axis.setRange(x_min, x_max)
        chart.addAxis(x_axis, Qt.AlignmentFlag.AlignBottom)
        series.attachAxis(x_axis)

        y_axis = QValueAxis()
        y_axis.setTitleText(axis_title)
        y_axis.setLabelFormat("%.4f")
        y_min = min(y_values)
        y_max = max(y_values)
        if y_min == y_max:
            y_max = y_min + 1.0
        y_axis.setRange(y_min, y_max)
        chart.addAxis(y_axis, Qt.AlignmentFlag.AlignLeft)
        series.attachAxis(y_axis)

        self.spectrum_chart.setChart(chart)

    def _populate_table(
        self,
        table: QTableWidget,
        headers: list[str],
        rows: list[list[str]],
    ) -> None:
        table.clear()
        if not headers:
            table.setColumnCount(1)
            table.setRowCount(1)
            table.setHorizontalHeaderLabels(["Status"])
            table.setItem(0, 0, QTableWidgetItem("No summary data found."))
            return

        table.setColumnCount(len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.setRowCount(len(rows))

        for row_index, row in enumerate(rows):
            for column_index, value in enumerate(row):
                table.setItem(row_index, column_index, QTableWidgetItem(value))

        table.resizeColumnsToContents()

    def _populate_key_value_table(
        self,
        table: QTableWidget,
        pairs: list[tuple[str, str]],
    ) -> None:
        table.clear()
        table.setColumnCount(2)
        table.setHorizontalHeaderLabels(["Metric", "Value"])

        if not pairs:
            table.setRowCount(1)
            table.setItem(0, 0, QTableWidgetItem("Status"))
            table.setItem(0, 1, QTableWidgetItem("No summary data found."))
            table.resizeColumnsToContents()
            return

        table.setRowCount(len(pairs))
        for row_index, (metric, value) in enumerate(pairs):
            table.setItem(row_index, 0, QTableWidgetItem(metric))
            table.setItem(row_index, 1, QTableWidgetItem(value))
        table.resizeColumnsToContents()

    def _handle_file_selection(
        self,
        current: QListWidgetItem | None,
        previous: QListWidgetItem | None,
    ) -> None:
        del previous
        if current is None:
            self._selected_file = None
            self._selected_pixmap = None
            self._auto_selected_file = False
            self.file_selection_hint.setText("Select an output file to preview it here.")
            self.file_path_label.setText("Select an output file to preview it here.")
            self.file_text_preview.clear()
            self.file_image_label.setText("Image preview will appear here.")
            self.file_image_label.setPixmap(QPixmap())
            self.open_file_button.setEnabled(False)
            return

        key = current.data(Qt.ItemDataRole.UserRole)
        self._selected_file = Path(key)
        self._auto_selected_file = self._programmatic_file_selection
        self._programmatic_file_selection = False
        if self._auto_selected_file:
            self.file_selection_hint.setText("Showing the file that best matches the current analysis state.")
        else:
            self.file_selection_hint.setText("Manual file selection. Metric changes will not override this choice.")
        self.file_path_label.setText(str(self._selected_file))
        self.open_file_button.setEnabled(True)

        if self._selected_file.suffix.lower() in {".jpg", ".jpeg", ".png"}:
            pixmap = QPixmap(str(self._selected_file))
            self._selected_pixmap = pixmap if not pixmap.isNull() else None
            self._update_image_preview()
            self.file_preview_stack.setCurrentWidget(self.file_image_label)
            return

        self._selected_pixmap = None
        self.file_text_preview.setPlainText(self._file_previews.get(key, "No preview available."))
        self.file_preview_stack.setCurrentWidget(self.file_text_preview)

    def _sync_preferred_file_with_metric(self) -> None:
        if not self._preferred_metric_files:
            return
        if self._selected_file is None or self._auto_selected_file:
            self._select_preferred_file()

    def _select_preferred_file(self) -> None:
        if not self.files_list.count():
            return

        desired = self._preferred_file_for_current_metric() or self._preferred_file
        if desired is None:
            self._programmatic_file_selection = True
            self.files_list.setCurrentRow(0)
            return

        for index in range(self.files_list.count()):
            item = self.files_list.item(index)
            if item.data(Qt.ItemDataRole.UserRole) == str(desired):
                self._programmatic_file_selection = True
                self.files_list.setCurrentRow(index)
                return

        self._programmatic_file_selection = True
        self.files_list.setCurrentRow(0)

    def _preferred_file_for_current_metric(self) -> Path | None:
        metric_key = self.spectrum_metric_combo.currentData()
        if metric_key in self._preferred_metric_files:
            return self._preferred_metric_files[metric_key]
        return None

    def _open_selected_file(self) -> None:
        if self._selected_file is None:
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self._selected_file)))

    def _export_plot(self) -> None:
        if not self._spectrum_points:
            return

        default_name = "spectrum_plot.png"
        if self._spectrum_file is not None:
            default_name = f"{self._spectrum_file.stem}_{self.spectrum_metric_combo.currentData()}.png"

        path, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Export spectrum plot",
            default_name,
            "PNG Images (*.png);;JPEG Images (*.jpg *.jpeg)",
        )
        if not path:
            return

        chosen = Path(path)
        if not chosen.suffix:
            if "JPEG" in selected_filter:
                chosen = chosen.with_suffix(".jpg")
            else:
                chosen = chosen.with_suffix(".png")

        self.spectrum_chart.grab().save(str(chosen))

    def _show_planned_outputs_placeholder(self) -> None:
        self.file_selection_hint.setText("Generated files will appear here after a run or when you point to an output folder.")
        self.file_path_label.setText("No generated files found yet.")
        if self._planned_outputs:
            self.file_text_preview.setPlainText(
                "Planned output files:\n\n" + "\n".join(str(path) for path in self._planned_outputs)
            )
        else:
            self.file_text_preview.setPlainText("No output files planned yet.")
        self.file_preview_stack.setCurrentWidget(self.file_text_preview)
        self.file_image_label.setPixmap(QPixmap())
        self.file_image_label.setText("Image preview will appear here.")

    def _update_image_preview(self) -> None:
        if self._selected_pixmap is None:
            self.file_image_label.setText("Unable to preview this image.")
            self.file_image_label.setPixmap(QPixmap())
            return

        scaled = self._selected_pixmap.scaled(
            self.file_image_label.size(),
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        self.file_image_label.setPixmap(scaled)
        self.file_image_label.setText("")

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if self.file_preview_stack.currentWidget() is self.file_image_label:
            self._update_image_preview()

    def _format_output_lines(self, heading: str, names: list[str], max_items: int = 5) -> list[str]:
        lines = [heading]
        visible = names[:max_items]
        lines.extend(f"- {name}" for name in visible)
        if len(names) > max_items:
            lines.append(f"... +{len(names) - max_items} more")
        return lines
