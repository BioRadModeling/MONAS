from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QObject, QProcess, Signal

from app.services.command_builder import CommandSpec


class ProcessRunner(QObject):
    commandStarted = Signal(str)
    outputReady = Signal(str)
    runFinished = Signal(bool)

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._process = QProcess(self)
        self._process.readyReadStandardOutput.connect(self._handle_stdout)
        self._process.readyReadStandardError.connect(self._handle_stderr)
        self._process.finished.connect(self._handle_finished)

        self._commands: list[CommandSpec] = []
        self._current_index = -1
        self._workdir = Path.cwd()
        self._running = False

    @property
    def is_running(self) -> bool:
        return self._running

    def run_commands(self, commands: list[CommandSpec], workdir: Path) -> None:
        if self._running:
            self.outputReady.emit("A run is already in progress.\n")
            return

        if not commands:
            self.outputReady.emit("No commands to run.\n")
            self.runFinished.emit(False)
            return

        self._commands = commands
        self._current_index = -1
        self._workdir = workdir
        self._running = True
        self._start_next()

    def cancel(self) -> None:
        if self._running:
            self._process.kill()

    def _start_next(self) -> None:
        self._current_index += 1
        if self._current_index >= len(self._commands):
            self._running = False
            self.runFinished.emit(True)
            return

        command = self._commands[self._current_index]
        self.commandStarted.emit(command.render())
        self._process.setWorkingDirectory(str(self._workdir))
        self._process.start(command.argv[0], command.argv[1:])

    def _handle_stdout(self) -> None:
        text = bytes(self._process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if text:
            self.outputReady.emit(text)

    def _handle_stderr(self) -> None:
        text = bytes(self._process.readAllStandardError()).decode("utf-8", errors="replace")
        if text:
            self.outputReady.emit(text)

    def _handle_finished(self, exit_code: int, exit_status: QProcess.ExitStatus) -> None:
        if exit_status != QProcess.NormalExit or exit_code != 0:
            self._running = False
            self.outputReady.emit(
                f"\nCommand failed with exit code {exit_code}.\n"
            )
            self.runFinished.emit(False)
            return

        self.outputReady.emit("\nCommand completed successfully.\n\n")
        self._start_next()
