# MONAS Desktop Wrapper

This folder contains the first `PySide6` scaffold for a local desktop wrapper
around the existing `microdosimetry_with_LUTs` executable.

## Scope of this scaffold

The current implementation focuses on:

- the main desktop shell
- the approved analysis-path selection flow
- path selection controls
- real command preview generation
- sequential CLI execution through `QProcess`
- a results area with an overview tab and run log

It does not yet implement:

- full input inspection from the phase-space file
- CSV loading into tables
- spectrum plotting
- output-file discovery and preview
- polished visual styling

## Install

```bash
python3 -m pip install -r requirements.txt
```

## Run

From this folder:

```bash
python3 main.py
```

## Notes

- The GUI delegates calculation work to the existing C++ executable in
  `../build/microdosimetry_with_LUTs`.
- If that executable is missing, the GUI will still open, but running commands
  will fail until the CLI is built.
