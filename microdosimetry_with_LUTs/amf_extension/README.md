# AMF TOPAS Extension

This folder is reserved for the TOPAS AMF extension source files:

- `AMFDose.cc` / `AMFDose.hh`
- `ScoreAMFSpectra.cc` / `ScoreAMFSpectra.hh`
- `ScoreAMF_yD.cc` / `ScoreAMF_yD.hh`
- `ScoreAMF_yS.cc` / `ScoreAMF_yS.hh`

Keep this folder as a thin staging area for the TOPAS extension itself. The
`microdosimetry_with_LUTs` executable should orchestrate TOPAS runs and parse
outputs rather than reimplementing the scorer as the reference path.

## Registering With OpenTOPAS

OpenTOPAS registers extensions during CMake configure by scanning the first line
of each `.cc` file. The AMF files use these registration lines:

```text
// Scorer for AMFDose
// Scorer for AMFSpectra
// Scorer for AMF_yD
// Scorer for AMF_yS
```

For the local TOPAS build used during validation, AMF was added while preserving
the existing OpenTOPAS-Microdosimetry extensions:

```bash
cmake -S /Applications/TOPAS/OpenTOPAS \
  -B /Applications/TOPAS/OpenTOPAS-build \
  -DTOPAS_EXTENSIONS_DIR="/Applications/OpenTOPAS-Microdosimetry;/Users/kylediamond/MONAS/microdosimetry_with_LUTs/amf_extension"

cmake --build /Applications/TOPAS/OpenTOPAS-build --target topas -j8
```

The AMF headers contain one small portability fix for AppleClang/OpenTOPAS 4.2:
`mparased` is declared as `static constexpr int` in `ScoreAMF_yD.hh` and
`ScoreAMF_yS.hh`, so it can be used as a compile-time array size.
