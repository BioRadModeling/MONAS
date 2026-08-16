#!/bin/sh
set -eu

MONAS_BIN=${MONAS_BIN:-../bin/monas}
OUTPUT_DIR=${MONAS_OUTPUT_DIR:-GSM2_from_spectra_smoke_output}

"$MONAS_BIN" \
  -MKMFromSpectra \
  -outputDir "$OUTPUT_DIR" \
  -spectrum DeCunha:poly:poly_spectrum_fixture.csv \
  -GSM2 \
  -fSetMultiEventStatistic 100 \
  -GSM2_rDomain 0.42 \
  -GSM2_rNucleus 6.0 \
  -GSM2_alphaX 0.19 \
  -GSM2_betaX 0.05 \
  -GSM2_a 0.1 \
  -GSM2_b 0.1 \
  -GSM2_r 0.1 \
  -Doses 0 1 1

test -f "$OUTPUT_DIR/MKM_from_spectra_summary.csv"
test -f "$OUTPUT_DIR/MKM_from_spectra_manifest.txt"
test -f "$OUTPUT_DIR/DeCunha_GSM2.csv"

echo "Binned-spectrum GSM2 smoke outputs written to $OUTPUT_DIR"
