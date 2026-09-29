#!/bin/sh
set -eu

MONAS_BIN=${MONAS_BIN:-../bin/monas}
OUTPUT_DIR=${MONAS_OUTPUT_DIR:-GSM2_from_spectra_smoke_output}

"$MONAS_BIN" \
  -MKMFromSpectra \
  -outputDir "$OUTPUT_DIR" \
  -spectrum DeCunha:poly:poly_spectrum_fixture.csv \
  -spectrum TOPAS:topas:topas_yspec_fixture.txt \
  -spectrum REAL:topas:ySpecfile.txt \
  -GSM2 \
  -H460 \
  -fSetMultiEventStatistic 5000 \
  -Doses 0 2 1

test -f "$OUTPUT_DIR/MKM_from_spectra_summary.csv"
test -f "$OUTPUT_DIR/MKM_from_spectra_manifest.txt"
test -f "$OUTPUT_DIR/DeCunha_GSM2.csv"
test -f "$OUTPUT_DIR/TOPAS_GSM2.csv"
test -f "$OUTPUT_DIR/REAL_GSM2.csv"
grep -q ",GSM2," "$OUTPUT_DIR/MKM_from_spectra_summary.csv"

echo "GSM2 binned-spectrum smoke outputs written to $OUTPUT_DIR"
