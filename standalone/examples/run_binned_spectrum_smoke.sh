#!/bin/sh
set -eu

MONAS_BIN=${MONAS_BIN:-../bin/monas}
OUTPUT_DIR=${MONAS_OUTPUT_DIR:-MKM_from_spectra_smoke_output}

"$MONAS_BIN" \
  -MKMFromSpectra \
  -outputDir "$OUTPUT_DIR" \
  -spectrum DeCunha:poly:poly_spectrum_fixture.csv \
  -spectrum AMF:amf:amf_spectra_fixture.csv:amf_moments_fixture.csv:0,0,0 \
  -spectrum TOPAS:topas:topas_yspec_fixture.txt \
  -MKMSatCorr -MKMnonPoiss -SMKM \
  -MKM_rDomain 0.44 \
  -MKM_rNucleus 8.0 \
  -MKM_alpha 0.19 \
  -MKM_beta 0.07 \
  -MKM_y0 150 \
  -MKM_alphaX 0.19 \
  -MKM_betaX 0.05 \
  -Doses 0 2 1

test -f "$OUTPUT_DIR/MKM_from_spectra_summary.csv"
test -f "$OUTPUT_DIR/MKM_from_spectra_manifest.txt"
test -f "$OUTPUT_DIR/DeCunha_MKM_SaturationCorrected.csv"
test -f "$OUTPUT_DIR/AMF_MKM_SaturationCorrected.csv"
test -f "$OUTPUT_DIR/TOPAS_MKM_SaturationCorrected.csv"

echo "Binned-spectrum smoke outputs written to $OUTPUT_DIR"
