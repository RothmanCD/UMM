#!/usr/bin/env bash
# Build all six UMM figures as vector PDF + 300 dpi PNG into figures/.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PROJ="$(cd "$ROOT/.." && pwd)"
BUILD_AUX="$PROJ/build/figures"
mkdir -p "$BUILD_AUX"
cd "$ROOT"

FIGS=(
  fig4_consensus_independent
  fig5_dirac_spectrum
  fig6_dimensional_reduction
  fig2_rotation_curves
  fig3_bullet_cluster
  fig1_logistic_regimes
)

for f in "${FIGS[@]}"; do
  echo "=== Compiling $f ==="
  pdflatex -interaction=nonstopmode "$f.tex" > "${f}_build.log" 2>&1
  pdflatex -interaction=nonstopmode "$f.tex" >> "${f}_build.log" 2>&1
  if [[ ! -f "$f.pdf" ]]; then
    echo "FAIL: missing $f.pdf" >&2
    tail -30 "${f}_build.log" >&2
    exit 1
  fi
  gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 \
     -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
     -sOutputFile="$f.png" "$f.pdf" >> "${f}_build.log" 2>&1
  # Park aux/logs under build/figures; keep PDF/PNG + sources here
  for aux in "$f.aux" "$f.log" "${f}_build.log"; do
    [[ -f "$aux" ]] && mv -f "$aux" "$BUILD_AUX/"
  done
  echo "OK $f.pdf + $f.png"
done

echo "=== Inventory ==="
ls -la "${FIGS[@]/%/.pdf}" "${FIGS[@]/%/.png}"
echo "Done."
