#!/bin/bash
# Build the element-level tangent assembly mex file for Octave.
# Run from the slope_stability/ directory.
#
# Usage:
#   bash mex/build_assemble_K_tangent.sh
#
# The compiled .mex file is placed in the current directory so that
# Octave can find it when running from slope_stability/.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC="$SCRIPT_DIR/assemble_K_tangent_vals.c"
OUT_DIR="$(pwd)"

echo "Compiling assemble_K_tangent_vals.c ..."

# Detect the right mkoctfile
if command -v mkoctfile &>/dev/null; then
    MKOCTFILE=mkoctfile
elif [ -x "${OCTAVE_BIN%/*}/mkoctfile" ]; then
    MKOCTFILE="${OCTAVE_BIN%/*}/mkoctfile"
else
    echo "ERROR: mkoctfile not found. Source activate_optimized_octave.sh first."
    exit 1
fi

$MKOCTFILE --mex -O2 \
    -DHAVE_OPENMP \
    "$SRC" \
    -o "$OUT_DIR/assemble_K_tangent_vals.mex" \
    -fopenmp -lgomp

echo "Done.  Output: $OUT_DIR/assemble_K_tangent_vals.mex"
