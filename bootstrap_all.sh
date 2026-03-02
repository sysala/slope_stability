#!/usr/bin/env bash
# =============================================================================
# bootstrap_all.sh — One-call build for the entire Octave + HYPRE + MEX stack
# =============================================================================
#
# This is the single entry-point to build everything from scratch:
#   1. (Optional) Clean previous local builds
#   2. Build optimized Octave with OpenBLAS, librsb, and sparsersb
#   3. Build HYPRE (OpenMP, no MPI) and the hypre_boomeramg MEX binding
#   4. Build the constitutive-problem and assembly MEX files
#   5. Set up a Jupyter venv with the Octave kernel
#   6. (Optional) Run verification benchmarks
#
# Usage:
#   ./bootstrap_all.sh [--no-clean] [--no-verify] [--skip-hypre]
#
# After a successful run, activate the environment with:
#   source setup/activate_optimized_octave.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SETUP_DIR="${SCRIPT_DIR}/setup"

RUN_CLEAN=1
RUN_VERIFY=1
BUILD_HYPRE=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-clean)
      RUN_CLEAN=0
      shift
      ;;
    --no-verify)
      RUN_VERIFY=0
      shift
      ;;
    --skip-hypre)
      BUILD_HYPRE=0
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      echo "Usage: $0 [--no-clean] [--no-verify] [--skip-hypre]" >&2
      exit 1
      ;;
  esac
done

# --- Step 1: Clean (optional) -----------------------------------------------
if [[ "${RUN_CLEAN}" == "1" ]]; then
  "${SETUP_DIR}/clean_local_builds.sh"
fi

# --- Step 2: Build Octave stack (OpenBLAS + Octave + librsb + sparsersb) -----
"${SETUP_DIR}/build_octave_stack.sh"

# --- Step 3: Build HYPRE + MEX binding (optional) ---------------------------
if [[ "${BUILD_HYPRE}" == "1" ]]; then
  echo ""
  echo "=== Building HYPRE and hypre_boomeramg MEX ==="
  "${SETUP_DIR}/setup_hypre_mex.sh" --skip-mex
  # Build the Octave MEX via mkoctfile (the .m builder handles Octave too)
  source "${SETUP_DIR}/activate_optimized_octave.sh"
  cd "${SCRIPT_DIR}/slope_stability"
  octave-cli --quiet --eval "LINEAR_SOLVERS.build_hypre_boomeramg_mex();"
  cd "${SCRIPT_DIR}"
fi

# --- Step 4: Build constitutive-problem and assembly MEX files ---------------
echo ""
echo "=== Building constitutive and assembly MEX files ==="
source "${SETUP_DIR}/activate_optimized_octave.sh"
cd "${SCRIPT_DIR}/slope_stability"
bash "./+CONSTITUTIVE_PROBLEM/mex/build_constitutive_3D_mex.sh"
bash "./+ASSEMBLY/mex/build_assemble_K_tangent.sh"
cd "${SCRIPT_DIR}"

# --- Step 5: Set up Jupyter venv + Octave kernel ----------------------------
"${SETUP_DIR}/setup_jupyter_octave_venv.sh"

# --- Step 6: Verify (optional) ----------------------------------------------
if [[ "${RUN_VERIFY}" == "1" ]]; then
  THREADS="${THREADS:-16}" "${SETUP_DIR}/verify_stack.sh"
fi

echo ""
echo "========================================"
echo " Bootstrap complete."
echo ""
echo " Activate environment:"
echo "   source setup/activate_optimized_octave.sh"
echo ""
echo " Then run demos from slope_stability/:"
echo "   cd slope_stability"
echo "   octave-cli slope_stability_2D_homo_SSR.m"
echo "========================================"

