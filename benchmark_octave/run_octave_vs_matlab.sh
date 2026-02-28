#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
PROFILE="${1:-medium}"
REPEATS="${2:-3}"
WARMUP="${3:-1}"
DEFAULT_LOCAL_OCTAVE="${SCRIPT_DIR}/local/install/octave-11.1.0-zen/bin/octave-cli"

if [[ -n "${OCTAVE_BIN:-}" ]]; then
  OCT_BIN="${OCTAVE_BIN}"
elif command -v octave >/dev/null 2>&1; then
  OCT_BIN="$(command -v octave)"
elif [[ -x "${DEFAULT_LOCAL_OCTAVE}" ]]; then
  OCT_BIN="${DEFAULT_LOCAL_OCTAVE}"
else
  echo "Could not find Octave. Set OCTAVE_BIN or install Octave."
  exit 1
fi

MAT_BIN="${MATLAB_BIN:-matlab}"

mkdir -p "${RESULTS_DIR}"

OCT_JSON="${RESULTS_DIR}/octave_${PROFILE}.json"
MAT_JSON="${RESULTS_DIR}/matlab_${PROFILE}.json"

echo "Running Octave benchmark..."
oct_rc=0
if ! OMP_NUM_THREADS="${OCTAVE_OMP_NUM_THREADS:-${OMP_NUM_THREADS:-}}" \
  "${OCT_BIN}" --quiet --eval "addpath('${SCRIPT_DIR}'); run_benchmarks('profile','${PROFILE}','repeats',${REPEATS},'warmup',${WARMUP},'out_file','${OCT_JSON}');"; then
  oct_rc=$?
  echo "Warning: Octave benchmark exited with code ${oct_rc}."
fi

echo "Running MATLAB benchmark..."
mat_rc=0
if ! "${MAT_BIN}" -batch "addpath('${SCRIPT_DIR}'); run_benchmarks('profile','${PROFILE}','repeats',${REPEATS},'warmup',${WARMUP},'out_file','${MAT_JSON}');"; then
  mat_rc=$?
  echo "Warning: MATLAB benchmark exited with code ${mat_rc}."
fi

if [[ -f "${OCT_JSON}" && -f "${MAT_JSON}" ]]; then
  echo "Comparing results..."
  python3 "${SCRIPT_DIR}/compare_results.py" "${OCT_JSON}" "${MAT_JSON}"
else
  echo "Could not compare because one or both result files are missing."
  echo "Octave result present: $( [[ -f "${OCT_JSON}" ]] && echo yes || echo no )"
  echo "MATLAB result present: $( [[ -f "${MAT_JSON}" ]] && echo yes || echo no )"
  exit 1
fi
