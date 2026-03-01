#!/usr/bin/env bash
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "Run this script with: source ./activate_optimized_octave.sh" >&2
  exit 1
fi

source "/home/beremi/repos/slope_stability/.octave_all/env.sh"
echo "Activated optimized local Octave:"
echo "  OCTAVE_BIN=${OCTAVE_BIN}"
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
