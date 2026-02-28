#!/usr/bin/env bash
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "Run this script with: source ./activate_optimized_octave.sh" >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${ROOT_DIR}/.octave_all/env.sh"

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "Missing ${ENV_FILE}. Run ./build_octave_stack.sh first." >&2
  return 1
fi

source "${ENV_FILE}"
echo "Activated optimized local Octave:"
echo "  OCTAVE_BIN=${OCTAVE_BIN}"
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
