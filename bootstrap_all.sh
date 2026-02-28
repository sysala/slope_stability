#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RUN_CLEAN=1
RUN_VERIFY=1

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
    *)
      echo "Unknown option: $1" >&2
      echo "Usage: $0 [--no-clean] [--no-verify]" >&2
      exit 1
      ;;
  esac
done

if [[ "${RUN_CLEAN}" == "1" ]]; then
  "${SCRIPT_DIR}/clean_local_builds.sh"
fi

"${SCRIPT_DIR}/build_octave_stack.sh"
"${SCRIPT_DIR}/setup_jupyter_octave_venv.sh"

if [[ "${RUN_VERIFY}" == "1" ]]; then
  THREADS="${THREADS:-16}" "${SCRIPT_DIR}/verify_stack.sh"
fi

echo "Bootstrap complete."

