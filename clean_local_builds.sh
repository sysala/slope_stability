#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

PURGE_SOURCES="${PURGE_SOURCES:-0}"
PURGE_RESULTS="${PURGE_RESULTS:-0}"

log "Cleaning local build artifacts under ${OCTAVE_ALL_DIR}"
rm -rf "${BUILD_DIR}" "${INSTALL_DIR}" "${OCTAVE_ALL_DIR}/bin" "${JUPYTER_VENV}" "${RUNTIME_ENV}"

if [[ "${PURGE_SOURCES}" == "1" ]]; then
  log "PURGE_SOURCES=1 -> removing downloaded source archives"
  rm -rf "${SRC_DIR}"
fi

if [[ "${PURGE_RESULTS}" == "1" ]]; then
  log "PURGE_RESULTS=1 -> removing benchmark result JSON files"
  rm -rf "${RESULTS_DIR}"
fi

ensure_dirs
log "Cleanup complete."
