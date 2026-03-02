#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BENCH_DIR="${ROOT_DIR}/benchmark_octave"
OCTAVE_ALL_DIR="${ROOT_DIR}/.octave_all"
SRC_DIR="${OCTAVE_ALL_DIR}/src"
BUILD_DIR="${OCTAVE_ALL_DIR}/build"
INSTALL_DIR="${OCTAVE_ALL_DIR}/install"
RESULTS_DIR="${BENCH_DIR}/results"

OPENBLAS_VERSION="${OPENBLAS_VERSION:-0.3.31}"
OPENBLAS_TARGET="${OPENBLAS_TARGET:-ZEN}"
OCTAVE_VERSION="${OCTAVE_VERSION:-11.1.0}"
LIBRSB_VERSION="${LIBRSB_VERSION:-1.3.0.2}"
SPARSERSB_VERSION="${SPARSERSB_VERSION:-1.0.9}"

OPENBLAS_TARBALL="${SRC_DIR}/OpenBLAS-${OPENBLAS_VERSION}.tar.gz"
OCTAVE_TARBALL="${SRC_DIR}/octave-${OCTAVE_VERSION}.tar.xz"
LIBRSB_TARBALL="${SRC_DIR}/librsb-${LIBRSB_VERSION}.tar.gz"
SPARSERSB_TARBALL="${SRC_DIR}/sparsersb-${SPARSERSB_VERSION}.tar.gz"
SPARSERSB_PATCHED_TARBALL="${SRC_DIR}/sparsersb-${SPARSERSB_VERSION}-oct11-auto.tar.gz"

OPENBLAS_PREFIX="${INSTALL_DIR}/openblas-${OPENBLAS_VERSION}-zen"
OCTAVE_PREFIX="${INSTALL_DIR}/octave-${OCTAVE_VERSION}-zen"
LIBRSB_PREFIX="${INSTALL_DIR}/librsb-${LIBRSB_VERSION}"
OCTAVE_BIN="${OCTAVE_PREFIX}/bin/octave-cli"
LOCAL_WRAPPER="${OCTAVE_ALL_DIR}/bin/octave-rsb"
RUNTIME_ENV="${OCTAVE_ALL_DIR}/env.sh"
ACTIVATE_SCRIPT="${SCRIPT_DIR}/activate_optimized_octave.sh"
JUPYTER_VENV="${ROOT_DIR}/.venv"

# HYPRE paths (all inside .octave_all/ to keep things self-contained)
HYPRE_SRC_DIR="${HYPRE_SRC_DIR:-${SRC_DIR}/hypre}"
HYPRE_BUILD_DIR="${HYPRE_BUILD_DIR:-${BUILD_DIR}/hypre-openmp}"
HYPRE_INSTALL_DIR="${HYPRE_INSTALL_DIR:-${INSTALL_DIR}/hypre-openmp}"

if command -v getconf >/dev/null 2>&1; then
  NPROC_DEFAULT="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"
fi
if [[ -z "${NPROC_DEFAULT:-}" ]] && command -v nproc >/dev/null 2>&1; then
  NPROC_DEFAULT="$(nproc)"
fi
NPROC="${NPROC:-${NPROC_DEFAULT:-4}}"

log() {
  printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

fail() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    fail "Missing required command: ${cmd}"
  fi
}

ensure_dirs() {
  mkdir -p "${SRC_DIR}" "${BUILD_DIR}" "${INSTALL_DIR}" "${RESULTS_DIR}" "${OCTAVE_ALL_DIR}/bin"
}

download_if_missing() {
  local out="$1"
  shift

  if [[ -f "${out}" ]]; then
    log "Using existing source archive: ${out}"
    return 0
  fi

  local url
  for url in "$@"; do
    log "Downloading ${url}"
    if command -v curl >/dev/null 2>&1; then
      if curl -fL --retry 3 --retry-delay 2 -o "${out}" "${url}"; then
        return 0
      fi
    elif command -v wget >/dev/null 2>&1; then
      if wget -O "${out}" "${url}"; then
        return 0
      fi
    else
      fail "Neither curl nor wget is available for downloads."
    fi
    log "Download failed from ${url}, trying next source."
  done

  fail "Could not download required archive: ${out}"
}

extract_fresh() {
  local tarball="$1"
  local dst_dir="$2"
  local parent
  parent="$(dirname "${dst_dir}")"
  rm -rf "${dst_dir}"
  mkdir -p "${parent}"
  tar -xf "${tarball}" -C "${parent}"
}

export_runtime_env() {
  export PATH="${OCTAVE_PREFIX}/bin:${LIBRSB_PREFIX}/bin:${PATH}"
  export LD_LIBRARY_PATH="${LIBRSB_PREFIX}/lib:${OPENBLAS_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
}
