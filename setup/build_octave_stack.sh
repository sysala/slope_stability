#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

COMMON_OPT_FLAGS="${COMMON_OPT_FLAGS:--O3 -march=native -mtune=native}"
FORCE_REBUILD="${FORCE_REBUILD:-0}"
FORCE_SPARSERSB_REINSTALL="${FORCE_SPARSERSB_REINSTALL:-0}"

prepare_sources() {
  ensure_dirs

  download_if_missing "${OPENBLAS_TARBALL}" \
    "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${OPENBLAS_VERSION}/OpenBLAS-${OPENBLAS_VERSION}.tar.gz"

  download_if_missing "${OCTAVE_TARBALL}" \
    "https://ftp.gnu.org/gnu/octave/octave-${OCTAVE_VERSION}.tar.xz"

  download_if_missing "${LIBRSB_TARBALL}" \
    "https://downloads.sourceforge.net/project/librsb/librsb-${LIBRSB_VERSION}.tar.gz" \
    "https://sourceforge.net/projects/librsb/files/librsb-${LIBRSB_VERSION}.tar.gz/download" \
    "https://downloads.sourceforge.net/project/librsb/librsb/librsb-${LIBRSB_VERSION}.tar.gz" \
    "https://sourceforge.net/projects/librsb/files/librsb/librsb-${LIBRSB_VERSION}.tar.gz/download"

  download_if_missing "${SPARSERSB_TARBALL}" \
    "https://packages.octave.org/download/sparsersb-${SPARSERSB_VERSION}.tar.gz" \
    "https://downloads.sourceforge.net/project/octave/Octave%20Forge%20Packages/Individual%20Package%20Releases/sparsersb-${SPARSERSB_VERSION}.tar.gz" \
    "https://sourceforge.net/projects/octave/files/Octave%20Forge%20Packages/Individual%20Package%20Releases/sparsersb-${SPARSERSB_VERSION}.tar.gz/download" \
    "https://downloads.sourceforge.net/project/librsb/sparsersb/sparsersb-${SPARSERSB_VERSION}.tar.gz" \
    "https://sourceforge.net/projects/librsb/files/sparsersb/sparsersb-${SPARSERSB_VERSION}.tar.gz/download"
}

patch_librsb_for_modern_gcc() {
  local header="$1"
  if grep -q "RSB_PATCH_OMP_EXTERN_C" "${header}"; then
    return 0
  fi

  perl -0777 -i -pe '
    s|
      \#if\ RSB_WANT_OMP_RECURSIVE_KERNELS\s*\n
      \#include\ <omp\.h>\s*\n
      \#endif\ /\*\ RSB_WANT_OMP_RECURSIVE_KERNELS\ \*/\s*\n
    |
#ifdef __cplusplus
}
#endif /* __cplusplus */
#if RSB_WANT_OMP_RECURSIVE_KERNELS
#include <omp.h>
#endif /* RSB_WANT_OMP_RECURSIVE_KERNELS */
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
/* RSB_PATCH_OMP_EXTERN_C */
|xs;
  ' "${header}"

  if ! grep -q "RSB_PATCH_OMP_EXTERN_C" "${header}"; then
    fail "Failed to patch librsb header for omp.h C++ compatibility: ${header}"
  fi
}

patch_sparsersb_for_octave11() {
  local src_dir="$1"
  local cfg="${src_dir}/src/configure"
  local cc="${src_dir}/src/sparsersb.cc"

  if [[ ! -f "${cfg}" || ! -f "${cc}" ]]; then
    fail "Unexpected sparsersb source layout in ${src_dir}"
  fi

  if ! grep -q "SPARSERSB_PATCH_CXX11_OVERRIDE" "${cfg}"; then
    perl -0777 -i -pe '
      s{
        SPARSERSB_CXX11="-std=gnu\+\+11"
      }{
        if test "x\$SPARSERSB_CXX11" = "x"; then\nSPARSERSB_CXX11="-std=gnu++11"\nfi\n# SPARSERSB_PATCH_CXX11_OVERRIDE
      }x;
    ' "${cfg}"
  fi

  perl -i -pe 's/\.mex_get_ir\s*\(\s*\)/.ridx()/g; s/\.mex_get_jc\s*\(\s*\)/.cidx()/g;' "${cc}"
}

build_openblas() {
  if [[ "${FORCE_REBUILD}" != "1" ]] && [[ -f "${OPENBLAS_PREFIX}/lib/libopenblas.so" ]]; then
    log "Skipping OpenBLAS build (already installed at ${OPENBLAS_PREFIX})"
    return 0
  fi

  local src_unpack="${BUILD_DIR}/OpenBLAS-${OPENBLAS_VERSION}"
  extract_fresh "${OPENBLAS_TARBALL}" "${src_unpack}"

  log "Building OpenBLAS ${OPENBLAS_VERSION} (${OPENBLAS_TARGET})"
  make -C "${src_unpack}" clean >/dev/null 2>&1 || true
  make -C "${src_unpack}" -j"${NPROC}" \
    TARGET="${OPENBLAS_TARGET}" \
    DYNAMIC_ARCH=0 \
    USE_OPENMP=1 \
    NO_AFFINITY=1 \
    BINARY=64
  make -C "${src_unpack}" PREFIX="${OPENBLAS_PREFIX}" install
}

build_octave() {
  if [[ "${FORCE_REBUILD}" != "1" ]] && [[ -x "${OCTAVE_BIN}" ]]; then
    log "Skipping Octave build (already installed at ${OCTAVE_PREFIX})"
    return 0
  fi

  local src_unpack="${BUILD_DIR}/octave-${OCTAVE_VERSION}"
  local build_dir="${BUILD_DIR}/octave-${OCTAVE_VERSION}-build"
  extract_fresh "${OCTAVE_TARBALL}" "${src_unpack}"

  log "Building Octave ${OCTAVE_VERSION}"
  rm -rf "${build_dir}"
  mkdir -p "${build_dir}"

  pushd "${build_dir}" >/dev/null
  CPPFLAGS="-I${OPENBLAS_PREFIX}/include" \
  LDFLAGS="-L${OPENBLAS_PREFIX}/lib -Wl,-rpath,${OPENBLAS_PREFIX}/lib" \
  CFLAGS="${COMMON_OPT_FLAGS}" \
  CXXFLAGS="${COMMON_OPT_FLAGS}" \
  FFLAGS="${COMMON_OPT_FLAGS}" \
  FCFLAGS="${COMMON_OPT_FLAGS}" \
  "${src_unpack}/configure" \
    --prefix="${OCTAVE_PREFIX}" \
    --without-fltk \
    --disable-docs \
    --with-blas="-L${OPENBLAS_PREFIX}/lib -lopenblas" \
    --with-lapack="-L${OPENBLAS_PREFIX}/lib -lopenblas"
  make -j"${NPROC}"
  make install
  popd >/dev/null
}

build_librsb() {
  if [[ "${FORCE_REBUILD}" != "1" ]] && [[ -f "${LIBRSB_PREFIX}/lib/librsb.so.0" ]]; then
    log "Skipping librsb build (already installed at ${LIBRSB_PREFIX})"
    return 0
  fi

  local src_unpack="${BUILD_DIR}/librsb-${LIBRSB_VERSION}"
  extract_fresh "${LIBRSB_TARBALL}" "${src_unpack}"
  patch_librsb_for_modern_gcc "${src_unpack}/rsb_common.h"

  log "Building librsb ${LIBRSB_VERSION}"
  pushd "${src_unpack}" >/dev/null
  CFLAGS="${COMMON_OPT_FLAGS}" \
  CXXFLAGS="${COMMON_OPT_FLAGS}" \
  FFLAGS="${COMMON_OPT_FLAGS}" \
  FCFLAGS="${COMMON_OPT_FLAGS}" \
  ./configure --prefix="${LIBRSB_PREFIX}" --enable-openmp
  make -j"${NPROC}"
  make install
  popd >/dev/null
}

build_and_install_sparsersb() {
  local pkg_prefix="${OCTAVE_PREFIX}/share/octave/packages"
  local pkg_archprefix="${OCTAVE_PREFIX}/lib/octave/packages"
  local pkg_list="${OCTAVE_PREFIX}/share/octave/octave_packages"

  if [[ "${FORCE_REBUILD}" != "1" ]] && [[ "${FORCE_SPARSERSB_REINSTALL}" != "1" ]] && [[ -x "${OCTAVE_BIN}" ]]; then
    if export_runtime_env && "${OCTAVE_BIN}" --quiet --eval "\
        pkg('local_list','${pkg_list}'); \
        pkg('prefix','${pkg_prefix}','${pkg_archprefix}'); \
        pkg load sparsersb; disp('sparsersb-present');" >/dev/null 2>&1; then
      log "Skipping sparsersb install (package already loadable from local prefix)"
      return 0
    fi
  fi

  local src_unpack="${BUILD_DIR}/sparsersb-${SPARSERSB_VERSION}"
  extract_fresh "${SPARSERSB_TARBALL}" "${src_unpack}"
  patch_sparsersb_for_octave11 "${src_unpack}"

  log "Creating patched sparsersb source archive"
  tar -czf "${SPARSERSB_PATCHED_TARBALL}" -C "${BUILD_DIR}" "sparsersb-${SPARSERSB_VERSION}"

  log "Installing sparsersb into local Octave prefix (not user home)"
  export_runtime_env
  export SPARSERSB_CXX11="-std=gnu++17"

  # Set pkg prefix to keep packages inside the Octave install tree,
  # avoiding writes to ~/.local/share/octave/ or ~/.config/octave/.
  mkdir -p "${pkg_prefix}" "${pkg_archprefix}"

  "${OCTAVE_BIN}" --quiet --eval "\
    pkg ('prefix', '${pkg_prefix}', '${pkg_archprefix}'); \
    pkg ('local_list', '${pkg_list}'); \
    pkg install -verbose '${SPARSERSB_PATCHED_TARBALL}'; \
    pkg load sparsersb; \
    A = sparsersb([1;2],[1;2],[1;1],2,2); \
    disp(full(A));"
}

write_runtime_wrappers() {
  mkdir -p "${OCTAVE_ALL_DIR}/bin"

  # Write site-level octaverc so pkg system uses local paths (not ~/.local/)
  local site_startup="${OCTAVE_PREFIX}/share/octave/site/m/startup"
  mkdir -p "${site_startup}"
  cat > "${site_startup}/octaverc" <<'OCTAVERC_EOF'
## Site-level octaverc — redirect package paths into this Octave prefix
## (keeps everything local to the repo, nothing written to ~/.local or ~/.config)
_prefix = OCTAVE_HOME();
_pkg_dir = fullfile(_prefix, 'share', 'octave', 'packages');
_pkg_arch = fullfile(_prefix, 'lib', 'octave', 'packages');
_pkg_list = fullfile(_prefix, 'share', 'octave', 'octave_packages');
pkg('prefix', _pkg_dir, _pkg_arch);
pkg('local_list', _pkg_list);
clear _prefix _pkg_dir _pkg_arch _pkg_list;
OCTAVERC_EOF
  log "Wrote site octaverc: ${site_startup}/octaverc"

  cat > "${RUNTIME_ENV}" <<EOF
#!/usr/bin/env bash
export OPENBLAS_PREFIX="${OPENBLAS_PREFIX}"
export OCTAVE_PREFIX="${OCTAVE_PREFIX}"
export LIBRSB_PREFIX="${LIBRSB_PREFIX}"
export OCTAVE_BIN="${OCTAVE_BIN}"
export PATH="${OCTAVE_PREFIX}/bin:${LIBRSB_PREFIX}/bin:\${PATH}"
export LD_LIBRARY_PATH="${LIBRSB_PREFIX}/lib:${OPENBLAS_PREFIX}/lib:\${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="\${OMP_NUM_THREADS:-16}"
EOF
  chmod +x "${RUNTIME_ENV}"

  cat > "${LOCAL_WRAPPER}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export LD_LIBRARY_PATH="${LIBRSB_PREFIX}/lib:${OPENBLAS_PREFIX}/lib:\${LD_LIBRARY_PATH:-}"
exec "${OCTAVE_BIN}" --no-gui "\$@"
EOF
  chmod +x "${LOCAL_WRAPPER}"

  cat > "${ACTIVATE_SCRIPT}" <<EOF
#!/usr/bin/env bash
if [[ "\${BASH_SOURCE[0]}" == "\${0}" ]]; then
  echo "Run this script with: source ./activate_optimized_octave.sh" >&2
  exit 1
fi

source "${RUNTIME_ENV}"
echo "Activated optimized local Octave:"
echo "  OCTAVE_BIN=\${OCTAVE_BIN}"
echo "  OMP_NUM_THREADS=\${OMP_NUM_THREADS}"
EOF
  chmod +x "${ACTIVATE_SCRIPT}"
}

main() {
  require_cmd tar
  require_cmd make
  require_cmd perl
  require_cmd gcc
  require_cmd g++
  require_cmd gfortran

  prepare_sources
  build_openblas
  build_octave
  build_librsb
  build_and_install_sparsersb
  write_runtime_wrappers

  export_runtime_env
  log "Build complete."
  log "Octave binary: ${OCTAVE_BIN}"
  log "Octave wrapper: ${LOCAL_WRAPPER}"
}

main "$@"
