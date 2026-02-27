#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: setup_hypre_mex.sh [options]

One-call preparation:
1) Clone or update HYPRE sources.
2) Build/install HYPRE (OpenMP ON, MPI OFF).
3) Build MATLAB MEX binding for BoomerAMG.

Options:
  --ref <git-ref>        HYPRE git ref to checkout (default: master)
  --jobs <N>             Parallel build jobs (default: number of CPUs)
  --no-pull              Do not fetch/pull updates when source already exists
  --skip-mex             Build HYPRE only, skip MATLAB MEX build
  --matlab <command>     MATLAB executable (default: matlab)
  -h, --help             Show this help

Environment overrides:
  HYPRE_GIT_URL
  HYPRE_SRC_DIR
  HYPRE_BUILD_DIR
  HYPRE_INSTALL_DIR
  MATLAB_BIN
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HYPRE_GIT_URL="${HYPRE_GIT_URL:-https://github.com/hypre-space/hypre.git}"
HYPRE_REF="master"
HYPRE_SRC_DIR="${HYPRE_SRC_DIR:-$ROOT_DIR/third_party/hypre}"
HYPRE_BUILD_DIR="${HYPRE_BUILD_DIR:-$ROOT_DIR/third_party/hypre-build-openmp}"
HYPRE_INSTALL_DIR="${HYPRE_INSTALL_DIR:-$ROOT_DIR/third_party/hypre-openmp}"
MATLAB_BIN="${MATLAB_BIN:-matlab}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
DO_PULL=1
DO_MEX=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref)
      HYPRE_REF="${2:?missing value for --ref}"
      shift 2
      ;;
    --jobs)
      JOBS="${2:?missing value for --jobs}"
      shift 2
      ;;
    --no-pull)
      DO_PULL=0
      shift
      ;;
    --skip-mex)
      DO_MEX=0
      shift
      ;;
    --matlab)
      MATLAB_BIN="${2:?missing value for --matlab}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

echo "[1/3] Preparing HYPRE source"
if [[ ! -d "$HYPRE_SRC_DIR/.git" ]]; then
  mkdir -p "$(dirname "$HYPRE_SRC_DIR")"
  git clone "$HYPRE_GIT_URL" "$HYPRE_SRC_DIR"
fi

if [[ $DO_PULL -eq 1 ]]; then
  git -C "$HYPRE_SRC_DIR" fetch --tags --prune origin
  git -C "$HYPRE_SRC_DIR" checkout "$HYPRE_REF"
  CURRENT_BRANCH="$(git -C "$HYPRE_SRC_DIR" rev-parse --abbrev-ref HEAD)"
  if [[ "$CURRENT_BRANCH" != "HEAD" ]]; then
    git -C "$HYPRE_SRC_DIR" pull --ff-only origin "$CURRENT_BRANCH"
  fi
else
  git -C "$HYPRE_SRC_DIR" checkout "$HYPRE_REF"
fi

echo "[2/3] Configuring and building HYPRE"
cmake -S "$HYPRE_SRC_DIR/src" -B "$HYPRE_BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HYPRE_INSTALL_DIR" \
  -DHYPRE_ENABLE_MPI=OFF \
  -DHYPRE_ENABLE_OPENMP=ON \
  -DHYPRE_BUILD_EXAMPLES=OFF \
  -DHYPRE_BUILD_TESTS=OFF

cmake --build "$HYPRE_BUILD_DIR" --parallel "$JOBS"
cmake --install "$HYPRE_BUILD_DIR"

if [[ $DO_MEX -eq 1 ]]; then
  echo "[3/3] Building MATLAB MEX binding"
  "$MATLAB_BIN" -batch "cd('$ROOT_DIR/slope_stability'); LINEAR_SOLVERS.build_hypre_boomeramg_mex('$HYPRE_INSTALL_DIR');"
else
  echo "[3/3] Skipped MATLAB MEX build (--skip-mex)"
fi

echo
echo "Done."
echo "HYPRE install: $HYPRE_INSTALL_DIR"
if [[ $DO_MEX -eq 1 ]]; then
  echo "MEX source built via: slope_stability/+LINEAR_SOLVERS/build_hypre_boomeramg_mex.m"
fi
