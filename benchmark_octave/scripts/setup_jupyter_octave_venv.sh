#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

KERNEL_NAME="${KERNEL_NAME:-octave-local-rsb}"
KERNEL_DISPLAY_NAME="${KERNEL_DISPLAY_NAME:-Octave (local-rsb)}"
DEFAULT_THREADS="${DEFAULT_THREADS:-16}"

require_cmd python3

if [[ ! -x "${LOCAL_WRAPPER}" ]]; then
  fail "Missing local octave wrapper: ${LOCAL_WRAPPER}. Run build_octave_stack.sh first."
fi

log "Creating Python venv at ${JUPYTER_VENV}"
python3 -m venv "${JUPYTER_VENV}"

log "Installing Jupyter + octave kernel packages"
"${JUPYTER_VENV}/bin/python" -m pip install --upgrade pip setuptools wheel
"${JUPYTER_VENV}/bin/python" -m pip install jupyterlab notebook ipykernel octave_kernel

KERNEL_DIR="${JUPYTER_VENV}/share/jupyter/kernels/${KERNEL_NAME}"
mkdir -p "${KERNEL_DIR}"

cat > "${KERNEL_DIR}/kernel.json" <<EOF
{
  "argv": [
    "${JUPYTER_VENV}/bin/python",
    "-m",
    "octave_kernel",
    "-f",
    "{connection_file}"
  ],
  "display_name": "${KERNEL_DISPLAY_NAME}",
  "language": "octave",
  "env": {
    "OCTAVE_EXECUTABLE": "${LOCAL_WRAPPER}",
    "OMP_NUM_THREADS": "${DEFAULT_THREADS}",
    "LD_LIBRARY_PATH": "${LIBRSB_PREFIX}/lib:${OPENBLAS_PREFIX}/lib"
  }
}
EOF

log "Kernel installed at ${KERNEL_DIR}"
log "Registering kernelspec into user Jupyter path"
"${JUPYTER_VENV}/bin/jupyter" kernelspec install --user --name "${KERNEL_NAME}" --replace "${KERNEL_DIR}"

log "Testing kernel registration"
"${JUPYTER_VENV}/bin/python" -m jupyter kernelspec list

log "Testing local octave wrapper"
"${LOCAL_WRAPPER}" --quiet --eval "pkg('list'); disp(version());"

log "Jupyter launch command:"
log "  source ${JUPYTER_VENV}/bin/activate && jupyter lab"
