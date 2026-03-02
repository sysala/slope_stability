#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

THREADS="${THREADS:-16}"
PROFILE="${PROFILE:-medium}"
REPEATS="${REPEATS:-5}"
WARMUP="${WARMUP:-2}"
ND="${ND:-12}"

FULL_JSON="${RESULTS_DIR}/octave_btDb_prebuilt_line_breakdown_${PROFILE}_verify_t${THREADS}_full.json"
SYM_JSON="${RESULTS_DIR}/octave_btDb_prebuilt_line_breakdown_${PROFILE}_verify_t${THREADS}_sym.json"
STACK_JSON="${RESULTS_DIR}/octave_${PROFILE}_verify_sparsersb_t${THREADS}.json"

if [[ ! -x "${OCTAVE_BIN}" ]]; then
  fail "Missing Octave binary: ${OCTAVE_BIN}. Run build_octave_stack.sh first."
fi

require_cmd python3
export_runtime_env
export OMP_NUM_THREADS="${THREADS}"

log "Functional check: Octave + sparsersb"
"${OCTAVE_BIN}" --quiet --eval "pkg load sparsersb; A=sparsersb([1;2],[1;2],[2;3],2,2); disp(full(A));"

log "Timing prebuilt assembly (full constructor path)"
"${OCTAVE_BIN}" --quiet --eval "addpath('${BENCH_DIR}'); run_btDb_prebuilt_line_breakdown('profile','${PROFILE}','constructor','full','repeats',${REPEATS},'warmup',${WARMUP},'n_d',${ND},'out_file','${FULL_JSON}');"

log "Timing prebuilt assembly (sym_unique constructor path)"
"${OCTAVE_BIN}" --quiet --eval "addpath('${BENCH_DIR}'); run_btDb_prebuilt_line_breakdown('profile','${PROFILE}','constructor','sym_unique','repeats',${REPEATS},'warmup',${WARMUP},'n_d',${ND},'out_file','${SYM_JSON}');"

log "Timing core benchmark suite with sparsersb backend"
"${OCTAVE_BIN}" --quiet --eval "addpath('${BENCH_DIR}'); run_benchmarks('profile','${PROFILE}','repeats',3,'warmup',1,'sparse_backend','sparsersb','out_file','${STACK_JSON}');"

log "Checking performance envelope against known-good ranges"
python3 - <<PY
import json
import math
import pathlib
import sys

full_path = pathlib.Path("${FULL_JSON}")
sym_path = pathlib.Path("${SYM_JSON}")
stack_path = pathlib.Path("${STACK_JSON}")

full = json.loads(full_path.read_text())
sym = json.loads(sym_path.read_text())
stack = json.loads(stack_path.read_text())

full_total = float(full["total_median_s"])
sym_total = float(sym["total_median_s"])
speedup = full_total / sym_total if sym_total > 0 else math.inf

full_cons = float(full["median_totals_s"][-1])
sym_cons = float(sym["median_totals_s"][-1])
cons_speedup = full_cons / sym_cons if sym_cons > 0 else math.inf

cases = {c["kernel"]: float(c["median_s"]) for c in stack["cases"]}

print("Verification summary")
print(f"  full_total_s       = {full_total:.6f}")
print(f"  sym_total_s        = {sym_total:.6f}")
print(f"  total_speedup      = {speedup:.3f}x")
print(f"  full_constructor_s = {full_cons:.6f}")
print(f"  sym_constructor_s  = {sym_cons:.6f}")
print(f"  ctor_speedup       = {cons_speedup:.3f}x")
for k, v in cases.items():
    print(f"  {k:24s} = {v:.6f}")

checks = []
checks.append((full_total < 0.35, f"full_total_s too high: {full_total:.6f}"))
checks.append((sym_total < 0.20, f"sym_total_s too high: {sym_total:.6f}"))
checks.append((speedup > 1.70, f"total speedup too low: {speedup:.3f}x"))
checks.append((cons_speedup > 2.20, f"constructor speedup too low: {cons_speedup:.3f}x"))

# Envelope for medium profile from prior runs in this repo.
if "B' * D_p * B" in cases:
    bt = cases["B' * D_p * B"]
    checks.append((bt < 0.80, f"B' * D_p * B too slow: {bt:.6f}"))
if "A * x (sparse)" in cases:
    ax = cases["A * x (sparse)"]
    checks.append((ax < 0.01, f"A*x too slow: {ax:.6f}"))
if "X' * X (dense Gram)" in cases:
    xx = cases["X' * X (dense Gram)"]
    checks.append((xx < 0.02, f"X'X too slow: {xx:.6f}"))

errors = [msg for ok, msg in checks if not ok]
if errors:
    print("\nVerification FAILED:")
    for msg in errors:
        print(" -", msg)
    sys.exit(1)

print("\nVerification PASSED.")
PY

log "All checks passed."
