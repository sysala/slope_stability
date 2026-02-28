#!/usr/bin/env python3
import json
import sys
from pathlib import Path


def load(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def case_map(data):
    return {item["kernel"]: item for item in data["cases"]}


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: compare_results.py <octave.json> <matlab.json>", file=sys.stderr)
        return 2

    oct_path = Path(sys.argv[1])
    mat_path = Path(sys.argv[2])
    oct_data = load(oct_path)
    mat_data = load(mat_path)

    oct_cases = case_map(oct_data)
    mat_cases = case_map(mat_data)
    kernels = [
        "B' * D_p * B",
        "A * x (sparse)",
        "X' * X (dense Gram)",
    ]

    print("\n=== Octave vs MATLAB (median seconds) ===")
    print(f"{'Kernel':24s} {'Octave':>10s} {'MATLAB':>10s} {'Oct/MAT':>10s}")
    for kernel in kernels:
        o = oct_cases[kernel]["median_s"]
        m = mat_cases[kernel]["median_s"]
        ratio = o / m if m != 0 else float("inf")
        print(f"{kernel:24s} {o:10.6f} {m:10.6f} {ratio:10.3f}")

    print(f"\nOctave result file: {oct_path}")
    print(f"MATLAB result file: {mat_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
