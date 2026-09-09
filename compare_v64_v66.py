#!/usr/bin/env python3
"""Compare short V6-4-1 and V6-6-1 correctness-run outputs."""

import argparse
import math
from pathlib import Path

import numpy as np


OUTPUT_FILES = (
    "bssn_BH.dat",
    "bssn_ADMQs.dat",
    "bssn_constraint.dat",
    "puncture_parameters_new.txt",
)


def locate(root: Path, name: str):
    candidates = (root / name, root / "binary_output" / name)
    for path in candidates:
        if path.is_file():
            return path
    return None


def load_numeric(path: Path):
    data = np.loadtxt(path, comments="#", ndmin=2)
    return np.asarray(data, dtype=float)


def error_stats(reference, candidate):
    delta = np.abs(candidate - reference)
    scale = np.maximum(np.abs(reference), 1.0e-300)
    relative = delta / scale
    return {
        "max_abs": float(np.max(delta)),
        "max_rel": float(np.max(relative)),
        "rms": float(math.sqrt(np.mean(np.square(delta)))),
        "mean_abs": float(np.mean(delta)),
        "mean_rel": float(np.mean(relative)),
    }


def compare_outputs(v64: Path, v66: Path):
    failed = False
    for name in OUTPUT_FILES:
        left = locate(v64, name)
        right = locate(v66, name)
        if left is None or right is None:
            print(f"{name}: skipped (V6-4-1={left}, V6-6-1={right})")
            continue
        a = load_numeric(left)
        b = load_numeric(right)
        if a.shape != b.shape:
            print(f"{name}: SHAPE MISMATCH {a.shape} vs {b.shape}")
            failed = True
            continue
        stats = error_stats(a, b)
        finite = np.isfinite(a).all() and np.isfinite(b).all()
        print(f"{name}: shape={a.shape} finite={finite} "
              f"max_abs={stats['max_abs']:.17e} "
              f"max_rel={stats['max_rel']:.17e} rms={stats['rms']:.17e}")
        if not finite:
            failed = True
    return failed


def summarize_samples(directory: Path):
    paths = sorted(directory.glob("batch17_verify_rank*.dat"))
    if not paths:
        print("batch17 samples: none found")
        return False
    rows = []
    for path in paths:
        with path.open(encoding="utf-8") as stream:
            for line in stream:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split()
                rows.append((int(fields[7]), fields[8], int(fields[6]),
                             float(fields[14]), float(fields[15])))
    if not rows:
        print("batch17 samples: files found, but no data rows")
        return True
    print(f"batch17 samples: {len(rows) // 17} points, "
          f"reflection_masks={sorted({row[2] for row in rows})}")
    for index in range(17):
        selected = [row for row in rows if row[0] == index]
        absolute = np.asarray([row[3] for row in selected])
        relative = np.asarray([row[4] for row in selected])
        print(f"  {index + 1:2d} {selected[0][1]:8s} n={len(selected):3d} "
              f"max_abs={absolute.max():.17e} mean_abs={absolute.mean():.17e} "
              f"max_rel={relative.max():.17e} mean_rel={relative.mean():.17e}")
    all_abs = np.asarray([row[3] for row in rows])
    all_rel = np.asarray([row[4] for row in rows])
    print(f"  ALL n={len(rows)} max_abs={all_abs.max():.17e} "
          f"mean_abs={all_abs.mean():.17e} max_rel={all_rel.max():.17e} "
          f"mean_rel={all_rel.mean():.17e}")
    if not np.isfinite(all_abs).all() or not np.isfinite(all_rel).all():
        print("  ALERT: NaN or Inf detected")
        return True
    if all_rel.max() >= 1.0e-6:
        print("  ALERT: sampled relative error is at least 1e-6")
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("v64", type=Path, help="V6-4-1 output directory")
    parser.add_argument("v66", type=Path, help="V6-6-1 output directory")
    parser.add_argument("--samples", type=Path,
                        help="directory containing batch17_verify_rank*.dat")
    args = parser.parse_args()
    failed = compare_outputs(args.v64, args.v66)
    if args.samples is not None:
        failed = summarize_samples(args.samples) or failed
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
