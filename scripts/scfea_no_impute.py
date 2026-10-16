#!/usr/bin/env python3
"""Run untouched upstream scFEA with its optional MAGIC import disabled."""
from __future__ import annotations

import runpy
import sys
import types
from pathlib import Path

import pandas as pd


def main():
    if "--sc_imputation" not in sys.argv:
        raise SystemExit("scFEA wrapper requires --sc_imputation False")
    position = sys.argv.index("--sc_imputation")
    if position + 1 >= len(sys.argv) or sys.argv[position + 1] != "False":
        raise SystemExit("this wrapper only supports sc_imputation=False")
    source = Path(__file__).resolve().parents[1] / "upstream" / "scFEA-master" / "src" / "scFEA.py"
    if not source.is_file():
        raise SystemExit(f"upstream scFEA source missing: {source}")
    # Upstream imports magic at module load even when imputation is disabled.
    # The original algorithm never accesses this placeholder in this mode.
    sys.modules["magic"] = types.ModuleType("magic")
    if not hasattr(pd.DataFrame, "append"):
        # scFEA uses the pandas 1.x method removed in pandas 2.x.
        def append(self, other, ignore_index=False, sort=False):
            return pd.concat([self, other], ignore_index=ignore_index, sort=sort).infer_objects()
        pd.DataFrame.append = append
    sys.path.insert(0, str(source.parent))
    sys.argv[0] = str(source)
    runpy.run_path(str(source), run_name="__main__")


if __name__ == "__main__":
    main()
