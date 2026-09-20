# run demo

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Run dep_main.py on the demo data")
    parser.add_argument("--data-dir", type=Path, default=Path("demo_data"))
    parser.add_argument("--dep-main", type=Path, default=Path("dep_main.py"))
    parser.add_argument("--n-bootstrap", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    args, _ = parser.parse_known_args()

    for path in (
        args.dep_main,
        args.data_dir / "mediation_info.csv",
        args.data_dir / "Data_dep_brain_icd.csv",
    ):
        if not path.is_file():
            raise FileNotFoundError(
                f"Not found: {path}\n"
                "Run make_demo.py first, and run this script from the folder "
                "that contains dep_main.py (or pass --dep-main /path/to/dep_main.py)."
            )

    command = [
        sys.executable,
        str(args.dep_main),
        "--data-dir", str(args.data_dir),
        "--n-bootstrap", str(args.n_bootstrap),
        "--seed", str(args.seed),
    ]
    print("Running:", " ".join(command), "\n")
    proc = subprocess.run(command, capture_output=True, text=True)
    if proc.stdout:
        print(proc.stdout.rstrip())
    if proc.returncode != 0:
        print(proc.stderr)
        raise RuntimeError("dep_main.py failed")

    results_path = args.data_dir / "results_mediation_chunk_0.csv"
    results = pd.read_csv(results_path)

    print("\n=== Mediation result ===")
    # one model -> print as a vertical table
    print(results.iloc[0].to_string())
    print(f"\nSaved: {results_path}")


if __name__ == "__main__":
    main()