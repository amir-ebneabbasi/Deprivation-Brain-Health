# make demo data

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

N_PCS = 10
N_OTHER_F = 20
SITE_VALUES = [1, 2, 3, 4]
F00_PREVALENCE = 0.08
MIN_F00_CASES = 60  # must exceed dep_main.py --min-cases (default 50)


def make_demo_data(n: int, seed: int, with_signal: bool) -> pd.DataFrame:
    """Create the demo analysis dataset."""
    rng = np.random.default_rng(seed)

    df = pd.DataFrame({"id_col": [f"sub-{i:05d}" for i in range(1, n + 1)]})

    # Ancestry PCs (random z-scores)
    for i in range(1, N_PCS + 1):
        df[f"PC{i}"] = rng.standard_normal(n)

    # Binary and categorical covariates
    df["site"] = rng.choice(SITE_VALUES, size=n)
    df["sex"] = rng.integers(0, 2, size=n)

    # Age (z-score) and derived terms
    df["age"] = rng.standard_normal(n)
    df["age2"] = df["age"] ** 2
    df["sex_age"] = df["sex"] * df["age"]
    df["age2_sex"] = df["age2"] * df["sex"]

    # Exposure, imaging covariates, mediator (random z-scores)
    df["IMD"] = rng.standard_normal(n)
    df["SurfaceHoles"] = rng.standard_normal(n)
    df["FD"] = rng.standard_normal(n)
    df["FD_max"] = rng.standard_normal(n)
    df["vol_bankssts"] = rng.standard_normal(n)

    if with_signal:
        # Higher deprivation -> lower volume
        df["vol_bankssts"] = -0.3 * df["IMD"] + 0.95 * rng.standard_normal(n)

    # F00: outcome with more than 50 cases
    if with_signal:
        base_logit = np.log(F00_PREVALENCE / (1 - F00_PREVALENCE))
        logit = base_logit + 0.15 * df["IMD"] - 0.4 * df["vol_bankssts"]
        p_f00 = 1 / (1 + np.exp(-logit))
        f00 = (rng.random(n) < p_f00).astype(int)
    else:
        f00 = (rng.random(n) < F00_PREVALENCE).astype(int)

    n_cases = int(f00.sum())
    if n_cases < MIN_F00_CASES:
        controls_idx = np.flatnonzero(f00 == 0)
        extra = rng.choice(controls_idx, size=MIN_F00_CASES - n_cases, replace=False)
        f00[extra] = 1
    df["F00"] = f00

    # 20 other random F diagnoses (names must match F## for dep_main.py)
    candidates = [f"F{i:02d}" for i in range(1, 100)]
    other_f = sorted(rng.choice(candidates, size=N_OTHER_F, replace=False))
    for code in other_f:
        prevalence = rng.uniform(0.02, 0.06)
        df[code] = (rng.random(n) < prevalence).astype(int)

    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Create demo data for dep_main.py")
    parser.add_argument("--out-dir", type=Path, default=Path("demo_data"))
    parser.add_argument("--n", type=int, default=2000, help="Number of participants")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--with-signal",
        action="store_true",
        help="Add a weak IMD -> vol_bankssts -> F00 relationship",
    )
    args, _ = parser.parse_known_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = make_demo_data(args.n, args.seed, args.with_signal)
    info = pd.DataFrame(
        {"X": ["IMD"], "Mediator": ["vol_bankssts"], "Y": ["F00"]}
    )

    data_path = args.out_dir / "Data_dep_brain_icd.csv"
    info_path = args.out_dir / "mediation_info.csv"
    df.to_csv(data_path, index=False)
    info.to_csv(info_path, index=False)

    # Summary using the same control definition as dep_main.py
    fg_cols = [c for c in df.columns if re.fullmatch(r"[FG]\d{2}", c)]
    n_controls = int(df[fg_cols].eq(0).all(axis=1).sum())
    print(f"Saved: {data_path}  ({df.shape[0]} rows x {df.shape[1]} columns)")
    print(f"Saved: {info_path}")
    print(f"F00 cases: {int(df['F00'].sum())} (need > 50)")
    print(f"F/G disease columns: {len(fg_cols)}; control pool: {n_controls}")
    print(f"Sites: {sorted(df['site'].unique().tolist())}")
    print()
    print("Run the mediation demo with:")
    print(
        f"  python dep_main.py --data-dir {args.out_dir} "
        "--n-bootstrap 200 --seed 42"
    )


if __name__ == "__main__":
    main()