#!/usr/bin/env python3
"""
Bootstrap mediation analysis for binary clinical outcomes.

This script estimates:

    a: X → Mediator
    b: Mediator → Y, adjusted for X
    direct: X → Y, adjusted for Mediator
    indirect: a × b

The mediator model is fitted using ordinary least squares, while the binary
outcome model is fitted using binomial logistic regression.

The script supports SLURM array jobs by processing a separate chunk of the
mediation specification file for each SLURM_ARRAY_TASK_ID.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm


DEFAULT_DATA_DIR = Path("path/to/working/dir")
DEFAULT_INFO_FILENAME = "mediation_info.csv"
DEFAULT_DATA_FILENAME = "Data_dep_brain_icd.csv"
DEFAULT_CHUNK_SIZE = 10
DEFAULT_N_BOOTSTRAP = 5000
DEFAULT_MIN_CASES = 50

RESULT_COLUMNS = [
    "X",
    "Mediator",
    "Y",
    "a",
    "b",
    "direct",
    "indirect",
    "a_ci_low",
    "a_ci_high",
    "a_p",
    "b_ci_low",
    "b_ci_high",
    "b_p",
    "direct_ci_low",
    "direct_ci_high",
    "direct_p",
    "indirect_ci_low",
    "indirect_ci_high",
    "indirect_p",
    "N",
    "N_CASES",
    "N_CONTROLS",
    "error",
]


# =============================================================================
# Configuration
# =============================================================================
def get_imaging_covariates() -> list[str]:
    """Return the covariates used in the imaging mediation models."""
    return [
        "PC1",
        "PC2",
        "PC3",
        "PC4",
        "PC5",
        "PC6",
        "PC7",
        "PC8",
        "PC9",
        "PC10",
        "site",
        "sex",
        "age",
        "age2",
        "sex_age",
        "age2_sex",
        "SurfaceHoles",
    ]


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run bootstrap mediation analyses in SLURM-compatible chunks."
    )

    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Directory containing the input CSV files.",
    )
    parser.add_argument(
        "--info-file",
        default=DEFAULT_INFO_FILENAME,
        help="CSV file containing X, Mediator, and Y columns.",
    )
    parser.add_argument(
        "--data-file",
        default=DEFAULT_DATA_FILENAME,
        help="CSV file containing the analysis dataset.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to --data-dir.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Number of mediation models processed by each array task.",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=DEFAULT_N_BOOTSTRAP,
        help="Number of bootstrap repetitions.",
    )
    parser.add_argument(
        "--min-cases",
        type=int,
        default=DEFAULT_MIN_CASES,
        help="Minimum required number of cases.",
    )
    parser.add_argument(
        "--task-id",
        type=int,
        default=None,
        help=(
            "Chunk task ID. When omitted, SLURM_ARRAY_TASK_ID is used. "
            "If unavailable, task ID 0 is used."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for reproducible bootstrap sampling.",
    )

    return parser.parse_args()


# =============================================================================
# Input and validation
# =============================================================================
def load_data(
    data_dir: Path,
    info_filename: str,
    data_filename: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the mediation specification and main analysis datasets."""
    info_path = data_dir / info_filename
    data_path = data_dir / data_filename

    if not info_path.is_file():
        raise FileNotFoundError(
            f"Mediation information file not found: {info_path}"
        )

    if not data_path.is_file():
        raise FileNotFoundError(
            f"Main data file not found: {data_path}"
        )

    print(f"Loading mediation information: {info_path}")
    info = pd.read_csv(info_path)

    print(f"Loading analysis data: {data_path}")
    df_main = pd.read_csv(data_path)

    validate_info_dataframe(info)

    return info, df_main


def validate_info_dataframe(info: pd.DataFrame) -> None:
    """Validate the mediation specification dataframe."""
    required_columns = {"X", "Mediator", "Y"}
    missing_columns = required_columns.difference(info.columns)

    if missing_columns:
        raise ValueError(
            "The mediation information file is missing required columns: "
            f"{sorted(missing_columns)}"
        )


def validate_model_columns(
    data: pd.DataFrame,
    x: str,
    mediator: str,
    y: str,
    covariates: Sequence[str],
) -> None:
    """Ensure that all required variables exist in the dataset."""
    required_columns = {x, mediator, y, *covariates}
    missing_columns = required_columns.difference(data.columns)

    if missing_columns:
        raise KeyError(
            "Required analysis columns are missing: "
            f"{sorted(missing_columns)}"
        )


# =============================================================================
# Control definition
# =============================================================================
def identify_fg_columns(data: pd.DataFrame) -> list[str]:
    """Identify disease indicator columns matching F## or G##."""
    return [
        column
        for column in data.columns
        if re.fullmatch(r"[FG]\d{2}", str(column))
    ]


def create_control_pool(data: pd.DataFrame) -> pd.DataFrame:
    """Create the control pool using all F## and G## disease variables."""
    fg_columns = identify_fg_columns(data)

    if not fg_columns:
        raise ValueError(
            "No disease indicator columns matching F## or G## were found."
        )

    control_mask = data[fg_columns].eq(0).all(axis=1)
    controls = data.loc[control_mask].copy()

    print(
        f"Control pool: {len(controls)} participants "
        f"using {len(fg_columns)} disease columns"
    )

    return controls


# =============================================================================
# Regression
# =============================================================================
def fit_mediation_models(
    data: pd.DataFrame,
    x: str,
    mediator: str,
    y: str,
    covariates: Sequence[str],
) -> tuple[float, float, float, float]:
    """Fit the mediator and binary outcome regression models."""
    mediator_predictors = sm.add_constant(
        data[[x, *covariates]],
        has_constant="add",
    )

    mediator_model = sm.OLS(
        data[mediator],
        mediator_predictors,
    ).fit()

    a_path = float(mediator_model.params[x])

    outcome_predictors = sm.add_constant(
        data[[x, mediator, *covariates]],
        has_constant="add",
    )

    outcome_model = sm.GLM(
        data[y],
        outcome_predictors,
        family=sm.families.Binomial(),
    ).fit()

    b_path = float(outcome_model.params[mediator])
    direct_effect = float(outcome_model.params[x])
    indirect_effect = a_path * b_path

    return a_path, b_path, direct_effect, indirect_effect


# =============================================================================
# Bootstrap statistics
# =============================================================================
def calculate_percentile_ci(
    distribution: Sequence[float],
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Calculate a percentile-based confidence interval."""
    values = np.asarray(distribution, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.nan, np.nan

    lower, upper = np.percentile(
        values,
        [100 * alpha / 2, 100 * (1 - alpha / 2)],
    )

    return float(lower), float(upper)


def calculate_bootstrap_pvalue(
    distribution: Sequence[float],
) -> float:
    """Calculate a two-sided sign-based bootstrap p-value."""
    values = np.asarray(distribution, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.nan

    n_positive = np.sum(values > 0)
    n_negative = np.sum(values < 0)
    n_nonzero = n_positive + n_negative

    if n_nonzero == 0:
        return np.nan

    return float(2 * min(n_positive, n_negative) / n_nonzero)


def summarize_distribution(
    distribution: Sequence[float],
) -> tuple[float, float, float, float]:
    """Return mean, lower CI, upper CI, and bootstrap p-value."""
    values = np.asarray(distribution, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.nan, np.nan, np.nan, np.nan

    ci_low, ci_high = calculate_percentile_ci(values)
    p_value = calculate_bootstrap_pvalue(values)

    return (
        float(np.mean(values)),
        ci_low,
        ci_high,
        p_value,
    )


# =============================================================================
# Result construction
# =============================================================================
def create_error_result(
    x: str,
    mediator: str,
    y: str,
    *,
    n: float = np.nan,
    n_cases: float = np.nan,
    n_controls: float = np.nan,
    error: str | float = np.nan,
) -> pd.DataFrame:
    """Create a one-row result dataframe for a failed analysis."""
    result = {
        "X": x,
        "Mediator": mediator,
        "Y": y,
        "a": np.nan,
        "b": np.nan,
        "direct": np.nan,
        "indirect": np.nan,
        "a_ci_low": np.nan,
        "a_ci_high": np.nan,
        "a_p": np.nan,
        "b_ci_low": np.nan,
        "b_ci_high": np.nan,
        "b_p": np.nan,
        "direct_ci_low": np.nan,
        "direct_ci_high": np.nan,
        "direct_p": np.nan,
        "indirect_ci_low": np.nan,
        "indirect_ci_high": np.nan,
        "indirect_p": np.nan,
        "N": n,
        "N_CASES": n_cases,
        "N_CONTROLS": n_controls,
        "error": error,
    }

    return pd.DataFrame([result], columns=RESULT_COLUMNS)


def create_success_result(
    x: str,
    mediator: str,
    y: str,
    a_distribution: Sequence[float],
    b_distribution: Sequence[float],
    direct_distribution: Sequence[float],
    indirect_distribution: Sequence[float],
    *,
    n: int,
    n_cases: int,
    n_controls: int,
) -> pd.DataFrame:
    """Summarize bootstrap distributions in a one-row dataframe."""
    a_mean, a_low, a_high, a_p = summarize_distribution(a_distribution)
    b_mean, b_low, b_high, b_p = summarize_distribution(b_distribution)

    direct_mean, direct_low, direct_high, direct_p = summarize_distribution(
        direct_distribution
    )

    indirect_mean, indirect_low, indirect_high, indirect_p = (
        summarize_distribution(indirect_distribution)
    )

    result = {
        "X": x,
        "Mediator": mediator,
        "Y": y,
        "a": a_mean,
        "b": b_mean,
        "direct": direct_mean,
        "indirect": indirect_mean,
        "a_ci_low": a_low,
        "a_ci_high": a_high,
        "a_p": a_p,
        "b_ci_low": b_low,
        "b_ci_high": b_high,
        "b_p": b_p,
        "direct_ci_low": direct_low,
        "direct_ci_high": direct_high,
        "direct_p": direct_p,
        "indirect_ci_low": indirect_low,
        "indirect_ci_high": indirect_high,
        "indirect_p": indirect_p,
        "N": n,
        "N_CASES": n_cases,
        "N_CONTROLS": n_controls,
        "error": np.nan,
    }

    return pd.DataFrame([result], columns=RESULT_COLUMNS)


# =============================================================================
# Bootstrap analysis
# =============================================================================
def run_bootstrap_mediation(
    df_main: pd.DataFrame,
    controls_pool: pd.DataFrame,
    x: str,
    mediator: str,
    y: str,
    covariates: Sequence[str],
    n_bootstrap: int,
    min_cases: int = DEFAULT_MIN_CASES,
    random_seed: int | None = None,
) -> pd.DataFrame:
    """Run bootstrap mediation analysis for one variable combination."""
    n_cases = 0
    n_controls = len(controls_pool)
    total_n = n_controls

    try:
        validate_model_columns(
            data=df_main,
            x=x,
            mediator=mediator,
            y=y,
            covariates=covariates,
        )

        cases = df_main.loc[
            df_main[y].eq(1) & df_main[x].notna()
        ].copy()

        n_cases = len(cases)
        total_n = n_cases + n_controls

        if n_cases <= min_cases:
            return create_error_result(
                x=x,
                mediator=mediator,
                y=y,
                n=total_n,
                n_cases=n_cases,
                n_controls=n_controls,
                error=f"too_few_cases: required > {min_cases}",
            )

        required_columns = list(
            dict.fromkeys([x, mediator, y, *covariates])
        )

        rng = np.random.default_rng(random_seed)

        a_distribution = np.full(n_bootstrap, np.nan)
        b_distribution = np.full(n_bootstrap, np.nan)
        direct_distribution = np.full(n_bootstrap, np.nan)
        indirect_distribution = np.full(n_bootstrap, np.nan)

        for bootstrap_index in range(n_bootstrap):
            case_seed = int(rng.integers(0, np.iinfo(np.int32).max))
            control_seed = int(rng.integers(0, np.iinfo(np.int32).max))

            sampled_cases = cases.sample(
                n=n_cases,
                replace=True,
                random_state=case_seed,
            )

            sampled_controls = controls_pool.sample(
                n=n_controls,
                replace=True,
                random_state=control_seed,
            )

            sampled_data = pd.concat(
                [sampled_cases, sampled_controls],
                ignore_index=True,
            )

            complete_data = sampled_data[required_columns].dropna()

            if complete_data.empty or complete_data[y].nunique() < 2:
                continue

            try:
                a_path, b_path, direct_effect, indirect_effect = (
                    fit_mediation_models(
                        data=complete_data,
                        x=x,
                        mediator=mediator,
                        y=y,
                        covariates=covariates,
                    )
                )
            except Exception:
                continue

            a_distribution[bootstrap_index] = a_path
            b_distribution[bootstrap_index] = b_path
            direct_distribution[bootstrap_index] = direct_effect
            indirect_distribution[bootstrap_index] = indirect_effect

        successful_bootstraps = np.isfinite(indirect_distribution).sum()

        if successful_bootstraps == 0:
            return create_error_result(
                x=x,
                mediator=mediator,
                y=y,
                n=total_n,
                n_cases=n_cases,
                n_controls=n_controls,
                error="all_bootstrap_models_failed",
            )

        print(
            f"{x} -> {mediator} -> {y}: "
            f"{successful_bootstraps}/{n_bootstrap} successful bootstraps"
        )

        return create_success_result(
            x=x,
            mediator=mediator,
            y=y,
            a_distribution=a_distribution,
            b_distribution=b_distribution,
            direct_distribution=direct_distribution,
            indirect_distribution=indirect_distribution,
            n=total_n,
            n_cases=n_cases,
            n_controls=n_controls,
        )

    except Exception as error:
        return create_error_result(
            x=x,
            mediator=mediator,
            y=y,
            n=total_n,
            n_cases=n_cases,
            n_controls=n_controls,
            error=f"{type(error).__name__}: {error}",
        )


# =============================================================================
# SLURM chunking
# =============================================================================
def get_task_id(task_id: int | None = None) -> int:
    """Return the supplied task ID or the SLURM array task ID."""
    if task_id is not None:
        resolved_task_id = task_id
    else:
        resolved_task_id = int(
            os.environ.get("SLURM_ARRAY_TASK_ID", "0")
        )

    if resolved_task_id < 0:
        raise ValueError("Task ID must be zero or greater.")

    return resolved_task_id


def get_analysis_chunk(
    info: pd.DataFrame,
    task_id: int,
    chunk_size: int,
) -> tuple[pd.DataFrame, int, int]:
    """Return the mediation specification rows assigned to one task."""
    if chunk_size <= 0:
        raise ValueError("Chunk size must be greater than zero.")

    start_index = task_id * chunk_size
    end_index = min(start_index + chunk_size, len(info))

    info_chunk = info.iloc[start_index:end_index].copy()

    print(
        f"Task {task_id}: processing rows "
        f"{start_index} to {end_index}"
    )

    if info_chunk.empty:
        print(
            f"No rows assigned to task {task_id}. "
            f"The information file contains {len(info)} rows."
        )

    return info_chunk, start_index, end_index


# =============================================================================
# Analysis workflow
# =============================================================================
def run_analysis_chunk(
    info_chunk: pd.DataFrame,
    df_main: pd.DataFrame,
    controls_pool: pd.DataFrame,
    covariates: Sequence[str],
    n_bootstrap: int,
    min_cases: int,
    random_seed: int | None = None,
) -> pd.DataFrame:
    """Run all mediation analyses assigned to a chunk."""
    results: list[pd.DataFrame] = []

    for model_number, (_, row) in enumerate(info_chunk.iterrows()):
        x = str(row["X"])
        mediator = str(row["Mediator"])
        y = str(row["Y"])

        print(
            f"Running model {model_number + 1}/{len(info_chunk)}: "
            f"{x} -> {mediator} -> {y}"
        )

        model_seed = (
            None
            if random_seed is None
            else random_seed + model_number
        )

        result = run_bootstrap_mediation(
            df_main=df_main,
            controls_pool=controls_pool,
            x=x,
            mediator=mediator,
            y=y,
            covariates=covariates,
            n_bootstrap=n_bootstrap,
            min_cases=min_cases,
            random_seed=model_seed,
        )

        results.append(result)

    if not results:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    return pd.concat(results, ignore_index=True)


def save_results(
    results: pd.DataFrame,
    output_dir: Path,
    task_id: int,
) -> Path:
    """Save one task's results to a CSV file."""
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"results_mediation_chunk_{task_id}.csv"
    results.to_csv(output_path, index=False)

    print(f"Saved results: {output_path}")

    return output_path


def run_pipeline(
    data_dir: Path,
    info_filename: str,
    data_filename: str,
    output_dir: Path,
    chunk_size: int,
    n_bootstrap: int,
    min_cases: int,
    task_id: int | None = None,
    random_seed: int | None = None,
) -> Path:
    """Run the complete mediation-analysis pipeline."""
    resolved_task_id = get_task_id(task_id)

    info, df_main = load_data(
        data_dir=data_dir,
        info_filename=info_filename,
        data_filename=data_filename,
    )

    covariates = get_imaging_covariates()

    missing_covariates = set(covariates).difference(df_main.columns)

    if missing_covariates:
        raise ValueError(
            "The main dataset is missing covariates: "
            f"{sorted(missing_covariates)}"
        )

    controls_pool = create_control_pool(df_main)

    info_chunk, _, _ = get_analysis_chunk(
        info=info,
        task_id=resolved_task_id,
        chunk_size=chunk_size,
    )

    results = run_analysis_chunk(
        info_chunk=info_chunk,
        df_main=df_main,
        controls_pool=controls_pool,
        covariates=covariates,
        n_bootstrap=n_bootstrap,
        min_cases=min_cases,
        random_seed=random_seed,
    )

    return save_results(
        results=results,
        output_dir=output_dir,
        task_id=resolved_task_id,
    )


def main() -> None:
    """Command-line entry point."""
    args = parse_arguments()

    if args.n_bootstrap <= 0:
        raise ValueError("--n-bootstrap must be greater than zero.")

    output_dir = args.output_dir or args.data_dir

    run_pipeline(
        data_dir=args.data_dir,
        info_filename=args.info_file,
        data_filename=args.data_file,
        output_dir=output_dir,
        chunk_size=args.chunk_size,
        n_bootstrap=args.n_bootstrap,
        min_cases=args.min_cases,
        task_id=args.task_id,
        random_seed=args.seed,
    )


if __name__ == "__main__":
    main()
