from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Tuple


# Fallback defaults (original CALVADOS2)
DEFAULT_CALVADOS_LAMBDA_REF: Dict[str, float] = {
    "R": 0.730762476752,
    "D": 0.041604048061,
    "N": 0.425585900979,
    "E": 0.000693546096,
    "K": 0.179021173899,
    "H": 0.466366729056,
    "Q": 0.393431855106,
    "S": 0.462541681161,
    "C": 0.561543509914,
    "G": 0.705884373367,
    "T": 0.371316297627,
    "A": 0.274329796904,
    "M": 0.530848113434,
    "Y": 0.977461144934,
    "V": 0.208376960817,
    "W": 0.989376474037,
    "L": 0.644000500778,
    "I": 0.542362361067,
    "P": 0.359312657636,
    "F": 0.867235898206,
}

DEFAULT_CALVADOS_SIGMA_REF: Dict[str, float] = {
    "R": 0.656,
    "D": 0.558,
    "N": 0.568,
    "E": 0.592,
    "K": 0.636,
    "H": 0.608,
    "Q": 0.602,
    "S": 0.518,
    "C": 0.548,
    "G": 0.450,
    "T": 0.562,
    "A": 0.504,
    "M": 0.618,
    "Y": 0.646,
    "V": 0.586,
    "W": 0.678,
    "L": 0.618,
    "I": 0.618,
    "P": 0.556,
    "F": 0.636,
}

DEFAULT_CALVADOS_CHARGE_REF: Dict[str, int] = {
    "R": 1,
    "D": -1,
    "N": 0,
    "E": -1,
    "K": 1,
    "H": 0,
    "Q": 0,
    "S": 0,
    "C": 0,
    "G": 0,
    "T": 0,
    "A": 0,
    "M": 0,
    "Y": 0,
    "V": 0,
    "W": 0,
    "L": 0,
    "I": 0,
    "P": 0,
    "F": 0,
}


def _detect_lambda_column(fieldnames: list[str]) -> str:
    candidates = ["lambdas", "lambda_", "lambda"]
    for col in candidates:
        if col in fieldnames:
            return col
    raise KeyError(
        f"Could not find a lambda column. Expected one of {candidates}. "
        f"Found columns: {fieldnames}"
    )


def load_calvados_refs_from_template(
    template_csv: str | Path,
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, int]]:
    """
    Read residue lambda/sigma/charge reference values from a CALVADOS residue CSV.

    Expected columns:
      - one
      - lambdas (or lambda_ / lambda)
      - sigmas
      - q
    """
    path = Path(template_csv)
    if not path.exists():
        raise FileNotFoundError(f"Template residue file not found: {path}")

    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

        required_base = {"one", "sigmas", "q"}
        missing = required_base - set(fieldnames)
        if missing:
            raise KeyError(
                f"Template is missing required columns: {sorted(missing)}. "
                f"Found columns: {fieldnames}"
            )

        lambda_col = _detect_lambda_column(fieldnames)

        lambda_ref: Dict[str, float] = {}
        sigma_ref: Dict[str, float] = {}
        charge_ref: Dict[str, int] = {}

        for row in reader:
            aa = row["one"].strip().upper()
            if len(aa) != 1:
                raise ValueError(f"Invalid residue code in template: {row['one']!r}")

            lambda_ref[aa] = float(row[lambda_col])
            sigma_ref[aa] = float(row["sigmas"])
            charge_ref[aa] = int(float(row["q"]))

    return lambda_ref, sigma_ref, charge_ref


def get_default_refs() -> Tuple[Dict[str, float], Dict[str, float], Dict[str, int]]:
    return (
        dict(DEFAULT_CALVADOS_LAMBDA_REF),
        dict(DEFAULT_CALVADOS_SIGMA_REF),
        dict(DEFAULT_CALVADOS_CHARGE_REF),
    )


def get_refs(
    template_csv: str | Path | None = None,
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, int]]:
    """
    If template_csv is given, read refs from that file.
    Otherwise return the built-in CALVADOS2 defaults.
    """
    if template_csv is None:
        return get_default_refs()
    return load_calvados_refs_from_template(template_csv)


if __name__ == "__main__":
    # Change this path if you want to inspect a specific template file
    template = None

    lambda_ref, sigma_ref, charge_ref = get_refs(template)

    for aa in sorted(lambda_ref):
        print(
            f"{aa}: "
            f"lambda={lambda_ref[aa]:.12f}, "
            f"sigma={sigma_ref[aa]:.3f}, "
            f"q={charge_ref[aa]}"
        )