from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

INTERESTING = ["A", "V", "P", "G", "I", "L", "F", "M", "Y", "W", "K", "R", "D", "E", "H", "S", "T"]


def read_csv(path: Path) -> list[dict]:
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))


def print_subset(path: Path, residues: Iterable[str] = INTERESTING) -> None:
    residue_set = set(residues)
    rows = [r for r in read_csv(path) if r["residue"] in residue_set]
    rows.sort(key=lambda r: list(residue_set).index(r["residue"]) if r["residue"] in residue_set else 999)

    print(f"\n=== {path.name} ===")
    print(
        f"{'res':>3}  {'lambda_ref':>10}  {'ddu_h':>10}  {'kappa':>8}  "
        f"{'ddG_ex':>10}  {'dLambda':>10}  {'lambda_T':>10}"
    )

    for r in rows:
        print(
            f"{r['residue']:>3}  "
            f"{float(r['lambda_ref_300K']):10.6f}  "
            f"{float(r['delta_delta_mu_h_kcal_per_mol']):10.6f}  "
            f"{float(r['kappa']):8.4f}  "
            f"{float(r['delta_delta_g_ex_kcal_per_mol']):10.6f}  "
            f"{float(r['delta_lambda']):10.6f}  "
            f"{float(r['lambda_T']):10.6f}"
        )
