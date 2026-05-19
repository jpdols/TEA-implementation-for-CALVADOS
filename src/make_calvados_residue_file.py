from __future__ import annotations

import csv
from pathlib import Path


def read_csv_rows(path: Path) -> list[dict]:
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))


def write_csv_rows(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_lambda_lookup(tea_csv: Path) -> dict[str, float]:
    rows = read_csv_rows(tea_csv)
    return {row["residue"]: float(row["lambda_T"]) for row in rows}


def write_calvados_residues(tea_csv: Path, output_csv: Path, template_csv: Path) -> None:
    lambda_lookup = build_lambda_lookup(tea_csv)
    template_rows = read_csv_rows(template_csv)

    if not template_rows:
        raise ValueError(f"Template file is empty: {template_csv}")

    fieldnames = list(template_rows[0].keys())
    required = {"one", "three", "MW", "lambdas", "sigmas", "q", "bondlength"}
    missing = required - set(fieldnames)
    if missing:
        raise KeyError(
            f"Template is missing required columns: {sorted(missing)}. Found columns: {fieldnames}"
        )

    out_rows = []
    for row in template_rows:
        aa = row["one"]
        if aa not in lambda_lookup:
            raise KeyError(f"Residue {aa!r} from template not found in TEA table.")
        row_new = dict(row)
        row_new["lambdas"] = f"{lambda_lookup[aa]:.12f}"
        out_rows.append(row_new)

    write_csv_rows(output_csv, out_rows, fieldnames)
    print(f"Wrote CALVADOS residues file: {output_csv}")
