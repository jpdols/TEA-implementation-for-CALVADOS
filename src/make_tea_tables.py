from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable

from calvados_lambda_ref import get_default_refs
from tea_lambda import DEFAULT_LAMBDA_MAPPING_SIGN, lambda_table_for_residues
from tea_thermo import T_REF

DEFAULT_TEMPERATURES = [293.0, 300.0, 310.0, 320.0]
DEFAULT_GAMMAS = [2.0, 2.5, 3.0]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_lambda_table_csv(
    output_csv: Path,
    temperature_k: float,
    gamma: float,
    lambda_ref_table: Dict[str, float],
    sigma_ref_table: Dict[str, float],
    charge_ref_table: Dict[str, int],
    use_arg_charmm36: bool = False,
    t_ref: float = T_REF,
    mapping_sign: float = DEFAULT_LAMBDA_MAPPING_SIGN,
) -> None:
    table = lambda_table_for_residues(
        lambda_ref_table=lambda_ref_table,
        temperature_k=temperature_k,
        gamma=gamma,
        t_ref=t_ref,
        use_arg_charmm36=use_arg_charmm36,
        mapping_sign=mapping_sign,
    )

    fieldnames = [
        "residue",
        "temperature_K",
        "gamma",
        "mapping_sign",
        "use_arg_charmm36",
        "lambda_ref_300K",
        "sigma_ref",
        "charge_ref",
        "delta_delta_mu_h_kcal_per_mol",
        "kappa",
        "delta_delta_g_ex_kcal_per_mol",
        "delta_lambda",
        "lambda_T",
    ]

    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for aa in sorted(table):
            vals = table[aa]
            writer.writerow(
                {
                    "residue": aa,
                    "temperature_K": float(temperature_k),
                    "gamma": float(gamma),
                    "mapping_sign": float(mapping_sign),
                    "use_arg_charmm36": bool(use_arg_charmm36),
                    "lambda_ref_300K": lambda_ref_table[aa],
                    "sigma_ref": sigma_ref_table[aa],
                    "charge_ref": charge_ref_table[aa],
                    "delta_delta_mu_h_kcal_per_mol": vals["delta_delta_mu_h"],
                    "kappa": vals["kappa"],
                    "delta_delta_g_ex_kcal_per_mol": vals["delta_delta_g_ex"],
                    "delta_lambda": vals["delta_lambda"],
                    "lambda_T": vals["lambda_t"],
                }
            )


def write_all_tables(
    output_dir: Path,
    temperatures_k: Iterable[float] = DEFAULT_TEMPERATURES,
    gammas: Iterable[float] = DEFAULT_GAMMAS,
    lambda_ref_table: Dict[str, float] | None = None,
    sigma_ref_table: Dict[str, float] | None = None,
    charge_ref_table: Dict[str, int] | None = None,
    use_arg_charmm36: bool = False,
    t_ref: float = T_REF,
    mapping_sign: float = DEFAULT_LAMBDA_MAPPING_SIGN,
) -> None:
    ensure_dir(output_dir)

    if lambda_ref_table is None or sigma_ref_table is None or charge_ref_table is None:
        lambda_ref_table, sigma_ref_table, charge_ref_table = get_default_refs()

    for gamma in gammas:
        gamma_tag = str(gamma).replace(".", "p")
        sign_tag = "minus" if mapping_sign < 0 else "plus"
        for T in temperatures_k:
            temp_tag = f"{int(round(T))}K"
            arg_tag = "argcharmm36" if use_arg_charmm36 else "argdefault"
            filename = f"tea_calvados_gamma_{gamma_tag}_{temp_tag}_{arg_tag}_{sign_tag}.csv"
            out_csv = output_dir / filename

            write_lambda_table_csv(
                output_csv=out_csv,
                temperature_k=T,
                gamma=gamma,
                lambda_ref_table=lambda_ref_table,
                sigma_ref_table=sigma_ref_table,
                charge_ref_table=charge_ref_table,
                use_arg_charmm36=use_arg_charmm36,
                t_ref=t_ref,
                mapping_sign=mapping_sign,
            )
            print(f"Wrote {out_csv}")


if __name__ == "__main__":
    base_dir = Path(__file__).resolve().parent
    output_dir = base_dir / "output" / "tea_tables"

    lambda_ref_table, sigma_ref_table, charge_ref_table = get_default_refs()

    write_all_tables(
        output_dir=output_dir,
        temperatures_k=DEFAULT_TEMPERATURES,
        gammas=DEFAULT_GAMMAS,
        lambda_ref_table=lambda_ref_table,
        sigma_ref_table=sigma_ref_table,
        charge_ref_table=charge_ref_table,
        use_arg_charmm36=False,
        t_ref=T_REF,
        mapping_sign=DEFAULT_LAMBDA_MAPPING_SIGN,
    )