from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

from calvados_lambda_ref import get_refs
from inspect_tea_tables import print_subset
from make_calvados_residue_file import write_calvados_residues
from make_tea_tables import write_all_tables
from tea_lambda import DEFAULT_LAMBDA_MAPPING_SIGN, delta_delta_g_ex, get_kappa, lambda_t
from tea_thermo import T_REF, delta_delta_mu_h, delta_mu_h

AUDIT_RESIDUES_DEFAULT = ["I", "L", "F", "M", "V", "P", "Y", "W", "A", "K", "R", "D", "E", "S", "H"]


def format_sign(sign: float) -> str:
    return "minus" if sign < 0 else "plus"


def read_csv_rows(path: Path) -> list[dict]:
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))


def compare_residue_files(original_csv: Path, new_csv: Path, residues: Iterable[str]) -> None:
    orig_rows = {row["one"]: row for row in read_csv_rows(original_csv)}
    new_rows = {row["one"]: row for row in read_csv_rows(new_csv)}

    print("\n=== Original CALVADOS vs TEA-adjusted residues file ===")
    print(f"{'res':>3}  {'lambda_orig':>12}  {'lambda_new':>12}  {'delta':>12}")
    for aa in residues:
        if aa not in orig_rows or aa not in new_rows:
            continue
        lam_orig = float(orig_rows[aa]["lambdas"])
        lam_new = float(new_rows[aa]["lambdas"])
        print(f"{aa:>3}  {lam_orig:12.6f}  {lam_new:12.6f}  {lam_new-lam_orig:12.6f}")


def audit_residues(
    temperature_k: float,
    gamma: float,
    residues: Iterable[str],
    mapping_sign: float,
    lambda_ref_table: dict[str, float],
) -> None:
    print("\n=== TEA audit (direct formulas) ===")
    print(
        f"T_ref = {T_REF:.1f} K, T = {temperature_k:.1f} K, gamma = {gamma:.3f}, "
        f"mapping_sign = {mapping_sign:+.0f}"
    )
    print(
        f"{'res':>3}  {'mu_300':>10}  {'mu_T':>10}  {'ddu_h':>10}  {'kappa':>8}  "
        f"{'ddG_ex':>10}  {'lambda_ref':>10}  {'lambda_T':>10}"
    )
    for aa in residues:
        mu_300 = delta_mu_h(aa, T_REF)
        mu_t = delta_mu_h(aa, temperature_k)
        ddu = delta_delta_mu_h(aa, temperature_k)
        kappa = get_kappa(aa)
        ddg = delta_delta_g_ex(aa, temperature_k)
        lam_ref = lambda_ref_table[aa]
        lam_t = lambda_t(
            residue=aa,
            temperature_k=temperature_k,
            lambda_ref=lam_ref,
            gamma=gamma,
            mapping_sign=mapping_sign,
        )
        print(
            f"{aa:>3}  {mu_300:10.4f}  {mu_t:10.4f}  {ddu:10.4f}  {kappa:8.4f}  "
            f"{ddg:10.4f}  {lam_ref:10.6f}  {lam_t:10.6f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the TEA -> CALVADOS residue generation pipeline with diagnostics.")
    parser.add_argument("--temperature", type=float, default=310.0, help="Target temperature in K for the main audit/output.")
    parser.add_argument("--gamma", type=float, default=3.0, help="TEA gamma scaling factor.")
    parser.add_argument("--mapping-sign", type=float, default=DEFAULT_LAMBDA_MAPPING_SIGN, help="Mapping sign from ΔΔG_E to lambda; use -1 for the corrected CALVADOS convention.")
    parser.add_argument(
        "--template",
        type=Path,
        required=True,
        help="Path to the input CALVADOS residues CSV template, e.g. data/residues_CALVADOS2.csv.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("output"),
        help="Output directory. Default: ./output",
    )
    parser.add_argument("--arg-charmm36", action="store_true", help="Use the alternative Arg kappa from the paper.")
    args = parser.parse_args()

    if not args.template.exists():
        raise FileNotFoundError(f"Template residue file not found: {args.template}")

    # Load lambda/sigma/charge refs from the template file the user passes
    lambda_ref_table, sigma_ref_table, charge_ref_table = get_refs(args.template)

    output_tables = args.outdir / "tea_tables"
    output_residues = args.outdir / "residues"
    output_tables.mkdir(parents=True, exist_ok=True)
    output_residues.mkdir(parents=True, exist_ok=True)

    audit_residues(
        args.temperature,
        args.gamma,
        AUDIT_RESIDUES_DEFAULT,
        args.mapping_sign,
        lambda_ref_table=lambda_ref_table,
    )

    write_all_tables(
        output_dir=output_tables,
        temperatures_k=[250.0, 260.0, 270.0, 278.0, 280.0, 288.0, 290.0, 293.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0],
        gammas=[args.gamma],
        use_arg_charmm36=args.arg_charmm36,
        mapping_sign=args.mapping_sign,
        lambda_ref_table=lambda_ref_table,
        sigma_ref_table=sigma_ref_table,
        charge_ref_table=charge_ref_table,
    )

    sign_tag = format_sign(args.mapping_sign)
    arg_tag = "argcharmm36" if args.arg_charmm36 else "argdefault"
    tea_csv = output_tables / f"tea_calvados_gamma_{str(args.gamma).replace('.', 'p')}_{int(round(args.temperature))}K_{arg_tag}_{sign_tag}.csv"
    residue_csv = output_residues / f"residues_gamma_{str(args.gamma).replace('.', 'p')}_{int(round(args.temperature))}K_{arg_tag}_{sign_tag}.csv"

    print_subset(tea_csv)
    write_calvados_residues(tea_csv=tea_csv, output_csv=residue_csv, template_csv=args.template)
    compare_residue_files(args.template, residue_csv, AUDIT_RESIDUES_DEFAULT)

    print("\nDone.")
    print(f"TEA table:    {tea_csv}")
    print(f"Residues CSV: {residue_csv}")


if __name__ == "__main__":
    main()