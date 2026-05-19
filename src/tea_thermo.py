from __future__ import annotations

from math import log
from typing import Dict, Iterable

from tea_hydration_data import HydrationParams, TEA_HYDRATION_PARAMS

T_REF = 300.0  # K


def _validate_temperature(temperature_k: float) -> float:
    """Validate and return temperature as float."""
    T = float(temperature_k)
    if T <= 0.0:
        raise ValueError(f"Temperature must be > 0 K, got {temperature_k}.")
    return T


def _validate_residue(residue: str) -> str:
    """Validate one-letter residue code against the TEA hydration table."""
    if not residue or len(residue) != 1:
        raise ValueError(f"Residue must be a one-letter code, got {residue!r}.")
    aa = residue.upper()
    if aa not in TEA_HYDRATION_PARAMS:
        supported = ", ".join(sorted(TEA_HYDRATION_PARAMS))
        raise KeyError(f"Unsupported residue {aa!r}. Supported: {supported}")
    return aa


def get_hydration_params(residue: str) -> HydrationParams:
    """Return TEA hydration parameters for a residue."""
    aa = _validate_residue(residue)
    return TEA_HYDRATION_PARAMS[aa]


def delta_mu_h(
    residue: str,
    temperature_k: float,
    t_ref: float = T_REF,
) -> float:
    """
    Compute hydration free energy Δu_h(T) in kcal/mol.

    Uses the integrated Gibbs–Helmholtz form with the fitted parameters
    Δu_h(T0), Δh, and Δc_p from TEA Table S1.

    Formula implemented:
        Δu_h(T) =
            Δu_h(T0) * (T / T0)
            + Δh * (1 - T / T0)
            + Δc_p * ((T - T0) - T * ln(T / T0))

    Notes
    -----
    - Δu_h(T0) and Δh are in kcal/mol
    - Δc_p is converted from cal/mol/K to kcal/mol/K
    - T0 defaults to 300 K, matching the paper
    """
    aa = _validate_residue(residue)
    T = _validate_temperature(temperature_k)
    T0 = _validate_temperature(t_ref)

    p = TEA_HYDRATION_PARAMS[aa]
    cp = p.dc_p_kcal  # kcal/mol/K

    return (
        p.du_h_300 * (T / T0)
        + p.dh * (1.0 - T / T0)
        + cp * ((T - T0) - T * log(T / T0))
    )


def delta_delta_mu_h(
    residue: str,
    temperature_k: float,
    t_ref: float = T_REF,
) -> float:
    """
    Compute Δu_h(T) - Δu_h(T0) in kcal/mol.
    """
    aa = _validate_residue(residue)
    T = _validate_temperature(temperature_k)
    T0 = _validate_temperature(t_ref)

    return delta_mu_h(aa, T, T0) - delta_mu_h(aa, T0, T0)


def hydration_table_row(
    residue: str,
    temperatures_k: Iterable[float],
    t_ref: float = T_REF,
) -> list[dict]:
    """
    Return a small table of Δu_h(T) and ΔΔu_h(T) values for one residue.
    """
    aa = _validate_residue(residue)
    rows: list[dict] = []

    for T in temperatures_k:
        T = _validate_temperature(T)
        rows.append(
            {
                "residue": aa,
                "T_K": T,
                "delta_mu_h_kcal_per_mol": delta_mu_h(aa, T, t_ref),
                "delta_delta_mu_h_kcal_per_mol": delta_delta_mu_h(aa, T, t_ref),
            }
        )
    return rows


def sequence_hydration_profile(
    sequence: str,
    temperature_k: float,
    t_ref: float = T_REF,
) -> Dict[str, float]:
    """
    Summarize residue-level hydration quantities for a sequence.

    Returns
    -------
    dict with:
        mean_delta_mu_h
        mean_delta_delta_mu_h
        sum_delta_mu_h
        sum_delta_delta_mu_h

    Notes
    -----
    This is only a simple residue-sum diagnostic.
    It is not yet the TEA interaction model.
    """
    if not sequence:
        raise ValueError("Sequence must not be empty.")

    seq = sequence.upper()
    values = [delta_mu_h(aa, temperature_k, t_ref) for aa in seq]
    dvalues = [delta_delta_mu_h(aa, temperature_k, t_ref) for aa in seq]

    n = len(seq)
    return {
        "length": float(n),
        "mean_delta_mu_h": sum(values) / n,
        "mean_delta_delta_mu_h": sum(dvalues) / n,
        "sum_delta_mu_h": sum(values),
        "sum_delta_delta_mu_h": sum(dvalues),
    }


if __name__ == "__main__":
    test_residues = ["V", "P", "G", "I", "K", "R", "D", "E", "H"]
    test_temps = [280.0, 293.0, 300.0, 310.0, 320.0, 330.0]

    for aa in test_residues:
        print(f"\nResidue {aa}")
        for row in hydration_table_row(aa, test_temps):
            print(
                f"T={row['T_K']:6.1f} K  "
                f"Δu_h={row['delta_mu_h_kcal_per_mol']:9.4f} kcal/mol  "
                f"ΔΔu_h={row['delta_delta_mu_h_kcal_per_mol']:9.4f} kcal/mol"
            )