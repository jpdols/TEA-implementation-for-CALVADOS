from __future__ import annotations

from typing import Dict, Iterable

from tea_kappa_data import ARG_KAPPA_CHARMM36, TEA_KAPPA_MAIN
from tea_thermo import T_REF, delta_delta_mu_h


DEFAULT_LAMBDA_MAPPING_SIGN = -1.0


def _validate_residue(residue: str) -> str:
    if not residue or len(residue) != 1:
        raise ValueError(f"Residue must be a one-letter code, got {residue!r}.")
    return residue.upper()


def _validate_temperature(temperature_k: float) -> float:
    T = float(temperature_k)
    if T <= 0.0:
        raise ValueError(f"Temperature must be > 0 K, got {temperature_k}.")
    return T


def get_kappa(
    residue: str,
    use_arg_charmm36: bool = False,
) -> float:
    aa = _validate_residue(residue)

    if aa == "R" and use_arg_charmm36:
        return ARG_KAPPA_CHARMM36

    if aa not in TEA_KAPPA_MAIN:
        supported = ", ".join(sorted(TEA_KAPPA_MAIN))
        raise KeyError(f"Unsupported residue {aa!r}. Supported: {supported}")

    return TEA_KAPPA_MAIN[aa]


def delta_delta_g_ex(
    residue: str,
    temperature_k: float,
    t_ref: float = T_REF,
    use_arg_charmm36: bool = False,
) -> float:
    """
    Compute ΔΔG_E(T) = ΔG_E(T) - ΔG_E(T0) in kcal/mol.

    TEA Eq. 3:
        ΔΔG_E(T) = kappa_i * [Δu_h(T) - Δu_h(T0)]
    """
    aa = _validate_residue(residue)
    T = _validate_temperature(temperature_k)
    T0 = _validate_temperature(t_ref)

    kappa_i = get_kappa(aa, use_arg_charmm36=use_arg_charmm36)
    ddu = delta_delta_mu_h(aa, T, t_ref=T0)
    return kappa_i * ddu


def lambda_delta_from_ddg(
    delta_delta_g_ex_value: float,
    gamma: float,
    mapping_sign: float = DEFAULT_LAMBDA_MAPPING_SIGN,
) -> float:
    """
    Map TEA ΔΔG_E onto the CALVADOS lambda convention.

    Practical note
    --------------
    In CALVADOS, larger lambda means stronger attraction. For hydrophobics,
    ΔG_E becomes more negative at higher temperature, so to make hydrophobics
    more attractive with temperature we use a negative mapping sign by default:

        Δlambda = mapping_sign * gamma * ΔΔG_E

    With mapping_sign = -1, a negative ΔΔG_E increases lambda.
    """
    return float(mapping_sign) * float(gamma) * float(delta_delta_g_ex_value)


def lambda_t(
    residue: str,
    temperature_k: float,
    lambda_ref: float,
    gamma: float,
    t_ref: float = T_REF,
    use_arg_charmm36: bool = False,
    mapping_sign: float = DEFAULT_LAMBDA_MAPPING_SIGN,
) -> float:
    aa = _validate_residue(residue)
    T = _validate_temperature(temperature_k)
    _validate_temperature(t_ref)

    ddg = delta_delta_g_ex(
        aa,
        T,
        t_ref=t_ref,
        use_arg_charmm36=use_arg_charmm36,
    )
    return float(lambda_ref) + lambda_delta_from_ddg(ddg, gamma, mapping_sign=mapping_sign)


def lambda_table_for_residues(
    lambda_ref_table: Dict[str, float],
    temperature_k: float,
    gamma: float,
    t_ref: float = T_REF,
    use_arg_charmm36: bool = False,
    mapping_sign: float = DEFAULT_LAMBDA_MAPPING_SIGN,
) -> Dict[str, Dict[str, float]]:
    T = _validate_temperature(temperature_k)
    T0 = _validate_temperature(t_ref)

    out: Dict[str, Dict[str, float]] = {}
    for residue, lam_ref in lambda_ref_table.items():
        aa = _validate_residue(residue)
        ddu = delta_delta_mu_h(aa, T, t_ref=T0)
        kappa_i = get_kappa(aa, use_arg_charmm36=use_arg_charmm36)
        ddg = kappa_i * ddu
        dlam = lambda_delta_from_ddg(ddg, gamma, mapping_sign=mapping_sign)
        lam_t = float(lam_ref) + dlam

        out[aa] = {
            "lambda_ref": float(lam_ref),
            "delta_delta_mu_h": float(ddu),
            "kappa": float(kappa_i),
            "delta_delta_g_ex": float(ddg),
            "delta_lambda": float(dlam),
            "lambda_t": float(lam_t),
        }
    return out


def lambda_table_over_temperatures(
    lambda_ref_table: Dict[str, float],
    temperatures_k: Iterable[float],
    gamma: float,
    t_ref: float = T_REF,
    use_arg_charmm36: bool = False,
    mapping_sign: float = DEFAULT_LAMBDA_MAPPING_SIGN,
) -> Dict[float, Dict[str, Dict[str, float]]]:
    out: Dict[float, Dict[str, Dict[str, float]]] = {}
    for T in temperatures_k:
        T = _validate_temperature(T)
        out[T] = lambda_table_for_residues(
            lambda_ref_table=lambda_ref_table,
            temperature_k=T,
            gamma=gamma,
            t_ref=t_ref,
            use_arg_charmm36=use_arg_charmm36,
            mapping_sign=mapping_sign,
        )
    return out


def sequence_mean_lambda(
    sequence: str,
    lambda_ref_table: Dict[str, float],
    temperature_k: float,
    gamma: float,
    t_ref: float = T_REF,
    use_arg_charmm36: bool = False,
    mapping_sign: float = DEFAULT_LAMBDA_MAPPING_SIGN,
) -> float:
    if not sequence:
        raise ValueError("Sequence must not be empty.")
    seq = sequence.upper()
    lambdas = [
        lambda_t(
            residue=aa,
            temperature_k=temperature_k,
            lambda_ref=lambda_ref_table[aa],
            gamma=gamma,
            t_ref=t_ref,
            use_arg_charmm36=use_arg_charmm36,
            mapping_sign=mapping_sign,
        )
        for aa in seq
    ]
    return sum(lambdas) / len(lambdas)


if __name__ == "__main__":
    lambda_ref_demo = {
        "V": 1.00,
        "P": 1.00,
        "G": 0.50,
        "I": 1.10,
        "K": 0.20,
        "R": 0.20,
        "D": 0.10,
        "E": 0.10,
        "H": 0.30,
    }

    temps = [280.0, 293.0, 300.0, 310.0, 320.0, 330.0]
    gamma = 2.5

    for T in temps:
        print(f"\n=== T = {T:.1f} K, gamma = {gamma:.2f}, mapping_sign = {DEFAULT_LAMBDA_MAPPING_SIGN:+.0f} ===")
        table = lambda_table_for_residues(
            lambda_ref_table=lambda_ref_demo,
            temperature_k=T,
            gamma=gamma,
            mapping_sign=DEFAULT_LAMBDA_MAPPING_SIGN,
        )
        for aa, vals in sorted(table.items()):
            print(
                f"{aa}: "
                f"λ_ref={vals['lambda_ref']:7.4f}  "
                f"ΔΔμ_h={vals['delta_delta_mu_h']:8.4f}  "
                f"κ={vals['kappa']:7.4f}  "
                f"ΔΔG_E={vals['delta_delta_g_ex']:8.4f}  "
                f"Δλ={vals['delta_lambda']:8.4f}  "
                f"λ(T)={vals['lambda_t']:8.4f}"
            )
