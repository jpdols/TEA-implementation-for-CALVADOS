from __future__ import annotations

from typing import Dict


# Table 1 from the main TEA paper:
# residue-specific kappa_i values used in Eq. 3
#
# Important:
# - use these kappa_i values, not the kappa_i' values from Table S2,
#   because Table S2 corresponds to the alternative molecular-specific-cutoff definition.
# - for Arg, the main table shows both "Arg" and "Arg (charmm36)".
#   For now we expose both explicitly.
TEA_KAPPA_MAIN: Dict[str, float] = {
    "A": -0.103,
    "V": -0.119,   # Val/Pro
    "P": -0.119,   # Val/Pro
    "L": -0.178,
    "I": -0.170,
    "M": -0.120,
    "F": -0.100,
    "C": -0.054,
    "Y": -0.039,
    "W": -0.039,
    "S":  0.003,
    "T": -0.009,
    "N": -0.002,
    "Q": -0.009,
    "H":  0.001,
    "G": -0.004,   # Backbone/Gly
    "R": -0.014,   # default Arg row
    "K": -0.007,
    "D": -0.002,
    "E":  0.001,
}


# Alternative Arg value shown separately in Table 1 / discussion.
# Keep this separate so we can decide later which one to use in CALVADOS+TEA.
ARG_KAPPA_CHARMM36 = 0.059


# Optional readable labels
TEA_KAPPA_LABELS: Dict[str, str] = {
    "A": "Ala",
    "V": "Val/Pro",
    "P": "Val/Pro",
    "L": "Leu",
    "I": "Ile",
    "M": "Met",
    "F": "Phe",
    "C": "Cys",
    "Y": "Tyr",
    "W": "Trp",
    "S": "Ser",
    "T": "Thr",
    "N": "Asn",
    "Q": "Gln",
    "H": "His",
    "G": "Backbone/Gly",
    "R": "Arg",
    "K": "Lys",
    "D": "Asp",
    "E": "Glu",
}