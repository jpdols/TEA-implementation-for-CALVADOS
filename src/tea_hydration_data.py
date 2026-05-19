from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class HydrationParams:
    """
    Thermodynamic parameters for hydration free energies from TEA Table S1.

    Units:
        du_h_300   : kcal/mol
        dh         : kcal/mol
        dc_p_cal   : cal/mol/K
    """
    du_h_300: float
    dh: float
    dc_p_cal: float

    @property
    def dc_p_kcal(self) -> float:
        """Return Δc_p in kcal/mol/K."""
        return self.dc_p_cal / 1000.0


# Table S1 from the TEA supplement.
# Note:
# - "Val/Pro" is shared by V and P
# - "Backbone/Gly" is used for G in the one-letter mapping
TEA_HYDRATION_PARAMS: Dict[str, HydrationParams] = {
    "A": HydrationParams(du_h_300=2.45,   dh=0.02,    dc_p_cal=35.30),
    "V": HydrationParams(du_h_300=2.54,   dh=-1.13,   dc_p_cal=40.50),
    "P": HydrationParams(du_h_300=2.54,   dh=-1.13,   dc_p_cal=40.50),
    "L": HydrationParams(du_h_300=2.90,   dh=-1.17,   dc_p_cal=50.49),
    "I": HydrationParams(du_h_300=2.74,   dh=-1.93,   dc_p_cal=62.68),
    "M": HydrationParams(du_h_300=1.03,   dh=-3.75,   dc_p_cal=43.66),
    "F": HydrationParams(du_h_300=0.37,   dh=-5.71,   dc_p_cal=63.67),
    "C": HydrationParams(du_h_300=-0.02,  dh=-3.69,   dc_p_cal=39.10),
    "Y": HydrationParams(du_h_300=-4.00,  dh=-10.60,  dc_p_cal=41.91),
    "W": HydrationParams(du_h_300=-5.25,  dh=-12.94,  dc_p_cal=50.57),
    "S": HydrationParams(du_h_300=-4.38,  dh=-8.83,   dc_p_cal=21.47),
    "T": HydrationParams(du_h_300=-4.24,  dh=-10.10,  dc_p_cal=37.74),
    "N": HydrationParams(du_h_300=-7.53,  dh=-13.01,  dc_p_cal=24.54),
    "Q": HydrationParams(du_h_300=-7.27,  dh=-13.80,  dc_p_cal=36.41),
    "H": HydrationParams(du_h_300=-10.55, dh=-18.72,  dc_p_cal=41.86),
    "G": HydrationParams(du_h_300=-7.15,  dh=-14.06,  dc_p_cal=41.81),
    "R": HydrationParams(du_h_300=-45.74, dh=-57.99,  dc_p_cal=58.02),
    "K": HydrationParams(du_h_300=-53.54, dh=-62.55,  dc_p_cal=41.93),
    "D": HydrationParams(du_h_300=-97.11, dh=-111.95, dc_p_cal=3.49),
    "E": HydrationParams(du_h_300=-97.05, dh=-112.90, dc_p_cal=4.68),
}


# Optional: readable analogue names, useful for debugging / output tables
TEA_ANALOGUE_LABELS: Dict[str, str] = {
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