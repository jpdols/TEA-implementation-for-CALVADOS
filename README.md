# CALVADOS–TEA Temperature Mapping for Intrinsically Disordered Proteins

Implementation of a temperature-dependent hydrophobicity mapping for the CALVADOS coarse-grained model using thermodynamic quantities from the TEA framework.

This repository generates temperature-adjusted CALVADOS residue parameter tables (`residues.csv`) by mapping residue-specific hydration thermodynamics onto CALVADOS λ (“stickiness”) parameters.

---

# Overview

This implementation combines:

* **CALVADOS** residue-specific λ interaction parameters
* **TEA** residue-specific hydration thermodynamics
* Temperature-dependent mapping:

\lambda(T)=\lambda_{ref}+s\gamma\Delta\Delta G_E(T)

where:

* λ(T) = temperature-adjusted CALVADOS lambda
* γ = scaling factor
* s = mapping sign (default = −1)
* ΔΔG_E(T) = temperature-dependent excess free energy contribution from TEA

The implementation is designed for studying:

* LCST behavior
* ELP collapse
* temperature-dependent compaction
* sequence-specific IDP behavior

# Repository Structure

```text
calvados-tea/
│
├── src/
│   ├── generate_tea_residues.py
│   ├── tea_thermo.py
│   ├── tea_lambda.py
│   ├── tea_hydration_data.py
│   ├── tea_kappa_data.py
│   ├── make_tea_tables.py
│   ├── make_calvados_residue_file.py
│   ├── inspect_tea_tables.py
│   └── calvados_lambda_ref.py
│
├── data/
│   └── residues_CALVADOS2.csv
│
├── output/
│
├── README.md
├── requirements.txt
└── .gitignore
```

---

# Installation

## Clone repository

```bash
git clone https://github.com/jpdolsE/TEA-implementation-for-CALVADOS.git
cd TEA-implementation-for-CALVADOS
```

## Create conda environment

```bash
conda create -n calvados-tea python=3.11
conda activate calvados-tea
```

---

## Install dependencies

```bash
pip install -r requirements.txt
```

---

# Requirements

Minimal requirements:

```text
python >= 3.10
```

Currently only standard-library modules are used.

---

# Usage

## Generate temperature-adjusted CALVADOS residue tables

Example:

```bash
python src/generate_tea_residues.py \
    --temperature 310 \
    --gamma 6.0 \
    --template data/residues_CALVADOS2.csv \
    --outdir output
```
run using:
python src/generate_tea_residues.py --temperature 310 --gamma 6.0 --template data/residues_CALVADOS2.csv --outdir output

---

# Command-Line Arguments

| Argument         | Description                           |
| ---------------- | ------------------------------------- |
| `--temperature`  | Target temperature in Kelvin          |
| `--gamma`        | TEA scaling factor (I use 6)          |
| `--template`     | Reference CALVADOS residue CSV        |
| `--outdir`       | Output directory                      |
| `--mapping-sign` | Sign used in λ mapping (default = -1) |
| `--arg-charmm36` | Use alternative Arg κ value           |

---

# Outputs

The pipeline generates:

## 1. TEA diagnostic tables

```text
output/tea_tables/
```

Contains:

* ΔΔμ_h
* κ
* ΔΔG_E
* Δλ
* λ(T)

for every residue.

---

## 2. CALVADOS residue tables

```text
output/residues/
```

Contains:

* temperature-adjusted CALVADOS `residues.csv`

These files can be used directly in CALVADOS simulations.

---

# Physical Interpretation

In CALVADOS:

* larger λ → stronger attraction
* stronger attraction → greater collapse propensity

Hydrophobic residues become effectively more attractive at elevated temperatures.

Therefore the implementation uses:

\Delta\lambda=-\gamma\Delta\Delta G_E

such that:

* more negative ΔΔG_E
* produces larger λ
* leading to stronger hydrophobic collapse

---

# Validation

A useful sanity check:

Hydrophobic residues (e.g. I, L, F, V):

* should generally obtain larger λ at higher temperature

Charged residues:

* should change less strongly

---

# Example Workflow

1. Generate temperature-adjusted residue table:

```bash
python src/generate_tea_residues.py \
    --temperature 350 \
    --gamma 6.0 \
    --template data/residues_CALVADOS2.csv \
    --outdir output
```

2. Use generated `residues.csv` in CALVADOS

3. Run simulation

4. Analyze collapse / Rg / LCST behavior

---

# Notes

* Reference temperature defaults to:

```text
300 K
```

* The implementation currently assumes:

```text
mapping_sign = -1
```

which is required to map TEA free-energy changes onto the CALVADOS attraction convention.

---

# Citation / References

Fangke Chen, Xiangze Zeng doi: https://doi.org/10.64898/2026.01.30.702805

---

# Author

Joep Dols
TU Delft – MSc Nanobiology
