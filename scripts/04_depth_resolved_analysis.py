"""
04_depth_resolved_analysis.py

Creates depth-resolved comparisons from the complete SRIM vacancy and
E2RECOIL profile CSV files.

Run from the project root:
    py scripts\\04_depth_resolved_analysis.py

Inputs:
    data/input/vacancy_depth_profiles_complete.csv
    data/input/recoil_energy_depth_profiles_complete.csv

Outputs:
    data/processed/depth_profile_totals_long.csv
    data/processed/bilayer_layer_depth_profiles.csv
    figures/vacancy_depth_profiles_by_energy.png
    figures/recoil_depth_profiles_by_energy.png
    figures/bilayer_layer_resolved_vacancies.png
    figures/bilayer_layer_resolved_recoil_energy.png

Source units:
    Element-resolved SRIM values: per Angstrom per incident ion.
    Convert to per nm by multiplying by 10.
    Profiles are plotted as depth-bin values (not cumulative totals).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "input"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
PROCESSED.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

VAC_FILE = INPUT / "vacancy_depth_profiles_complete.csv"
REC_FILE = INPUT / "recoil_energy_depth_profiles_complete.csv"

A_TO_NM = 10.0
ENERGIES = [10, 30, 50, 100]
MATERIAL_ORDER = [
    "MAPbI3", "Lead_reduced_Sn_Pb",
    "Triple_cation", "MAPbI3_FAPbI3_bilayer"
]
ELEMENTS = ["C", "H", "N", "Sn", "Cs", "Pb", "I", "Br"]
BILAYER_ELEMENTS = ["C", "H", "N", "Pb", "I"]

def normalize_material(value):
    s = str(value).strip().lower()
    if "bilayer" in s:
        return "MAPbI3_FAPbI3_bilayer"
    if "lead" in s and "sn" in s:
        return "Lead_reduced_Sn_Pb"
    if "triple" in s:
        return "Triple_cation"
    if "mapbi3" in s:
        return "MAPbI3"
    return str(value).strip()

def read_profiles(path, kind):
    if not path.exists():
        raise FileNotFoundError(f"Cannot find input file: {path}")
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]
    needed = ["Material", "Energy_keV", "Depth_nm"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name} missing required columns: {missing}")
    df["Material_key"] = df["Material"].map(normalize_material)
    df["Energy_keV"] = pd.to_numeric(df["Energy_keV"], errors="coerce")
    df["Depth_nm"] = pd.to_numeric(df["Depth_nm"], errors="coerce")
    if df[["Energy_keV", "Depth_nm"]].isna().any().any():
        raise ValueError(f"Non-numeric energy/depth values in {path.name}")
    return df

def numeric_sum(row, columns):
    vals = []
    for c in columns:
        if c in row.index:
            v = pd.to_numeric(pd.Series([row[c]]), errors="coerce").iloc[0]
            if pd.notna(v):
                vals.append(float(v))
    return sum(vals) if vals else np.nan

def build_long(df, kind):
    is_vac = kind == "vacancy"
    value_name = "vacancies_per_nm_per_ion" if is_vac else "recoil_energy_eV_per_nm_per_ion"
    rows = []
    for _, row in df.iterrows():
        mat = row["Material_key"]
        energy = int(row["Energy_keV"])
        depth = float(row["Depth_nm"])
        if mat == "MAPbI3_FAPbI3_bilayer":
            layer_columns = {
                "MAPbI3": [f"{e}_layer1" for e in BILAYER_ELEMENTS],
                "FAPbI3": [f"{e}_layer2" for e in BILAYER_ELEMENTS],
            }
        else:
            layer_columns = {"Bulk": ELEMENTS}

        layer_values = {}
        for layer, cols in layer_columns.items():
            val = numeric_sum(row, cols)
            if pd.notna(val):
                layer_values[layer] = val * A_TO_NM
                rows.append({
                    "material": mat,
                    "energy_keV": energy,
                    "depth_nm": depth,
                    "layer": layer,
                    value_name: val * A_TO_NM
                })
        if layer_values:
            rows.append({
                "material": mat,
                "energy_keV": energy,
                "depth_nm": depth,
                "layer": "Total",
                value_name: sum(layer_values.values())
            })
    return pd.DataFrame(rows), value_name

vac = read_profiles(VAC_FILE, "vacancy")
rec = read_profiles(REC_FILE, "recoil")

vac_long, vac_col = build_long(vac, "vacancy")
rec_long, rec_col = build_long(rec, "recoil")

# Merge on material, energy, depth and layer to preserve bilayer layer identity.
combined = vac_long.merge(
    rec_long,
    on=["material", "energy_keV", "depth_nm", "layer"],
    how="outer",
    validate="one_to_one"
)
combined.to_csv(PROCESSED / "depth_profile_totals_long.csv", index=False)

bilayer = combined[combined["material"] == "MAPbI3_FAPbI3_bilayer"].copy()
bilayer.to_csv(PROCESSED / "bilayer_layer_depth_profiles.csv", index=False)

# General material comparison: each panel is one incident proton energy.
short_names = {
    "MAPbI3": "MAPbI3",
    "Lead_reduced_Sn_Pb": "Mixed Sn-Pb",
    "Triple_cation": "Triple-cation",
    "MAPbI3_FAPbI3_bilayer": "MAPbI3/FAPbI3 bilayer"
}

def plot_all_materials(df, value_col, ylabel, title, filename):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    axes = axes.flatten()
    for ax, energy in zip(axes, ENERGIES):
        sub = df[(df["energy_keV"] == energy) & (df["layer"] == "Total")]
        for material in MATERIAL_ORDER:
            g = sub[sub["material"] == material].sort_values("depth_nm")
            if not g.empty:
                ax.plot(g["depth_nm"], g[value_col], linewidth=2,
                        label=short_names.get(material, material))
        ax.set_title(f"{energy} keV protons")
        ax.set_xlabel("Depth (nm)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.35)
        ax.legend(fontsize=8)
    fig.suptitle(title, fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)

plot_all_materials(
    combined, vac_col, "Vacancy density (vacancies nm$^{-1}$ proton$^{-1}$)",
    "Depth-resolved vacancy production across perovskite compositions",
    "vacancy_depth_profiles_by_energy.png"
)
plot_all_materials(
    combined, rec_col, "Recoil-energy deposition (eV nm$^{-1}$ proton$^{-1}$)",
    "Depth-resolved target recoil energy across perovskite compositions",
    "recoil_depth_profiles_by_energy.png"
)

# Bilayer plots: show MAPbI3 and FAPbI3 layer contributions separately.
def plot_bilayer(df, value_col, ylabel, title, filename):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    axes = axes.flatten()
    for ax, energy in zip(axes, ENERGIES):
        sub = df[(df["energy_keV"] == energy) &
                 (df["layer"].isin(["MAPbI3", "FAPbI3"]))]
        for layer in ["MAPbI3", "FAPbI3"]:
            g = sub[sub["layer"] == layer].sort_values("depth_nm")
            if not g.empty:
                ax.plot(g["depth_nm"], g[value_col], linewidth=2, label=layer)
        ax.axvline(500, linestyle="--", linewidth=1.5, label="Layer interface (500 nm)")
        ax.set_title(f"{energy} keV protons")
        ax.set_xlabel("Depth from incident surface (nm)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.35)
        ax.legend(fontsize=8)
    fig.suptitle(title, fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)

plot_bilayer(
    combined, vac_col,
    "Vacancy density (vacancies nm$^{-1}$ proton$^{-1}$)",
    "Layer-resolved vacancy production in MAPbI3/FAPbI3 bilayer",
    "bilayer_layer_resolved_vacancies.png"
)
plot_bilayer(
    combined, rec_col,
    "Recoil-energy deposition (eV nm$^{-1}$ proton$^{-1}$)",
    "Layer-resolved target recoil energy in MAPbI3/FAPbI3 bilayer",
    "bilayer_layer_resolved_recoil_energy.png"
)

# Coverage and profile checks.
print("\nDEPTH-RESOLVED ANALYSIS COMPLETE")
print("Vacancy profiles:", vac[["Material_key", "Energy_keV"]].drop_duplicates().shape[0])
print("Recoil profiles:", rec[["Material_key", "Energy_keV"]].drop_duplicates().shape[0])
print("Combined long rows:", len(combined))
print("Bilayer rows:", len(bilayer))
print("Expected depth bins per simulation: 100")
print("\nFiles written:")
for p in [
    PROCESSED / "depth_profile_totals_long.csv",
    PROCESSED / "bilayer_layer_depth_profiles.csv",
    FIGURES / "vacancy_depth_profiles_by_energy.png",
    FIGURES / "recoil_depth_profiles_by_energy.png",
    FIGURES / "bilayer_layer_resolved_vacancies.png",
    FIGURES / "bilayer_layer_resolved_recoil_energy.png",
]:
    print(" -", p)
