"""
02_process_recoil_energy.py

Processes the consolidated SRIM E2RECOIL depth-profile CSV.

Input:
  data/input/recoil_energy_depth_profiles_complete.csv

Outputs:
  data/processed/recoil_energy_depth_profiles_long.csv
  data/processed/recoil_energy_element_summary.csv
  data/processed/recoil_energy_integral_validation.csv
  data/processed/recoil_energy_material_summary.csv
  figures/recoil_total_profile_<material>.png
  figures/recoil_element_profiles_<material>_<energy>keV.png
  figures/recoil_energy_vs_incident_energy.png
  figures/recoil_element_contributions_<energy>keV.png

SRIM E2RECOIL elemental columns are in eV/(Angstrom * incident ion).
Depth bins are 100 Angstrom = 10 nm.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------- PATHS -----------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_DIR / "data" / "input" / "recoil_energy_depth_profiles_complete.csv"
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
FIGURES_DIR = PROJECT_DIR / "figures"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

BIN_WIDTH_A = 100.0
A_PER_NM = 10.0
BIN_WIDTH_NM = BIN_WIDTH_A / A_PER_NM

# Columns used by homogeneous-material profiles.
ELEMENTS = ["C", "H", "N", "Sn", "Cs", "Pb", "I", "Br"]

# Layer-specific columns used by the bilayer profile.
BILAYER_COLUMNS = {
    "MAPbI3": ["C_layer1", "H_layer1", "N_layer1", "Pb_layer1", "I_layer1"],
    "FAPbI3": ["C_layer2", "H_layer2", "N_layer2", "Pb_layer2", "I_layer2"],
}

# --------------------------- HELPERS --------------------------------
def find_column(columns, candidates):
    """Return first matching column, case-insensitive."""
    lookup = {str(c).strip().lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def clean_name(value):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_")


# --------------------------- LOAD ------------------------------------
if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Missing input file: {INPUT_FILE}\n"
        "Place recoil_energy_depth_profiles_complete.csv in data/input/."
    )

df = pd.read_csv(INPUT_FILE)
df.columns = [str(c).strip() for c in df.columns]

material_col = find_column(df.columns, ["Material", "material"])
energy_col = find_column(df.columns, ["Energy_keV", "energy_keV"])
depth_col = find_column(df.columns, ["Depth_nm", "depth_nm"])

if not all([material_col, energy_col, depth_col]):
    raise ValueError(
        "Required columns not found. Need Material, Energy_keV and Depth_nm. "
        f"Found: {df.columns.tolist()}"
    )

df[energy_col] = pd.to_numeric(df[energy_col], errors="coerce")
df[depth_col] = pd.to_numeric(df[depth_col], errors="coerce")

if df[[material_col, energy_col, depth_col]].isna().any().any():
    raise ValueError("Missing or non-numeric material, energy or depth values.")

if df.duplicated([material_col, energy_col, depth_col]).any():
    raise ValueError("Duplicate material-energy-depth rows detected.")

print("INPUT CHECK")
print(f"Rows: {len(df)}")
print(f"Material-energy profiles: {df.groupby([material_col, energy_col]).ngroups}")
print(df.groupby([material_col, energy_col]).size().to_string())

# ---------------------- WIDE TO LONG ---------------------------------
records = []

for _, row in df.iterrows():
    material = str(row[material_col])
    energy = int(row[energy_col])
    depth_nm = float(row[depth_col])

    is_bilayer = "bilayer" in material.lower()

    if is_bilayer:
        for layer, columns in BILAYER_COLUMNS.items():
            for col in columns:
                if col not in df.columns:
                    continue
                value = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
                if pd.notna(value):
                    element = col.split("_")[0]
                    # Raw elemental input is eV/(A ion), convert to eV/(nm ion).
                    records.append({
                        "material": material,
                        "energy_keV": energy,
                        "layer": layer,
                        "depth_nm": depth_nm,
                        "element": element,
                        "recoil_energy_eV_per_nm_ion": float(value) * A_PER_NM,
                    })
    else:
        for element in ELEMENTS:
            if element not in df.columns:
                continue
            value = pd.to_numeric(pd.Series([row[element]]), errors="coerce").iloc[0]
            if pd.notna(value):
                records.append({
                    "material": material,
                    "energy_keV": energy,
                    "layer": "Single",
                    "depth_nm": depth_nm,
                    "element": element,
                    "recoil_energy_eV_per_nm_ion": float(value) * A_PER_NM,
                })

long_df = pd.DataFrame(records)
if long_df.empty:
    raise ValueError("No elemental recoil columns were found. Check input column names.")

long_df = long_df.sort_values(
    ["material", "energy_keV", "layer", "element", "depth_nm"]
).reset_index(drop=True)

long_path = PROCESSED_DIR / "recoil_energy_depth_profiles_long.csv"
long_df.to_csv(long_path, index=False)

# --------------------- ELEMENT INTEGRALS -----------------------------
# Converted profile is eV/(nm ion); integrate over 10 nm bins.
element_summary = (
    long_df.groupby(["material", "energy_keV", "layer", "element"], as_index=False)
    .agg(
        recoil_energy_eV_per_ion=(
            "recoil_energy_eV_per_nm_ion",
            lambda x: x.sum() * BIN_WIDTH_NM,
        )
    )
)

material_totals = (
    element_summary.groupby(["material", "energy_keV"], as_index=False)
    .agg(total_target_recoil_energy_eV_per_ion=("recoil_energy_eV_per_ion", "sum"))
)

element_summary = element_summary.merge(
    material_totals, on=["material", "energy_keV"], how="left"
)
element_summary["element_fraction"] = (
    element_summary["recoil_energy_eV_per_ion"]
    / element_summary["total_target_recoil_energy_eV_per_ion"]
)
element_summary["element_fraction_pct"] = element_summary["element_fraction"] * 100

element_path = PROCESSED_DIR / "recoil_energy_element_summary.csv"
element_summary.to_csv(element_path, index=False)

# --------------------- PROFILE VALIDATION ----------------------------
# Check that the per-bin elemental sum equals the sum of elemental columns.
# If a precomputed total column exists, compare against it but do not trust
# its label alone: the earlier consolidated file had a factor-of-10 mismatch.
total_col = find_column(
    df.columns,
    [
        "Total_target_recoil_energy_eV_per_nm_per_proton",
        "Total_target_recoil_energy_eV_per_nm_ion",
        "Total_recoil_energy_eV_per_nm_ion",
    ],
)

validation_rows = []
for (material, energy), group in long_df.groupby(["material", "energy_keV"]):
    per_bin = (
        group.groupby("depth_nm")["recoil_energy_eV_per_nm_ion"].sum()
        .sort_index()
    )
    integrated_from_elements = float(per_bin.sum() * BIN_WIDTH_NM)

    comparison = np.nan
    difference_pct = np.nan
    if total_col:
        source_group = df[
            (df[material_col].astype(str) == str(material))
            & (df[energy_col].astype(int) == int(energy))
        ].sort_values(depth_col)
        source_total = pd.to_numeric(source_group[total_col], errors="coerce")
        # Column is expected to be per nm; integrate using 10 nm bins.
        comparison = float(source_total.sum() * BIN_WIDTH_NM)
        if comparison != 0:
            difference_pct = abs(integrated_from_elements - comparison) / comparison * 100

    validation_rows.append({
        "material": material,
        "energy_keV": int(energy),
        "integrated_element_sum_eV_per_ion": integrated_from_elements,
        "integrated_source_total_eV_per_ion": comparison,
        "difference_percent": difference_pct,
        "depth_bins": int(group["depth_nm"].nunique()),
    })

validation = pd.DataFrame(validation_rows)
validation_path = PROCESSED_DIR / "recoil_energy_integral_validation.csv"
validation.to_csv(validation_path, index=False)

# --------------------- MATERIAL SUMMARY ------------------------------
summary = (
    element_summary.groupby(["material", "energy_keV"], as_index=False)
    .agg(
        total_target_recoil_energy_eV_per_ion=("recoil_energy_eV_per_ion", "sum")
    )
)
summary_path = PROCESSED_DIR / "recoil_energy_material_summary.csv"
summary.to_csv(summary_path, index=False)

# --------------------- PLOTS -----------------------------------------
plt.rcParams.update({
    "figure.figsize": (10, 6),
    "axes.grid": True,
    "grid.alpha": 0.25,
    "font.size": 11,
})

# Total recoil energy depth profiles, one figure per material.
for material in sorted(long_df["material"].unique()):
    subset = long_df[long_df["material"] == material]
    profile = (
        subset.groupby(["energy_keV", "depth_nm"], as_index=False)
        ["recoil_energy_eV_per_nm_ion"].sum()
    )

    fig, ax = plt.subplots()
    for energy in sorted(profile["energy_keV"].unique()):
        curve = profile[profile["energy_keV"] == energy]
        ax.plot(
            curve["depth_nm"],
            curve["recoil_energy_eV_per_nm_ion"],
            label=f"{energy} keV",
            linewidth=2,
        )
    ax.set_xlabel("Depth (nm)")
    ax.set_ylabel("Target recoil energy (eV nm$^{-1}$ ion$^{-1}$)")
    ax.set_title(f"Target recoil-energy profile: {material}")
    ax.legend(title="Incident proton energy")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"recoil_total_profile_{clean_name(material)}.png", dpi=300)
    plt.close(fig)

# Element-specific profiles for each material-energy pair.
for (material, energy), subset in long_df.groupby(["material", "energy_keV"]):
    fig, ax = plt.subplots()
    for element in sorted(subset["element"].unique()):
        curve = subset[subset["element"] == element].sort_values("depth_nm")
        ax.plot(
            curve["depth_nm"],
            curve["recoil_energy_eV_per_nm_ion"],
            label=element,
            linewidth=1.8,
        )
    ax.set_xlabel("Depth (nm)")
    ax.set_ylabel("Recoil energy (eV nm$^{-1}$ ion$^{-1}$)")
    ax.set_title(f"Element-resolved recoil energy: {material}, {energy} keV")
    ax.legend(title="Element")
    fig.tight_layout()
    fig.savefig(
        FIGURES_DIR / f"recoil_element_profiles_{clean_name(material)}_{energy}keV.png",
        dpi=300,
    )
    plt.close(fig)

# Total recoil energy versus incident energy.
fig, ax = plt.subplots()
for material in sorted(summary["material"].unique()):
    curve = summary[summary["material"] == material].sort_values("energy_keV")
    ax.plot(
        curve["energy_keV"],
        curve["total_target_recoil_energy_eV_per_ion"],
        marker="o",
        linewidth=2,
        label=material,
    )
ax.set_xlabel("Incident proton energy (keV)")
ax.set_ylabel("Integrated target recoil energy (eV/ion)")
ax.set_title("Integrated target recoil energy versus proton energy")
ax.legend()
fig.tight_layout()
fig.savefig(FIGURES_DIR / "recoil_energy_vs_incident_energy.png", dpi=300)
plt.close(fig)

# Element contributions, grouped by material and energy.
for energy in sorted(element_summary["energy_keV"].unique()):
    subset = element_summary[element_summary["energy_keV"] == energy].copy()
    subset["sample"] = subset["material"] + " | " + subset["layer"]
    pivot = subset.pivot_table(
        index="sample",
        columns="element",
        values="element_fraction_pct",
        aggfunc="sum",
        fill_value=0,
    )
    ax = pivot.plot(kind="bar", stacked=True, figsize=(12, 7), width=0.8)
    ax.set_ylabel("Contribution to target recoil energy (%)")
    ax.set_xlabel("")
    ax.set_title(f"Elemental recoil-energy contributions at {energy} keV")
    ax.legend(title="Element", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f"recoil_element_contributions_{energy}keV.png", dpi=300)
    plt.close()

# --------------------- REPORT ----------------------------------------
print("\nPROCESSING COMPLETE")
print(f"Long-format profiles: {long_path}")
print(f"Element summary: {element_path}")
print(f"Validation: {validation_path}")
print(f"Material summary: {summary_path}")
print(f"Figures: {FIGURES_DIR}")
print("\nIntegrated recoil-energy summary:")
print(summary.to_string(index=False))
print("\nValidation report:")
print(validation.to_string(index=False))
