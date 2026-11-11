"""
01_process_vacancy_profiles.py

Processes SRIM/TRIM VACANCY depth profiles from the consolidated wide-format CSV.

Input:
  data/input/vacancy_depth_profiles_complete.csv

Outputs:
  data/processed/vacancy_depth_profiles_long.csv
  data/processed/element_resolved_vacancies.csv
  data/processed/vacancy_profile_validation.csv
  figures/vacancy_profiles_by_material.png
  figures/elemental_vacancy_contributions.png

Units:
  Input vacancy columns: vacancies / (Angstrom * incident ion)
  Output profile: vacancies / (nm * incident ion)
  Integrated elemental vacancies: vacancies / incident ion

For SRIM's 100 depth bins across 10,000 Angstrom, each bin is 100 Angstrom.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------
# 1. PATHS AND SETTINGS
# ---------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_DIR / "data" / "input" / "vacancy_depth_profiles_complete.csv"
MASTER_FILE = PROJECT_DIR / "data" / "input" / "master_simulation_summary.csv"

PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
FIGURES_DIR = PROJECT_DIR / "figures"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

BIN_WIDTH_ANGSTROM = 100.0
ANGSTROM_PER_NM = 10.0

# Element columns for homogeneous compositions.
ELEMENT_COLUMNS = ["C", "H", "N", "Sn", "Cs", "Pb", "I", "Br"]

# Bilayer columns are layer-specific in the source CSV.
BILAYER_LAYER_COLUMNS = {
    "MAPbI3": ["C_layer1", "H_layer1", "N_layer1", "Pb_layer1", "I_layer1"],
    "FAPbI3": ["C_layer2", "H_layer2", "N_layer2", "Pb_layer2", "I_layer2"],
}


# ---------------------------------------------------------------------
# 2. LOAD AND CHECK INPUT
# ---------------------------------------------------------------------
if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Could not find {INPUT_FILE}\n"
        "Place vacancy_depth_profiles_complete.csv in data/input/."
    )

df = pd.read_csv(INPUT_FILE)

required = {"Material", "Energy_keV", "Depth_nm"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

expected_energies = {10, 30, 50, 100}
simulation_counts = df.groupby("Material")["Energy_keV"].nunique()

if df.duplicated(["Material", "Energy_keV", "Depth_nm"]).any():
    raise ValueError("Duplicate material-energy-depth records found.")

print("\nINPUT CHECK")
print(f"Rows: {len(df)}")
print(f"Materials: {df['Material'].nunique()}")
print(f"Material-energy combinations: {df.groupby(['Material', 'Energy_keV']).ngroups}")
print("\nDepth bins per simulation:")
print(df.groupby(["Material", "Energy_keV"]).size().to_string())

if not set(df["Energy_keV"].unique()).issubset(expected_energies):
    raise ValueError("Unexpected proton energy found in input.")


# ---------------------------------------------------------------------
# 3. CONVERT WIDE FORMAT TO LONG FORMAT
# ---------------------------------------------------------------------
records = []

for _, row in df.iterrows():
    material = row["Material"]
    energy = int(row["Energy_keV"])
    depth_nm = float(row["Depth_nm"])

    if material == "MAPbI3_FAPbI3_bilayer":
        for layer_name, columns in BILAYER_LAYER_COLUMNS.items():
            for col in columns:
                value = row.get(col, np.nan)
                if pd.notna(value):
                    element = col.split("_")[0]
                    records.append({
                        "material": material,
                        "energy_keV": energy,
                        "layer": layer_name,
                        "depth_nm": depth_nm,
                        "element": element,
                        "vacancy_density_per_nm_ion": float(value) * ANGSTROM_PER_NM,
                    })
    else:
        for element in ELEMENT_COLUMNS:
            if element in df.columns:
                value = row.get(element, np.nan)
                if pd.notna(value):
                    records.append({
                        "material": material,
                        "energy_keV": energy,
                        "layer": "Single",
                        "depth_nm": depth_nm,
                        "element": element,
                        "vacancy_density_per_nm_ion": float(value) * ANGSTROM_PER_NM,
                    })

long_df = pd.DataFrame(records)
long_df = long_df.sort_values(
    ["material", "energy_keV", "layer", "element", "depth_nm"]
).reset_index(drop=True)

long_file = PROCESSED_DIR / "vacancy_depth_profiles_long.csv"
long_df.to_csv(long_file, index=False)


# ---------------------------------------------------------------------
# 4. INTEGRATE EACH ELEMENTAL PROFILE
# ---------------------------------------------------------------------
# Original SRIM values are per Angstrom per ion.
# Therefore integrate with the original bin width in Angstrom:
# V_element/ion = sum(profile_per_A_ion) * 100 A.
# In the converted long data, profile is per nm per ion and bin width
# is 10 nm, giving the same result.
element_summary = (
    long_df.groupby(["material", "energy_keV", "layer", "element"], as_index=False)
    .agg(
        vacancies_per_ion=(
            "vacancy_density_per_nm_ion",
            lambda x: x.sum() * (BIN_WIDTH_ANGSTROM / ANGSTROM_PER_NM),
        )
    )
)

group_totals = (
    element_summary.groupby(["material", "energy_keV"], as_index=False)
    .agg(total_elemental_vacancies_per_ion=("vacancies_per_ion", "sum"))
)

element_summary = element_summary.merge(
    group_totals, on=["material", "energy_keV"], how="left"
)
element_summary["vacancy_fraction"] = (
    element_summary["vacancies_per_ion"]
    / element_summary["total_elemental_vacancies_per_ion"]
)
element_summary["vacancy_fraction_pct"] = element_summary["vacancy_fraction"] * 100

element_file = PROCESSED_DIR / "element_resolved_vacancies.csv"
element_summary.to_csv(element_file, index=False)


# ---------------------------------------------------------------------
# 5. VALIDATE PROFILE INTEGRALS
# ---------------------------------------------------------------------
# The source file's Total_vacancies_per_A_per_ion column is integrated
# independently as a profile-level check.
validation_records = []

for (material, energy), group in df.groupby(["Material", "Energy_keV"]):
    group = group.sort_values("Depth_nm")
    if "Total_vacancies_per_A_per_ion" in group.columns:
        total_from_profile = (
            pd.to_numeric(group["Total_vacancies_per_A_per_ion"], errors="coerce")
            .sum()
            * BIN_WIDTH_ANGSTROM
        )
    else:
        total_from_profile = np.nan

    elemental_total = group_totals.loc[
        (group_totals["material"] == material)
        & (group_totals["energy_keV"] == energy),
        "total_elemental_vacancies_per_ion",
    ].iloc[0]

    difference_pct = (
        abs(elemental_total - total_from_profile) / total_from_profile * 100
        if pd.notna(total_from_profile) and total_from_profile != 0
        else np.nan
    )

    validation_records.append({
        "material": material,
        "energy_keV": int(energy),
        "elemental_integral_vacancies_per_ion": elemental_total,
        "total_profile_integral_vacancies_per_ion": total_from_profile,
        "difference_percent": difference_pct,
        "depth_bins": len(group),
    })

validation = pd.DataFrame(validation_records)
validation_file = PROCESSED_DIR / "vacancy_profile_validation.csv"
validation.to_csv(validation_file, index=False)


# ---------------------------------------------------------------------
# 6. OPTIONAL CHECK AGAINST MASTER TDATA SUMMARY
# ---------------------------------------------------------------------
if MASTER_FILE.exists():
    master = pd.read_csv(MASTER_FILE)
    print("\nMaster summary found. Columns:")
    print(master.columns.tolist())
    print(
        "The master file is retained for a later TDATA cross-check. "
        "Check/standardise its column names before joining automatically."
    )
else:
    print(
        "\nMaster summary not found at the expected path. "
        "Profile integration validation has still been completed."
    )


# ---------------------------------------------------------------------
# 7. PLOT VACANCY PROFILES
# ---------------------------------------------------------------------
sns.set_theme(context="talk", style="whitegrid")

# One figure per material, with energy curves for total vacancy profile.
for material in sorted(long_df["material"].unique()):
    subset = long_df[long_df["material"] == material].copy()
    profile = (
        subset.groupby(["energy_keV", "depth_nm"], as_index=False)
        ["vacancy_density_per_nm_ion"].sum()
    )

    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=profile,
        x="depth_nm",
        y="vacancy_density_per_nm_ion",
        hue="energy_keV",
        palette="viridis",
        linewidth=2,
    )
    plt.xlabel("Depth (nm)")
    plt.ylabel("Vacancy density (vacancies nm$^{-1}$ ion$^{-1}$)")
    plt.title(f"Depth-dependent vacancy production: {material}")
    plt.tight_layout()
    safe_name = material.replace("/", "_").replace(" ", "_")
    plt.savefig(FIGURES_DIR / f"vacancy_profile_{safe_name}.png", dpi=300)
    plt.close()


# ---------------------------------------------------------------------
# 8. PLOT ELEMENTAL CONTRIBUTIONS
# ---------------------------------------------------------------------
for energy in sorted(element_summary["energy_keV"].unique()):
    subset = element_summary[element_summary["energy_keV"] == energy].copy()
    subset["sample"] = subset["material"] + " | " + subset["layer"]

    pivot = subset.pivot_table(
        index="sample",
        columns="element",
        values="vacancy_fraction_pct",
        aggfunc="sum",
        fill_value=0,
    )

    ax = pivot.plot(
        kind="bar",
        stacked=True,
        figsize=(12, 7),
        colormap="tab20",
        width=0.8,
    )
    ax.set_ylabel("Contribution to vacancies (%)")
    ax.set_xlabel("")
    ax.set_title(f"Elemental vacancy contributions at {energy} keV")
    ax.legend(title="Element", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f"elemental_contributions_{energy}keV.png", dpi=300)
    plt.close()


# ---------------------------------------------------------------------
# 9. PRINT COMPLETION SUMMARY
# ---------------------------------------------------------------------
print("\nPROCESSING COMPLETE")
print(f"Long-format depth profiles: {long_file}")
print(f"Element-resolved vacancies: {element_file}")
print(f"Profile validation: {validation_file}")
print(f"Figures saved in: {FIGURES_DIR}")
print("\nValidation preview:")
print(validation.to_string(index=False))
