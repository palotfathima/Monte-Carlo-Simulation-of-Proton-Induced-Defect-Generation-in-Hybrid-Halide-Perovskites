"""
03_combined_analysis.py

Combines SRIM vacancy and E2RECOIL depth profiles for 4 materials x 4 energies.

Expected project structure:
  data/input/vacancy_depth_profiles_complete.csv
  data/input/recoil_energy_depth_profiles_complete.csv
  data/input/master_simulation_data.csv

Outputs:
  data/processed/combined_material_energy_summary.csv
  data/processed/combined_element_summary.csv
  data/processed/combined_depth_profiles.csv
  data/processed/combined_validation_report.csv
  figures/combined_vacancies_vs_recoil_energy.png
  figures/combined_vacancies_by_energy.png
  figures/combined_recoil_energy_by_energy.png
  figures/combined_element_vacancy_contributions.png
  figures/combined_element_recoil_contributions.png

Units:
  Vacancy source: vacancies / (Angstrom * ion)
  E2RECOIL source: eV / (Angstrom * ion)
  Both are integrated over 100 Angstrom (10 nm) bins.
  Integrated outputs are per incident proton.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = PROJECT_DIR / "data" / "input"
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
FIGURES_DIR = PROJECT_DIR / "figures"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

VACANCY_FILE = INPUT_DIR / "vacancy_depth_profiles_complete.csv"
RECOIL_FILE = INPUT_DIR / "recoil_energy_depth_profiles_complete.csv"
MASTER_FILE = INPUT_DIR / "master_simulation_data.csv"

BIN_WIDTH_A = 100.0
A_PER_NM = 10.0
BIN_WIDTH_NM = BIN_WIDTH_A / A_PER_NM
EXPECTED_ENERGIES = {10, 30, 50, 100}
ELEMENTS = ["C", "H", "N", "Sn", "Cs", "Pb", "I", "Br"]
BILAYER_ELEMENTS = ["C", "H", "N", "Pb", "I"]

def canonical_material(value):
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

def check_profile_file(path, label):
    if not path.exists():
        raise FileNotFoundError(f"Missing {label} input:\n{path}")
    d = pd.read_csv(path)
    d.columns = [str(c).strip() for c in d.columns]
    required = {"Material", "Energy_keV", "Depth_nm"}
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"{label} file is missing columns: {sorted(missing)}")
    d["Material_key"] = d["Material"].map(canonical_material)
    d["Energy_keV"] = pd.to_numeric(d["Energy_keV"], errors="coerce")
    d["Depth_nm"] = pd.to_numeric(d["Depth_nm"], errors="coerce")
    if d[["Material_key", "Energy_keV", "Depth_nm"]].isna().any().any():
        raise ValueError(f"{label} file contains missing/non-numeric keys.")
    if d.duplicated(["Material_key", "Energy_keV", "Depth_nm"]).any():
        raise ValueError(f"{label} file has duplicate material-energy-depth rows.")
    return d

def extract_integrals(df, kind):
    """Return elemental integrals and per-bin long profiles."""
    is_vacancy = kind == "vacancy"
    value_col = (
        "vacancy_density_per_nm_ion" if is_vacancy
        else "recoil_energy_eV_per_nm_ion"
    )
    records = []
    depth_records = []
    for _, row in df.iterrows():
        material = row["Material_key"]
        energy = int(row["Energy_keV"])
        depth = float(row["Depth_nm"])
        bilayer = material == "MAPbI3_FAPbI3_bilayer"
        if bilayer:
            layer_map = {
                "MAPbI3": [f"{e}_layer1" for e in BILAYER_ELEMENTS],
                "FAPbI3": [f"{e}_layer2" for e in BILAYER_ELEMENTS],
            }
        else:
            layer_map = {"Single": ELEMENTS}

        for layer, cols in layer_map.items():
            for col in cols:
                if col not in df.columns:
                    continue
                val = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
                if pd.isna(val):
                    continue
                element = col.split("_")[0]
                # SRIM raw elemental values are per Angstrom per incident ion.
                per_nm = float(val) * A_PER_NM
                records.append({
                    "material": material, "energy_keV": energy,
                    "layer": layer, "element": element,
                    "integrated_value": per_nm * BIN_WIDTH_NM,
                })
                depth_records.append({
                    "material": material, "energy_keV": energy,
                    "layer": layer, "element": element,
                    "depth_nm": depth, value_col: per_nm,
                })

    detail = pd.DataFrame(records)
    profiles = pd.DataFrame(depth_records)
    if detail.empty or profiles.empty:
        raise ValueError(f"No elemental columns could be parsed from {kind} data.")
    return detail, profiles, value_col

print("Loading vacancy, recoil-energy and master summary files...")
vac = check_profile_file(VACANCY_FILE, "Vacancy")
rec = check_profile_file(RECOIL_FILE, "Recoil-energy")

# Confirm expected simulation coverage before calculating comparisons.
def coverage(d):
    return set((r.Material_key, int(r.Energy_keV)) for r in d[["Material_key", "Energy_keV"]].drop_duplicates().itertuples(index=False))
vac_cov, rec_cov = coverage(vac), coverage(rec)
if vac_cov != rec_cov:
    missing_recoil = sorted(vac_cov - rec_cov)
    missing_vacancy = sorted(rec_cov - vac_cov)
    raise ValueError(
        f"Vacancy/recoil simulation coverage differs.\n"
        f"Missing recoil profiles: {missing_recoil}\n"
        f"Missing vacancy profiles: {missing_vacancy}"
    )
if len(vac_cov) != 16:
    print(f"WARNING: expected 16 profiles, found {len(vac_cov)}.")

vac_detail, vac_profiles, vac_profile_col = extract_integrals(vac, "vacancy")
rec_detail, rec_profiles, rec_profile_col = extract_integrals(rec, "recoil")

# Sum the 100 depth-bin contributions for each element before merging.
# The raw detail tables contain one row per depth bin, not one row per element.
vac_detail = (
    vac_detail.groupby(
        ["material", "energy_keV", "layer", "element"],
        as_index=False
    )["integrated_value"].sum()
    .rename(columns={"integrated_value": "vacancies_per_ion"})
)
rec_detail = (
    rec_detail.groupby(
        ["material", "energy_keV", "layer", "element"],
        as_index=False
    )["integrated_value"].sum()
    .rename(columns={"integrated_value": "recoil_energy_eV_per_ion"})
)

# Combine elemental values. A layer is retained so the bilayer's two materials
# are not incorrectly merged into a single elemental contribution.
elements = vac_detail.merge(
    rec_detail,
    on=["material", "energy_keV", "layer", "element"],
    how="outer",
    validate="one_to_one"
)
elements["vacancies_per_ion"] = pd.to_numeric(elements["vacancies_per_ion"], errors="coerce")
elements["recoil_energy_eV_per_ion"] = pd.to_numeric(elements["recoil_energy_eV_per_ion"], errors="coerce")
elements["vacancy_fraction_pct"] = (
    elements["vacancies_per_ion"] /
    elements.groupby(["material", "energy_keV"])["vacancies_per_ion"].transform("sum") * 100
)
elements["recoil_energy_fraction_pct"] = (
    elements["recoil_energy_eV_per_ion"] /
    elements.groupby(["material", "energy_keV"])["recoil_energy_eV_per_ion"].transform("sum") * 100
)
elements.to_csv(PROCESSED_DIR / "combined_element_summary.csv", index=False)

# Integrate full depth profiles and independently compare their totals.
def total_from_profile(d, total_column, label):
    rows = []
    for (material, energy), g in d.groupby(["Material_key", "Energy_keV"]):
        g = g.sort_values("Depth_nm")
        if total_column not in g.columns:
            rows.append({
                "material": material, "energy_keV": int(energy),
                f"{label}_profile_total": np.nan,
                f"{label}_depth_bins": len(g)
            })
            continue
        vals = pd.to_numeric(g[total_column], errors="coerce")
        total = vals.sum() * BIN_WIDTH_A
        rows.append({
            "material": material, "energy_keV": int(energy),
            f"{label}_profile_total": total,
            f"{label}_depth_bins": len(g)
        })
    return pd.DataFrame(rows)

vac_element_totals = (
    vac_detail.groupby(["material", "energy_keV"], as_index=False)
    .agg(vacancies_from_elements_per_ion=("vacancies_per_ion", "sum"))
)
rec_element_totals = (
    rec_detail.groupby(["material", "energy_keV"], as_index=False)
    .agg(recoil_energy_from_elements_eV_per_ion=("recoil_energy_eV_per_ion", "sum"))
)

vac_profile_totals = total_from_profile(
    vac, "Total_vacancies_per_A_per_ion", "vacancy"
)
rec_total_candidates = [
    c for c in rec.columns
    if c.lower() in {
        "total_target_recoil_energy_ev_per_nm_per_proton",
        "total_target_recoil_energy_ev_per_nm_ion",
        "total_recoil_energy_ev_per_nm_ion",
    }
]
# The consolidated recoil total is per nm, so integrate over 10 nm bins.
if rec_total_candidates:
    rec_profile_totals = []
    col = rec_total_candidates[0]
    for (material, energy), g in rec.groupby(["Material_key", "Energy_keV"]):
        val = pd.to_numeric(g.sort_values("Depth_nm")[col], errors="coerce").sum() * BIN_WIDTH_NM
        rec_profile_totals.append({
            "material": material, "energy_keV": int(energy),
            "recoil_profile_total": val, "recoil_depth_bins": len(g)
        })
    rec_profile_totals = pd.DataFrame(rec_profile_totals)
else:
    rec_profile_totals = pd.DataFrame(columns=[
        "material", "energy_keV", "recoil_profile_total", "recoil_depth_bins"
    ])

summary = vac_element_totals.merge(rec_element_totals, on=["material", "energy_keV"], validate="one_to_one")
summary = summary.merge(vac_profile_totals, on=["material", "energy_keV"], how="left")
summary = summary.merge(rec_profile_totals, on=["material", "energy_keV"], how="left")

# Bring in independent TDATA/master values if the file is available.
if MASTER_FILE.exists():
    master = pd.read_csv(MASTER_FILE)
    master.columns = [str(c).strip() for c in master.columns]
    if {"Material", "Energy_keV"}.issubset(master.columns):
        master["material"] = master["Material"].map(canonical_material)
        master["energy_keV"] = pd.to_numeric(master["Energy_keV"], errors="coerce").astype(int)
        keep = ["material", "energy_keV"] + [
            c for c in [
                "Vacancies_per_ion", "Mean_Range_nm",
                "Longitudinal_Straggling_nm", "Backscattered_fraction_pct",
                "Transmitted_fraction_pct", "Retained_fraction_pct"
            ] if c in master.columns
        ]
        summary = summary.merge(master[keep], on=["material", "energy_keV"], how="left", validate="one_to_one")

# Ratios are descriptive only, not a material property or radiation-tolerance metric.
summary["recoil_energy_per_vacancy_eV"] = (
    summary["recoil_energy_from_elements_eV_per_ion"] /
    summary["vacancies_from_elements_per_ion"]
)
summary.to_csv(PROCESSED_DIR / "combined_material_energy_summary.csv", index=False)

# Validation report: elemental sums vs independent total profiles and TDATA.
checks = []
for _, r in summary.iterrows():
    for name, a, b in [
        ("vacancy_elements_vs_total_profile", r["vacancies_from_elements_per_ion"], r["vacancy_profile_total"]),
        ("recoil_elements_vs_total_profile", r["recoil_energy_from_elements_eV_per_ion"], r["recoil_profile_total"]),
    ]:
        if pd.notna(a) and pd.notna(b) and b != 0:
            diff = abs(a-b)/abs(b)*100
            status = "PASS" if diff < 0.5 else "CHECK"
        else:
            diff, status = np.nan, "NOT_AVAILABLE"
        checks.append({
            "material": r["material"], "energy_keV": int(r["energy_keV"]),
            "check": name, "difference_percent": diff, "status": status
        })
    if "Vacancies_per_ion" in summary.columns and pd.notna(r.get("Vacancies_per_ion")):
        a = r["vacancies_from_elements_per_ion"]
        b = r["Vacancies_per_ion"]
        diff = abs(a-b)/abs(b)*100 if b else np.nan
        checks.append({
            "material": r["material"], "energy_keV": int(r["energy_keV"]),
            "check": "vacancy_elements_vs_TDATA", "difference_percent": diff,
            "status": "PASS" if pd.notna(diff) and diff < 0.5 else "CHECK"
        })
validation = pd.DataFrame(checks)
validation.to_csv(PROCESSED_DIR / "combined_validation_report.csv", index=False)

# Combine depth profiles in a common long table, preserving layer and element.
vp = vac_profiles.rename(columns={"vacancy_density_per_nm_ion": "vacancies_per_nm_ion"})
rp = rec_profiles.rename(columns={"recoil_energy_eV_per_nm_ion": "recoil_energy_eV_per_nm_ion"})
depth = vp.merge(
    rp,
    on=["material", "energy_keV", "layer", "element", "depth_nm"],
    how="outer",
    validate="one_to_one"
)
depth.to_csv(PROCESSED_DIR / "combined_depth_profiles.csv", index=False)

# ----------------------------- FIGURES -------------------------------
plt.rcParams.update({"figure.figsize": (10, 6), "axes.grid": True})

# 1. Scatter: total vacancies vs total recoil energy, with one point per run.
fig, ax = plt.subplots()
for material, g in summary.groupby("material"):
    ax.scatter(
        g["vacancies_from_elements_per_ion"],
        g["recoil_energy_from_elements_eV_per_ion"],
        label=material, s=75
    )
    for _, r in g.iterrows():
        ax.annotate(f'{int(r["energy_keV"])} keV',
                    (r["vacancies_from_elements_per_ion"], r["recoil_energy_from_elements_eV_per_ion"]),
                    xytext=(4, 4), textcoords="offset points", fontsize=8)
ax.set_xlabel("Total vacancies per incident proton")
ax.set_ylabel("Target recoil energy (eV per incident proton)")
ax.set_title("Vacancy production vs target recoil energy")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES_DIR / "combined_vacancies_vs_recoil_energy.png", dpi=300)
plt.close(fig)

# 2 and 3. Cross-material bars for each energy.
for metric, ylabel, prefix in [
    ("vacancies_from_elements_per_ion", "Vacancies per incident proton", "combined_vacancies"),
    ("recoil_energy_from_elements_eV_per_ion", "Target recoil energy (eV per incident proton)", "combined_recoil_energy"),
]:
    fig, ax = plt.subplots()
    pivot = summary.pivot(index="energy_keV", columns="material", values=metric).sort_index()
    pivot.plot(kind="bar", ax=ax)
    ax.set_xlabel("Incident proton energy (keV)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"Cross-material comparison: {ylabel}")
    ax.legend(title="Material", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"{prefix}_by_energy.png", dpi=300)
    plt.close(fig)

# 4 and 5. Element contribution plots, grouped by material and energy.
for metric, ylabel, filename in [
    ("vacancies_per_ion", "Vacancies per incident proton", "combined_element_vacancy_contributions.png"),
    ("recoil_energy_eV_per_ion", "Recoil energy (eV per incident proton)", "combined_element_recoil_contributions.png"),
]:
    plot_data = elements.copy()
    plot_data["species"] = plot_data["element"] + np.where(
        plot_data["layer"].eq("Single"), "", " (" + plot_data["layer"] + ")"
    )
    pivot = plot_data.pivot_table(
        index=["material", "energy_keV"], columns="species",
        values=metric, aggfunc="sum", fill_value=0
    )
    ax = pivot.plot(kind="bar", stacked=True, figsize=(14, 7))
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Material and incident proton energy")
    ax.set_title(f"Element-resolved {ylabel.lower()}")
    ax.legend(title="Element / layer", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches="tight")
    plt.close()

print("\nCOMBINED ANALYSIS COMPLETE")
print(f"Matched simulations: {len(summary)}")
print(f"Element rows: {len(elements)}")
print(f"Depth-profile rows: {len(depth)}")
print("\nCombined summary:")
print(summary[[
    "material", "energy_keV", "vacancies_from_elements_per_ion",
    "recoil_energy_from_elements_eV_per_ion", "recoil_energy_per_vacancy_eV"
]].sort_values(["material", "energy_keV"]).to_string(index=False))
print("\nValidation status:")
print(validation.groupby(["check", "status"]).size().to_string())
print(f"\nProcessed files: {PROCESSED_DIR}")
print(f"Figures: {FIGURES_DIR}")
