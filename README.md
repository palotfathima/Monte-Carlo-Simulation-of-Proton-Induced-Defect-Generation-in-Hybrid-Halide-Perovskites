# Monte Carlo Simulation of Proton-Induced Defect Generation in Hybrid Halide Perovskites

A computational study using **SRIM/TRIM-2013** and Python to examine proton transport, target recoil energy, and displacement-vacancy production in four hybrid halide perovskite absorber compositions.

> **Scientific scope:** This repository reports SRIM/TRIM simulation outputs and derived comparisons. Vacancy production is used as a displacement-damage proxy. These simulations do not directly predict photovoltaic efficiency, device lifetime, defect recovery, or experimentally measured radiation tolerance.

## Materials and simulation matrix

The initial comparison contains four absorber configurations and four incident proton energies, for 16 simulations.

| Material/configuration | Model composition | Geometry in baseline comparison |
|---|---|---|
| MAPbI3 | CH6N·PbI3 (elemental stoichiometry C1H6N1Pb1I3) | 1 µm |
| MAPbI3/FAPbI3 bilayer | MAPbI3 followed by FAPbI3 | 500 nm + 500 nm |
| Lead-reduced mixed Sn-Pb | (FASn)0.6(MAPb)0.4I3 | 1 µm |
| Triple-cation | Cs0.05MA0.17FA0.83Pb(I0.83Br0.17)3 | 1 µm |

Incident proton energies: **10, 30, 50, and 100 keV**.

Baseline run settings: 10,000 incident protons, 0° incidence, full damage cascades, and a 1 µm total absorber thickness. The bilayer is irradiated from the MAPbI3 side.

The triple-cation formula is retained as reported in the source composition; its printed A-site coefficients sum to 1.05. Do not silently renormalize this formula.

## Repository structure

```text
.
├── data/
│   ├── input/
│   │   ├── vacancy_depth_profiles_complete.csv
│   │   └── recoil_energy_depth_profiles_complete.csv
│   └── processed/
├── docs/
│   └── DATA_DICTIONARY.md
├── figures/
├── scripts/
│   ├── 01_process_vacancy_profiles.py
│   ├── 02_process_recoil_energy.py
│   ├── 03_combined_analysis.py
│   └── 04_depth_resolved_analysis.py
├── srim_raw/
├── .gitignore
├── requirements.txt
└── README.md
```

`data/input/` contains the consolidated depth-profile datasets used by the Python scripts. `data/processed/` contains derived tables. `figures/` contains generated plots. `srim_raw/` is reserved for original SRIM output files and input setup files; add these only after checking file naming and provenance.

## Software requirements

- Windows 10/11 (the commands below use Windows Command Prompt)
- Python 3.10 or later
- SRIM/TRIM-2013 for reproducing the Monte Carlo simulations
- Python packages listed in `requirements.txt`

Install dependencies:

```bash
py -m pip install -r requirements.txt
```

## Reproducing the Python analysis

Run commands from the repository root directory.

```bash
py scripts\01_process_vacancy_profiles.py
py scripts\02_process_recoil_energy.py
py scripts\03_combined_analysis.py
py scripts\04_depth_resolved_analysis.py
```

The scripts expect the two consolidated CSV files to be in `data/input/`. Script 03 can optionally use `data/input/master_simulation_data.csv` to add independent TDATA summary fields. That master file is not included in this initial repository package; add it if you want that cross-check.

## Data processing and units

### Vacancy profiles

The SRIM vacancy profile is expressed as vacancies/(Å·incident ion). The analysis converts the depth-density values to vacancies/(nm·incident proton) by multiplying by 10. With 100 Å (10 nm) bins, the per-proton integrated value is calculated as:

\[
N_\mathrm{vac/ion} = \sum_i \rho_{\mathrm{vac},i}\,(100\ \mathrm{\AA})
\]

### Target recoil energy

The elemental E2RECOIL values are expressed as eV/(Å·incident ion). They are converted to eV/(nm·incident proton) by multiplying by 10, then integrated over 100 Å bins:

\[
E_\mathrm{recoil/ion} = \sum_i \rho_{\mathrm{recoil},i}\,(100\ \mathrm{\AA})
\]

The `ENERGY from IONS` column is a separate SRIM quantity and is not added to the elemental target-recoil total. Target recoil energy is **not** the same as total electronic plus nuclear stopping power.

### Bilayer

The MAPbI3 and FAPbI3 contributions are kept separate for layer-resolved analysis. The interface is at 500 nm from the incident surface. The combined bilayer profile is the sum of the two layer contributions at each depth bin.

## Current preliminary observations

These are descriptive summaries of the baseline 1 µm simulation dataset, not final claims of radiation tolerance.

- Vacancy production generally increases from 10 to 50 keV in all four configurations.
- At 100 keV, the triple-cation model has the largest integrated vacancy count within the 1 µm domain (about 21.98 vacancies per incident proton).
- At 100 keV, the triple-cation model also has the largest integrated target recoil energy within the 1 µm domain (about 503.70 eV per incident proton).
- In the bilayer, the fraction of integrated damage assigned to FAPbI3 increases with incident energy. The current layer-resolved dataset gives approximately 72% of vacancies and 72% of target recoil energy in FAPbI3 at 100 keV.

Use the processed CSV files for exact values. These observations should be revisited after the finite-thickness and displacement-energy sensitivity checks below.

## Important limitations and validation tasks

1. **Finite target thickness at 100 keV:** the original 1 µm simulations show substantial proton transmission. Integrated values therefore describe energy transfer and vacancy production inside the 1 µm domain, not necessarily the complete proton trajectory. A thickness-convergence study (e.g., 2 µm, then 3 µm if transmission remains appreciable) is still required.
2. **Displacement energies:** results depend on the element-specific displacement-energy values supplied to TRIM. Treat these as model inputs. Document their sources and perform sensitivity analysis before making strong comparisons.
3. **Material density and phase:** document the source and phase/film basis for each density. The FAPbI3 density and phase assumption should be checked and cited.
4. **Bilayer interface:** verify layer ordering, depth-bin assignment, and interface treatment against the original TRIM setup and raw outputs.
5. **Triple-cation stoichiometry:** the printed A-site coefficients sum to 1.05. Preserve the source formula and explicitly document how the elemental composition was entered in SRIM.
6. **Statistical uncertainty:** a single 10,000-ion run is not a set of independent experimental replicates. Do not plot error bars unless they are derived from an appropriate statistical method or repeat simulations.
7. **Interpretation:** SRIM vacancy counts are displacement-damage proxies. They do not establish defect survival, recombination, self-healing, carrier transport, PCE degradation, or device lifetime.

## Suggested citation

If you use SRIM/TRIM, cite the SRIM software and the relevant SRIM/TRIM documentation according to the citation instructions supplied with your installation. Also cite the primary literature used for material compositions, densities, and displacement-energy assumptions. Add full bibliographic entries in `docs/REFERENCES.md` before public release.

## License

No license has been selected. Until one is added, the repository is publicly viewable if uploaded to GitHub, but others do not automatically receive permission to reuse, modify, or redistribute the contents. Add a license only after confirming the rights to all included data, scripts, and figures.
