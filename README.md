# Monte Carlo Simulation of Proton-Induced Defect Generation in Hybrid Halide Perovskites

A computational study using **SRIM/TRIM-2013** and Python to examine proton transport, target recoil energy, and displacement-vacancy production in four hybrid halide perovskite absorber compositions.

## Materials and simulation matrix

| Material/configuration | Model composition | Geometry in baseline comparison | Reference |
|---|---|---|---|
| MAPbI3 | CH6N·PbI3 (elemental stoichiometry C1H6N1Pb1I3) | 1 µm | https://doi.org/10.1002/adma.201603326 |
| MAPbI3/FAPbI3 bilayer | MAPbI3 followed by FAPbI3 | 500 nm + 500 nm | https://doi.org/10.1016/j.jpcs.2025.113302 |
| Lead-reduced mixed Sn-Pb | (FASn)0.6(MAPb)0.4I3 | 1 µm | https://doi.org/10.1016/j.solmat.2021.111232 |
| Triple-cation | Cs0.05MA0.17FA0.83Pb(I0.83Br0.17)3 | 1 µm | https://doi.org/10.1039/c9ee00077a |

Incident proton energies: **10, 30, 50, and 100 keV**.

Baseline run settings: 10,000 incident protons, 0° incidence, full damage cascades, and a 1 µm total absorber thickness. The bilayer is irradiated from the MAPbI3 side.

## Observations

- Vacancy production generally increases from 10 to 50 keV in all four configurations.
- At 100 keV, the triple-cation model has the largest integrated vacancy count within the 1 µm domain (about 21.98 vacancies per incident proton).
- At 100 keV, the triple-cation model also has the largest integrated target recoil energy within the 1 µm domain (about 503.70 eV per incident proton).
- In the bilayer, the fraction of integrated damage assigned to FAPbI3 increases with incident energy. The current layer-resolved dataset gives approximately 72% of vacancies and 72% of target recoil energy in FAPbI3 at 100 keV.
