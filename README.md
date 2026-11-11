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

- Vacancy production per proton rises from 10 to 50 keV in all compositions, reflecting the growing share of the proton’s path that lies inside the film. At 100 keV, MAPbI3, the bilayer and the Sn-Pb absorber fall back (18.5, 17.5 and 18.7 vacancies per proton). This fall-off is a thickness effect, not a change in material response. At 100 keV, 65% of protons are transmitted through MAPbI3, and the damage profiles are still rising at the rear surface. The triple-cation film is the exception (22.0 vacancies per proton), because its higher density (4.59 g/cm³) stops more ions (29% transmitted). The 100 keV ranking therefore depends on film thickness and would change for a different absorber thickness.

- For ≤50 keV, composition has a small effect. At 10 keV, where all ions stop in the film, the integrated target recoil energy is nearly the same for all four materials (~266–268 eV per proton). Composition changes where the damage sits rather than how much there is: the denser triple-cation film has a higher, shallower peak (~1.77 vs ~1.5 eV/nm) and a damage profile that ends ~100 nm earlier at 50 keV. Across all materials and energies, vacancy yield scales with recoil energy at ~0.04–0.046 vacancies per eV, so vacancy production is governed mainly by the nuclear energy deposited. The Sn-Pb absorber shows 6–8% fewer vacancies than MAPbI3 at 10 and 30 keV, which suggests modestly lower damage production.

- In MAPbI3 at 10 keV, iodine accounts for ~43% of vacancies and hydrogen ~27%, together ~70%, while Pb contributes ~16% (~20% at 50 keV). The dominant inorganic defect channel is therefore iodine displacement, consistent with the importance of iodide vacancies in perovskite ion migration.

- At 10 keV the bilayer is identical to MAPbI3, because all damage lies in the first layer. At 30 and 50 keV the total vacancy production agrees with MAPbI3 within 0.3%, so the interface redistributes damage without changing the total. At 100 keV the bilayer yields 5% fewer vacancies. However, the vacancy profile shows a sharp drop just past the interface that is absent from the recoil-energy profile, which points to a layer-boundary artefact. 
