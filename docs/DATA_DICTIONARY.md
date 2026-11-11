# Data dictionary

## Input files

### `vacancy_depth_profiles_complete.csv`
Consolidated SRIM vacancy depth profiles for 16 material-energy combinations. The raw elemental vacancy-density columns retain SRIM units of vacancies/(Å·incident ion), unless otherwise stated in the source header. Depth is provided in nm.

### `recoil_energy_depth_profiles_complete.csv`
Consolidated E2RECOIL depth profiles for the same 16 combinations. Elemental recoil-energy density columns retain SRIM units of eV/(Å·incident ion). The `ENERGY from IONS` quantity, where present, is separate from elemental recoil energy. Check column headers before reprocessing.

## Processed files

- `combined_material_energy_summary.csv`: integrated vacancy and target recoil-energy totals by material and incident energy.
- `combined_element_summary.csv`: integrated elemental contributions and fractions.
- `combined_depth_profiles.csv`: combined element/layer/depth records.
- `combined_validation_report.csv`: comparisons between elemental sums and total profiles.
- `continuous_profile_integrated_comparison.csv`: total integrated values for each material and energy.
- `continuous_total_depth_profiles.csv`: total vacancy and target recoil-energy profiles by depth.
- `bilayer_layer_damage_partition.csv`: layer-resolved integrated values and percentage shares.
- `bilayer_layer_depth_profiles.csv`: separate MAPbI3 and FAPbI3 layer profiles.

All integrated vacancy and recoil-energy values are normalized per incident proton. They are not totals over all 10,000 incident protons.
