# Lineage-specific protection and immune imprinting shape the age distributions of influenza B cases

Code and data for Vieira et al., 2020 (Lineage-specific protection and immune imprinting shape the age distributions of influenza B cases).

## Data
- `results/processed_data/`: processed data in the format directly used by the model. `case_data_nz_all_surveillance_untyped_assigned.csv` is the default dataset to which the model was fitted and includes sentinel and non-sentinel cases, with cases that had missing lineage information assigned to the dominant lineage in particular seasons. `case_data_nz_all_surveillance.csv` excludes cases with missing lineage information altogether. `lineage_frequencies_gisaid-genbank_noVicin1990s.csv` contains the default lineage frequencies, while `lineage_frequencies_gisaid-genbank.csv` contains the frequencies used in the sensitivity analysis of B/Yamagata's frequency in the 1990s. `netherlands_seroprevalence.csv` contains the cross-sectional serology data from children in the Netherlands.

- `data/gisaid_metadata/` and `data/genbank_data/`: Data from sequence databases used to estimate lineage frequencies.

- `data/sequence_divergence/`: Sequence data used to estimate divergence between influenza B lineages and influenza A subtypes.

- `data/surveillance_intensity/`: Data from surveillance reports used to estimate the intensity of influenza B circulation in different seasons.

- `data/surveillance_report_frequencies`: Historical frequencies of B/Victoria and B/Yamagata from surveillance reports.

- `data/demographic_data/`: Data on the age distribution of the general population in New Zealand and Australia.


## Code






## Dependencies
The analyses were run using R version 3.4.3 and the following packages
- dplyr
- tidyr
- stringr
- optimParallel
- lubridate
- FlexParamCurve
- ggplot2
- cowplot
- reshape2 
- viridis
- ggridges

