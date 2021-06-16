# Lineage-specific protection and immune imprinting shape the age distributions of influenza B cases

Code and data for Vieira et al., 2021 (Lineage-specific protection and immune imprinting shape the age distributions of influenza B cases).

## Data
- `results/processed_data/`: processed data in the format directly used by the model. `case_data_nz_all_surveillance_untyped_assigned.csv` is the default dataset to which the model was fitted and includes sentinel and non-sentinel cases, with cases that had missing lineage information assigned to the dominant lineage in particular seasons. `case_data_nz_all_surveillance.csv` excludes cases with missing lineage information altogether. `lineage_frequencies_gisaid-genbank_noVicin1990s.csv` contains the default lineage frequencies, while `lineage_frequencies_gisaid-genbank.csv` contains the frequencies used in the sensitivity analysis of B/Yamagata's frequency in the 1990s. `netherlands_seroprevalence.csv` contains the cross-sectional serology data from children in the Netherlands.

- `data/gisaid_metadata/` and `data/genbank_data/`: Data from sequence databases used to estimate lineage frequencies.

- `data/sequence_divergence/`: Sequence data used to estimate divergence between influenza B lineages and influenza A subtypes.

- `data/surveillance_intensity/`: Data from surveillance reports used to estimate the intensity of influenza B circulation in different seasons.

- `data/surveillance_report_frequencies`: Historical frequencies of B/Victoria and B/Yamagata from surveillance reports.

- `data/demographic_data/`: Data on the age distribution of the general population in New Zealand and Australia.


## Code
The core analysis is executed by `likelihood_profile.R`, which takes paths to the fitting data and covariates along with the names of one or two parameters to profile over the specified intervals at the specified increments. For each parameter value (or combination of values) to profile, `likelihood_profile.R` then generates a specification file which is passed on to `fixed_pars_likelihood.R` via a sbatch file created by the script assuming a cluster with SLURM (detailed configurations are user-specific). `fixed_pars_likelihood.R` then uses package `optimParallel` to find the maximum likelihood combination of the remaining parameters.

Given a path to a directory, `combine_likprofile_csvs.R` looks for a `likelihood_profiles/` folder containing subfolders with likelihood profiles of different parameters (one csv file for each parameter value/combination) and combines them into `likelihood_profiles/combined_likelihood_profiles.csv`. 

`plot_model_fits.R` searches for the global maximum likelihood parameter estimates and plots the model predictions under the MLE. It requires that most of the inputs used by `likelihood_profile.R` be collected into an input file assumed to be in the directory containing `combined_likelihood_profiles.csv`. An example input file is given in the comments.

### Note: 
While the 'VY' subscript indicates protection from B/Yamagata against B/Victoria in the paper, variables with the `_VY` in the code and its output refer to protection from B/Victoria against B/Yamagata (and similarly for other cross-lineage protection subscripts).

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

