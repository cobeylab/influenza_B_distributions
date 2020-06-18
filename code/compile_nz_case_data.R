library(dplyr)
library(tidyr)
library(stringr)

# Read nz data
nz_data <- lapply(as.list(list.files('../results/processed_data/case_data_batches/', full.names = T)),
                  FUN = read.csv)
names(nz_data) <- str_replace(list.files('../results/processed_data/case_data_batches/', full.names = F),'\\.csv','')

# Compile data from both surveillance types combined, separately with and without assignment of unidentified cases
bind_rows(nz_data$`case_data_nz_2001-2012_all_surveillance`,
          nz_data$`case_data_nz_2012-2019_all_surveillance`) %>%
  write.csv('../results/processed_data/case_data_nz_all_surveillance.csv', row.names = F)

bind_rows(nz_data$`case_data_nz_2001-2012_all_surveillance_untyped_assigned`,
          nz_data$`case_data_nz_2012-2019_all_surveillance_untyped_assigned`) %>%
  write.csv('../results/processed_data/case_data_nz_all_surveillance_untyped_assigned.csv', row.names = F)

# Data from sentinel cases only (meaning GP, since hospitalizations post 2012 are labelled sentinel by ESR)
# No overlap since 2012 data from second batch are all hospitalizations...
bind_rows(nz_data$`case_data_nz_2001-2012_sentinel_only_untyped_assigned`,
          nz_data$`case_data_nz_2012-2019_sentinel_only_untyped_assigned`) %>%
  write.csv('../results/processed_data/case_data_nz_sentinel_only_untyped_assigned.csv', row.names = F)

# Data from non-sentinel cases only (meaning hospitalizations, since hospitalizations post 2012 are labelled sentinel by ESR)
# 2012 non-sentinel cases from first batch removed to avoid double counting with 2012 from second batch.
bind_rows(nz_data$`case_data_nz_2001-2012_nonsentinel_only_untyped_assigned`,
          nz_data$`case_data_nz_2012-2019_nonsentinel_only_untyped_assigned`) %>%
  write.csv('../results/processed_data/case_data_nz_nonsentinel_only_untyped_assigned.csv', row.names = F)
