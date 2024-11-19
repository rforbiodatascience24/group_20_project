
### METADATA - Clean-up

# Clean environment
rm(list=ls())


# Load library
library("tidyverse")


# Read data
data <- read.delim('data/_raw/meta_data/E-MTAB-12807.sdrf.txt', sep = "\t")


# Check for differences in columns - we only want the columns that shows differences between patients
filtered_data <- data |> 
  select(where(~ n_distinct(.) > 1))


# Remove duplicates in rows based on source name
filtered_data <- filtered_data |> 
  distinct(Source.Name, .keep_all = TRUE) 


# Remove assay name & extract name as they are the same as source name
# Furthermore, we remove columns we do not intend to use
filtered_data <- filtered_data |> 
  select(-Extract.Name, 
         -Assay.Name,
         -Scan.Name,
         -Factor.Value.disease.,
         -starts_with("Comment."))


#Rename the columns to more meaningful names
filtered_data <- filtered_data |> 
  rename(
    source_name = Source.Name,
    individual = Characteristics.individual.,
    age = Characteristics.age.,
    sex = Characteristics.sex.,
    disease = Characteristics.disease.,
    data_file = Derived.Array.Data.File,
    clinical_information = Factor.Value.clinical.information.
  )


# Change datafile .sf to .csv
filtered_data <- filtered_data |> 
  mutate(data_file = sub("\\.sf$", ".csv", data_file))
  












