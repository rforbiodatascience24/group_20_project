# Group 20 Project

## Project Contributors

Lea Eschen Skytthe (leasn-bot), s203531, contact: s203531@dtu.dk

Trine Søgaard (trinesoe), s204655, contact: triso@dtu.dk

Amy Surry (AmySurry), s201902, contact: s201902@dtu.dk

Olivier Gaufres (nollive), s243252, contact: s243252@dtu.dk

Antoine Andreoletti (apollis44), s242830, contact: s242830@dtu.dk

## Describtion
This project investigates differences in hepatic gene expression between healthy individuals, patients with non-alcoholic fatty liver disease (NAFLD), and patients with cirrhosis. Using RNA-sequence of coding RNA, we analyze differences between these groups.

## Project Structure
- **Data Setup and Download:** `01_load` creates a `data/` directory, with the subdirectories `data/_raw/count_data/`, and `data/_raw/metadata/`. It downloads the raw data from the specified URL (found under `Data`) and saves it into the subfolders.

- **Data Cleaning:** `02_clean`

- **Data Augmentation:** `03_augment`

- **Data Description:** `04_describe`: Generates descriptive plots of data 

- **Data Analysis:**
    `05_analysis_1`: Performs PCA
    `06_analysis_2`: Makes heatmap
    `07_analysis_3`: Performs Differential Expression Analysis (DEA) using DESeq2. Does also perform PCA and generate heatmap & volcano plots
    `08_analysis_4`: Performs Gene Set Enrichment Analysis (GSEA) and plot the 5 most significantly enriched pathways

- **All:** `00_all`: Runs the project in one go

## Data
Data is available here: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12807

## Link to presentation
XXXXX
