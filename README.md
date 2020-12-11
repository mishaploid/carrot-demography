# Carrot demography

This repo is a snakemake workflow to run demographic analysis using the SMC++ software for carrot.

## Key dependencies
* Snakemake (see setup)
* Singularity or docker (note that workflow is set up to access a docker image through singularity)

## Project organization

### data/
Data is not uploaded to git, but includes the following files:
* `kevinG.d3.vcf.gz` - VCF file for carrot accessions, filtered for quality but not for MAF  
* `kevinG.d3.vcf.gz.tbi` - tabix index for the VCF file (created using `tabix -p vcf kevinG.d3.vcf.gz`)  
* `Resequencing_population_assignments.xlsx` - Excel spreadsheet with population assignments based on ADMIXTURE
* `sample_ids.txt` - tab-delimited text file containing two columns with population and sample ID, generated from the Excel file above using the code described below  
* `DCARv3.4.all.out.gff.gz` - results from RepeatMasker to mask regions that are highly repetitive/ambiguously mapped  
* `masked_regions.bed.gz` - results from RepeatMasker reformatted into BED format (tab-delimited, chromosome, start position, stop position) using the code below.

### Snakefile
This is the workhorse file for Snakemake and tracks the expected output files in addition to including some python code to format input variables.

### rules/  
Folder containing rules files (`.smk`) that are called by snakemake. For this project, `rules/run_smc.smk` contains the code necessary to run SMC++ to estimate effective population size.  

### config.yaml
This file can be edited to adjust file paths for input files and parameters to run SMC++

### submit.json
Specifies cluster configuration, log files, and CPU/mem usage for each rule - needs to be adjusted depending on HPC environment. 

## Setup

### Installing snakemake
This project uses snakemake for workflow management. Follow [instructions to install via conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). If running outside of the farm2 cluster at UC-Davis, will need to edit the `submit.json` file accordingly.

```
conda create --name carrot-demography
source activate carrot-demography
conda install -c conda-forge mamba
mamba install -c bioconda snakemake
```

### A few data wrangling steps to get things formatted for SMC++:

1. **Format population and sample ids to create a python dictionary object**  
This step allows easier formatting to paste input strings of population/samples for SMC++. The input is an excel file with ADMIXTURE results and two sheets, and the desired format is a text file with two columns for Population and Sample_ID.

```
import pandas as pd

# Read in ADMIXTURE results
samp_info = pd.read_excel(r'data/Resequencing_population_assignments.xlsx', sheet_name = 'Ultralow_admixture')

# select population group and sample id & convert to string
samp_ids = samp_info[['Population', 'Sample_ID']].astype(str)

# replace spaces with '_'
samp_ids = samp_ids.replace(' ', '_', regex = True)

# export formatted population and sample ids
samp_ids.to_csv('data/sample_ids.txt', header = True, index = False, sep = '\t')

```

2. **Reformat RepeatMasker results as a BED mask file**  
The mask file indicates regions that should be excluded from analysis (e.g. highly repetitive, poorly mapped regions) and is necessary to distinguish from long runs of homozygosity. It should be BED format (tab delimited) with no headers and three columns: chromosome, start position, stop position.

```
import pandas as pd
import csv

mask = pd.read_csv('data/DCARv3.4.all.out.gff.gz', compression = 'gzip', skiprows = 3, sep = "\t", header = None)

# select chromosome and start/stop positions
mask = mask[[0,3,4]]

# sort by position
sorted_mask = mask.sort_values(by = [0,3,4])

# export formatted mask file in bed format
sorted_mask.to_csv('data/masked_regions.bed', header = False, index = False, sep = '\t', quoting = csv.QUOTE_NONE)

```

In the terminal, index the mask file:
```
bgzip data/masked_regions.bed
tabix -p bed data/masked_regions.bed.gz
```
