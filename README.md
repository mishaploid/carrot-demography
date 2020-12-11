# Carrot demography

This repo is a snakemake workflow to run demographic analysis using the SMC++ software for carrot.

## Setup

### Installing snakemake
Follow [instructions to install via conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

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
