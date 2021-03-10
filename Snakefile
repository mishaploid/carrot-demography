################################################################################
## Demographic analysis for carrot
## Sarah Turner-Hissong
################################################################################

# Setup
import pandas as pd

# specify config file (contains file paths, etc.)
configfile: "config.yaml"

# submit.sh - bash script to run workflow
# submit.json - configuration file for cluster settings

################################################################################
## List of chromosome names
################################################################################
CHR = ['DCARv3_Chr' + str(n) for n in range(1,10)]

################################################################################
## Distinguished individuals to use for SMC++ input files
## These are 10 randomly sampled individuals from each population
## Can be used to estimate a composite likelihood and uncertainty in the recent past
################################################################################
dist_ind = pd.read_table(config['distinguished_inds'])

dist_list = dist_ind.groupby('Population')['Sample_ID'].apply(lambda x: x.tolist())

dist_dict = dist_list.to_dict()

# create a list of filenames for smc++ input files using different distinguished individuals
smc_input_files = [expand('models/smc/input/{population}.{distinguished_ind}.{chr}.smc.gz',
                    population = key, distinguished_ind = value, chr = CHR)
                    for key, value in dist_dict.items()]

# create a list of filenames to feed into smc++ cv command
# want to include all .smc.gz files for each population (includes all chromosomes and distinguished individuals)
def smc_cv_input(wildcards):
    files = expand('models/smc/input/{population}.{distinguished_ind}.{chr}.smc.gz', chr=CHR, population=wildcards.population, distinguished_ind=dist_dict[wildcards.population])
    return files

################################################################################
## Create a dictionary of sample names and population ids for SMC++ estimate
## https://stackoverflow.com/questions/35029731/convert-two-columns-pandas-data-frame-to-dictionary-of-list-with-first-column-as
################################################################################

# read in list of sample and population ids
samp_ids = pd.read_table(config['samp_ids'])

# create a list of samples for each population
samp_list = samp_ids.groupby('Population')['Sample_ID'].apply(lambda x: x.tolist())

# convert sample list to a python dictionary object
popdict = samp_list.to_dict()

# join sample names and separate them with a ','
from collections import defaultdict
for key in popdict.keys():
    popdict[key] = ','.join(popdict[key])

# format input for SMC++
# Population:sample1,sample2,sample3...samplen
for key, value in popdict.items() :
    popdict[key] = key + ':' + value

# function to print formatted list of population: sample ids using wildcards
def pop_choose(wildcards):
	list = popdict[wildcards.population]
	return list

################################################################################
## Create a dictionary of paired populations for SMC++ split
################################################################################

# # dictionary of population pairs to estimate divergence
pop_pair_dict = {'Wild_CentralAsianCultivated':['Wild', 'Central_Asian_Cultivated'],
'Wild_EasternCultivated':['Wild', 'Eastern_Cultivated'],
'Wild_WesternCultivated':['Wild', 'Western_Cultivated'],
'Wild_Cultivar':['Wild', 'Cultivar'],
'CentralAsian_Eastern':['Central_Asian_Cultivated', 'Eastern_Cultivated'],
'CentralAsian_Western':['Central_Asian_Cultivated', 'Western_Cultivated'],
'CentralAsian_Cultivar':['Central_Asian_Cultivated', 'Cultivar'],
'Eastern_Western':['Eastern_Cultivated', 'Western_Cultivated'],
'Eastern_Cultivar':['Eastern_Cultivated', 'Cultivar'],
'Western_Cultivar':['Western_Cultivated', 'Cultivar']}

# 'Wild_Cultivar':['Wild', 'Cultivar'],


# def pop_pair_choose(wildcards):
# 	list = popdict[wildcards.pop_pair]
# 	return list

# for split time:
# based on current wildcard, find pop pair in pop_pair_dict
# for each population in current pop pair, find pop in popdict
# combine strings to create a super long string
# return the super long string

def pair_string_choose12(wildcards):
    pops = pop_pair_dict[wildcards.pop_pair]
    pop1 = popdict[pops[0]]
    pop2 = popdict[pops[1]]
    pop_pair_string12 = pop1 + " " + pop2
    return pop_pair_string12

def pair_string_choose21(wildcards):
    pops = pop_pair_dict[wildcards.pop_pair]
    pop1 = popdict[pops[0]]
    pop2 = popdict[pops[1]]
    pop_pair_string21 = pop2 + " " + pop1
    return pop_pair_string21

# # create string: location of model.final.json for each pop

def smc_split_input(wildcards):
    pops = pop_pair_dict[wildcards.pop_pair]
    models = expand("models/smc_cv/{population}/model.final.json", population = pops)
    pop1_input = expand("models/smc/input/{population}.{distinguished_ind}.{chr}.smc.gz", chr=CHR, population=pops[0], distinguished_ind=dist_dict[pops[0]])
    pop2_input = expand("models/smc/input/{population}.{distinguished_ind}.{chr}.smc.gz", chr=CHR, population=pops[1], distinguished_ind=dist_dict[pops[1]])
    joint_input = expand("models/smc_split/input/{pop_pair}{order}.{chr}.smc.gz", pop_pair = wildcards.pop_pair, order = [12, 21], chr = CHR)
    return models + pop1_input + pop2_input + joint_input


################################################################################
## Define rule all
## This is a pseudo-rule that collects the target files (expected outputs)
################################################################################

rule all:
	input:
		vcf2smc = smc_input_files,
		smc_cv = expand("models/smc_cv/{population}/model.final.json", population = popdict.keys()),
		plot_estimate = "reports/smc_cv_results.png",
        joint_vcf2smc = expand("models/smc_split/input/{pop_pair}{order}.{chr}.smc.gz", pop_pair = pop_pair_dict.keys(), order = [12, 21], chr = CHR),
        smc_split = expand("models/smc_split/{pop_pair}/model.final.json", pop_pair = pop_pair_dict.keys()),
        plot_split = expand("reports/figures/{pop_pair}.split.png", pop_pair = pop_pair_dict.keys())

################################################################################
## Rule files to include
################################################################################

include: "rules/01-population_size_history.smk"
include: "rules/02-population_splits.smk"
