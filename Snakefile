################################################################################
## Demographic analysis for carrot
## Sarah Turner-Hissong
################################################################################

# Setup
import pandas as pd
import click

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
# import pandas as pd
# samp_ids = pd.read_table("data/samp_ids.txt")
# dist_ind = samp_ids.groupby('Population', as_index = False).apply(pd.DataFrame.sample, n = 10, replace = False)
# dist_ind.to_csv('data/distinguished_individuals.txt', header = True, index = False, sep = '\t')
################################################################################
dist_ind = pd.read_table(config['distinguished_inds'])

dist_ind = dist_ind[dist_ind.Population != 'LandraceAWild']
# dist_ind = dist_ind[(dist_ind.Population == 'Landrace_A') | (dist_ind.Population == 'Landrace_B')]

dist_list = dist_ind.groupby('Population')['Sample_ID'].apply(lambda x: x.tolist())

dist_dict = dist_list.to_dict()


# create a list of filenames for smc++ input files using different distinguished individuals (output of vcf2smc)
smc_input_files = [expand('models/smc_estimate_input/{population}.{distinguished_ind}.{chr}.smc.gz',
                    population = key, distinguished_ind = value, chr = CHR)
                    for key, value in dist_dict.items()]

# for bootstrapped inputs
smc_bootstrap_input = [expand('models/smc_estimate_bootstrap_input/{population}_{distinguished_ind}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz',
                    population = key, distinguished_ind = value, n_bootstrap = range(1,11), boot_chr = range(1,10))
                    for key, value in dist_dict.items()]

# create a list of filenames to feed into smc++ cv command
# want to include all .smc.gz files for each population (includes all chromosomes and distinguished individuals)
def smc_estimate_input(wildcards):
    files = expand('models/smc_estimate_input/{population}.{distinguished_ind}.{chr}.smc.gz', chr=CHR, population=wildcards.population, distinguished_ind=dist_dict[wildcards.population])
    return files

def smc_estimate_boot_input(wildcards):
    files = expand('models/smc_estimate_bootstrap_input/{population}_{distinguished_ind}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz',  population=wildcards.population, distinguished_ind=dist_dict[wildcards.population],
    n_bootstrap = range(1,11),
    boot_chr = range(1,10))
    return files

################################################################################
## Create a dictionary of sample names and population ids for SMC++ estimate
## https://stackoverflow.com/questions/35029731/convert-two-columns-pandas-data-frame-to-dictionary-of-list-with-first-column-as
################################################################################

# read in list of sample and population ids
samp_ids = pd.read_table(config['samp_ids'])

samp_ids = samp_ids[samp_ids.Population != 'LandraceAWild']

# create a list of samples for each population
samp_list =  samp_ids.groupby('Population')['Sample_ID'].apply(lambda x: x.tolist())

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

pop_pair_dict = {'LandraceA_LandraceB':['Landrace_A', 'Landrace_B'],
'LandraceA_Wild':['Landrace_A', 'Wild'],
'LandraceA_EarlyCultivar':['Landrace_A', 'Early_Cultivar'],
'LandraceA_ImprovedCultivar':['Landrace_A', 'Improved_Cultivar'],
'EarlyCultivar_ImprovedCultivar':['Early_Cultivar', 'Improved_Cultivar'],
'EarlyCultivar_Wild':['Early_Cultivar', 'Wild'],
'ImprovedCultivar_Wild':['Improved_Cultivar', 'Wild'],
'LandraceB_EarlyCultivar':['Landrace_B', 'Early_Cultivar'],
'LandraceB_ImprovedCultivar':['Landrace_B', 'Improved_Cultivar'],
'LandraceB_Wild':['Landrace_B', 'Wild']}
# 'LandraceA_LAwild':['Landrace_A', 'LandraceAWild'],

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

# # create a list of filenames for smc++ input files using different distinguished individuals
smc_split_input_files12 = [expand('models/smc_split_input/{pop_pair}_12.{distinguished_ind1}.{chr}.smc.gz',
                    pop_pair = key, distinguished_ind1 = dist_dict[value[0]], chr = CHR)
                    for key, value in pop_pair_dict.items()]

smc_split_input_files21 = [expand('models/smc_split_input/{pop_pair}_21.{distinguished_ind2}.{chr}.smc.gz',
                    pop_pair = key, distinguished_ind2 = dist_dict[value[1]], chr = CHR)
                    for key, value in pop_pair_dict.items()]


# list of filenames for split bootstrapping
smc_split_bootstrap_input_12 = [expand('models/smc_split_bootstrap_input/{pop_pair}_12.{distinguished_ind1}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz',
                    pop_pair = key, distinguished_ind1 = dist_dict[value[0]], n_bootstrap = range(1,11), boot_chr = range(1,10))
                    for key, value in pop_pair_dict.items()]

smc_split_bootstrap_input_21 = [expand('models/smc_split_bootstrap_input/{pop_pair}_21.{distinguished_ind2}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz',
                    pop_pair = key, distinguished_ind2 = dist_dict[value[1]], n_bootstrap = range(1,11), boot_chr = range(1,10))
                    for key, value in pop_pair_dict.items()]


# # create string: location of model.final.json for each pop

def smc_split_input(wildcards):
    pops = pop_pair_dict[wildcards.pop_pair]
    models = expand("models/smc_estimate_no_timepoints/{population}/model.final.json", population = pops)
    pop1_input = expand("models/smc_estimate_input/{population}.{distinguished_ind}.{chr}.smc.gz", population=pops[0], distinguished_ind=dist_dict[pops[0]], chr=CHR)
    pop2_input = expand("models/smc_estimate_input/{population}.{distinguished_ind}.{chr}.smc.gz", population=pops[1], distinguished_ind=dist_dict[pops[1]], chr=CHR)
    joint_input12 = expand("models/smc_split_input/{pop_pair}_12.{distinguished_ind}.{chr}.smc.gz", pop_pair = wildcards.pop_pair, distinguished_ind=dist_dict[pops[0]], chr = CHR)
    joint_input21 = expand("models/smc_split_input/{pop_pair}_21.{distinguished_ind}.{chr}.smc.gz", pop_pair = wildcards.pop_pair, distinguished_ind=dist_dict[pops[1]], chr = CHR)
    return models + pop1_input + pop2_input + joint_input12 + joint_input21

def smc_split_bootstrap_input(wildcards):
    pops = pop_pair_dict[wildcards.pop_pair]
    models = expand("models/smc_estimate_no_timepoints_bootstrap/{population}_{n_bootstrap}/model.final.json", population = pops, n_bootstrap = wildcards.n_bootstrap)
    pop1_input = expand('models/smc_estimate_bootstrap_input/{population}_{distinguished_ind}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz',  population=pops[0], distinguished_ind=dist_dict[pops[0]],
    n_bootstrap = wildcards.n_bootstrap, boot_chr = range(1,10))
    pop2_input = expand('models/smc_estimate_bootstrap_input/{population}_{distinguished_ind}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz',  population=pops[1], distinguished_ind=dist_dict[pops[1]],
    n_bootstrap = wildcards.n_bootstrap, boot_chr = range(1,10))
    joint_input12 = expand("models/smc_split_bootstrap_input/{pop_pair}_12.{distinguished_ind}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz", pop_pair = wildcards.pop_pair, distinguished_ind=dist_dict[pops[0]], n_bootstrap = wildcards.n_bootstrap, boot_chr = range(1,10))
    joint_input21 = expand("models/smc_split_bootstrap_input/{pop_pair}_21.{distinguished_ind}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz", pop_pair = wildcards.pop_pair, distinguished_ind=dist_dict[pops[1]], n_bootstrap = wildcards.n_bootstrap, boot_chr = range(1,10))
    return models + pop1_input + pop2_input + joint_input12 + joint_input21


################################################################################
## Define rule all
## This is a pseudo-rule that collects the target files (expected outputs)
################################################################################

rule all:
	input:
		vcf2smc = smc_input_files,
		# smc_cv = expand("models/smc_cv_no_timepoints/{population}/model.final.json", population = popdict.keys()),
		smc_estimate = expand("models/smc_estimate_no_timepoints/{population}/model.final.json", population = dist_dict.keys()),
		plot_estimate = "reports/carrot_all_pops_smc_estimate_no_timepoints.png",
        smc_bootstrap = smc_bootstrap_input,
        # smc_cv_bootstrap = expand("models/smc_cv_bootstrap/{population}_{n_bootstrap}/model.final.json", population = popdict.keys(), n_bootstrap = range(1,11)),
        joint_vcf2smc12 = smc_split_input_files12,
        joint_vcf2smc21 = smc_split_input_files21,
        smc_split = expand("models/smc_split/{pop_pair}/model.final.json", pop_pair = pop_pair_dict.keys()),
        joint_bootstrap_vcf2smc12 = smc_split_bootstrap_input_12,
        joint_bootstrap_vcf2smc21 = smc_split_bootstrap_input_21,
        smc_split_bootstrap = expand("models/smc_split_bootstrap/{pop_pair}_{n_bootstrap}/model.final.json", pop_pair = pop_pair_dict.keys(), n_bootstrap = range(1,11))
        # plot_split = expand("reports/figures/{pop_pair}.split.png", pop_pair = pop_pair_dict.keys())

################################################################################
## Rule files to include
################################################################################

include: "rules/01-population_size_history.smk"
include: "rules/02-population_splits.smk"
