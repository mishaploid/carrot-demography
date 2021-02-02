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
## Create a dictionary of sample names and population ids for SMC++ estimate
## https://stackoverflow.com/questions/35029731/convert-two-columns-pandas-data-frame-to-dictionary-of-list-with-first-column-as
################################################################################

# read in list of sample and population ids
samp_ids = pd.read_table(config['samp_ids'])

# create a list of samples for each population
samp_list = samp_ids.groupby('Population')['Sample_ID'].apply(lambda x: x.tolist())

# convert sample list to a python dictionary object
popdict = samp_list.to_dict()
print(popdict)

# join sample names and separate them with a ','
from collections import defaultdict
for key in popdict.keys():
    popdict[key] = ','.join(popdict[key])

# format input for SMC++
# Population:sample1,sample2,sample3...samplen
for key, value in popdict.items() :
    popdict[key] = key + ':' + value

# function to print formatted list of population: sample ids using wildcards (WC)
def pop_choose(WC):
	list = popdict[WC.population]
	return list

################################################################################
## Create a dictionary of paired populations for SMC++ split
################################################################################

# # dictionary of population pairs to estimate divergence
pop_pair_dict = {'Wild_CentralAsianCultivated':['Wild', 'Central_Asian_Cultivated']}
# pop_pair_dict = {'botrytis_italica':['botrytis', 'italica'],
# 'italica_botrytis':['italica', 'botrytis'],
# 'italica_wild':['italica', 'wild'],
# 'wild_italica':['wild', 'italica']}

def pop_pair_choose(WC):
	list = popdict[WC.pop_pair]
	return list

# # for split time:
# # based on current wildcard, find pop pair in pop_pair_dict
# # for each population in current pop pair, find pop in popdict
# # combine strings to create a super long string
# # return the super long string
#
def pair_string_choose(WC):
	pops = pop_pair_dict[WC.pop_pair]
	pop1 = popdict[pops[0]]
	pop2 = popdict[pops[1]]
	pop_pair_string = pop1 + " " + pop2
	return pop_pair_string
#
# # create string: location of model.final.json for each pop
#
# def model_chooser(WC):
# 	pops = pop_pair_dict[WC.pop_pair]
# 	pop1 = pops[0]
# 	pop2 = pops[1]
# 	model_string = "models/smc/cv_1e3_1e6/" + pop1 + "/model.final.json " \
# 	+ "models/smc/cv_1e3_1e6/" + pop2 + "/model.final.json"
# 	return model_string


################################################################################
## Define rule all
## This is a pseudo-rule that collects the target files (expected outputs)
################################################################################

rule all:
	input:
		vcf2smc = expand("models/smc/input/{population}.{chr}.smc.gz", population = popdict.keys(), chr = CHR),
		smc_cv = expand("models/smc/{population}/fold{fold}/model.final.json", population = popdict.keys(), fold = ['0','1']),
		plot_estimate = "reports/smc_cv_results.png",
        joint_vcf2smc = expand("models/smc/split_input/{pop_pair}.{chr}.smc.gz", pop_pair = pop_pair_dict.keys(), chr = CHR)

################################################################################
## Rule files to include
################################################################################

include: "rules/01-population_size_history.smk"
include: "rules/02-population_splits.smk"
