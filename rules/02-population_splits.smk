# Estimating divergence times with SMC++ 
# Check out SMC++ git repo for additional documentation:
#   https://github.com/popgenmethods/smcpp

################################################################################
# STEP 1
# Create SMC++ input files for the joint frequency spectrum of two populations
#   --mask = poorly mapped regions to exclude from analysis
#            otherwise assumes missing data = regions of homozygosity
#   -d = distinguished individual for first pop listed
#   argument order:
#       mask, distinguished individual, input vcf, output directory, chromosome, 
#       population A: list of individuals, population B: list of individuals
################################################################################

# PopA_PopB
rule joint_vcf2smc12:
    input:
        vcf = config['vcf'],
        index = config['vcf_index'],
        mask = config['mask']
    output:
        out12 = "models/smc_split_input/{pop_pair}_12.{distinguished_ind1}.{chr}.smc.gz",
    params:
    	chrom = "{chr}",
        distind1 = "{distinguished_ind1}",
    	pop_pair_string12 = pair_string_choose12
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        """
        smc++ vcf2smc \
        --mask {input.mask} \
        -d {params.distind1} {params.distind1} \
        {input.vcf} {output.out12} {params.chrom} {params.pop_pair_string12}
        """

# PopB_PopA
rule joint_vcf2smc21:
    input:
        vcf = config['vcf'],
        index = config['vcf_index'],
        mask = config['mask']
    output:
        out21 = "models/smc_split_input/{pop_pair}_21.{distinguished_ind2}.{chr}.smc.gz"
    params:
    	chrom = "{chr}",
        distind2 = "{distinguished_ind2}",
        pop_pair_string21 = pair_string_choose21
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        """
        smc++ vcf2smc \
        --mask {input.mask} \
        -d {params.distind2} {params.distind2} \
        {input.vcf} {output.out21} {params.chrom} {params.pop_pair_string21}
        """

################################################################################
# STEP 2
# Estimate a joint demography for two populations
# inputs include marginal estimates of demography from each population,
# and smc.gz input files for each individual population and for the joint 
# frequency spectrum 
################################################################################

rule smc_split:
    input:
        smc_split_input
    params:
        model_out_dir = "models/smc_split/{pop_pair}/"
    output:
        model_out = "models/smc_split/{pop_pair}/model.final.json"
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        "smc++ split \
        --cores 8 \
        -o {params.model_out_dir} \
        {input}"

################################################################################
# STEP 3
# Bootstrapping
#   uses smc.gz inputs split into chunks and resampled 
################################################################################

rule joint_bootstrap_vcf2smc12:
    input:
        expand("models/smc_split_input/{{pop_pair}}_12.{{distinguished_ind1}}.{chr}.smc.gz", chr = CHR)
    output:
        expand('models/smc_split_bootstrap_input/{{pop_pair}}_12.{{distinguished_ind1}}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz', n_bootstrap = range(1,11), boot_chr = range(1,10)),
    params:
    	pop_pair = "{pop_pair}",
        distind1 = "{distinguished_ind1}",
        nr_bootstraps = 10,
        chunk_size = 5000000,
        nr_chr = 9,
        input_dir = "models/smc_split_input/{pop_pair}_12.{distinguished_ind1}*"
    shell:
        """
        python3 scripts/smc_bootstrap.py \
        --nr_bootstraps {params.nr_bootstraps} \
        --chunk_size {params.chunk_size} \
        --chunks_per_chromosome 10 \
        --nr_chromosomes {params.nr_chr} \
        models/smc_split_bootstrap_input/{params.pop_pair}_12.{params.distind1}_rep \
        {params.input_dir}
        """

rule joint_bootstrap_vcf2smc21:
    input:
        expand("models/smc_split_input/{{pop_pair}}_21.{{distinguished_ind1}}.{chr}.smc.gz", chr = CHR)
    output:
        expand('models/smc_split_bootstrap_input/{{pop_pair}}_21.{{distinguished_ind1}}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz', n_bootstrap = range(1,11), boot_chr = range(1,10)),
    params:
    	pop_pair = "{pop_pair}",
        distind1 = "{distinguished_ind1}",
        nr_bootstraps = 10,
        chunk_size = 5000000,
        nr_chr = 9,
        input_dir = "models/smc_split_input/{pop_pair}_21.{distinguished_ind1}*"
    shell:
        """
        python3 scripts/smc_bootstrap.py \
        --nr_bootstraps {params.nr_bootstraps} \
        --chunk_size {params.chunk_size} \
        --chunks_per_chromosome 10 \
        --nr_chromosomes {params.nr_chr} \
        models/smc_split_bootstrap_input/{params.pop_pair}_21.{params.distind1}_rep \
        {params.input_dir}
        """

################################################################################
# STEP 4
# Estimate a joint demography for two populations using bootstrapped samples 
################################################################################
rule smc_split_bootstrap:
    input:
        smc_split_bootstrap_input
    params:
        model_out_dir = "models/smc_split_bootstrap/{pop_pair}_{n_bootstrap}"
    output:
        model_out = "models/smc_split_bootstrap/{pop_pair}_{n_bootstrap}/model.final.json"
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        "smc++ split \
        --cores 2 \
        -o {params.model_out_dir} \
        {input}"

################################################################################
# STEP 5
# Plot results for divergence time estimates 
#   --csv exports a csv formatted file with results
#   -g specifies the number of years per generation
################################################################################

rule plot_split:
    input:
        smc_split = expand("models/smc_split/{pop_pair}/model.final.json", pop_pair = pop_pair_dict.keys()),
        smc_split_bootstrap = expand("models/smc_split_bootstrap/{pop_pair}_{n_bootstrap}/model.final.json", pop_pair = pop_pair_dict.keys(), n_bootstrap = range(1,11))
    output:
        split = "results/carrot_all_pops_smc_split.png",
        bootstrap = "results/carrot_all_pops_smc_split_bootstrap.png"
    params:
        gen = config['gen']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        """
        smc++ plot \
        --csv \
        -g {params.gen} \
        {output.split} \
        {input.smc_split}
        smc++ plot \
        --csv \
        -g {params.gen} \
        {output.bootstrap} \
        {input.smc_split_bootstrap}
        """
