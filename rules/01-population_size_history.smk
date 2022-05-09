# Demographic history with SMC++
# adapted from lovely script by cattlefriends Harly Durbin & Troy Rowan
# Check out SMC++ git repo for additional documentation:
#   https://github.com/popgenmethods/smcpp

################################################################################
# STEP 1
# Create input files for SMC++ (.smc.gz format)
#   --mask = poorly mapped regions to exclude from analysis
#            otherwise assumes missing data = regions of homozygosity
################################################################################

rule vcf2smc:
    input:
        vcf = config['vcf'],
        index = config['vcf_index'],
        mask = config['mask']
    output:
        "models/smc_estimate_input/{population}.{distinguished_ind}.{chr}.smc.gz"
    params:
        chrom = "{chr}",
    	# pop_choose exports a string of all individuals in current population
        list = pop_choose,
        # specify a distinguished individual
        # iterates over 10 randomly sampled individuals per population
        distind = "{distinguished_ind}"
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
    # argument order:
    # mask, input vcf, output directory, chromosome, population: list of individuals
        "smc++ vcf2smc \
        --mask {input.mask} \
        -d {params.distind} {params.distind} \
        {input.vcf} {output} {params.chrom} {params.list}"

################################################################################
# STEP 2
# Fit population size history to data
#   cv method uses cross-validation to obtain model parameters
################################################################################

# rule smc_cv:
#     input:
#     	smc_cv_input
#     output:
#     	cv_folds = expand("models/smc_cv_no_timepoints/{{population}}/fold{fold}/model.final.json", fold = ['0','1']),
#         final_model = "models/smc_cv_no_timepoints/{population}/model.final.json"
#     threads: 25
#     params:
#     	model_in = "models/smc/input/{population}.*",
#     	model_out_dir = "models/smc_cv_no_timepoints/{population}",
#     	mu = config['mu']
#     singularity:
#         "docker://terhorst/smcpp:latest"
#     shell:
#     	"smc++ cv \
#         --cores {threads} \
#         --spline cubic \
#         --ftol 0.001 \
#     	-o {params.model_out_dir} {params.mu} {params.model_in}"

rule smc_estimate:
    input:
    	smc_estimate_input
    output:
        final_model = "models/smc_estimate_no_timepoints/{population}/model.final.json"
    threads: 25
    params:
    	model_in = "models/smc_estimate_input/{population}.*",
    	model_out_dir = "models/smc_estimate_no_timepoints/{population}",
    	mu = config['mu']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
    	"smc++ estimate \
        --cores {threads} \
        --spline cubic \
    	-o {params.model_out_dir} {params.mu} {params.model_in}"


################################################################################
# STEP 3
# Bootstrapping
#   split smc input files into chunks and resample
################################################################################

rule smc_bootstrap_input:
    input:
        expand("models/smc_estimate_input/{{population}}.{{distinguished_ind}}.{chr}.smc.gz", chr = CHR)
    output:
        expand("models/smc_estimate_bootstrap_input/{{population}}_{{distinguished_ind}}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz", n_bootstrap = range(1,11), boot_chr = range(1,10))
    params:
        population = "{population}",
        distind = "{distinguished_ind}", # lambda wildcards: dist_dict[wildcards.population],
        nr_bootstraps = 10,
        chunk_size = 5000000,
        nr_chr = 9,
        input_dir = "models/smc_estimate_input/{population}.{distinguished_ind}*"
    shell:
        "python3 scripts/smc_bootstrap.py \
        --nr_bootstraps {params.nr_bootstraps} \
        --chunk_size {params.chunk_size} \
        --chunks_per_chromosome 10 \
        --nr_chromosomes {params.nr_chr} \
        models/smc_estimate_bootstrap_input/{params.population}_{params.distind}_rep \
        {params.input_dir}"

# rule smc_cv_bootstrap:
#     input:
#         smc_cv_boot_input
#     output:
#     	cv_folds = expand("models/smc_cv_bootstrap/{{population}}_{{n_bootstrap}}/fold{fold}/model.final.json", fold = ['0','1']),
#         final_model = "models/smc_cv_bootstrap/{population}_{n_bootstrap}/model.final.json"
#     threads: 16
#     params:
#     	model_in = "models/smc/bootstrap_input/{population}_*rep_{n_bootstrap}/*",
#     	model_out_dir = "models/smc_cv_bootstrap/{population}_{n_bootstrap}",
#     	mu = config['mu']
#     singularity:
#         "docker://terhorst/smcpp:latest"
#     shell:
#     	"smc++ cv \
#         --cores {threads} \
#         --spline cubic \
#     	-o {params.model_out_dir} {params.mu} {params.model_in}"

rule smc_estimate_bootstrap:
    input:
    	smc_estimate_boot_input
    output:
        final_model = "models/smc_estimate_no_timepoints_bootstrap/{population}_{n_bootstrap}/model.final.json"
    threads: 16
    params:
    	model_in = "models/smc_estimate_bootstrap_input/{population}_*rep_{n_bootstrap}/*",
    	model_out_dir = "models/smc_estimate_no_timepoints_bootstrap/{population}_{n_bootstrap}",
    	mu = config['mu']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
    	"smc++ estimate \
        --cores {threads} \
        --spline cubic \
    	-o {params.model_out_dir} {params.mu} {params.model_in}"

################################################################################
# STEP 4
# Plot results for population size history estimates
#   --csv exports a csv formatted file with results
#   -g specifies the number of years per generation
################################################################################

rule plot_estimate:
    input:
        smc_out = expand("models/smc_estimate_no_timepoints/{population}/model.final.json", population = popdict.keys()),
        smc_bootstrap_out = expand("models/smc_estimate_no_timepoints_bootstrap/{population}_{n_bootstrap}/model.final.json", population = popdict.keys(), n_bootstrap = range(1,11))
    output:
        smc = "reports/carrot_all_pops_smc_estimate_no_timepoints.png",
        bootstrap = "reports/carrot_all_pops_smc_estimate_no_timepoints_bootstrap.png"
    params:
        gen = config['gen']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        """
        smc++ plot \
        --csv \
        -g {params.gen} \
        {output.smc} \
        {input.smc_out}
        smc++ plot \
        --csv \
        -g {params.gen} \
        {output.bootstrap} \
        {input.smc_bootstrap_out}
        """
