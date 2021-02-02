# Demographic history with SMC++
# adapted from lovely script by cattlefriends Harly Durbin & Troy Rowan
# Check out SMC++ git repo for additional documentation:
#   https://github.com/popgenmethods/smcpp

################################################################################
# STEP 1
# Create input files for SMC++ (.smc.gz format)
#   --mask = poorly mapped regions to exclude from analysis (otherwise assumes missing data = regions of homozygosity)
################################################################################

rule vcf2smc:
    input:
        vcf = config['vcf'],
        index = config['vcf_index'],
        mask = config['mask']
    output:
        "models/smc/input/{population}.{chr}.smc.gz"
    params:
        chrom = "{chr}",
    	# pop_choose exports a string of all individuals in current population
        list = pop_choose
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
    # argument order:
    # mask, input vcf, output directory, chromosome, population: list of individuals
        "smc++ vcf2smc \
        --mask {input.mask} {input.vcf} {output} {params.chrom} {params.list}"

################################################################################
# STEP 2
# Fit population size history to data
#   cv method uses cross-validation to obtain model parameters
################################################################################

rule smc_cv:
    input:
    	smc_chr = expand("models/smc/input/{{population}}.{chr}.smc.gz", chr = CHR)
    output:
    	cv_folds = expand("models/smc/{{population}}/fold{fold}/model.final.json", fold = ['0','1']),
        final_model = "models/smc/{population}/model.final.json"
    threads: 16
    params:
    	model_in = "models/smc/input/{population}.*",
    	model_out_dir = "models/smc/{population}",
    	mu = config['mu']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
    	"smc++ cv \
        --cores {threads} \
        --spline cubic \
    	-o {params.model_out_dir} {params.mu} {params.model_in}"

################################################################################
# STEP 3
# Plot results for population size history estimates
#   --csv exports a csv formatted file with results
#   -g specifies the number of years per generation
################################################################################

rule plot_estimate:
    input:
        smc_out = expand("models/smc/{population}/model.final.json", population = popdict.keys())
    output:
        "reports/smc_cv_results.png"
    params:
        gen = 2
    shell:
        "smc++ plot \
        --csv \
        -g {params.gen} \
        {output} \
        {input.smc_out}"