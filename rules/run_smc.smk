# Demographic history with SMC++
# adapted from lovely script by cattlefriends Harly Durbin & Troy Rowan

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


### START HERE
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
#         #         --timepoints 1e3 1e6 \

smc_plot:
    input:
        smc_out = expand("models/smc/{population}/model.final.json", population = popdict.keys()) 
    output:
        "reports/smc/{population}_smc_cv.png"
    params:
        gen = 2
    shell:
        "smc++ plot --c \
        -g {params.gen} \
        {output} \
        {input.smc_out}"

# # #Generate vcf2smc files containing the joint frequency spectrum for both populations
# rule joint_vcf2smc:
#     input:
#         vcf = "data/processed/filtered_snps/{chr}.filtered.snps.vcf.gz",
#         index = "data/processed/filtered_snps/{chr}.filtered.snps.vcf.gz.tbi"
#     output:
#         pop_pair_out = "models/smc/split/{pop_pair}.{chr}.smc.gz"
#     threads: 12
#     params:
#     	chrom = "{chr}",
#         mask = "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz",
#     	pop_pair_string = pair_string_choose
#     shell:
#     	"smc++ vcf2smc \
#         --cores {threads} \
#         -m {params.mask} \
#         {input.vcf} \
#         {output.pop_pair_out} {params.chrom} {params.pop_pair_string}"
#
# rule split:
#     input:
#         expand("models/smc/split/{pop_pair}.{chr}.smc.gz", pop_pair = ['botrytis_italica', 'italica_botrytis'], chr = chr),
#         expand("models/smc/input/{pop}.{chr}.smc.gz", pop = pops, chr = chr)
#     threads: 16
#     params:
#         model_out_dir = "models/smc/split/test",
#         marginal_models = model_chooser
#     output:
#         model_out = "models/smc/split/model.final.json"
#     shell:
#         "smc++ split \
#         -o {params.model_out_dir} \
#         --cores {threads} \
#         {params.marginal_models} \
#         {input}"

# # rule plot_split:
# #     input:
# #         model = "models/split/{rundate}.{dataset}/{pop_pair}/model.final.json"
# #     output:
# #         plot = "reports/figures/{rundate}.{dataset}.{pop_pair}.split.png"
# #     shell:
# #         "smc++ plot -c -g 2 {output.plot} {input.model}"
