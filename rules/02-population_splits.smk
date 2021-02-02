# # #Generate vcf2smc files containing the joint frequency spectrum for both populations

rule joint_vcf2smc:
    input:
        vcf = config['vcf'],
        index = config['vcf_index'],
        mask = config['mask']
    output:
        "models/smc/split_input/{pop_pair}.{chr}.smc.gz"
    params:
    	chrom = "{chr}",
    	pop_pair_string = pair_string_choose
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
    	"smc++ vcf2smc \
        --cores {threads} \
        --mask {input.mask} \
        {input.vcf} {output} {params.chrom} {params.pop_pair_string}"
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
