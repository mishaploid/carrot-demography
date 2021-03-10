# # #Generate vcf2smc files containing the joint frequency spectrum for both populations

rule joint_vcf2smc:
    input:
        vcf = config['vcf'],
        index = config['vcf_index'],
        mask = config['mask']
    output:
        out12 = "models/smc_split/input/{pop_pair}12.{chr}.smc.gz",
        out21 = "models/smc_split/input/{pop_pair}21.{chr}.smc.gz"
    params:
    	chrom = "{chr}",
    	pop_pair_string12 = pair_string_choose12,
        pop_pair_string21 = pair_string_choose21
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        """
        smc++ vcf2smc \
        --mask {input.mask} \
        {input.vcf} {output.out12} {params.chrom} {params.pop_pair_string12}
        smc++ vcf2smc \
        --mask {input.mask} \
        {input.vcf} {output.out21} {params.chrom} {params.pop_pair_string21}
        """

rule smc_split:
    input:
        smc_split_input
    threads: 20
    params:
        model_out_dir = "models/smc_split/{pop_pair}/"
    output:
        model_out = "models/smc_split/{pop_pair}/model.final.json"
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        "smc++ split \
        -o {params.model_out_dir} \
        --cores {threads} \
        {input}"

rule plot_split:
    input:
        "models/smc_split/{pop_pair}/model.final.json"
    output:
        plot = "reports/figures/{pop_pair}.split.png"
    params:
        gen = config['gen']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        "smc++ plot \
        --csv \
        -g {params.gen} \
        {output} \
        {input}"
