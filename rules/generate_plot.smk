rule generate_plot:
    input:
        nb_SNPs_in_pangenome = Path(config["output_folder"]) / "nb_SNPs_in_pangenome",
        nb_SNPs_in_each_genome = Path(config["output_folder"]) / "nb_SNPs_in_each_genome"
    output:
        plot = Path(config["output_folder"]) / "perfect_caller_sensitivity_in_pangenome.pdf"
    params:
        dataframe = lambda wildcards, output: Path(output.plot).with_suffix('.dataframe')
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/generate_plot.log"
    shell:
        """
        python scripts/get_perfect_caller_sensitivity_on_pangenome.py {input.nb_SNPs_in_pangenome} {input.nb_SNPs_in_each_genome} > {params.dataframe} &&
        Rscript scripts/create_plot.R  {params.dataframe} {output.plot}
        """