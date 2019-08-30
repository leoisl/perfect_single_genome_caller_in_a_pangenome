rule get_SNPs_using_mummer:
    input:
        ref_genome = Path(config["input_folder"]) / "{genome_1}.fna",
        all_other_genomes = expand( str(Path(config["input_folder"]) / "{genomes_2}.fna"), genomes_2 = genomes_names)
    output:
        touch(Path(config["output_folder"]) / f"{{genome_1}}.all_SNPs_from_mummer.done")
    params:
        probe_length = config["probe_length"], #TODO: vary several probe lengths?
        genomes = genomes,
        genomes_names = genomes_names
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{genome_1}_get_SNPs_using_mummer.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/get_SNPs_using_mummer.py"


rule get_length_of_longest_contig:
    input:
        all_genomes = expand(str(Path(config["input_folder"]) / "{genome}.fna"), genome = genomes_names)
    output:
        length_of_longest_contig_file = Path(config["output_folder"]) / "length_of_longest_contig"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/get_length_of_longest_contig.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/get_length_of_longest_contig.py"


rule output_SNP_clusters:
    input:
        all_snp_probes_files_done_flag_file = expand( str(Path(config["output_folder"]) / f"{{genomes}}.all_SNPs_from_mummer.done"), genomes = genomes_names),
        length_of_longest_contig_file = rules.get_length_of_longest_contig.output.length_of_longest_contig_file
    output:
        SNP_clusters = Path(config["output_folder"]) / "SNP_clusters"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    params:
        genomes = genomes,
        genomes_names = genomes_names
    log:
        "logs/output_SNP_clusters"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/output_SNP_clusters.py"


rule get_perfect_genotyper_recall:
    input:
        SNP_clusters = rules.output_SNP_clusters.output.SNP_clusters
    output:
        perfect_genotyper_recall_df = Path(config["output_folder"]) / "perfect_genotyper.dataframe"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    params:
        genomes = genomes,
        genomes_names = genomes_names
    log:
        "logs/get_perfect_genotyper_recall"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/get_perfect_genotyper_recall.py"


rule generate_plot:
    input:
        perfect_genotyper_recall_df = rules.get_perfect_genotyper_recall.output.perfect_genotyper_recall_df
    output:
        plot = Path(config["output_folder"]) / "perfect_caller_sensitivity_in_pangenome.pdf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/generate_plot.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/generate_plot.R"
