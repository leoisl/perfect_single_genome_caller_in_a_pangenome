rule get_SNPs_using_mummer:
    input:
        ref_genome = "{input_folder}/{genome_1}.fna"
    output:
        touch("{input_folder}_out/{genome_1}.all_SNPs_from_mummer.done")
    params:
        probe_length = config["probe_length"], #TODO: vary several probe lengths?
        genomes = lambda wildcards: species_to_genomes[wildcards.input_folder],
        genomes_names = lambda wildcards: species_to_genomes_names[wildcards.input_folder]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{input_folder}_{genome_1}_get_SNPs_using_mummer.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/get_SNPs_using_mummer.py"


rule get_length_of_longest_contig:
    input:
        all_genomes = lambda wildcards: species_to_genomes[wildcards.input_folder]
    output:
        length_of_longest_contig_file = "{input_folder}_out/length_of_longest_contig"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{input_folder}_get_length_of_longest_contig.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/get_length_of_longest_contig.py"


rule output_SNP_clusters:
    input:
        all_snp_probes_files_done_flag_file = lambda wildcards: expand(wildcards.input_folder+"_out/{genomes}.all_SNPs_from_mummer.done", genomes = species_to_genomes_names[wildcards.input_folder]),
        length_of_longest_contig_file = rules.get_length_of_longest_contig.output.length_of_longest_contig_file
    output:
        SNP_clusters = "{input_folder}_out/SNP_clusters"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    params:
        genomes = lambda wildcards: species_to_genomes[wildcards.input_folder],
        genomes_names = lambda wildcards: species_to_genomes_names[wildcards.input_folder]
    log:
        "logs/{input_folder}_output_SNP_clusters.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/output_SNP_clusters.py"


rule get_perfect_genotyper_recall:
    input:
        SNP_clusters = rules.output_SNP_clusters.output.SNP_clusters
    output:
        perfect_genotyper_recall_df = "{input_folder}_out/perfect_genotyper.dataframe"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    params:
        genomes = lambda wildcards: species_to_genomes[wildcards.input_folder],
        genomes_names = lambda wildcards: species_to_genomes_names[wildcards.input_folder]
    log:
        "logs/{input_folder}_get_perfect_genotyper_recall.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/get_perfect_genotyper_recall.py"


rule generate_geom_line_plot:
    input:
        perfect_genotyper_recall_df = rules.get_perfect_genotyper_recall.output.perfect_genotyper_recall_df
    output:
        plot = "{input_folder}_out/perfect_caller_sensitivity_in_pangenome.geom_line.pdf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{input_folder}_generate_geom_line_plot.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/generate_geom_line_plot.R"


rule concatenate_dataframes:
    input:
        perfect_genotyper_recall_dfs = expand(rules.get_perfect_genotyper_recall.output.perfect_genotyper_recall_df, input_folder = config["input_folders"])
    output:
        concatenated_df = "perfect_genotyper.dataframe"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/concatenate_dataframes.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/concatenate_dataframes.py"


rule generate_violin_plot:
    input:
        concatenated_df = rules.concatenate_dataframes.output.concatenated_df
    output:
        plot = "perfect_caller_sensitivity_in_pangenome.violin.pdf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/generate_violin_plot.log"
    singularity:
        config["singularity_image"]
    script:
        "../scripts/generate_violin_plot.R"
