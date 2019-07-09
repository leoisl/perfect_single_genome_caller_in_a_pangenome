rule transform_SNPs_into_coordinate_SNPs:
    input:
        show_snps = Path(config["output_folder"]) / f"{{genome_1}}{SEPARATOR}{{genome_2}}.show_snps"
    output:
        pos_and_base_change = Path(config["output_folder"]) / f"{{genome_1}}{SEPARATOR}{{genome_2}}.pos_and_base_change"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    conda:
        "../envs/global.yaml"
    log:
        "logs/{genome_1}_{genome_2}transform_SNPs_into_coordinate_SNPs.log"
    shell:
        "awk '{{print $1, $2, $3}}' {input} | sort | uniq > {output} 2> log"


rule get_all_unique_coordinate_SNPs_from_a_genome:
    input:
        pos_and_base_change_from_a_genome = lambda wildcard: expand( str(Path(config["output_folder"]) / f"{{genome_1}}{SEPARATOR}{{genomes_2}}.pos_and_base_change"), genome_1 = wildcard.genome_1, genomes_2 = genomes_names)
    output:
        all_coordinate_SNPs_from_a_genome = Path(config["output_folder"]) / f"{{genome_1}}.all_coordinate_SNPs_from_a_genome"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    conda:
        "../envs/global.yaml"
    log:
        "logs/{genome_1}_get_all_unique_coordinate_SNPs_from_a_genome.log"
    shell:
        "sort {input} --parallel={threads} | uniq > {output} 2> log"

rule get_nb_unique_coordinate_SNPs_from_all_genomes:
    input:
        all_coordinate_SNPs_from_a_genome = expand(str(Path(config["output_folder"]) / f"{{genome}}.all_coordinate_SNPs_from_a_genome"), genome=genomes_names)
    output:
        nb_coordinate_SNPs_from_all_genomes = Path(config["output_folder"]) / "nb_coordinate_SNPs_from_all_genomes"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    conda:
        "../envs/global.yaml"
    log:
        "logs/get_nb_unique_coordinate_SNPs_from_all_genomes.log"
    shell:
        "bash scripts/get_nb_unique_coordinate_SNPs_from_all_genomes.sh {input} > {output} 2> log"