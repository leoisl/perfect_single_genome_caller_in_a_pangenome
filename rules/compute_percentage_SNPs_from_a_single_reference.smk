#transforms all_unique_canonical_snps_for_a_single_genome in a SNP panel fasta file
rule build_SNP_panel_fasta_file_for_a_single_genome:
    input:
        all_unique_canonical_snps_for_a_single_genome = Path(config["output_folder"]) / "{genome_1}.all_unique_canonical_snps"
    output:
        SNP_panel_fasta_file_for_a_single_genome = Path(config["output_folder"]) / "{genome_1}.SNP_panel.fa"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{genome_1}_build_SNP_panel_fasta_file_for_a_single_genome.log"
    run:
        index = 0
        with open(input.all_unique_canonical_snps_for_a_single_genome) as all_unique_canonical_snps_file_for_a_single_genome, \
             open(output.SNP_panel_fasta_file_for_a_single_genome, "w") as SNP_panel_fasta_file_for_a_single_genome:
            for line in all_unique_canonical_snps_file_for_a_single_genome:
                first_allele, second_allele = line.strip().split("$")
                print(f">SNP_{index}_allele_1\n{first_allele}", file=SNP_panel_fasta_file_for_a_single_genome)
                print(f">SNP_{index}_allele_2\n{second_allele}", file=SNP_panel_fasta_file_for_a_single_genome)
                index+=1




#run bwa mem to get which SNPs from SNP_refined_panel.fa we can find with this genome
rule get_SNPs_mapping_to_SNP_refined_panel:
    input:
        SNP_panel_fasta_file_for_a_single_genome = Path(config["output_folder"]) / "{genome_1}.SNP_panel.fa",
        SNP_refined_panel = Path(config["output_folder"]) / "SNP_refined_panel.fa",
        amb_file = Path(config["output_folder"]) / "SNP_refined_panel.fa.amb",
        ann_file = Path(config["output_folder"]) / "SNP_refined_panel.fa.ann",
        bwt_file = Path(config["output_folder"]) / "SNP_refined_panel.fa.bwt",
        pac_file = Path(config["output_folder"]) / "SNP_refined_panel.fa.pac",
        sa_file = Path(config["output_folder"]) / "SNP_refined_panel.fa.sa"
    output:
        SNPs_found_in_the_pangenome_if_ref_is_this_genome = Path(config["output_folder"]) / "SNPs_found_in_the_pangenome_if_ref_is_{genome_1}"
    params:
        minimum_score_to_output = int((float(config["probe_length"])*2+1) * float(config["proportion_of_match_in_probes_to_say_SNPs_are_the_same"]))
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    log:
        "logs/{genome_1}_get_SNPs_mapping_to_SNP_refined_panel.log"
    shell:
        """
        bwa mem -t {threads} -A 1 -B 0 -O [6,6] -E [1,1] -L [5,5] -U 0 -T {params.minimum_score_to_output} {input.SNP_refined_panel} {input.SNP_panel_fasta_file_for_a_single_genome} |
        grep -v '^@' | awk '{{print $3}}' | awk -F '_' '{{if($2 ~  /^[0-9]+$/)print $2}}' | sort -u > {output}
        """



rule count_nb_SNPs_in_each_genome:
    input:
        SNPs_found_in_the_pangenome_if_ref_is_this_genome = expand(str(Path(config["output_folder"]) / "SNPs_found_in_the_pangenome_if_ref_is_{genomes}"), genomes=genomes_names)
    output:
        nb_SNPs_in_this_genome = Path(config["output_folder"]) / "nb_SNPs_in_each_genome",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/count_nb_SNPs_in_each_genome.log"
    shell:
        "bash scripts/get_nb_SNPs_in_each_genome.sh {input} > {output} 2> {log}"

