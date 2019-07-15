import os
import subprocess

rule run_dnadiff:
    input:
        ref_genome = Path(config["input_folder"]) / "{genome_1}.fna",
        all_other_genomes = expand( str(Path(config["input_folder"]) / "{genomes_2}.fna"), genomes_2 = genomes_names)
    output:
        all_delta_files_done_flag_file = Path(config["output_folder"]) / f"{{genome_1}}.all_delta_files_done"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    log:
        "logs/{genome_1}_run_dnadiff.log"
    run:
        for genome_path, genome_name in zip(genomes, genomes_names):
            out_dir = Path(config['output_folder']) / wildcards.genome_1
            os.makedirs(out_dir, exist_ok=True)
            prefix = f"{out_dir / wildcards.genome_1}{SEPARATOR}{genome_name}"
            subprocess.check_output(f"dnadiff {input.ref_genome} {genome_path} -p {prefix} ", shell=True)
        with open(output.all_delta_files_done_flag_file, "w") as fout: pass

rule run_show_snps:
    input:
        all_delta_files_done_flag_file = Path(config["output_folder"]) / f"{{genome_1}}.all_delta_files_done"
    output:
        all_snp_probes_files_done_flag_file = Path(config["output_folder"]) / f"{{genome_1}}.all_snp_probes_files_done"
    params:
        probe_length = config["probe_length"] #TODO: vary several probe lengths?
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    log:
        "logs/{genome_1}_run_show_snps.log"
    run:
        #output file header: [P1]    [SUB]   [SUB]   [P2]    [BUFF]  [DIST]  [CTX R] [CTX Q] [FRM]   [TAGS]
        for genome_path, genome_name in zip(genomes, genomes_names):
            out_dir = Path(config['output_folder']) / wildcards.genome_1
            prefix = f"{out_dir / wildcards.genome_1}{SEPARATOR}{genome_name}"
            delta_file = f"{prefix}.delta"
            show_snps_file = f"{prefix}.show_snps"
            subprocess.check_output(f"show-snps -C -H -I -r -T -x {params.probe_length} {delta_file} > {show_snps_file}", shell=True)
        with open(output.all_snp_probes_files_done_flag_file, "w") as fout: pass


#this is done so that SNPs that are the same are not represented duplicatedly
rule transform_SNPs_into_canonical_SNPs:
    input:
        all_snp_probes_files_done_flag_file = Path(config["output_folder"]) / f"{{genome_1}}.all_snp_probes_files_done"
    output:
        all_canonical_snps_files_done_flag_file = Path(config["output_folder"]) / f"{{genome_1}}.all_canonical_snps_files_done"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{genome_1}_transform_SNPs_into_canonical_SNPs.log"
    run:
        #output file header: [P1]    [SUB]   [SUB]   [P2]    [BUFF]  [DIST]  [CTX R] [CTX Q] [FRM]   [TAGS]
        for genome_path, genome_name in zip(genomes, genomes_names):
            out_dir = Path(config['output_folder']) / wildcards.genome_1
            prefix = f"{out_dir / wildcards.genome_1}{SEPARATOR}{genome_name}"
            show_snps_file = f"{prefix}.show_snps"
            canonical_snps_file = f"{prefix}.canonical_snps"
            subprocess.check_output(f"python scripts/transform_snps_into_canonical.py < {show_snps_file} > {canonical_snps_file}", shell=True)
        with open(output.all_canonical_snps_files_done_flag_file, "w") as fout: pass


rule get_unique_canonical_SNPs_for_a_single_genome:
    input:
        all_canonical_snps_files_done_flag_file = Path(config["output_folder"]) / f"{{genome_1}}.all_canonical_snps_files_done"
    output:
        all_unique_canonical_snps_for_a_single_genome = Path(config["output_folder"]) / "{genome_1}.all_unique_canonical_snps"
    params:
        all_canonical_snps = lambda wildcard: expand( str(Path(config["output_folder"]) / f"{wildcard.genome_1}/{wildcard.genome_1}{SEPARATOR}{{genomes_2}}.canonical_snps"), genomes_2=genomes_names)
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    log:
        "logs/{genome_1}_get_unique_canonical_SNPs_for_a_single_genome.log"
    shell:
        "sort {params.all_canonical_snps} --parallel={threads} | uniq > {output} 2> log"



rule get_unique_canonical_SNPs_from_the_pangenome:
    input:
         all_unique_canonical_snps_for_a_single_genome = expand( str(Path(config["output_folder"]) / "{genomes_1}.all_unique_canonical_snps"), genomes_1=genomes_names)
    output:
        all_unique_canonical_snps = Path(config["output_folder"]) / "all_unique_canonical_snps"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    log:
        "logs/get_unique_canonical_SNPs.log"
    shell:
        "sort {input} --parallel={threads} | uniq > {output} 2> log"

#transforms all_unique_canonical_snps in a SNP panel fasta file
rule build_SNP_panel_fasta_file:
    input:
        all_unique_canonical_snps = Path(config["output_folder"]) / "all_unique_canonical_snps"
    output:
        SNP_panel_fasta_file = Path(config["output_folder"]) / "SNP_panel.fa"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/build_SNP_panel_fasta_file.log"
    run:
        index = 0
        with open(input.all_unique_canonical_snps) as all_unique_canonical_snps_file, \
             open(output.SNP_panel_fasta_file, "w") as SNP_panel_fasta_file:
            for line in all_unique_canonical_snps_file:
                first_allele, second_allele = line.strip().split("$")
                print(f">SNP_{index}_allele_1\n{first_allele}", file=SNP_panel_fasta_file)
                print(f">SNP_{index}_allele_2\n{second_allele}", file=SNP_panel_fasta_file)
                index+=1

#run bwa mem to get the unrefined clusters
rule get_unrefined_clusters_using_bwa_mem:
    input:
        SNP_panel_fasta_file = Path(config["output_folder"]) / "SNP_panel.fa"
    output:
        unrefined_clusters = Path(config["output_folder"]) / "unrefined_clusters"
    params:
        minimum_score_to_output = int((float(config["probe_length"])*2+1) * float(config["proportion_of_match_in_probes_to_say_SNPs_are_the_same"])),
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    log:
        "logs/get_unrefined_clusters_using_bwa_mem.log"
    shell:
        """
        bwa index {input} &&
        bwa mem -t {threads} -A 1 -B 0 -O [6,6] -E [1,1] -L [5,5] -U 0 -T {params.minimum_score_to_output} -a {input} {input} |
        grep -v '^@' | awk '{{print $1, $3}}' | awk -F '_' '{{print $2, $5}}' | sort | uniq > {output}
        """

# refine clusters:
#   * bwa mem give us some clues of the clusters;
#   * given two SNPs that bwa mem told us that are similar, we have to ensure that they are indeed similar;
#   * we do this by aligning them (each path of each SNP should align to one path) and ensure the middle bases are the same;

# we represent the SNPs and their relationships as an undirected graph:
#   * SNPs are nodes in the graph;
#   * we have an edge between two SNPs if they are similar enough;
#   * we get the connected components of this graph (each connected component is a SNP cluster);
#   * the representative SNP of each SNP cluster is a node with highest degree;

# the SNP refined panel are the representative SNPs
rule refine_clusters_and_output_SNP_refined_panel:
    input:
        unrefined_clusters = Path(config["output_folder"]) / "unrefined_clusters"
    output:
        SNP_refined_panel = Path(config["output_folder"]) / "SNP_refined_panel.fa"
    params:
        SNP_panel = Path(config["output_folder"]) / "SNP_panel.fa"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    log:
        "logs/refine_clusters_and_output_SNP_refined_panel.log"
    shell:
         "python scripts/refine_clusters_and_output_a_representative.py {input} {params.SNP_panel} > {output}"

rule index_refined_panel:
    input:
        SNP_refined_panel = Path(config["output_folder"]) / "SNP_refined_panel.fa"
    output:
        Path(config["output_folder"]) / "SNP_refined_panel.fa.amb",
        Path(config["output_folder"]) / "SNP_refined_panel.fa.ann",
        Path(config["output_folder"]) / "SNP_refined_panel.fa.bwt",
        Path(config["output_folder"]) / "SNP_refined_panel.fa.pac",
        Path(config["output_folder"]) / "SNP_refined_panel.fa.sa"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb_heavy_jobs"][attempt-1]
    log:
        "logs/index_refined_panel.log"
    shell:
         "bwa index {input}"


rule count_nb_SNPs_in_pangenome:
    input:
        SNP_refined_panel = Path(config["output_folder"]) / "SNP_refined_panel.fa"
    output:
        nb_SNPs_in_pangenome = Path(config["output_folder"]) / "nb_SNPs_in_pangenome",
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/count_nb_SNPs_in_pangenome.log"
    shell:
         " grep '>' {input} | echo $((`wc -l`/2)) > {output}"