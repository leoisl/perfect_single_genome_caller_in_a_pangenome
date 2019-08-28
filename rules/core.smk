import os
from pandora1_paper.evaluate.evaluation import *
from pandora1_paper.evaluate.Positioned_SNPs import PositionedSNPsIndex
import itertools
import pandas as pd
import pysam

rule get_SNPs_using_mummer:
    input:
        ref_genome = Path(config["input_folder"]) / "{genome_1}.fna",
        all_other_genomes = expand( str(Path(config["input_folder"]) / "{genomes_2}.fna"), genomes_2 = genomes_names)
    output:
        all_SNPs_from_mummer_done_flag_file = Path(config["output_folder"]) / f"{{genome_1}}.all_SNPs_from_mummer.done"
    params:
        probe_length = config["probe_length"] #TODO: vary several probe lengths?
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{genome_1}_get_SNPs_using_mummer.log"
    run:
        for genome_path, genome_name in zip(genomes, genomes_names):
            if wildcards.genome_1 < genome_name:
                out_dir = Path(config['output_folder']) / wildcards.genome_1
                os.makedirs(out_dir, exist_ok=True)
                prefix = f"{out_dir / wildcards.genome_1}{SEPARATOR}{genome_name}.mummer"
                snps_as_StringIO = generate_mummer_snps(reference=Path(input.ref_genome), query=Path(genome_path), prefix=Path(prefix), flank_width=params.probe_length)
                snps_dataframe = ShowSnps.to_dataframe(snps_as_StringIO)
                snps_dataframe = snps_dataframe.translate_to_FWD_strand()
                snps_dataframe.to_csv(prefix+".csv", sep="\t")
        shell("touch {output.all_SNPs_from_mummer_done_flag_file}")


def get_length_of_longest_contig(path):
    length_of_longest_contig = 0
    with pysam.FastxFile(path) as fh:
        for entry in fh:
            length_of_longest_contig = max(length_of_longest_contig, len(entry.sequence))
    return length_of_longest_contig

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
    run:
        length_of_longest_contig = 0
        for genome in input.all_genomes:
            length_of_longest_contig = max(length_of_longest_contig, get_length_of_longest_contig(genome))
        with open(output.length_of_longest_contig_file, "w") as fout:
            print(length_of_longest_contig, file=fout)


rule output_SNP_clusters:
    input:
        all_snp_probes_files_done_flag_file = expand( str(Path(config["output_folder"]) / f"{{genomes}}.all_SNPs_from_mummer.done"), genomes = genomes_names),
        length_of_longest_contig_file = rules.get_length_of_longest_contig.output.length_of_longest_contig_file
    output:
        SNP_clusters = Path(config["output_folder"]) / "SNP_clusters"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/output_SNP_clusters"
    run:
        with open(input.length_of_longest_contig_file) as fin:
            length_of_longest_contig = int(fin.readlines()[0])
        positionedSNPsIndex = PositionedSNPsIndex(length_of_longest_contig)
        for genome_name_1, genome_name_2 in itertools.product(genomes_names, genomes_names):
            if genome_name_1 < genome_name_2:
                out_dir = Path(config['output_folder']) / genome_name_1
                prefix = f"{out_dir / genome_name_1}{SEPARATOR}{genome_name_2}.mummer"
                positionedSNPsIndex.add_SNPs_from_csv(prefix + ".csv", genome_name_1, genome_name_2)
        positionedSNPsIndex.save(output.SNP_clusters)

rule get_perfect_genotyper_recall:
    input:
        SNP_clusters = rules.output_SNP_clusters.output.SNP_clusters
    output:
        perfect_genotyper_recall_df = Path(config["output_folder"]) / "perfect_genotyper.dataframe"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/get_perfect_genotyper_recall"
    run:
        positionedSNPsIndex = PositionedSNPsIndex.load(input.SNP_clusters)
        genomes_to_nb_of_SNPs = positionedSNPsIndex.get_nb_SNPs_that_can_be_found_with_the_given_genomes(genomes_names)
        nb_SNPs_in_pangenome = genomes_to_nb_of_SNPs["all"]

        df = pd.DataFrame(columns=["genome_index", "sensitivity"])
        for index, genome_name in enumerate(genomes_names):
            nb_SNPs_in_genome = genomes_to_nb_of_SNPs[genome_name]
            df.loc[index] = [index, nb_SNPs_in_genome/nb_SNPs_in_pangenome]
        df.to_csv(output.perfect_genotyper_recall_df, sep="\t")

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
    shell: "Rscript scripts/create_plot.R {input.perfect_genotyper_recall_df} {output.plot}"
