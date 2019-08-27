import os
from pandora1_paper.evaluate.evaluation import *
from pandora1_paper.evaluate.Positioned_SNPs import PositionedSNPsIndex
import itertools
import pandas as pd

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

rule output_SNP_clusters:
    input:
        all_snp_probes_files_done_flag_file = expand( str(Path(config["output_folder"]) / f"{{genomes}}.all_SNPs_from_mummer.done"), genomes = genomes_names)
    output:
        SNP_clusters = Path(config["output_folder"]) / "SNP_clusters"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/output_SNP_clusters"
    run:
        positionedSNPsIndex = PositionedSNPsIndex()
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
        nb_SNPs_in_pangenome = positionedSNPsIndex.get_nb_SNPs_in_pangenome()
        df = pd.DataFrame(columns=["genome_index", "sensitivity"])
        for index, genome_name in enumerate(genomes_names):
            nb_SNPs_in_genome = positionedSNPsIndex.get_nb_SNPs_that_can_be_found_with_a_given_genome(genome_name)
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
