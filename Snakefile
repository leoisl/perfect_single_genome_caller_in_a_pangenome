from scripts.utils import *


# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

# ======================================================
# Rules
# ======================================================
genomes = get_all_genomes_files(config["input_folder"])
genomes_names = get_all_genomes_names(genomes)

rule all:
    #input: Path(config["output_folder"]) / "perfect_caller_sensitivity_in_pangenome.pdf"
    input:
        get_rule_build_SNP_panel_fasta_file_for_a_single_genome_final_files(config["output_folder"]),
        get_rule_count_nb_SNPs_in_pangenome_final_files(config["output_folder"])

rules_dir = Path("rules/")
include: str(rules_dir / "count_nb_SNPs_in_pangenome.smk")
include: str(rules_dir / "compute_percentage_SNPs_from_a_single_reference.smk")
#include: str(rules_dir / "generate_plot.smk")
