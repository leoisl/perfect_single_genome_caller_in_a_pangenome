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
    input: Path(config["output_folder"]) / "perfect_caller_sensitivity_in_pangenome.pdf"

rules_dir = Path("rules/")
include: str(rules_dir / "count_nb_SNPs_in_pangenome.smk")
include: str(rules_dir / "compute_percentage_SNPs_from_a_single_reference.smk")
include: str(rules_dir / "generate_plot.smk")
