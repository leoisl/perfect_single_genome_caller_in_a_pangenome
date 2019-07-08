from scripts.utils import *


# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

# ======================================================
# Rules
# ======================================================
genomes = get_all_genomes_files(config["input_folder"])
#all_final_files = get_all_final_files(genomes, config["output_folder"])
all_final_files = get_rule_run_dnadiff_final_files(config['output_folder'])

rule all:
    input: all_final_files

rules_dir = Path("rules/")
include: str(rules_dir / "count_nb_SNPs_in_pangenome.smk")
