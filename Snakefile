from scripts.utils import *
import logging

# ======================================================
# Logging
# ======================================================
logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s]:%(levelname)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
)

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"
genomes = get_all_genomes_files(config["input_folder"])
genomes_names = get_all_genomes_names(genomes)

# ======================================================
# Rules
# ======================================================
rules_dir = Path("rules/")
include: str(rules_dir / "core.smk")

rule all:
    input: rules.generate_plot.output
