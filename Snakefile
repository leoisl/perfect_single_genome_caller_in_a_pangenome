from scripts.utils import *
import logging

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"


# ======================================================
# Logging
# ======================================================
logging_level = logging.getLevelName(config["logging_level"])
logging.basicConfig(
        level=logging_level,
        format="[%(asctime)s]:%(levelname)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
)


# ======================================================
# Getting genomes to be processed
# ======================================================
genomes = get_all_genomes_files(config["input_folder"])
genomes_names = get_all_genomes_names(genomes)


# ======================================================
# Rules
# ======================================================
rules_dir = Path("rules/")
include: str(rules_dir / "core.smk")


# ======================================================
# Running all rules
# ======================================================
rule all:
    input: rules.generate_plot.output
