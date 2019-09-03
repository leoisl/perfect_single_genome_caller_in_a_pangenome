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
species_to_genomes = get_all_genomes_files(config["input_folders"])
species_to_output_folder = get_all_output_folders(config["input_folders"])
species_to_genomes_names = get_all_genomes_names(species_to_genomes)


# ======================================================
# Rules
# ======================================================
rules_dir = Path("rules/")
include: str(rules_dir / "core.smk")


# ======================================================
# Running all rules
# ======================================================
rule all:
    input:
         rules.generate_violin_plot.output.plot,
         expand(rules.generate_geom_line_plot.output.plot, input_folder = config["input_folders"])
