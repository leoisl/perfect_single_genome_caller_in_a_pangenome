from pathlib import Path

SEPARATOR = "-SEP-"

def get_all_genomes_files(input_folder):
    return list(Path().rglob("*.fna"))

def get_all_genomes_names(genomes):
    return [genome.with_suffix("").name for genome in genomes]

def get_rule_run_dnadiff_final_files(output_folder):
    nb_of_genomes = 4
    return [f"{output_folder}/genome.{i}{SEPARATOR}genome.{j}.delta" for i in range(nb_of_genomes) for j in range(nb_of_genomes)]

def get_rule_run_dnadiff_final_files(output_folder):
    nb_of_genomes = 4
    return [f"{output_folder}/genome.{i}{SEPARATOR}genome.{j}.show_snps" for i in range(nb_of_genomes) for j in range(nb_of_genomes)]

def get_rule_run_transform_SNPs_into_canonical_SNPs_final_files(output_folder):
    nb_of_genomes = 4
    return [f"{output_folder}/genome.{i}{SEPARATOR}genome.{j}.canonical_snps" for i in range(nb_of_genomes) for j in range(nb_of_genomes)]

def get_rule_run_transform_SNPs_into_coordinate_SNPs_final_files(output_folder):
    nb_of_genomes = 4
    return [f"{output_folder}/genome.{i}{SEPARATOR}genome.{j}.pos_and_base_change" for i in range(nb_of_genomes) for j in range(nb_of_genomes)]

def get_all_unique_coordinate_SNPs_from_a_genome_final_files(output_folder):
    nb_of_genomes = 4
    return [f"{output_folder}/genome.{i}.all_coordinate_SNPs_from_a_genome" for i in range(nb_of_genomes)]

def get_all_final_files(genomes, output_folder, suffix=".percentage"):
    return [ (Path(output_folder) / genome.name).with_suffix(".percentage") for genome in genomes]

if __name__ == "__main__":
    from yaml import load, dump
    try:
        from yaml import CLoader as Loader, CDumper as Dumper
    except ImportError:
        from yaml import Loader, Dumper

    with open("config.yaml") as config_file:
        config = load(config_file, Loader=Loader)

    genomes = get_all_genomes_files(config["input_folder"])
    print(f"genomes: {genomes}")
    final_files = get_all_final_files(genomes, config["output_folder"])
    print(f"final_files: {final_files}")
    print(f"get_rule_run_dnadiff_final_files(): {get_rule_run_dnadiff_final_files(config['output_folder'])}")
    print(f"get_all_genomes_names(): {get_all_genomes_names(genomes)}")
