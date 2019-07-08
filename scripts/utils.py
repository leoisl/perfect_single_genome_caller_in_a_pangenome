from pathlib import Path

def get_all_genomes_files(input_folder):
    return list(Path().rglob("*.fna"))

def get_rule_run_dnadiff_final_files(output_folder):
    nb_of_genomes = 4
    return [f"{output_folder}/genome_{i}.genome_{j}.delta" for i in range(nb_of_genomes) for j in range(nb_of_genomes)]

def get_rule_run_dnadiff_final_files(output_folder):
    nb_of_genomes = 4
    return [f"{output_folder}/genome_{i}.genome_{j}.show_snps" for i in range(nb_of_genomes) for j in range(nb_of_genomes)]

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
