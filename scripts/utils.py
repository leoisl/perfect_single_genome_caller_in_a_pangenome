from pathlib import Path
import copy

SEPARATOR = "-SEP-"

def get_all_genomes_files(input_folders):
    species_to_genomes = {}
    for species in input_folders:
        species_to_genomes[species] = [str(x) for x in list(Path(species).glob("*.fna"))]
    return species_to_genomes

def get_all_output_folders(input_folders):
    return {input_folder: f"{input_folder}_out" for input_folder in input_folders}

def get_all_genomes_names(species_to_genomes):
    species_to_genomes_names = copy.deepcopy(species_to_genomes)
    for genomes in species_to_genomes_names.values():
        for i, genome in enumerate(genomes):
            genomes[i] = Path(genome).with_suffix("").name
    return species_to_genomes_names


def get_SNPs_using_mummer_final_files(species_to_output_folder, species_to_genomes_names):
    final_files=[]
    for species, output_folder in species_to_output_folder.items():
        for genome_name in species_to_genomes_names[species]:
            final_files.append(f"{output_folder}/{genome_name}.all_SNPs_from_mummer.done")
    return final_files


# if __name__ == "__main__":
#     all_genomes_files = get_all_genomes_files(["5_random_genomes_ecoli", "5_random_genomes_saureus", "5_random_genomes_paeru"])
#     all_genomes_names = get_all_genomes_names(all_genomes_files)
#     print(all_genomes_files)
#     print(all_genomes_names)
