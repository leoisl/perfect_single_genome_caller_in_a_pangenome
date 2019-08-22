from pathlib import Path
from Bio.Seq import Seq
from Bio import pairwise2

def rev_comp(s):
    return str(Seq(s).reverse_complement())

def similarity(s1, s2):
    return pairwise2.align.globalxx(s1, s2)[0][2]

SEPARATOR = "-SEP-"
nb_of_genomes = 2

def get_all_genomes_files(input_folder):
    return list(Path(input_folder).glob("*.fna"))

def get_all_genomes_names(genomes):
    return [genome.with_suffix("").name for genome in genomes]

def get_rule_run_dnadiff_final_files(output_folder):
    return [f"{output_folder}/genome.{i}.all_delta_files_done" for i in range(nb_of_genomes)]

def get_rule_run_show_snps_final_files(output_folder):
    return [f"{output_folder}/genome.{i}.all_snp_probes_files_done" for i in range(nb_of_genomes)]

def get_rule_get_SNPs_using_mummer_final_files(output_folder):
    return [f"{output_folder}/genome.{i}.all_SNPs_from_mummer.done" for i in range(nb_of_genomes)]

def get_rule_run_transform_SNPs_into_canonical_SNPs_final_files(output_folder):
    return [f"{output_folder}/genome.{i}.all_canonical_snps_files_done" for i in range(nb_of_genomes)]

def get_rule_get_unique_canonical_SNPs_from_the_pangenome_final_files(output_folder):
    return Path(output_folder) / "all_unique_canonical_snps"

def get_rule_build_SNP_panel_fasta_file_final_files(output_folder):
    return Path(output_folder) / "SNP_panel.fa"

def get_rule_get_unrefined_clusters_using_bwa_mem_final_files(output_folder):
    return Path(output_folder) / "unrefined_clusters"

def get_rule_refine_clusters_and_output_SNP_refined_panel_final_files(output_folder):
    return Path(output_folder) / "SNP_refined_panel.fa"

def get_rule_count_nb_SNPs_in_pangenome_final_files(output_folder):
    return Path(output_folder) / "nb_SNPs_in_pangenome"

def get_rule_get_unique_canonical_SNPs_for_a_single_genome_final_files(output_folder):
    return [f"{output_folder}/genome.{i}.all_unique_canonical_snps" for i in range(nb_of_genomes)]

def get_rule_build_SNP_panel_fasta_file_for_a_single_genome_final_files(output_folder):
    return [f"{output_folder}/genome.{i}.SNP_panel.fa" for i in range(nb_of_genomes)]

def get_rule_get_SNPs_mapping_to_SNP_refined_panel_final_files(output_folder):
    return [f"{output_folder}/SNPs_found_in_the_pangenome_if_ref_is_genome.{i}" for i in range(nb_of_genomes)]

def get_rule_count_nb_SNPs_in_each_genome_final_files(output_folder):
    return f"{output_folder}/nb_SNPs_in_each_genome"

def get_all_unique_coordinate_SNPs_from_a_genome_final_files(output_folder):
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
