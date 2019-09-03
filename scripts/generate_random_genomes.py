import random
import copy
ACGT = list("ACGT")

def generate_random_genome(length):
    random_genome = []
    for i in range(length):
        random_genome.append(random.choice(ACGT))
    genome = "".join(random_genome)
    assert len(genome) == length
    return genome

def nb_differences_between_two_strings(str1, str2):
    assert len(str1) == len(str2)
    return len([i for i in range(len(str1)) if str1[i] != str2[i]])

def mutate_genome(genome, nb_positions_to_mutate):
    genome_as_list = list(genome)
    positions_to_mutate = random.sample(list(range(len(genome))), k=nb_positions_to_mutate)

    for position_to_mutate in positions_to_mutate:
        ACGT_without_the_original_base = copy.deepcopy(ACGT)
        ACGT_without_the_original_base.remove(genome[position_to_mutate])
        genome_as_list[position_to_mutate] = random.choice(ACGT_without_the_original_base)

    mutated_genome = "".join(genome_as_list)

    assert nb_differences_between_two_strings(genome, mutated_genome) == nb_positions_to_mutate, str(locals())
    return mutated_genome

def submutate_genome(original_genome, mutated_genome, positions_to_mutate):
    submutated_genome_as_list = list(original_genome)
    for i in range(len(submutated_genome_as_list)):
        if i in positions_to_mutate:
            submutated_genome_as_list[i] = mutated_genome[i]

    submutated_genome = "".join(submutated_genome_as_list)
    assert nb_differences_between_two_strings(original_genome, submutated_genome) == len(positions_to_mutate)
    return submutated_genome

def create_random_genome_and_mutated(genome_length, nb_bases_to_mutate_in_each_genome, basename):
    original_genome = generate_random_genome(genome_length)
    with open(f"{basename}.fna", "w") as fout:
        print(">random_original_genome", file=fout)
        print(original_genome, file=fout)

    nb_mutated_genomes = len(nb_bases_to_mutate_in_each_genome)
    subgenome_length = int(genome_length/nb_mutated_genomes)
    for mutation_index, nb_bases_to_mutate in enumerate(nb_bases_to_mutate_in_each_genome):
        original_subgenome = original_genome[mutation_index*subgenome_length : mutation_index*subgenome_length+subgenome_length]
        assert len(original_subgenome) == subgenome_length
        mutated_subgenome = mutate_genome(original_subgenome, nb_bases_to_mutate)
        with open(f"{basename}.mutation_{mutation_index}.fna", "w") as fout:
            print(f">mutated_genome_{mutation_index}", file=fout)
            print(mutated_subgenome, file=fout)


def sample_run():
    create_random_genome_and_mutated(200000, [1000, 500, 250, 125, 75, 50, 25, 10, 5], "random_genome")

sample_run()