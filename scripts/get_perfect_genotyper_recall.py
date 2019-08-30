import pandas as pd
from pandora1_paper.evaluate.Positioned_SNPs import PositionedSNPsIndex

positionedSNPsIndex = PositionedSNPsIndex.load(snakemake.input.SNP_clusters)
genomes_to_nb_of_SNPs = positionedSNPsIndex.get_nb_SNPs_that_can_be_found_with_the_given_genomes(snakemake.params.genomes_names)
nb_SNPs_in_pangenome = genomes_to_nb_of_SNPs["all"]

df = pd.DataFrame(columns=["genome_index", "sensitivity"])
for index, genome_name in enumerate(snakemake.params.genomes_names):
    nb_SNPs_in_genome = genomes_to_nb_of_SNPs[genome_name]
    df.loc[index] = [index, nb_SNPs_in_genome / nb_SNPs_in_pangenome]
df.to_csv(snakemake.output.perfect_genotyper_recall_df, sep="\t")
