from pandora1_paper.evaluate.Positioned_SNPs import PositionedSNPsIndex
from utils import *
import itertools

with open(snakemake.input.length_of_longest_contig_file) as fin:
    length_of_longest_contig = int(fin.readlines()[0])

positionedSNPsIndex = PositionedSNPsIndex(length_of_longest_contig)
for genome_name_1, genome_name_2 in itertools.product(snakemake.params.genomes_names, snakemake.params.genomes_names):
    if genome_name_1 < genome_name_2:
        out_dir = Path(snakemake.config['output_folder']) / genome_name_1
        prefix = f"{out_dir / genome_name_1}{SEPARATOR}{genome_name_2}.mummer"
        positionedSNPsIndex.add_SNPs_from_csv(prefix + ".csv", genome_name_1, genome_name_2)
positionedSNPsIndex.save(snakemake.output.SNP_clusters)
