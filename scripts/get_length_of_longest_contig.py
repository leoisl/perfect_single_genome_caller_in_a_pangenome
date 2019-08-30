import pysam
def get_length_of_longest_contig(path):
    length_of_longest_contig = 0
    with pysam.FastxFile(path) as fh:
        for entry in fh:
            length_of_longest_contig = max(length_of_longest_contig, len(entry.sequence))
    return length_of_longest_contig


length_of_longest_contig = 0
for genome in snakemake.input.all_genomes:
    length_of_longest_contig = max(length_of_longest_contig, get_length_of_longest_contig(genome))
with open(snakemake.output.length_of_longest_contig_file, "w") as fout:
    print(length_of_longest_contig, file=fout)
