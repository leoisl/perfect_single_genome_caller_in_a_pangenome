import os
from utils import *
from pandora1_paper.evaluate.evaluation import *

for genome_path, genome_name in zip(snakemake.params.genomes, snakemake.params.genomes_names):
    if snakemake.wildcards.genome_1 < genome_name:
        out_dir = Path(snakemake.config['output_folder']) / snakemake.wildcards.genome_1
        os.makedirs(out_dir, exist_ok=True)
        prefix = f"{out_dir / snakemake.wildcards.genome_1}{SEPARATOR}{genome_name}.mummer"
        snps_as_StringIO = generate_mummer_snps(reference=Path(snakemake.input.ref_genome), query=Path(genome_path),
                                                prefix=Path(prefix), flank_width=snakemake.params.probe_length)
        snps_dataframe = ShowSnps.to_dataframe(snps_as_StringIO)
        snps_dataframe = snps_dataframe.translate_to_FWD_strand()
        snps_dataframe.to_csv(prefix + ".csv", sep="\t")
