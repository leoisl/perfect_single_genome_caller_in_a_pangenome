#getting ecoli reference
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

#installing dsk
git clone --recursive https://github.com/GATB/dsk.git && cd dsk && bash INSTALL

#creating the plot
bash run_dsk_and_create_plot.sh
