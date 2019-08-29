# interactive build
sudo singularity build --sandbox ubuntu_18_04/ docker://ubuntu:18.04
sudo singularity shell -B /home/leandro/git/perfect_single_genome_caller_in_a_pangenome:/perfect_single_genome_caller_in_a_pangenome ubuntu_18_04/

# recipe build
sudo singularity build -B /home/leandro/git/perfect_single_genome_caller_in_a_pangenome:/perfect_single_genome_caller_in_a_pangenome leoisl-perfect_single_genome_caller_in_a_pangenome.sif Singularity.recipe
