# interactive build
sudo singularity build --sandbox ubuntu_18_04/ docker://ubuntu:18.04
sudo singularity shell ubuntu_18_04/

# recipe build
sudo singularity build leoisl-perfect_single_genome_caller_in_a_pangenome.sif Singularity.recipe
