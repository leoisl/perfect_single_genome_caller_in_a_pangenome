import sys
from pysam import FastxFile
import networkx as nx

###################################################################
# HELPERS
###################################################################
#RC function
def RC_base (base):
    assert base in ["A", "C", "G", "T"]
    if base=="A": return "T"
    if base == "T": return "A"
    if base == "C": return "G"
    if base == "G": return "C"


#class that will represent all versions of a canonical SNP
class SNP:
    def __init__(self, id, canonical_upper_path, canonical_lower_path):
        self.canonical_upper_path = canonical_upper_path
        self.canonical_lower_path = canonical_lower_path
        self.id = id

        #build the four SNPs this canonical SNP represents
        canonical_upper_path_SNP_base = canonical_upper_path.sequence[int(len(canonical_upper_path.sequence) / 2)]
        canonical_lower_path_SNP_base = canonical_lower_path.sequence[int(len(canonical_lower_path.sequence) / 2)]
        self.the_four_SNPs =\
            set([(canonical_upper_path_SNP_base, canonical_lower_path_SNP_base), \
             (canonical_lower_path_SNP_base, canonical_upper_path_SNP_base), \
             (RC_base(canonical_upper_path_SNP_base), RC_base(canonical_lower_path_SNP_base)), \
             (RC_base(canonical_lower_path_SNP_base), RC_base(canonical_upper_path_SNP_base))])

    def match(self, other):
        return len(self.the_four_SNPs.intersection(other.the_four_SNPs)) > 0

    def __str__(self):
        return f"{self.canonical_upper_path}\n{self.canonical_lower_path}"

###################################################################
# HELPERS
###################################################################

#load all SNPs
index=0
all_SNPs = []
with FastxFile(sys.argv[2]) as snps_fasta_file:
    first_allele = None
    for record in snps_fasta_file:
        if first_allele is None:
            first_allele = record
        else:
            second_allele = record

            all_SNPs.append(SNP(index, first_allele, second_allele))
            index += 1
            first_allele = None


#load unrefined_clusters
unrefined_clusters = {}
with open(sys.argv[1]) as unrefined_clusters_file:
    for line in unrefined_clusters_file:
        line_split = line.strip().split()
        cluster_1, cluster_2 = all_SNPs[int(line_split[0])], all_SNPs[int(line_split[1])]

        if cluster_1 not in unrefined_clusters:
            unrefined_clusters[cluster_1] = []

        unrefined_clusters[cluster_1].append(cluster_2)



#refine the clusters
refined_clusters = {}
for cluster_leader, cluster_group in unrefined_clusters.items():
    refined_clusters[cluster_leader] = []
    for cluster_candidate in cluster_group:
        if cluster_leader.match(cluster_candidate):
            refined_clusters[cluster_leader].append(cluster_candidate)


#build the graph on the refined clusters
G = nx.Graph()
for cluster_leader in refined_clusters:
    G.add_node(cluster_leader.id)

for cluster_leader, cluster_group in refined_clusters.items():
    for cluster_candidate in cluster_group:
        G.add_edge(cluster_leader.id, cluster_candidate.id)

for component in nx.connected_components(G):
    max_degree_of_nodes_in_this_component = max([G.degree(node) for node in component])
    nodes_with_max_degree = [node for node in component if G.degree(node) == max_degree_of_nodes_in_this_component]
    print(all_SNPs[nodes_with_max_degree[0]])
