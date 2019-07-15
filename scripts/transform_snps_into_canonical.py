import sys
from utils import rev_comp

def canonical_snp(snp_1, snp_2):
    snp_1_RC = rev_comp(snp_1)
    snp_2_RC = rev_comp(snp_2)
    all_possible_snps = [f"{snp_1}${snp_2}", f"{snp_2}${snp_1}", f"{snp_1_RC}${snp_2_RC}", f"{snp_2_RC}${snp_1_RC}"]
    return min(all_possible_snps)

for line in sys.stdin:
    line_split = line.strip().split()
    if len(line_split)==0: break
    snp_1,snp_2 = line_split[6], line_split[7]
    print(canonical_snp(snp_1, snp_2))
