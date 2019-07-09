import sys

with open(sys.argv[1]) as fin:
    nb_all_unique_canonical_snps = int(fin.readline())

print("genome sensitivity")

with open(sys.argv[2]) as fin:
    for line in fin:
        line_split = line.strip().split()
        print(f"{line_split[0][line_split[0].rindex('/')+1:]} {float(line_split[1])/nb_all_unique_canonical_snps}")
