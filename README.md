# About

Plots how good would be a 100%-recall and 100%-precision (perfect) single-genome caller to call variants in a pangenome in comparison to [pandora](https://github.com/rmcolq/pandora/).

The perfect single-genome caller is based on getting SNPs with [Mummer](http://mummer.sourceforge.net/):
1. All pairs of genomes are compared and SNPs are inferred using [Mummer](http://mummer.sourceforge.net/);
2. SNPs from different pairs of genomes are then clustered together if they represent the same SNP.
Such SNP clusters have thus a pangenomic coordinate (all the positions each SNP in the cluster appear).
Calling any SNP in a cluster translates into calling all SNPs (since they represent the same SNP).
Thus these clusters are the pangenomic SNPs to be found;
3. The number of clusters (pangenomic SNPs) found by each genome is then computed; 

For pandora, however, we don't have positioned SNPs. As such, we get probes from SNPs found by pandora and from the SNP clusters from Mummer to infer of pandora found or not the pangenomic SNPs. 

# Dependencies

```
python 3.6+
Snakemake 5
Singularity
```

# Rulegraph

![Rulegraph](rulegraph.png)
