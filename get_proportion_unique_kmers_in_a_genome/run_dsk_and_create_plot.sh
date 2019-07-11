nb_cores=8
file=GCF_000005845.2_ASM584v2_genomic.fna
tsv_file=kmersize_percentageuniquekmers_${file}.tsv
pdf_file=kmersize_percentageuniquekmers_${file}.pdf

echo "kmer_size percentage_unique_kmers" > $tsv_file

for kmer_size in {5..51..1}
do
  outfile=dsk_out_${file}_k_${kmer_size}
  dsk/build/bin/dsk -nb-cores $nb_cores -kmer-size $kmer_size -abundance-min 0 -histo 1 -out ${outfile} -file ${file}
  awk -v kmer_size="$kmer_size" 'BEGIN{unique=0; all_kmers=0} {all_kmers+=$2; if ($1==1){unique+=$2}} END{print kmer_size, unique/all_kmers*100}' \
  ${outfile}.histo >> $tsv_file
done

#create the plot
Rscript create_plot.R $tsv_file $pdf_file
