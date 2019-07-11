library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

pdf(args[2])
data = read.table(args[1], header=TRUE)
ggplot(data=data, aes(x=kmer_size, y=percentage_unique_kmers, group=1)) + geom_line() + geom_point() + scale_y_continuous(limits = c(0, 100))
dev.off()
