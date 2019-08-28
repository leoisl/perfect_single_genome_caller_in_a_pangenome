library(ggplot2)

#args = c("../assemblies_sample_out/perfect_genotyper.dataframe")
args = commandArgs(trailingOnly=TRUE)
pdf(args[2])

  data = read.table(args[1], header=TRUE)
  data$genome_index <- factor(data$genome_index, levels = data$genome_index[order(data$sensitivity)])
  ggplot(data=data, aes(x=genome_index, y=sensitivity, group=1)) + geom_line() + geom_point() + scale_y_continuous(limits = c(0, 1), breaks=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

dev.off()