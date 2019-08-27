library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
#args = c("../assemblies_sample_out/perfect_genotyper.dataframe")

pdf(args[2])

  data = read.table(args[1], header=TRUE)
  data$genome_index <- factor(data$genome_index, levels = data$genome_index[order(data$sensitivity)])
  ggplot(data=data, aes(x=genome_index, y=sensitivity, group=1)) + geom_line() + geom_point() + scale_y_continuous(limits = c(0, 1))

dev.off()