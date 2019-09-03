library(ggplot2)

pdf(snakemake@output[["plot"]])
  data = read.table(snakemake@input[["concatenated_df"]], header=TRUE)
  ggplot(data=data, aes(x=dataset, y=sensitivity)) + geom_violin() + geom_boxplot(width=0.1)
dev.off()