import pandas as pd

general_df = None
for df_filename in snakemake.input.perfect_genotyper_recall_dfs:
    df = pd.read_csv(df_filename, sep="\t")
    if general_df is None:
        general_df = df
    else:
        general_df = general_df.append(df, ignore_index=True)
general_df.to_csv(snakemake.output.concatenated_df, sep="\t", index=False)