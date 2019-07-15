sort pipeline_output/SNP_refined_panel.fa.refined_clusters.debug | uniq -c | awk '{print $1, $2}' | sort -nk2,2 | cat <(echo "freq cluster_size") -
