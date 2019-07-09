#!/usr/bin/env bash
CLUSTER_CMD=("bsub -n {threads} -R \"select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1]\" -M {resources.mem_mb} -o {cluster.output} -e {cluster.error} -J {cluster.name}")
JOB_NAME=snakemake_master_process
LOG_DIR=logs/
MEMORY=4000

    snakemake --use-conda \
    --cluster-config cluster.yaml \
    --jobs 2000 \
    --restart-times 3 \
    --cluster "${CLUSTER_CMD[@]}" \
    -p -n

exit 0
