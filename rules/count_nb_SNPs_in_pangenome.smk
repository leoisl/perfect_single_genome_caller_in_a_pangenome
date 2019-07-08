rule run_dnadiff:
    input:
        genome_1_path = Path(config["input_folder"]) / "{genome_1}.fna",
        genome_2_path = Path(config["input_folder"]) / "{genome_2}.fna"
    params:
        output_prefix = lambda wildcards, output: Path(output.delta_file).with_suffix('')
    output:
        delta_file = Path(config["output_folder"]) / "{genome_1}.{genome_2}.delta"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{genome_1}_{genome_2}_run_dnadiff.log"
    shell:
        "dnadiff {input} -p {params.output_prefix} 2> {log}"


rule run_show_snps:
    input:
        delta_file = Path(config["output_folder"]) / "{genome_1}.{genome_2}.delta"
    output:
        snp_probes = Path(config["output_folder"]) / "{genome_1}.{genome_2}.show_snps"
    params:
        probe_length = config["probe_length"] #TODO: make this as part of the input
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    log:
        "logs/{genome_1}_{genome_2}_run_show_snps.log"
    shell:
        "show-snps -C -H -I -r -T -x 10 {input} > {output} 2> log"

#header: [P1]    [SUB]   [SUB]   [P2]    [BUFF]  [DIST]  [CTX R] [CTX Q] [FRM]   [TAGS]


rule concatenate
