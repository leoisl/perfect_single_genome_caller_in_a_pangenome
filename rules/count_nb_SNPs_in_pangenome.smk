rule run_dnadiff:
    input:
        genome_1_path = Path(config["input_folder"]) / "{genome_1}.fna",
        genome_2_path = Path(config["input_folder"]) / "{genome_2}.fna"
    params:
        output_prefix = lambda wildcards, output: Path(output.delta_file).with_suffix('')
    output:
        delta_file = Path(config["output_folder"]) / f"{{genome_1}}{SEPARATOR}{{genome_2}}.delta"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    conda:
        "../envs/global.yaml"
    log:
        "logs/{genome_1}_{genome_2}_run_dnadiff.log"
    shell:
        "dnadiff {input} -p {params.output_prefix} 2> {log}"


rule run_show_snps:
    input:
        delta_file = Path(config["output_folder"]) / f"{{genome_1}}{SEPARATOR}{{genome_2}}.delta"
    output:
        snp_probes = Path(config["output_folder"]) / f"{{genome_1}}{SEPARATOR}{{genome_2}}.show_snps"
    params:
        probe_length = config["probe_length"] #TODO: vary several probe lengths?
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    conda:
        "../envs/global.yaml"
    log:
        "logs/{genome_1}_{genome_2}_run_show_snps.log"
    shell:
        "show-snps -C -H -I -r -T -x {params.probe_length} {input} > {output} 2> log"
        #output file header: [P1]    [SUB]   [SUB]   [P2]    [BUFF]  [DIST]  [CTX R] [CTX Q] [FRM]   [TAGS]



rule transform_SNPs_into_canonical_SNPs:
    input:
        snp_probes = Path(config["output_folder"]) / f"{{genome_1}}{SEPARATOR}{{genome_2}}.show_snps"
    output:
        canonical_snps = Path(config["output_folder"]) / f"{{genome_1}}{SEPARATOR}{{genome_2}}.canonical_snps"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    conda:
        "../envs/global.yaml"
    log:
        "logs/{genome_1}_{genome_2}_transform_SNPs_into_canonical_SNPs.log"
    shell:
        "python scripts/transform_snps_into_canonical.py < {input} > {output} 2> log"

rule get_unique_canonical_SNPs:
    input:
        all_canonical_snps = expand( str(Path(config["output_folder"]) / f"{{genomes_1}}{SEPARATOR}{{genomes_2}}.canonical_snps"), genomes_1=genomes_names, genomes_2=genomes_names)
    output:
        all_unique_canonical_snps = Path(config["output_folder"]) / "all_unique_canonical_snps"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    conda:
        "../envs/global.yaml"
    log:
        "logs/get_unique_canonical_SNPs.log"
    shell:
        "sort {input} --parallel={threads} | uniq > {output} 2> log"

rule get_number_unique_canonical_SNPs:
    input:
        all_unique_canonical_snps = Path(config["output_folder"]) / "all_unique_canonical_snps"
    output:
        nb_all_unique_canonical_snps = Path(config["output_folder"]) / "nb_all_unique_canonical_snps"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: config["mem_mb"][attempt-1]
    conda:
        "../envs/global.yaml"
    log:
        "logs/get_number_unique_canonical_SNPs.log"
    shell:
        "wc -l {input} | awk '{{print $1}}' > {output} 2> log"
