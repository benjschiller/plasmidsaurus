# Define your sample names
SAMPLES = ["SRR30810013"]
container: "plasmidsaurus/bacterial:latest"

# Define directories
READS_DIR = "."                      # Directory containing FASTQ files
ASSEMBLY_DIR = "assemblies"          # Output directory for SPAdes assemblies
QUAST_DIR = "quast_reports"          # Output directory for QUAST reports
CHECKM_DIR = "checkm_reports"        # Output directory for CheckM reports
PROKKA_DIR = "prokka_annotations"    # Output directory for Prokka annotations

# Define the final targets
rule all:
    input:
        expand(f"{PROKKA_DIR}/{{sample}}/PROKKA_{{sample}}.gff", sample=SAMPLES),
        expand(f"{QUAST_DIR}/{{sample}}/report.txt", sample=SAMPLES),
        expand(f"{CHECKM_DIR}/{{sample}}/storage/bin_stats_ext.tsv", sample=SAMPLES)

rule cutadapt:
    input:
        reads = lambda wildcards: f"{READS_DIR}/{wildcards.sample}.fastq.gz"
    output:
        trimmed_reads = f"{READS_DIR}/trimmed/{{sample}}.fastq.gz"
    log:
        f"logs/cutadapt/{{sample}}.log"
    shell:
        """
        cutadapt -o {output.trimmed_reads} {input.reads} > {log} 2>&1
        """

# SPAdes assembly rule
rule spades:
    input:
        reads = f"{READS_DIR}/trimmed/{{sample}}.fastq.gz"
    output:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta"
    params:
        outdir = f"{ASSEMBLY_DIR}/{{sample}}"
    threads: 8
    log:
        f"logs/spades/{{sample}}.log"
    shell:
        """
        spades.py --iontorrent -s {input.reads} -o {params.outdir} -t {threads} > {log} 2>&1
        """

# QUAST quality assessment rule
rule quast:
    input:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta"
    output:
        report = f"{QUAST_DIR}/{{sample}}/report.txt"
    params:
        outdir = f"{QUAST_DIR}/{{sample}}"
    threads: 4
    log:
        f"logs/quast/{{sample}}.log"
    shell:
        """
        quast {input.assembly} -o {params.outdir} -t {threads} > {log} 2>&1
        """

# CheckM completeness checking rule
rule checkm:
    input:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta"
    output:
        bin_stats = f"{CHECKM_DIR}/{{sample}}/storage/bin_stats_ext.tsv"
    params:
        outdir = f"{CHECKM_DIR}/{{sample}}"
    threads: 4
    log:
        f"logs/checkm/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        checkm lineage_wf -x fasta {ASSEMBLY_DIR}/{wildcards.sample} {params.outdir} -t {threads} > {log} 2>&1
        """

# Prokka annotation rule
rule prokka:
    input:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta"
    output:
        gff = f"{PROKKA_DIR}/{{sample}}/PROKKA_{{sample}}.gff"
    params:
        outdir = f"{PROKKA_DIR}/{{sample}}"
    threads: 4
    log:
        f"logs/prokka/{{sample}}.log"
    shell:
        """
        prokka --outdir {params.outdir} --prefix PROKKA_{wildcards.sample} --cpus {threads} {input.assembly} > {log} 2>&1
        """