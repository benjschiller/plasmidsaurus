# container: "docker://plasmidsaurus/bacterial:latest"

# Define your sample names
SAMPLES = ["SRR30810013"]

# Define directories
READS_DIR = "input"                      # Directory containing FASTQ files
FILTERED_READS_DIR = "filtered"      # Directory for filtered reads
ASSEMBLY_DIR = "assemblies"          # Directory for individual assemblies
QUAST_DIR = "quast_reports"          # Output directory for QUAST reports
CHECKM_DIR = "checkm_reports"        # Output directory for CheckM reports
PROKKA_DIR = "prokka_annotations"    # Output directory for Prokka annotations

# Define the final targets
rule all:
    input:
        # expand(f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta", sample=SAMPLES)#,
        expand(f"{PROKKA_DIR}/{{sample}}/PROKKA_{{sample}}.gff", sample=SAMPLES),
        expand(f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta", sample=SAMPLES),
        expand(f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta.fai", sample=SAMPLES),
        expand(f"{QUAST_DIR}/{{sample}}/report.txt", sample=SAMPLES),
        expand(f"{CHECKM_DIR}/{{sample}}/storage/bin_stats_ext.tsv", sample=SAMPLES)

rule filtlong:
    input:
        reads = lambda wildcards: f"{READS_DIR}/{wildcards.sample}.fastq.gz"
    output:
        filtered_reads = f"{FILTERED_READS_DIR}/{{sample}}.fastq.gz"
    params:
        outdir =  f"{FILTERED_READS_DIR}"
    threads: 8
    log:
        f"logs/filtlong/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 {input.reads} | bgzip -@{threads} > {output.filtered_reads} 2>&1
        """

# Rule to assemble with Flye
rule flye_assembly:
    input:
        reads = f"{FILTERED_READS_DIR}/{{sample}}.fastq.gz"
    output:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/assembly.fasta",
        assembly_info = f"{ASSEMBLY_DIR}/{{sample}}/assembly_info.txt"
    params:
        outdir = f"{ASSEMBLY_DIR}/{{sample}}"
    threads: 8
    log:
        f"logs/{{sample}}/flye.log"
    shell:
        """
        flye --nano-raw {input.reads} --out-dir {params.outdir} --threads {threads} --genome-size 5m > {log} 2>&1
        """

# Rule to run minimap2 and racon for polishing
rule racon_polish:
    input:
        reads = f"{FILTERED_READS_DIR}/{{sample}}.fastq.gz",
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/assembly.fasta"
    output:
        polished_assembly = f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta"
    params:
        outdir = f"{ASSEMBLY_DIR}/{{sample}}/polished"
    threads: 8
    log:
        f"logs/racon/{{sample}}.log"
    shell:
        """
        2>&1
        mkdir -p {params.outdir}
        minimap2 -t {threads} -x map-ont {input.assembly} {input.reads} > {params.outdir}/overlaps.paf
        racon -t {threads} {input.reads} {params.outdir}/overlaps.paf {input.assembly} > {output.polished_assembly}
        """

rule fa_idx:
    input:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta"
    output:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta.fai"
    log:
        "logs/faidx/{{sample}}.log"
    shell:
        """
        samtools faidx {input.assembly} 2>&1
        """

# QUAST quality assessment rule
rule quast:
    input:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta"
    output:
        report = f"{QUAST_DIR}/{{sample}}/report.txt"
    params:
        outdir = f"{QUAST_DIR}/{{sample}}"
    threads: 8
    log:
        f"logs/quast/{{sample}}.log"
    shell:
        """
        quast {input.assembly} -o {params.outdir} -t {threads} > {log} 2>&1
        """

# CheckM completeness checking rule
rule checkm:
    input:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta"
    output:
        bin_stats = f"{CHECKM_DIR}/{{sample}}/storage/bin_stats_ext.tsv"
    params:
        outdir = f"{CHECKM_DIR}/{{sample}}"
    threads: 8
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
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/polished/assembly.fasta"
    output:
        gff = f"{PROKKA_DIR}/{{sample}}/PROKKA_{{sample}}.gff"
    params:
        outdir = f"{PROKKA_DIR}/{{sample}}"
    threads: 8
    log:
        f"logs/prokka/{{sample}}.log"
    shell:
        """
        # --force because prokka complains incorrectly about the output dir already existing
        prokka --force --outdir {params.outdir} --prefix PROKKA_{wildcards.sample} --cpus {threads} {input.assembly} > {log} 2>&1
        """