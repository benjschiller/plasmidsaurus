# container: "docker://plasmidsaurus/bacterial:latest"

# Define your sample names
SAMPLES = ["SRR30810013"]

# Define directories
READS_DIR = "input"                      # Directory containing FASTQ files
FILTERED_READS_DIR = "filtered"      # Directory for filtered reads
ASSEMBLY_DIR = "assemblies"          # Directory for individual assemblies
TRYCYCLER_DIR = "trycycler"          # Directory for Trycycler outputs
QUAST_DIR = "quast_reports"          # Output directory for QUAST reports
CHECKM_DIR = "checkm_reports"        # Output directory for CheckM reports
PROKKA_DIR = "prokka_annotations"    # Output directory for Prokka annotations

# Define the final targets
rule all:
    input:
        # expand(f"{ASSEMBLY_DIR}/{{sample}}/flye/assembly.fasta", sample=SAMPLES)#,
        expand(f"{TRYCYCLER_DIR}/{{sample}}/clusters", sample=SAMPLES),
        expand(f"{TRYCYCLER_DIR}/{{sample}}/.reconcile_complete", sample=SAMPLES)
        # expand(f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/2_all_seqs.fasta", 
        #        sample=SAMPLES, 
        #        cluster=lambda wildcards: [path.split('_')[-1] for path in get_cluster_dirs(wildcards.sample)])
            #    cluster=[path.split('_')[-1] for path in get_cluster_dirs({'sample': sample})] for sample in SAMPLES),
        # expand(f"{PROKKA_DIR}/{{sample}}/PROKKA_{{sample}}.gff", sample=SAMPLES),
        # expand(f"{QUAST_DIR}/{{sample}}/report.txt", sample=SAMPLES),
        # expand(f"{CHECKM_DIR}/{{sample}}/storage/bin_stats_ext.tsv", sample=SAMPLES)

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
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/flye/assembly.fasta"
    params:
        outdir = f"{ASSEMBLY_DIR}/{{sample}}/flye"
    threads: 8
    log:
        f"logs/{{sample}}/flye.log"
    shell:
        """
        flye --nano-raw {input.reads} --out-dir {params.outdir} --threads {threads} --genome-size 5m > {log} 2>&1
        """

# Had to comment this out because of incompatibilities with the environment (seems to be Python3-related)
# Rule to assemble with Miniasm and Minipolish
# rule miniasm_minipolish_assembly:
#     input:
#         reads = f"{FILTERED_READS_DIR}/{{sample}}.fastq.gz"
#     output:
#         assembly = f"{ASSEMBLY_DIR}/{{sample}}/miniasm/assembly.fasta"
#     params:
#         outdir = f"{ASSEMBLY_DIR}/{{sample}}/miniasm"
#     threads: 8
#     log:
#         f"logs/{{sample}}/miniasm.log"
#     shell:
#         """
#         mkdir -p {params.outdir}
#         # Generate overlaps with minimap2
#         minimap2 -x ava-ont -t {threads} {input.reads} {input.reads} > {params.outdir}/overlaps.paf
#         # Assemble with miniasm
#         miniasm -f {input.reads} {params.outdir}/overlaps.paf > {params.outdir}/miniasm.gfa
#         # Convert GFA to FASTA
#         awk '/^S/{{print ">"$2"\\n"$3}}' {params.outdir}/miniasm.gfa > {params.outdir}/unpolished_assembly.fasta
#         # Polish assembly with minipolish
#         minipolish -t {threads} {params.outdir}/unpolished_assembly.fasta {input.reads} > {output.assembly} 2>> {log}
#         """

# Rule to assemble with Raven
rule raven_assembly:
    input:
        reads = f"{FILTERED_READS_DIR}/{{sample}}.fastq.gz"
    output:
        assembly = f"{ASSEMBLY_DIR}/{{sample}}/raven/assembly.fasta"
    params:
        outdir = f"{ASSEMBLY_DIR}/{{sample}}/raven"
    threads: 8
    log:
        f"logs/{{sample}}/raven.log"
    shell:
        """
        mkdir -p {params.outdir}
        raven --threads {threads} {input.reads} > {output.assembly} 2> {log}
        """

# Run Trycycler cluster
checkpoint generate_clusters:
    input:
        flye = f"{ASSEMBLY_DIR}/{{sample}}/flye/assembly.fasta",
        # miniasm = f"{ASSEMBLY_DIR}/{{sample}}/miniasm/assembly.fasta",
        raven = f"{ASSEMBLY_DIR}/{{sample}}/raven/assembly.fasta",
        reads = f"{FILTERED_READS_DIR}/{{sample}}.fastq.gz"
    output:
        clusters_dir = directory(f"{TRYCYCLER_DIR}/{{sample}}/clusters")
    threads: 8
    log:
        f"logs/trycycler/{{sample}}_cluster.log"
    shell:
        """
        mkdir -p {TRYCYCLER_DIR}/{wildcards.sample}
        trycycler cluster --threads {threads} --assemblies {input.flye} {input.raven} --reads {input.reads} --out_dir {TRYCYCLER_DIR}/{wildcards.sample}/clusters > {log} 2>&1
        """

def collect_reconcile_files(wildcards):
    sample = wildcards.sample

    # Wait for the checkpoint to finish for this sample
    ck = checkpoints.generate_clusters.get(sample=sample)
    # n.b. current checkpoint api returns a list of outputs, not named outputs
    clusters_dir = ck.output[0]

    # Now, list the clusters in that directory
    import glob
    cluster_dirs = glob.glob(f"{clusters_dir}/cluster_*")
    clusters = [os.path.basename(path).split('_')[-1] for path in cluster_dirs]
    reconcile_files = [
        f"{TRYCYCLER_DIR}/{sample}/clusters/cluster_{cluster}/2_all_seqs.fasta"
        for cluster in clusters
    ]
    return reconcile_files
  
# Run Trycycler reconcile on each cluster. Use checkpoint instead of rule to support dynamic file creation
rule trycycler_reconcile:
    input:
        cluster_dir = f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}",
        reads = f"{FILTERED_READS_DIR}/{{sample}}.fastq.gz"
    output:
        reconcile_fasta = f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/2_all_seqs.fasta"
    threads: 4
    log:
        f"logs/trycycler/{{sample}}_cluster_{{cluster}}_reconcile.log"
    shell:
        """
        trycycler reconcile --reads {input.reads} --cluster_dir {input.cluster_dir} > {log} 2>&1
        """

# temp
rule reconcile_sample:
    input:
        reconcile_files = collect_reconcile_files
    output:
        temp(f"{TRYCYCLER_DIR}/{{sample}}/.reconcile_complete")
    shell:
        "touch {output}"

# # Rule to run Trycycler msa on each cluster
# rule trycycler_msa:
#     input:
#         cluster_dir = lambda wildcards: f"{TRYCYCLER_DIR}/{wildcards.sample}/clusters/cluster_{wildcards.cluster}"
#     output:
#         msa_done = touch(f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/3_msa.fasta")
#     threads: 4
#     log:
#         f"logs/trycycler/{{sample}}_cluster_{{cluster}}_msa.log"
#     shell:
#         """
#         trycycler msa --cluster_dir {TRYCYCLER_DIR}/{wildcards.sample}/clusters/cluster_{wildcards.cluster} > {log} 2>&1
#         touch {output.msa_done}
#         """


# # Define a function to get cluster directories
# def get_cluster_dirs(wildcards):
#     import glob
#     cluster_path = f"{TRYCYCLER_DIR}/{wildcards.sample}/clusters"
#     cluster_dirs = glob.glob(f"{cluster_path}/cluster_*")
#     return cluster_dirs
#   return reconcile_files



# rule trycycler_partition:
#     input:
#         reads = f"{FILTERED_READS_DIR}/{{sample}}.fastq.gz",
#         cluster_dirs = lambda wildcards: get_cluster_dirs(wildcards)
#     output:
#         partitioned_reads = lambda wildcards: expand(
#             f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/4_reads.fastq",
#             sample=wildcards.sample,
#             cluster=[path.split('_')[-1] for path in get_cluster_dirs(wildcards)]
#         )
#     threads: 4
#     log:
#         f"logs/trycycler/{{sample}}_partition.log"
#     shell:
#         """
#         trycycler partition --reads {input.reads} --cluster_dirs {input.cluster_dirs} > {log} 2>&1
#         """

# # Rule to run Trycycler consensus on each cluster
# rule trycycler_consensus:
#     input:
#         all_seqs = f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/2_all_seqs.fasta",
#         msa = f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/3_msa.fasta",
#         reads = f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/4_reads.fastq"
#     output:
#         final_consensus = f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/7_final_consensus.fasta"
#     threads: 4
#     log:
#         f"logs/trycycler/{{sample}}_cluster_{{cluster}}_consensus.log"
#     shell:
#         """
#         trycycler consensus --cluster_dir {TRYCYCLER_DIR}/{wildcards.sample}/clusters/cluster_{wildcards.cluster} > {log} 2>&1
#         """

# # Rule to run Medaka consensus and combine results
# rule medaka_and_combine:
#     input:
#         final_consensus = lambda wildcards: expand(
#             f"{TRYCYCLER_DIR}/{{sample}}/clusters/cluster_{{cluster}}/7_final_consensus.fasta",
#             sample=wildcards.sample,
#             cluster=[path.split('_')[-1] for path in get_cluster_dirs(wildcards)]
#         )
#     output:
#         assembly = f"{ASSEMBLY_DIR}/{{sample}}/assembly.fasta"
#     threads: 8
#     log:
#         f"logs/medaka/{{sample}}_medaka_and_combine.log"
#     shell:
#         """
#         for c in {TRYCYCLER_DIR}/{wildcards.sample}/clusters/cluster_*; do
#             medaka_consensus -i "$c/4_reads.fastq" -d "$c/7_final_consensus.fasta" -o "$c/medaka" -m r941_min_sup_g507 -t {threads}
#             mv "$c/medaka/consensus.fasta" "$c/8_medaka.fasta"
#             rm -r "$c/medaka" "$c"/*.fai "$c"/*.mmi  # clean up
#         done
#         cat {TRYCYCLER_DIR}/{wildcards.sample}/clusters/cluster_*/8_medaka.fasta > {output.assembly}
#         """

# # QUAST quality assessment rule
# rule quast:
#     input:
#         assembly = f"{ASSEMBLY_DIR}/{{sample}}/assembly.fasta"
#     output:
#         report = f"{QUAST_DIR}/{{sample}}/report.txt"
#     params:
#         outdir = f"{QUAST_DIR}/{{sample}}"
#     threads: 4
#     log:
#         f"logs/quast/{{sample}}.log"
#     shell:
#         """
#         quast {input.assembly} -o {params.outdir} -t {threads} > {log} 2>&1
#         """

# # CheckM completeness checking rule
# rule checkm:
#     input:
#         assembly = f"{ASSEMBLY_DIR}/{{sample}}/assembly.fasta"
#     output:
#         bin_stats = f"{CHECKM_DIR}/{{sample}}/storage/bin_stats_ext.tsv"
#     params:
#         outdir = f"{CHECKM_DIR}/{{sample}}"
#     threads: 4
#     log:
#         f"logs/checkm/{{sample}}.log"
#     shell:
#         """
#         mkdir -p {params.outdir}
#         checkm lineage_wf -x fasta {ASSEMBLY_DIR}/{wildcards.sample} {params.outdir} -t {threads} > {log} 2>&1
#         """

# # Prokka annotation rule
# rule prokka:
#     input:
#         assembly = f"{ASSEMBLY_DIR}/{{sample}}/assembly.fasta"
#     output:
#         gff = f"{PROKKA_DIR}/{{sample}}/PROKKA_{{sample}}.gff"
#     params:
#         outdir = f"{PROKKA_DIR}/{{sample}}"
#     threads: 4
#     log:
#         f"logs/prokka/{{sample}}.log"
#     shell:
#         """
#         prokka --outdir {params.outdir} --prefix PROKKA_{wildcards.sample} --cpus {threads} {input.assembly} > {log} 2>&1
#         """