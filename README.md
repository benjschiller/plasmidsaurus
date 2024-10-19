# Plasmidsaurus Bacterial Assembly and Annotation Pipeline

This pipeline can be run inside the provided Docker container.

The following docker commands use `--platform linux/amd64` to maximize compatibility with other hosts (such as M1 MacBook Pros); if running on a native linux platform, it may be omitted.

To build the docker container, run something like

    docker build --platform linux/amd64 -t plasmidsaurus/bacterial .

The samples are explicitly listed in the Snakefile under `SAMPLES`. SRR30810013 is hardcoded, and the pipeline expects the sample reads to be present at "input/SRR30810013.fastq.gz"; it can be downloaded from [SRA](https://www.ncbi.nlm.nih.gov/sra)

For cross-platform compatibility, the pipeline should be run inside the container:

    docker run --platform linux/amd64 -v .:/data -w /data --rm plasmidsaurus/bacterial:latest time snakemake -c8 --verbose

This command must be run from inside the repo's root directory (same level as README.md). The number of CPUs is specified above as `-c8` (this assumes 8). Change this if you have a different number available.

The resulting GFF file will be available as prokka_annotations/SRR30810013/PROKKA_SRR30810013.gff (or named according to the sample, if using a different sample)

There are additional outputs that may be useful:

* The final polished genome and its index are at assemblies/SRR30810013/polished/assembly.fasta(.fai)
* reports for QUAST and CheckM at quast_reports and checkm_reports, respectively
