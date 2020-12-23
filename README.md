# quality_control_CellRanger

Snakemake workflow to run some basic quality control programs.

Programs run:

* Fastqc
* RSeQC:
  * read_GC.py
  * bam_stat.py
  * tin.py
  * geneBody_coverage.py
* MultiQC
* DropletUtils - short script is used to run emptyDrops, whose results might be used for subsequent cell filtering

## Usage

1. Install Anaconda if you do not have, install FastQC and MultiQC and make it available on your $PATH, or edit the right rules in workflow/rules/bam_quality_control.smk
2. Edit samples.yml to match you CellRanger output
3. Download the right BED reference files for RSeQC from [here](https://sourceforge.net/projects/rseqc/files/BED/), for eg.  Human_Homo_Sapiens/hg38_RefSeq.bed.gz. Unpack it, make sure that chromosome naming matches your bam files. If not - correct it. Edit config.yml 
3. Edit workflow/cluster.yml and run_snakemake.sh according to your needs.
