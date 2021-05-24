# ctg-wgs 
## nextflow pipeline for WGS analysis with Illumina Dragen server

### The following steps are performed by the pipeline:

* `Demultiplexing` (dragen bcl-conversion): Converts raw basecalls to fastq, and demultiplex samples based on index. Adapters sequences are trimmed if added to samplesheet [Settings].
* `Alignment` and `Variant calling` (dragen align + calling): Fastq from each sample is aligned to the reference genome. Variants (SNV+SV) are called (CNV will soon be supported).
* `Dragen metrics`: Compiling Dragen alignment and coverage metrics to table.
* `MultiQC`: Summarizes FastQC and Dragen metrics into one document (https://multiqc.info/).

### Output:
* CTG-output
    * `fastq`: Contains raw fastq files from demultiplexing.
    * `dragen`: Output from dragen alignment + variant calling. Here you find the .bam files for the aligned sample and .vcf files. The *hard-filtered.vcf contain the PASS variants only.
    * `qc`: Quality control output. 
      * `multiqc`: Summary of metrics 
      * `dragen`: text files with dragen metrics (compiled to the multiqc report)

### Specs
- Dragen version: 3.8 
    - User guide: (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf). 
    - Relase notes: (https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/dragen/Illumina-DRAGEN-Bio-IT-Platform-3.7-Release-Notes-1000000142362-v00.pdf)


### Requirements
- Slurm (recommended)
- Dragen Server
- multiQC (https://multiqc.info/)

- The ngs-tools container contain the multiqc among other NGS-related tools https://github.com/perllb/ctg-containers/tree/main/ngs-tools

### Reference genomes
- Reference genomes supported: hg38/hg37 / mm10




