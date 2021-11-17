# ctg-wgs 
## nextflow pipeline for WGS analysis with Illumina Dragen server

## Input files
The following files must be in the runfolder to start pipeline successfully.

1. Samplesheet (`CTG_SampleSheet.ctg-wgs.csv`)

### Samplesheet requirements

- The samplesheet format is standard IEM sheet, with the following modifications:
- It can also use just a simple csv containting only the [Data] section (In this case adapters will not be trimmed). 

#### I. Additional columns added in Data-table (table coming after [Data]):

| Column | Supported values |
| ------ | -------- |
| Sample_Ref | hg38 / hg19 / mm10 : hg38, hg19 and mm10 are currently set up for dragen. Note that annotation is not supported with mm10 |
| somatic | y / n : Set to "y" if it is somatic samples (e.g. tumor-only). If set to "n", VC will be germline mode. |

- Also note that no Sample_Name should be added. Leave that column blank!

### Samplesheet template (.csv)

Samplesheet name: `CTG_SampleSheet.ctg-wgs.csv`

```
[Header]
IEMFileVersion,5
Date,2021-04-29
Workflow,GenerateFASTQ
Application,NovaSeq FASTQ Only
Instrument Type,NovaSeq
"Index Adapters,""IDT- UD Indexes (96 Indexes)"""
Chemistry,Amplicon

[Reads]
151
151

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Sample_ref,somatic
1_NA24385_300ng,,,,A01,UDP0001,GAACTGAGCG,UDP0001,CGCTCCACGA,2021_060,hg38,y
2_NA24385_300ng,,,,B01,UDP0002,AGGTCAGATA,UDP0002,TATCTTGTAG,2021_060,hg38,y
3_NA24385_300ng,,,,C01,UDP0003,CGTCTCATAT,UDP0003,AGCTACTATA,2021_060,hg38,y
4_NA24385_300ng,,,,D01,UDP0004,ATTCCATAAG,UDP0004,CCACCAGGCA,2021_060,hg38,y
5_NA24385_300ng,,,,E01,UDP0005,GACGAGATTA,UDP0005,AGGATAATGT,2021_060,hg38,y
6_NA24385_1000ng,,,,F01,UDP0006,AACATCGCGC,UDP0006,ACAAGTGGAC,2021_060,hg38,y
7_NA24385_1000ng,,,,G01,UDP0007,CTAGTGCTCT,UDP0007,TACTGTTCCA,2021_060,hg38,y
8_NA24385_1000ng,,,,H01,UDP0008,GATCAAGGCA,UDP0008,ATTAACAAGG,2021_060,hg38,y
9_NA24385_1000ng,,,,A02,UDP0009,GACTGAGTAG,UDP0009,CACTATCAAC,2021_060,hg38,y
10_NA24385_1000ng,,,,B02,UDP0010,AGTCAGACGA,UDP0010,TGTCGCTGGT,2021_060,hg38,y 
```
## Running without demux (run with existing fastq file):

1. Run wgs-driver with -d flag (turns demux off) : `wgs-driver -d`
2. Samplesheet (Still called CTG_SampleSheet.ctg-wgs.csv) will have different format:
- It needs header with `metaid` (name of project/run that will be used to generate "project log folder" etc) and `fastqpath` (which should point to the directory that contain the fastq files).
- [Data] Section only needs Sample_ID,Sample_Project,Sample_ref,somatic. Same rules as described above applies here.

**Example samplesheet for demux-off analysis**
```
metaid,2020_test_giab_truseq
fastqpath,/projects/fs1/medpvb/proj/wgs/giab/2020_test_giab_CMD_truseq/fastq/
[Data]
Sample_ID,Sample_Project,Sample_ref,somatic
ALL503A1255_24-114140,2020_test_giab_truseq,hg38,y
```

## The following steps are performed by the pipeline:

* `Demultiplexing` (dragen bcl-conversion): Converts raw basecalls to fastq, and demultiplex samples based on index. Adapters sequences are trimmed if added to samplesheet [Settings].
* `Alignment` and `Variant calling` (dragen align + calling): Fastq from each sample is aligned to the reference genome. Variants (SNV+SV+CNV) are called.
* `Annotation` (nirvana): Clinical grade annotation of all variants passing basic filter in Dragen. 
* `Dragen metrics`: Compiling Dragen alignment and coverage metrics to table.
* `MultiQC`: Summarizes FastQC and Dragen metrics into one document (https://multiqc.info/).

## Output:
* CTG-output
    * `fastq`: Contains raw fastq files from demultiplexing.
    * `dragen`: Output from dragen alignment + variant calling. Here you find the .bam files for the aligned sample and .vcf files. The *hard-filtered.vcf contain the PASS variants only.
    * `qc`: Quality control output. 
      * `multiqc`: Summary of metrics 
      * `dragen`: text files with dragen metrics (compiled to the multiqc report)

## Specs
- Dragen version: 3.9
    - User guide: (https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/GPipelineIntro_fDG.htm). 
    - Relase notes: (https://jp.support.illumina.com/content/dam/illumina-support/documents/downloads/software/dragen/1000000157666_00_DRAGEN-3.8.4-Customer-Release-Notes.pdf)

## Container
- `ngs-tools` Singularity container contain NGS-related tools, embedded in the repo: 
https://github.com/perllb/ctg-wgs/tree/master/container 

## Requirements
- Slurm (recommended)
- Dragen Server

## Reference genomes
- Reference genomes supported: hg38/hg19 / mm10




