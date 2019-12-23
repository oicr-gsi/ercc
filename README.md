# ercc
workflow for analysis of samples with RNA spike-ins (ERCC)

ercc workflow, analysis of spike-ins in RNA-seq data. The workflow counts total reads in fastq files, evaluate expression of ERCC spike-ins using bwa aligner and creates a report in pdf format. Workflow provisions .json file with rpkm information for all ERCC transcripts, .pdf with a plot and rpkm table in .csv format.

![ercc flowchart](docs/ercc_flowchart.png)

## Usage

## Cromwell

``` 
 java -jar cromwell.jar run ercc.wdl --inputs inputs.json 
```

## Running Pipeline

```
 
 task countTotal:   procedure for counting reads in fastq file

 task countTranscripts: procedure for counting ERCC transcripts, counts bwa-aligned reads with chimeric reference

 task rpkmTable: this task produces rpkm counts using outputs from previous tasks

 task makeReport: making report and formatting results in json format

```

The workflow will accept paired or single-end sequencing data in fastq format. Assemblies for hg19 and hg38 are both supported. 

## Optional Assembly-specific Parameters:

hg19-specific data, for other assemblies these should be changed:

Paramter|Value
---|---
refGenome | String? (optional, default = "$HG19_ERCC_BWA_INDEX_ROOT/hg19_random_ercc.fa")
erccData | String? (optional, default = "$HG19_ERCC_ROOT/ERCC92.gtf")

## Other Parameters with default values:

Paramter|Value
---|---
cnv_file | String? (optional, default = "ercc_counts.csv")
jobMemory | Int? (optional, default = 8 or 10 for computationally-intensive tasks)
threads | Int? (optional, default = 8)
imagingScriptPath | String? (optional, default = "$ERCC_SCRIPTS_ROOT/ercc_plots.R")
ercc.makeReport.Rscript | String? (optional, default = "$RSTATS_ROOT/bin/Rscript")

## Required Inputs:

Paramter|Value
---|---
outputFileNamePrefix | String? (optional, will be a basename of the first input fastq if not customized)
mixId | String - need to specify Mix type as in MISO (ERCC Mix 1 or 2)
fastqR1 | File - first (or the only in case of single-end data) fastq file
sampleId | String - id of the input sample
fastqR2 | File? (optional) - second fastq file, optional in case of single-end data

## Outputs

```
  rpkmData  - rpkm numbers for a sample in csv format
  json  - rpkm data in json format
  image  - plot in .pdf format

```
