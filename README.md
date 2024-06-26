# ercc

Ercc 1.1, Niassa-wrapped Cromwell (widdle) workflow for running ERCC spike-in analysis on RNAseq data. The workflow counts total reads in fastq files, evaluate expression of ERCC spike-ins using bwa aligner and creates a report in .json format, along with plots in .png and .pdf formats. The .json file contains RPKMs for all ERCC transcripts as well as other information. RPKM table in .csv format is also provisioned.

## Overview

## Dependencies

* [rstats 4.0](https://www.r-project.org/)
* [bwa 0.7.17](http://bio-bwa.sourceforge.net/)
* [samtools 0.1.19](http://www.htslib.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run ercc.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastqR1`|File|File with reads for mate 1 or fastq file for single-read data
`mixId`|String|LIMS-approved identification of the Spike-In Mix
`sampleId`|String|The Id used to identify the data in the report
`reference`|String|hg19 or hg38 (assembly id)


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR2`|File?|None|File with reads for mate 2, is available
`outputFileNamePrefix`|String|basename(fastqR1,'.fastq.gz')|Prefix for the output file


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`countTotal.jobMemory`|Int|8|Memory allocated to this job
`countTotal.timeout`|Int|20|Timeout in hours for this task
`countTranscripts.timeout`|Int|20|Timeout in hours for this task
`countTranscripts.jobMemory`|Int|20|Memory allocated to this job
`countTranscripts.threads`|Int|8|Threads to use with bwa
`countTranscripts.cnv_file`|String|"ercc_counts.csv"|Output, contains ERCC ids and their respective number of reads
`rpkmTable.timeout`|Int|20|Timeout in hours for this task
`rpkmTable.jobMemory`|Int|10|Memory allocated to sort task
`makeReport.imagingScriptPath`|String|"$ERCC_SCRIPTS_ROOT/ercc_plots.R"|path to R script ercc_plots.R
`makeReport.rScript`|String|"$RSTATS_ROOT/bin/Rscript"|Path to Rscript command
`makeReport.jobMemory`|Int|10|Memory allocated to classify task


### Outputs

Output | Type | Description | Labels
---|---|---|---
`rpkmData`|File|ERCC readouts in RPKMs|vidarr_label: rpkmData
`image`|File|png with a plot (Dose Response supported at the moment)|vidarr_label: image
`pdf`|File|pdf with the same plot as png, a legacy output|vidarr_label: pdf
`json`|File|json file with ERCC numbers|vidarr_label: json


## Commands
This section lists command(s) run by ercc workflow
 
* Running ercc
 
ercc workflow checks out spike-in signal (reads) and plots the QC graph(s).
 
### Count reads in fastq file:
 
```
   zcat FASTQ_R1 | paste - - - - | wc -l  
 
```
### Align reads to chimeric reference (modified reference also containing all spike-in sequences)
 
```
   bwa mem -t THREADS -M REF_GENOME FASTQ_R1 FASTQ_R2 | samtools view -S - | 
       cut -f 3 | grep ERCC | sort | uniq -c | sed s/^\ *// > CNV_FILE
```
 
###  In this custom script we build the RPKM table

``` 
  python <<CODE
  import os
  ercc = os.path.expandvars("~{erccData}")
  counts = "~{erccCounts}"
  total = ~{totalReads}
  output_file = "~{basename(erccCounts, '.csv')}_rpkm.csv"
 
  # 1. Open file, read and collect length data
  ercc_data = open(ercc, 'r')
  length_hash = {}
  for nextLine in ercc_data.readlines():
      nextFields = nextLine.split("\t")
      length_hash[nextFields[0]] = int(nextFields[4])
  ercc_data.close()
 
  # 2. Calculate RPKMs, put them in a hash
  count_data = open(counts, 'r')
  rpkm_hash = {}
  count_hash = {}
 
  for nextLine in count_data.readlines():
      nextLine = nextLine.strip("\n")
      countChunk = nextLine.split(" ")
      if len(countChunk) == 2 and countChunk[1] in length_hash.keys():
          count_hash[countChunk[1]] = int(countChunk[0])
 
  pm = int(total)/1000000.0
  for e in sorted(length_hash.keys()):
     rpkm = 0.0
     if e in count_hash.keys():
         rpkm = count_hash[e]/pm/(length_hash[e]/1000.0)
     rpkm_hash[e] = "{0:.2f}".format(rpkm)
  count_data.close()
 
  # 3. print out RPKMS sorted by
  f = open(output_file, "w+")
  for e in sorted(rpkm_hash.keys()):
      f.write('\t'.join([e, rpkm_hash[e]]) + '\n')
  f.close()
  CODE
 
```
 
### Report Making using oputputs from the other steps:
 
```
  Rscript IMAGING_SCRIPT RPKM_TABLE CONTROL_DATA ERCC_DATA PREFIX DoseResponse SAMPLES
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
