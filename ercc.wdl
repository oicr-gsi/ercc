version 1.0

struct rccResources {
 String name
 String refBwaModule
 String refFasta
 String refGtf
 String refModule
 String controlData
}


workflow ercc {
input {
  File fastqR1
  File? fastqR2
  String outputFileNamePrefix = basename(fastqR1, '.fastq.gz')
  String mixId
  String sampleId
  String reference
}

Map[String, rccResources] resources = {
 "hg19": {
   "name": "hg19",
   "refBwaModule": "hg19-ercc-bwa-index/0.7.17",
   "refFasta": "$HG19_ERCC_BWA_INDEX_ROOT/hg19_random_ercc.fa",
   "refGtf": "$HG19_ERCC_ROOT/ERCC92.gtf",
   "refModule": "hg19-ercc/p13",
   "controlData": "$HG19_ERCC_ROOT/ERCC_Controls_Analysis_v2.txt"
 },
 "hg38": {
   "name": "hg38",
   "refBwaModule": "hg38-ercc-bwa-index/0.7.17",
   "refFasta": "$HG38_ERCC_BWA_INDEX_ROOT/hg38_random_ercc.fa",
   "refGtf": "$HG38_ERCC_ROOT/ERCC92.gtf",
   "refModule": "hg38-ercc/p12",
   "controlData": "$HG38_ERCC_ROOT/ERCC_Controls_Analysis_v2.txt"
 }
}

call countTotal {
  input:
    fastqR1 = fastqR1
}

call countTranscripts { 
  input:
    fastqR1 = fastqR1,
    fastqR2 = fastqR2,
    refGenome = resources[reference].refFasta,
    modules = "bwa/0.7.17 samtools/0.1.19 ~{resources[reference].refBwaModule}"
}

Array[File] mateFiles = select_all([fastqR1,fastqR2])

call rpkmTable {
  input:
    erccCounts = countTranscripts.erccCounts,
    totalReads = if length(mateFiles) == 1 then countTotal.firstMateCounts else countTotal.firstMateCounts*2,
    modules  = resources[reference].refModule,
    erccData = resources[reference].refGtf
}

call makeReport {
  input:
    rpkmTable = rpkmTable.rpkmTable,
    prefix = outputFileNamePrefix,
    samples = if mixId == "ERCC Mix 1" then "~{sampleId} NA" else "NA ~{sampleId}",
    controlData = resources[reference].controlData,
    erccData = resources[reference].refGtf,
    modules  = "rstats/4.0 ercc-scripts/1.1 ~{resources[reference].refModule}"
}

parameter_meta {
  fastqR1: "File with reads for mate 1 or fastq file for single-read data"
  fastqR2: "File with reads for mate 2, is available"
  reference: "hg19 or hg38 (assembly id)"
  outputFileNamePrefix: "Prefix for the output file"
  mixId: "LIMS-approved identification of the Spike-In Mix"
  sampleId: "The Id used to identify the data in the report"
}

output {
  File rpkmData = rpkmTable.rpkmTable
  File image = makeReport.reportImage
  File pdf   = makeReport.reportPdf
  File json  = makeReport.reportJson
}

meta {
    author: "Peter Ruzanov"
    email: "peter.ruzanov@oicr.on.ca"
    description: "Ercc 1.1, Niassa-wrapped Cromwell (widdle) workflow for running ERCC spike-in analysis on RNAseq data. The workflow counts total reads in fastq files, evaluate expression of ERCC spike-ins using bwa aligner and creates a report in .json format, along with plots in .png and .pdf formats. The .json file contains RPKMs for all ERCC transcripts as well as other information. RPKM table in .csv format is also provisioned."
    dependencies: [
      {
        name: "rstats/4.0",
        url: "https://www.r-project.org/"
      },
      {
        name: "bwa/0.7.17",
        url: "http://bio-bwa.sourceforge.net/"
      },
      {
        name: "samtools/0.1.19",
        url: "http://www.htslib.org/"
      }
    ]
    output_meta: {
    rpkmData: {
        description: "ERCC readouts in RPKMs",
        vidarr_label: "rpkmData"
    },
    image: {
        description: "png with a plot (Dose Response supported at the moment)",
        vidarr_label: "image"
    },
    pdf: {
        description: "pdf with the same plot as png, a legacy output",
        vidarr_label: "pdf"
    },
    json: {
        description: "json file with ERCC numbers",
        vidarr_label: "json"
    }
}
}

}

# ====================================================
# TASK 1 of 4: count number of total reads using fastq
# ====================================================
task countTotal {
input {
  File fastqR1
  Int jobMemory = 8
  Int timeout = 20
}

parameter_meta {
  fastqR1: "Read 1 fastq file, we report the number of reads in it"
  timeout: "Timeout in hours for this task"
  jobMemory: "Memory allocated to this job"
}

command <<<
  zcat ~{fastqR1} | paste - - - - | wc -l  
>>>

runtime {
  memory:  "~{jobMemory} GB"
  timeout: "~{timeout}"
}

output {
  Int firstMateCounts = read_int(stdout())
}
}


# ===================================================
#  TASK 2 of 4: calculate ERCC counts, make csv
#       use piping from alignment with bwa
# ===================================================
task countTranscripts {
input {
  File fastqR1
  File? fastqR2
  Int timeout = 20
  Int jobMemory = 20
  Int threads = 8
  String refGenome
  String cnv_file = "ercc_counts.csv"
  String modules
}

parameter_meta {
  fastqR1: "File with reads for mate 1 or fastq file for single-read data"
  fastqR2: "File with reads for mate 2, is available"
  refGenome: "path to fasta file for genome"
  jobMemory: "Memory allocated to this job"
  timeout: "Timeout in hours for this task"
  threads: "Threads to use with bwa"
  cnv_file: "Output, contains ERCC ids and their respective number of reads"
  modules: "Names and versions of modules needed for alignment"
}

command <<<
  bwa mem -t ~{threads} -M ~{refGenome} ~{fastqR1} ~{fastqR2} | samtools view -S - | \
      cut -f 3 | grep ERCC | sort | uniq -c | sed s/^\ *// > ~{cnv_file}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  cpu: "~{threads}"
  timeout: "~{timeout}"
  modules: "~{modules}"
}

output {
  File erccCounts = "${cnv_file}"
}
}

# ===============================================
#  TASK 3 of 4: count RPKMs using external data
# ===============================================
task rpkmTable {
input {
  File erccCounts
  Int totalReads
  Int timeout = 20
  Int jobMemory   = 10
  String modules 
  String erccData 
}

parameter_meta {
  erccCounts: ".csv file with counts for ERCC transcripts"
  totalReads: "Total reads fot current sample, used to calculate RPKMs" 
  jobMemory: "Memory allocated to sort task"
  timeout: "Timeout in hours for this task"
  modules: "Names and versions of modules needed for alignment"
  erccData: "Reference file from Agilent"
}

command <<<
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
>>>

runtime {
  memory:  "~{jobMemory} GB"
  timeout: "~{timeout}"
  modules: "~{modules}"
}

output {
  File rpkmTable = "~{basename(erccCounts, '.csv')}_rpkm.csv"
}
}

# ===========================================
#  TASK 4 of 4: make an image and json with R
# ===========================================
task makeReport {
input {
  File rpkmTable
  String imagingScriptPath = "$ERCC_SCRIPTS_ROOT/ercc_plots.R"
  String controlData
  String erccData 
  String prefix = "UnnamedReport"
  String samples
  String rScript = "$RSTATS_ROOT/bin/Rscript"
  Int jobMemory = 10
  String modules 
}

parameter_meta {
  rpkmTable:  "Input file with RPKMshost .bam file"
  controlData: "ERCC supporting data"
  erccData: "Reference file from Agilent"
  jobMemory: "Memory allocated to classify task"
  modules: "Names and versions of modules needed for making report"
  imagingScriptPath: "path to R script ercc_plots.R"
  rScript: "Path to Rscript command"
  prefix: "prefix to use with images"
  samples: "space-separated mix1 and mix2 samples"
}

command <<<
 ~{rScript} ~{imagingScriptPath} ~{rpkmTable} ~{controlData} ~{erccData} ~{prefix} DoseResponse ~{samples}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File reportPdf    = "~{prefix}_DoseResponse.pdf"
  File reportImage  = "~{prefix}_DoseResponse.png"
  File reportJson   = "~{prefix}_rpkm.json"
}
}
