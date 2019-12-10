version 1.0

workflow ercc {
input {
  File fastqR1
  File? fastqR2
  String? outputFileNamePrefix = basename(fastqR1, '.fastq.gz')
  String mixId
  String sampleId
}

call countTotal {
  input:
    fastqR1 = fastqR1
}

call countTranscripts { 
  input:
    fastqR1 = fastqR1,
    fastqR2 = fastqR2
}

Array[File] mateFiles = select_all([fastqR1,fastqR2])

call rpkmTable {
  input:
    erccCounts = countTranscripts.erccCounts,
    totalReads = if length(mateFiles) == 1 then countTotal.firstMateCounts else countTotal.firstMateCounts*2
}

call makeReport {
  input:
    rpkmTable = rpkmTable.rpkmTable,
    prefix = outputFileNamePrefix,
    samples = if mixId == "ERCC Mix 1" then "~{sampleId} NA" else "NA ~{sampleId}"
}


output {
  File rpkmData = rpkmTable.rpkmTable
  File image = makeReport.reportImage
  File json  = makeReport.reportJson
}

meta {
    author: "Peter Ruzanov"
    email: "peter.ruzanov@oicr.on.ca"
    description: "Ercc 1.0"
}

}

# ====================================================
# TASK 1 of 4: count number of total reads using fastq
# ====================================================
task countTotal {
input {
  File fastqR1
  Int? jobMemory = 8
}

parameter_meta {
  fastqR1: "Read 1 fastq file, we report the number of reads in it"
}

command <<<
  zcat ~{fastqR1} | paste - - - - | wc -l  
>>>

runtime {
  memory:  "~{jobMemory} GB"
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
  Int? jobMemory = 20
  Int? threads = 8
  String? refGenome = "$HG19_ERCC_BWA_INDEX_ROOT/hg19_random_ercc.fa"
  String? cnv_file = "ercc_counts.csv"
  String? modules = "bwa/0.7.17 samtools/0.1.19 hg19-ercc-bwa-index/0.7.17"
}

parameter_meta {
  fastqR1: "File with reads for mate 1 or fastq file for single-read data"
  fastqR2: "File with reads for mate 2, is available"
  refGenome: "path to fasta file for genome"
  jobMemory: "Memory allocated to this job"
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
  Int? jobMemory   = 10
  String? modules  = "hg19-ercc/p13"
  String? erccData = "$HG19_ERCC_ROOT/ERCC92.gtf"
}

command <<<
 python <<CODE
 ercc = "~{erccData}"
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

parameter_meta {
  erccData: "Reference file from Agilent"
  erccCounts: ".csv file with counts for ERCC transcripts"
  totalReads: "Total reads fot current sample, used to calculate RPKMs" 
  jobMemory: "Memory allocated to sort task"
}

runtime {
  memory:  "~{jobMemory} GB"
}

output {
  File rpkmTable = "~{basename(erccCounts, '.csv')}_rpkm.csv"
}
}

# =======================================
#  TASK 4 of 4: make an image with R
# =======================================
task makeReport {
input {
  File rpkmTable
  String? imagingScriptPath = "$ERCC_SCRIPTS_ROOT/ercc_plots.R"
  String? controlData = "$HG19_ERCC_ROOT/ERCC_Controls_Analysis_v2.txt"
  String? erccData = "$HG19_ERCC_ROOT/ERCC92.gtf"
  String? prefix = "UnnamedReport"
  String samples
  String? Rscript = "$RSTATS_ROOT/bin/Rscript"
  Int? jobMemory = 10
  String? modules  = "rstats/3.6 ercc-scripts/1.0 hg19-ercc/p13"
}

# TODO: test with header-less rpkm table, assign columnames using passed arguments
command <<<
 ~{Rscript} ~{imagingScriptPath} ~{rpkmTable} ~{controlData} ~{erccData} ~{prefix} DoseResponse ~{samples}
>>>

parameter_meta {
  rpkmTable:  "Input file with RPKMshost .bam file"
  controlData: "ERCC supporting data"
  erccData: "Reference file from Agilent"
  jobMemory: "Memory allocated to classify task"
  modules: "Names and versions of modules needed for making report"
  imagingScriptPath: "path to Rscript ercc_plots.R"
  prefix: "prefix to use with images"
  samples: "space-separated mix1 and mix2 samples"
}

runtime {
  memory:  "~{jobMemory} GB"
}

output {
  File reportImage  = "~{prefix}_DoseResponse.pdf"
  File reportJson   = "~{prefix}_rpkm.json"
}
}
