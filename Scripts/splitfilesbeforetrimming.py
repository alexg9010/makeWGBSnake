rule fastqc_after_trimming_pe:
  input:
  # DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz",
  # DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz"
  DIR_trimmed+"{sample}/{sample}_1.fq.part-{part}_val_1.fq.gz",
DIR_trimmed+"{sample}/{sample}_2.fq.part-{part}_val_2.fq.gz"
output:
  DIR_posttrim_QC+"{sample}/{sample}_1.fq.part-{part}_val_1_fastqc.html",
DIR_posttrim_QC+"{sample}/{sample}_1.fq.part-{part}_val_1_fastqc.zip",
DIR_posttrim_QC+"{sample}/{sample}_2.fq.part-{part}_val_2_fastqc.zip",
DIR_posttrim_QC+"{sample}/{sample}_2.fq.part-{part}_val_2_fastqc.html"
params:
  fastqc_args = config['args']['fastqc'],
outdir = "--outdir "+DIR_posttrim_QC + "{sample}/"
log:
  DIR_posttrim_QC+"{sample}/{sample}_{part}_trimmed_fastqc.log"
message:
  "Quality checking trimmmed paired-end data from {input}"
shell:
  "{tools}/fastqc {params} {input} > {log} 2> {log}.err"


rule trim_reads_pe:
  input:
  #qc    = [ DIR_rawqc+"{sample}_1_fastqc.html",
  #          DIR_rawqc+"{sample}_2_fastqc.html"],
  #files = [ inputdir+"{sample}_1.fq",
  #          inputdir+"{sample}_2.fq"]
  files = [ inputdir+"{sample}_1.fq.part-{part}",
            inputdir+"{sample}_2.fq.part-{part}"]
output:
  DIR_trimmed+"{sample}/{sample}_1.fq.part-{part}_val_1.fq.gz",
DIR_trimmed+"{sample}/{sample}_2.fq.part-{part}_val_2.fq.gz"
params:
  extra          = config['args']['trim_galore'],
outdir         = "--output_dir "+DIR_trimmed+"{sample}/",
phred          = "--phred33",
gz             = "--gzip",
cutadapt       = "--path_to_cutadapt " + tools +"cutadapt",
paired         = "--paired"
log:
  DIR_trimmed+"{sample}/{sample}.log"
message:
  "Trimming raw paired-end read data from {input}"
shell:
  "{tools}/trim_galore {params} {input.files}"


rule split_intput_into_pieces:
  input: inputdir+"{sample}_{ext}.fq"
output: [inputdir+"{sample}_{ext}.fq.part-"+i for i in ["1","2","3","4"]]
shell: "/fast/users/kwreczy_m/programs/bin/fastq-splitter.pl {input} --n-parts {NPARTS} --check"


rule fastqc_raw:
  input:
  inputdir+"{sample}_{ext}.fq"
output:
  DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.html",
DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.zip"
params:
  fastqc_args = config['args']['fastqc'],
outdir = "--outdir "+ DIR_rawqc
log:
  DIR_rawqc+"{sample}/{sample}_{ext}.log"
shell:
  config["tools"]+"fastqc {params} {input} > {log} 2> {log}.err"


rule gunzip:
  input:
  inputdir+"{sample}_{ext}.fq.gz"
output:
  inputdir+"{sample}_{ext}.fq"
shell: "gunzip {input}"



