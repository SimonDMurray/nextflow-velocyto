#!/usr/bin/env nextflow

//Using DSL1 as there is no need for workflows due to being only two processes
nextflow.enable.dsl=1

def helpMessage() {
    log.info"""
    =================
    velocyto pipeline
    =================
    This pipeline runs Velocyto. 
    The only parameters you need to input are:
      --SAMPLEFILE /full/path/to/sample/file
      --sampleID user99
    This file should be a tsv with 3 columns: SAMPLEID\t/PATH/TO/BAM/t/PATH/TO/BARCODES
    Each line should only contain information on a single sample.
    An example can be seen here: https://github.com/cellgeni/velocyto/blob/master/example-data/example.txt
    The default reference GTF used is: GRCh38_v32_filtered.gtf
    The default masked gtf used is: GRCh38_rmsk.gtf
    To change these defaults input:
      --GTF /path/to/reference/gtf
      --RMSK /path/to/mask/gtf
    """.stripIndent()
}

if (params.HELP) {
  helpMessage()
  exit 0
}

def errorMessage() {
    log.info"""
    ==============
    velocyto error
    ==============
    You failed to provide the --SAMPLEFILE input parameter
    Please provide this parameter as follows:
      --SAMPLEFILE /full/path/to/sample/file
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}

process email_startup {
  
  shell:
  '''
  contents=`cat !{params.SAMPLEFILE}`
  sendmail "!{params.sangerID}@sanger.ac.uk" <<EOF
  Subject: Launched pipeline
  From: noreply-cellgeni-pipeline@sanger.ac.uk

  Hi there, you've launched Cellular Genetics Informatics' Velocyto pipeline.
  Your parameters are:
  Samplefile: !{params.SAMPLEFILE}
  The Genome GTF file used is: !{params.GTF}
  The Mask file used is: !{params.RMSK}
  
  Your sample file looks like:
  $contents

  Thanks,
  Cellular Genetics Informatics
  EOF
  '''
}

//Puts samplefile into a channel unless it is null, if it is null then it displays error message and exits with status 1.
ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage() 

//Each line of the sample file is read and then emitted to its own set of channels, so each sample will be ran in parallel

ch_get_local_barcodes = Channel.create()
ch_get_irods_barcodes = Channel.create()
ch_get_local_bam = Channel.create()
ch_get_irods_bam = Channel.create()

ch_sample_list
  .flatMap{ it.readLines() }
  .tap( ch_get_local_barcodes )
  .tap( ch_get_irods_barcodes )
  .tap( ch_get_local_bam )
  .tap( ch_get_irods_bam )
  .set { ch_get_sample_id }

process get_local_barcodes {

  when:
  params.barcodes_on_irods == false

  input:
  val(sample) from ch_get_local_barcodes

  output:
  env(name) into ch_local_sample_id
  path('*barcodes.tsv') into ch_from_local_barcodes

  shell:
  '''
  name=`echo !{sample} | cut -f 1 -d " "`
  barcodes_path=`echo !{sample} | cut -f 3 -d " "`
  cp "${barcodes_path}" "${name}.barcodes.tsv.gz"
  gunzip "${name}.barcodes.tsv.gz"
  '''
}

process get_irods_barcodes {

  when:
  params.barcodes_on_irods == true

  input:
  val(sample) from ch_get_irods_barcodes

  output:
  env(name) into ch_irods_sample_id
  path('*barcodes.tsv') into ch_from_irods_barcodes

  shell:
  '''
  name=`echo !{sample} | cut -f 1 -d " "`
  barcodes_path=`echo !{sample} | cut -f 3 -d " "`
  iget -f -v -N 4 -K "${barcodes_path}" "${name}.barcodes.tsv.gz"
  gunzip "${name}.barcodes.tsv.gz"
  '''
}

process get_local_bam {

  when:
  params.bam_on_irods == false

  input:
  val(sample) from ch_get_local_bam

  output:
  path('*bam') into ch_from_local_bam
  path('*bam.bai') into ch_from_local_index

  shell:
  '''
  name=`echo !{sample} | cut -f 1 -d " "`
  bam_path=`echo !{sample} | cut -f 2 -d " "`
  cp "${bam_path}" "${name}.bam"
  cp "${bam_path}.bai" "${name}.bam.bai"
  '''
}

process get_irods_bam {

  when:
  params.bam_on_irods == true

  input:
  val(sample) from ch_get_irods_bam

  output:
  path('*bam') into ch_from_irods_bam
  path('*bam.bai') into ch_from_irods_index

  shell:
  '''
  name=`echo !{sample} | cut -f 1 -d " "`
  bam_path=`echo !{sample} | cut -f 2 -d " "`
  iget -f -v -N 4 -K "${bam_path}" "${name}.bam"
  iget -f -v -N 4 -K "${bam_path}.bai" "${name}.bam.bai" 
  '''
}

//ternary operators used in combination with process conditions to ensure files are grabbed from the right place
ch_run_velocyto_sample = params.barcodes_on_irods == true ? ch_irods_sample_id : ch_local_sample_id
ch_run_velocyto_barcodes = params.barcodes_on_irods == true ? ch_from_irods_barcodes : ch_from_local_barcodes
ch_run_velocyto_bam = params.bam_on_irods == true ? ch_from_irods_bam : ch_from_local_bam
ch_run_velocyto_index = params.bam_on_irods == true ? ch_from_irods_index : ch_from_local_index

process run_velocyto {

  //output velocyto files to results directory
  publishDir "/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/${params.sangerID}/${params.timestamp}/velocyto-results", mode: 'copy'

  input:
  val(name) from ch_run_velocyto_sample
  path(barcodes) from ch_run_velocyto_barcodes
  path(bam) from ch_run_velocyto_bam
  path(index) from ch_run_velocyto_index

  output:
  path('*.velocyto')
  val(name) into ch_collect

  shell:
  '''
  export LC_ALL=C.UTF-8
  export LANG=C.UTF-8
 
  mkdir !{name}.velocyto
  echo "velocyto run -t uint32 --samtools-threads !{params.THREADS} --samtools-memory !{params.MEM} -b !{barcodes} -o !{name}.velocyto -m !{params.RMSK} !{bam} !{params.GTF}" > "!{name}.velocyto/cmd.txt"

  velocyto run \
    -t uint32 \
    --samtools-threads !{params.THREADS} \
    --samtools-memory !{params.MEM} \
    -b !{barcodes} \
    -o !{name}.velocyto \
    -m !{params.RMSK} \
    !{bam} \
    !{params.GTF}
  '''
}

ch_collect
  .collect()
  .set{ ch_email_finish_sample }

process email_finish {
  
  input:
  val(name) from ch_email_finish_sample
  
  shell:
  '''
  sendmail "!{params.sangerID}@sanger.ac.uk" <<EOF
  Subject: Finished pipeline
  From: noreply-cellgeni-pipeline@sanger.ac.uk

  Hi there, your run of Cellular Genetics Informatics' Velocyto pipeline is complete.
  
  Results are available here: "/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/velocyto-results"
  
  The results will be deleted in a week so please copy your data to a sensible location, i.e.:
  cp -r "/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/velocyto-results" /path/to/sensible/location
  
  The general velocyto command run was:
  velocyto run \
    -t uint32 \
    --samtools-threads !{params.THREADS} \
    --samtools-memory !{params.MEM} \
    -b "sample-barcodes" \
    -o "sample-velocyto" \
    -m !{params.RMSK} \
    "sample-bam" \
    !{params.GTF}
  
  Each sample has the command run documented inside: "sampleID/cmd.txt"

  Thanks,
  Cellular Genetics Informatics
  EOF
  '''
}
