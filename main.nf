#!/usr/bin/env nextflow

//Using DSL1 as there is no need for workflows due to being only two processes
nextflow.enable.dsl=1

def helpMessage() {
    log.info"""
    =================
    velocyto pipeline
    =================
    This pipeline runs Velocyto. 
    The only parameter you need to input is:
      --SAMPLEFILE /full/path/to/sample/file
    This file should be a tsv with 2 columns: SAMPLEID\tIRODSPATH
    Each line should only contain information on a single sample.
    An example can be seen here: https://github.com/cellgeni/velocyto/blob/master/example-data/irods.txt
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
  publishDir "$params.outdir"

  input:
  val(name) from ch_run_velocyto_sample
  path(barcodes) from ch_run_velocyto_barcodes
  path(bam) from ch_run_velocyto_bam
  path(index) from ch_run_velocyto_index

  output:
  path('*.velocyto')

  shell:
  '''
  export LC_ALL=C.UTF-8
  export LANG=C.UTF-8
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
