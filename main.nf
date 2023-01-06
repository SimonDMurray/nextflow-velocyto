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
ch_sample_list = params.SAMPLEFILE   != null? Channel.fromPath(params.SAMPLEFILE)    : errorMessage() 

//Each line of the sample file is read and then emitted to its own channel which is used as input to the first process, so each sample will be ran in parallel
ch_sample_list
  .flatMap{ it.readLines() }
  .set { ch_samplelines_sf }

process get_velocyto {

  input:
  val(sample) from ch_samplelines_sf

  output:
  //outputs multiple objects: sample name (generated in the process shell so needs env) and a list of all files generated in processing
  set env(NAME), path('*') into ch_run_velocyto 

  shell:
  '''
  NAME=`echo !{sample} | cut -f 1 -d " "`
  IRODSPATH=`echo !{sample} | cut -f 2 -d " "`
  iget -f -v -N 4 -K $IRODSPATH/!{params.barcode_path}
  gunzip barcodes.tsv.gz
  mv barcodes.tsv $NAME.barcodes.tsv
  iget -f -v -N 4 -K $IRODSPATH/!{params.bam_file}
  mv !{params.bam_file} $NAME.bam
  iget -f -v -N 4 -K $IRODSPATH/!{params.index_file}
  mv !{params.index_file} $NAME.bam.bai
  '''
}

process run_velocyto {

  input:
  set NAME, val(files_list) from ch_run_velocyto

  output:
  path('*.velocyto') into ch_output

  shell:
  '''
  export LC_ALL=C.UTF-8
  export LANG=C.UTF-8
  bam_file=!{files_list[0]}
  index_file=!{files_list[1]}
  barcodes_file=!{files_list[2]}
  velocyto run \
    -t uint32 \
    --samtools-threads !{params.THREADS} \
    --samtools-memory !{params.MEM} \
    -b $barcodes_file \
    -o !{NAME}.velocyto \
    -m !{params.RMSK} \
    $bam_file \
    !{params.GTF}
  '''
}

//the directory containing the loom file is passed to this channel where it is copied to outdir
ch_output
  .subscribe {
      it.copyTo("${outdir}/")
  }
