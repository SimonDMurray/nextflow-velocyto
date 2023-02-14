#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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

process get_data {

  input:
  val(sample) 

  output:
  env(NAME), emit: name 
  path('*barcodes.tsv'), emit: barcodes  
  path('*bam'), emit: bam 
  path('*bam.bai'), emit: index 

  shell:
  '''
  NAME=`echo !{sample} | cut -f 1 -d " "`
  bam_path=`echo !{sample} | cut -f 2 -d " "`
  barcodes_path=`echo !{sample} | cut -f 3 -d " "`
  
  if [[ "!{params.bam_on_irods}" == "no" ]]; then
    cp "${bam_path}" "${NAME}.bam"
    cp "${bam_path}.bai" "${NAME}.bam.bai"
  elif [[ "!{params.bam_on_irods}" == "yes" ]]; then
    iget -f -v -K "${bam_path}" "${NAME}.bam"
    iget -f -v -K "${bam_path}.bai" "${NAME}.bam.bai"
  else
    echo "incorrect bam option"
    exit 1
  fi

  if [[ "!{params.barcodes_on_irods}" == "no" ]]; then
    cp "${barcodes_path}" "${NAME}.barcodes.tsv.gz"
    gunzip -f "${NAME}.barcodes.tsv.gz"
  elif [[ "!{params.barcodes_on_irods}" == "yes" ]]; then
    iget -f -v -K "${barcodes_path}" "${NAME}.barcodes.tsv.gz"
    gunzip -f "${NAME}.barcodes.tsv.gz"
  else
    echo "incorrect barcodes option"
    exit 1
  fi
  '''
}

process run_velocyto {

  //output velocyto files to results directory
  publishDir "/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/${params.sangerID}/${params.timestamp}/velocyto-results", mode: 'copy'

  input:
  val(NAME) 
  path(barcodes) 
  path(bam) 
  path(index) 

  output:
  path('*.velocyto')

  shell:
  '''
  export LC_ALL=C.UTF-8
  export LANG=C.UTF-8
 
  velocyto_cmd="velocyto run"
  if [[ "!{params.bam_has_umis}" == "no" ]]; then
    velocyto_cmd="velocyto run -U"
  fi

  mkdir -p !{NAME}.velocyto
  echo "${velocyto_cmd} -t uint32 --samtools-threads !{params.THREADS} --samtools-memory !{params.MEM} -b !{barcodes} -o !{NAME}.velocyto -m !{params.RMSK} !{bam} !{params.GTF}" > "!{NAME}.velocyto/cmd.txt"

  $velocyto_cmd \
    -t uint32 \
    --samtools-threads !{params.THREADS} \
    --samtools-memory !{params.MEM} \
    -b !{barcodes} \
    -o !{NAME}.velocyto \
    -m !{params.RMSK} \
    !{bam} \
    !{params.GTF}
  '''
}

process email_finish {
  
  input:
  val(NAME) 
  
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

workflow {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    //Puts samplefile into a channel unless it is null, if it is null then it displays error message and exits with status 1.
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    if (params.sangerID == null) {
      errorMessage()
    }
    else {
      email_startup()
      ch_sample_list | flatMap{ it.readLines() } | get_data 
      run_velocyto(get_data.out.name, get_data.out.barcodes, get_data.out.bam, get_data.out.index) \
      | collect \
      | email_finish
    }
  }
}
