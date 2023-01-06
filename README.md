# nextflow-velocyto
Nextflow version of cellgeni/velocyto for testing on Nextflow Tower

`examples/RESUME-velocyto` - script that executes velocyto nextflow pipeline, has 3 hardcoded arguments (for testing sake):
* `/path/to/sample/file`
* `/path/to/config/file`
* `/path/to/nextflow/script`

`nextflow.config` - the configuration script that allows the processes to be submittede to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path). Some global default parameters are also set in this file with absolute paths.

`main.nf` - the pipeline script that executes velocyto

`examples/irods.txt` - samplefile tsv containing 2 fields, sampleID followed by IRODS path to data
