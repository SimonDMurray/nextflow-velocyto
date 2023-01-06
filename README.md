# nextflow-velocyto
Nextflow version of cellgeni/velocyto for testing on Nextflow Tower

`RESUME-velocyto` - script that executes velocyto nextflow pipeline, has 3 hardcode arguments (for testing sake):
* `/path/to/sample/file`
* `/path/to/config/file`
* `/path/to/nextflow/script`

`velocyto.config` - the configuration script that allows the processes to be submittede to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path)

`velocyto.nf` - the pipeline script that executes velocyto

`irods.sh` - samplefile tsv containing 2 fields, sampleID followed by IRODS path to data
