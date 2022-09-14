# RNA-seq Nextflow pipeline and associated files
The pipeline can run in one of two modes: 
   - Operend - This will retrieve source files from and store the results in the Operend System. FastQ input files are retrieved from the Operend system and the pipeline results are posted back to it. The group of FastQ input files is identifid by the `params.seq_read_set` value in the nextflow config file. To activate reading from and writing to Operend, set `params.source_operend=true` in the nextflow config file.
   - File System - This operates exclusively on files on the file system. The FastQ input files and the `params.infile` .tsv input file must reside on the file system before running the pipeline and the results are written to the directory described by `params.output_dir` in the nextflow config file. To use the file system, set the `params.source_operend=false` in the nextflow config file.
## Setting up a new run

1. Create a run directory, change current directory to it, and retrieve files from GitHub using the command:

   `git clone https://github.com/Operend/RNA_Seq`

2. Rename the `RNA-seq_template.config` file to something more meaningful (e.g., the value of `params.prefix` with the extension `.config`)
3. Edit the renamed `.config` file:
   - Set `params.source_operend` to false if reading input files and writing output files to the file system. Set to true to source input files from and post results back to Operend.
   - Operend Specific Params:
      - Set `params.seq_read_set_name` to the name of the IlluminaSeqReadsSet that contain the FASTQ files as inputs. This is only applicable when sourcing the input FastQ files from Operend.
      - Set `params.opapi_base_url` to the url of the Operend Server - "http://operend.bu.edu/api". This is ignored if not using Operend and only using the File System to read and write files.
      - Set `params.opapi_verify_https` to false if the Operend service is not operating over HTTPS.  This is ignored if not using Operend and only using the File System to read and write files.
      - Set `params.job_type_name` to the name of the Operend Job Type for this nextflow run - "CBM_RNASeq" currently. This JobType MUST exist in Operend. There is no reason CMB_RNASeq JobType should not exist in Operend but the pipeline will not run against Operend if it does not.
      - Set `params.opapi_storage_location` to the name of the Operend Storage Location where the results should be stored. If not set, the User's default location will be used. 
   - Set `params.infile` to the full path to the tab-delimited file describing the FASTQ input files. Only applicable when reading input FastQ files from the file system and `params.source_operend` is false. This will be ignored if using Operend.
   - Set `params.output_dir` to the full path to the Nextflow run directory. If not using Operend, this will be where the final output of the pipeline resides. If using Operend, this will serve as a staging directory until the files are posted to the system and can be deleted after.
   - Set `params.prefix` to a meaningful name for the project.  This string will be used as a prefix to label many output files and act as a label/identifier for the Operend JobRun created by this pipeline run.
   - Set the fields of `params.genome` as needed depending on the species being analyzed.
   - Uncomment fields of `params.createSE.biomart_attributes` depending on the species being analyzed and the Ensembl version being used.
   - Change `params.read_length`, `params.paired_end`, and `params.stranded` if needed (rare).



 

4. Start the Nextflow run using the qsub file as follows:

   `qsub submit_RNA_Seq.qsub [config filename]`

This will run the Nextflow script `RNA_Seq.nf` using the input file and config files from steps 2-3.

## Description of files

### `RNA_Seq.nf`
This Nextflow pipeline contains the following processes:
#### FASTQ-level QC
- `runFastQC`: Use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to perform QC on each pair of FASTQ files
- `runMultiQCFastq`: Compile output from `runFastQC` into TSV tables and HTML report using [MultiQC](https://multiqc.info/)
#### Alignments
- `runSTAR1pass`: Perform a first-pass alignment to a specified genome using [STAR](https://github.com/alexdobin/STAR)
- `runSTARgenomeGenerate`: Create a new genome reference from splice junctions inferred from first-pass STAR alignment
- `runSTAR2pass`: Perform a second-pass alignment to the genome reference produced by `runSTARgenomeGenerate`
#### BAM QC, performed using [RSeQC](http://rseqc.sourceforge.net)
- `runRSeQC_bam_stat`
- `runRSeQC_clipping_profile`
- `runRSeQC_deletion_profile`
- `runRSeQC_geneBody_coverage` *Note: this is a very slow step*
- `runRSeQC_infer_experiment`
- `runRSeQC_inner_distance`
- `runRSeQC_insertion_profile`
- `runRSeQC_junction_annotation`
- `runRSeQC_junction_saturation`
- `runRSeQC_read_distribution`
- `runRSeQC_read_duplication`
- `runRSeQC_read_GC`
- `runRSeQC_read_NVC`
- `runRSeQC_read_quality`
- `runRSeQC_tin` *Note: this is a very slow step*
#### Other BAM-level QC
- `runMultiQCSample`: Compile output from RSeQC and STAR using MultiQC
#### Post-alignment processing
- `runRSEM`: Use [RSEM](https://deweylab.github.io/RSEM/) to estimate gene- and transcript (isoform)-level expression
- `createSE`: Create [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) R objects, one for gene-level expression and one for transcript-level expression; each object contains several estimates of expression from RSEM, as well as feature-level annotation

### `RNA_Seq_template.config`
This configuration file is intended to be used only with this Nextflow script.  It makes a number of assumptions about underlying directory structures and filenames.
The parameters that are typically changed are:

#### `params.infile`
Path to a TSV file containing the following columns (only necassary if working off of the file system and not Operend):
- `INDIVIDUAL_ID`: An ID for an individual from which one or more samples was obtained.
- `SAMPLE_ID`: An ID for each sample.  **This field cannot be left blank or set to `NA`.**
- `LIBRARY_ID`: An ID for each library prepared from a sample.
- `RG_ID`: Read Group ID: the flowcell ID, optionally followed by a lane-specific suffix (for instruments with independent lanes)
- `PLATFORM_UNIT`: The RG ID, followed by a suffix specific to a sample/library
- `PLATFORM`: Sequencing platform, e.g., "illumina" for Illumina instruments
- `PLATFORM_MODEL`: Instrument model, e.g., "NextSeq", "HiSeq", etc.)
- `RUN_DATE`: Optional run date
- `CENTER`: Optional name for center at which sequencing was performed
- `R1`: Full path to FASTQ file containing first paired-end read
- `R2`: Full path to FASTQ file containing second paired-end read

Note: some of these fields are discussed in more detail within the [GATK read groups documentation](https://software.broadinstitute.org/gatk/documentation/article.php?id=6472).

#### `params.output_dir`
Path where the Nextflow output should be written

#### `params.prefix`
Prefix for Nextflow output files

#### `params.paired_end`
Indicates whether paired-end sequencing was used (`true` or `false`)

#### `params.genome`
Parameters specific to the genome and annotation build used; these can be uncommented and edited as needed

### `submit_RNA_Seq.qsub`
This script kicks off the Nextflow process on the SGE using the .config file specified in its sole argument.

## Accessing Result Files from Operend using FUSE
Pipeline result files within Operend can be accessed as files on the file system via a FUSE mount using Operend's entityfs package. The entityfs package can be run from within a Singularity container or within a Conda environment. Note that although the Singularity container also allows for use /on any computer where Singularity is installed, the Environment Modules will not be available from within the container.

1. copy FUSE/entityfs_scc_template.ini to FUSE/entityfs.ini
2. edit FUSE/entityfs.ini
      - Set the api_token_secret field 
      - Set the url to the (if not set)
      - Set the rules (link to FPR documentation please)
         - If using the default rules (rna_seq_res_per_sample,rna_seq_res_per_job) please ensure the File Path Rules exist on the system
         - Use an array of JSON rules
         - Use names of exsiting File Path Rules
3. Ensure that a directory to mount the entityfs file system exists. Lets say `/home/username/mountpoint`. Please don't copy that verbatim.
#### FUSE with Singularity
   
   - Start a Singularity container instance using the Singularity `instance` command and the use the singularity `shell` to access that instance. Make sure to replace the paths with real paths. 
      - ```singularity instance start --fusemount "container:/opt/operend/entityfs/entityfs.py --config {PATH_TO_THIS_DIR}/RNA_Seq_Op/FUSE/entityfs.ini /home/username/mountpoint" library://operend/libs/entityfs.sif {NAME_OF_INSTANCE}```
      - ```singularity shell instance://{NAME_OF_INSTANCE}```
      - from inside the container (assuming Singularity mounts the user's home directory) - ```cd /home/username/mountpoint```
   
   - OR

   - Start a Singularity container instance and shell in to that container using one singularity `shell` command
      - ```singularity shell --fusemount "container:/opt/operend/entityfs/entityfs.py --config {PATH_TO_THIS_DIR}/RNA_Seq_Op/FUSE/entityfs.ini /home/username/mountpoint" library://operend/libs/opyrnd.sif {NAME_OF_INSTANCE}```
      - from inside the container (assuming Singularity mounts the user's home directory) - ```cd /home/username/mountpoint```

#### FUSE with entityfs and conda (from the FUSE directory)
1) Load the miniconda module
   - ```module load miniconda```
2) Create and activate conda envronment
   - ```conda env create -f entityfs.yml```
   - ```conda activate entityfs``` 
3) Copy entityfs from github (sorry its not deployed to conda-forget yet)
   - ```wget    https://raw.githubusercontent.com/Operend/entityfs/main/entityfs.py```
   - ```chmod +x entityfs.py```
4) Mount it
   - ```./entityfs.py --config entityfs.ini mountpoint &```
5) Enjoy your files!!
   - ```cd mountpoint``` 
6) Unmount (please do so when done or you'll have issues)   
   - ```fusermount3 -u mountpoint```

## Clean Up
When the pipeline saves the data to the Operend system it creates a JobRun. That JobRun generates many associated WorkFiles
and Entities. If unwanted JobRuns are left laying around, much storage space will be wasted. Deleting this JobRun and its 
associated data can be very tedious to do manually. Use the handy_scripts/job_run_cleanup.py script to accomplish this.
#### Delete the JobRun and associated data - use with caution as this will PERMANENTLY delete all data generated by the pipeline
Execute the Python script using Singularity:

```singularity exec library://operend/libs/opyrnd.sif:latest python3 job_run_cleanup.py -id {JOB_RUN_ID} -c {operend_config.ini}```

Where the JOB_RUN_ID is the system ID of the JobRun and operend_config.ini is Operend config file.

