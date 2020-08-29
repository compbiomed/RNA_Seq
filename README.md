# RNA-seq Nextflow pipeline and associated files

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
Path to a TSV file containing the following columns:
- `INDIVIDUAL_ID`: An ID for an individual from which one or more samples was obtained
- `SAMPLE_ID`: An ID for each sample
- `LIBRARY_ID`: An ID for each library prepared from a sample
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

## Setting up a new run

1. Create a run directory, change current directory to it, and retrieve files from GitHub using the command:

   `git clone https://github.com/compbiomed/RNA_Seq`

2. Generate a tab-delimited Nextflow input file following the format described above under `params.infile`.

3. Edit the `RNA-seq_template.config` file:
   - Set `params.infile` to the full path to the tab-delimited file describing the FASTQ input files.
   - Set `params.output_dir` to the full path to the Nextflow run directory.
   - Set `params.prefix` to a meaningful name for the project.  This string will be used as a prefix to label many output files.
   - Set the fields of `params.genome` as needed depending on the species being analyzed.
   - Uncomment fields of `params.createSE.biomart_attributes` depending on the species being analyzed and the Ensembl version being used.
   - Change `params.read_length`, `params.paired_end`, and `params.stranded` if needed (rare).

4. Rename the `RNA-seq_template.config` file to something more meaningful (e.g., the value of `params.prefix` with the extension `.config`)

5. Start the Nextflow run using the qsub file as follows:

   `qsub submit_RNA_Seq.qsub [config filename]`

This will run the Nextflow script `RNA_Seq.nf` using the input file and config files from steps 2-4.
