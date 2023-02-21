# RNA-seq Nextflow pipeline and associated files

## Setting up a new run

1. Create a run directory, change current directory to it, and retrieve files from GitHub to it using the command:

   `git clone https://github.com/compbiomed/RNA_Seq ./`

2. Generate a tab-delimited Nextflow input file following the format described below under `params.infile`.

3. Edit the `RNA-seq_template.config` file:
   - Set `params.infile` to the full path to the tab-delimited file describing the FASTQ input files.
   - Set `params.output_dir` to the full path to the Nextflow run directory.
   - Set `params.prefix` to a meaningful name for the project.  This string will be used as a prefix to label many output files.
   - Set the fields of `params.genome`:
      - `species`: The scientific name of the species being used (e.g., `"Homo sapiens"`, `"Mus musculus"`)
      - `ucsc`: The UCSC build corresponding to the FASTA reference that will be used (e.g., `"hg38"`, `"mm10"`), as generated, for example, by [make_ucsc_references.qsub](https://github.com/compbiomed/genome-reference-scripts/blob/master/make_ucsc_references.qsub)
      If an Ensembl genome reference is to be used instead, either remove this field or leave it set to `""`.
      - `assembly`: The corresponding genome assembly (e.g., `"GRCh38"`, `"GRCm38"`)
      - `set`: The subset of sequences that will be used: `"base"` (autosomes, sex chromosomes, and mitochondrial chromosome), `"base_random"` (base sequences plus random/unplaced contigs), or `"base_random_althap"` (base and random sequences plus alternative haplotype sequences); see [make_ucsc_references.qsub](https://github.com/compbiomed/genome-reference-scripts/blob/master/make_ucsc_references.qsub) for more details
      - `ensembl`: The Ensembl build number that will be used (e.g., `100`)
   - Change `params.read_length`, `params.paired_end`, and `params.stranded` if needed (rare).

4. Rename the `RNA-seq_template.config` file to something more meaningful (e.g., `[params.prefix].config`)

5. Start the Nextflow run using the qsub file as follows:

   `qsub -P [SGE project name] submit_RNA_Seq.qsub [config filename]`

This will run the Nextflow script `RNA_Seq.nf` using the input file and config files from steps 2-4.

## Description of files

### `RNA_Seq.nf`
This Nextflow pipeline contains the following processes:
#### Setup
- `generateGTF`:
   - Retrieve a GTF file from Ensembl for the species, genome assembly, and Ensembl build number specified in the .config file
   - Use the UCSC `chromAlias` table for the corresponding UCSC genome build to convert the sequence names (i.e., chromosomes) from Ensembl to UCSC nomenclature
   - Limit the result to the subset of sequences specified in `params.genome.set` (e.g., `"base_random"`)
- `generateBED`: Use the UCSC command-line utilities `gtfToGenePred` and `genePredToBed` to convert the GTF file to a BED file (for use with [RSeQC](http://rseqc.sourceforge.net) below)
- `runRSEMprepareReference`: Prepare a set of [RSEM](https://deweylab.github.io/RSEM/) reference files using the FASTA reference specified in the .config file and the Ensembl
#### FASTQ-level QC
- `runFastQC`: Use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to perform QC on each pair of FASTQ files
- `runMultiQCFastq`: Compile output from `runFastQC` into TSV tables and HTML report using [MultiQC](https://multiqc.info/)
#### Alignments
- `runSTAR1pass`: Perform a first-pass alignment to a specified genome using [STAR](https://github.com/alexdobin/STAR)
- `runSTARgenomeGenerate`: Create a new genome reference from splice junctions inferred from first-pass STAR alignment
- `runSTAR2pass`: Perform a second-pass alignment to the genome reference produced by `runSTARgenomeGenerate`
#### BAM QC, performed using [RSeQC](http://rseqc.sourceforge.net)
- `runRSeQCbamStat`
- `runRSeQCclippingProfile`
- `runRSeQCdeletionProfile`
- `runRSeQCgeneBodyCoverage` *Note: this is a very slow step*
- `runRSeQCinferExperiment`
- `runRSeQCinnerDistance`
- `runRSeQCinsertionProfile`
- `runRSeQCjunctionAnnotation`
- `runRSeQCjunctionSaturation`
- `runRSeQCreadDistribution`
- `runRSeQCreadDuplication`
- `runRSeQCreadGC`
- `runRSeQCreadNVC`
- `runRSeQCreadQuality`
- `runRSeQCtin` *Note: this is a very slow step*
#### Other BAM-level QC
- `runMultiQCSample`: Compile output from RSeQC and STAR using MultiQC
#### Post-alignment processing
- `runRSEMcalculateExpression`: Use [RSEM](https://deweylab.github.io/RSEM/) to estimate gene- and transcript (isoform)-level expression
- `createSE`: Create [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) R objects, one for gene-level expression and one for transcript-level expression; each object contains several estimates of expression from RSEM, as well as QC parameters and feature-level annotation

### `RNA_Seq_template.config`
This configuration file is intended to be used only with this Nextflow script.  It makes a number of assumptions about underlying directory structures and filenames.
The parameters that are typically changed are:

#### `params.infile`
Full path to a TSV file containing the following columns (those in bold are **not optional** and cannot be left blank):
- `INDIVIDUAL_ID`: An ID for an individual from which one or more samples was obtained (if only one sample was sequenced from each individual, this can be the same as `SAMPLE_ID`, or left blank)
- **`SAMPLE_ID`**: An ID for each sample
- `LIBRARY_ID`: An ID for each library prepared from a sample (if only one library was sequenced from each sample, this can be the same as `SAMPLE_ID`, or left blank)
- `RG_ID`: Read Group ID: the flowcell ID, optionally followed by a lane-specific suffix (for instruments with independent lanes), followed by a sample-specific identifier (e.g., `SAMPLE_ID`)
- `PLATFORM_UNIT`: `RG ID`, optionally followed by a suffix specific to a library (if more than one library was sequenced per sample)
- `PLATFORM`: Sequencing platform, e.g., "illumina" for Illumina instruments
- `PLATFORM_MODEL`: Instrument model, e.g., "NextSeq", "HiSeq", etc.)
- `RUN_DATE`: Optional run date
- `CENTER`: Optional name for center at which sequencing was performed
- **`R1`**: Full path to FASTQ file containing first paired-end read
- **`R2`**: Full path to FASTQ file containing second paired-end read

Note: some of these fields are discussed in more detail within the [GATK read groups documentation](https://software.broadinstitute.org/gatk/documentation/article.php?id=6472).

#### `params.output_dir`
Full path where the Nextflow output should be written

#### `params.prefix`
Prefix for Nextflow output files

#### `params.paired_end`
Indicates whether paired-end sequencing was used (`true` or `false`)

#### `params.genome`
Parameters specific to the genome and annotation build used; these can be uncommented and edited as needed

### `submit_RNA_Seq.qsub`
This script kicks off the Nextflow process on the SGE using the .config file specified in its sole argument.
