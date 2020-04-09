# RNA_Seq

### `RNA_Seq.nf`

This is a Nextflow pipeline that contains the following processes:
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

This is a template configuration file intended to be used only with this Nextflow script.  It makes a number of assumptions about underlying directory structures and filenames.
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
