#!/usr/bin/env nextflow

// Rewritten by Adam Gower, based on scripts by Yusuke Koga

// Global variables for required parameters
inputFile = file(params.infile)
inputFileHeader = params.infile_header

logParams(params, "nextflow_parameters.txt")

// Header log info
log.info ""
log.info "========================================="
log.info "Nextflow Version: $workflow.nextflow.version"
log.info "Command Line:     $workflow.commandLine"
log.info "========================================="
log.info ""

// Send FASTQ files to FastQC and STAR /////////////////////////////////////////

Channel.from(inputFile)
  .splitCsv(sep: '\t', header: inputFileHeader)
  .into {
    readInput_to_runFastQC
    readInput_to_runSTAR1pass
    readInput_to_runSTAR2pass
  }

// First-pass STAR alignment ///////////////////////////////////////////////////

process runSTAR1pass {
  tag "Running STAR first pass on ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/STAR_1Pass"

  module params.modules.star

  input:
  set indivID, sampleID, libraryID, rgID, platform_unit, platform,
    platform_model, run_date, center, R1, R2 \
    from readInput_to_runSTAR1pass

  output:
  file("${star_outFileNamePrefix}SJ.out.tab") \
    into runSTAR1pass_to_runSTARgenomeGenerate

  script:
  star_outFileNamePrefix = "${sampleID}.1pass."

  """
  module list

  STAR \
    --genomeDir ${params.STAR.genomeDir} \
    --readFilesIn ${R1} \
      \$([[ "${params.paired_end}" == "true" ]] && echo "${R2}") \
    --runThreadN \$NSLOTS \
    --outFileNamePrefix ${star_outFileNamePrefix} \
    --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax ${params.STAR.outFilterMultimapNmax} \
    --outFilterType BySJout \
    --readFilesCommand zcat
  """
}

process runSTARgenomeGenerate {
  tag "Generating STAR genome reference with splice junctions"
  publishDir "${params.output_dir}/Output/STAR"

  module params.modules.star

  input:
  // Method toSortedList() used to ensure task will be cached, not resubmitted
  val sjdb_files from runSTAR1pass_to_runSTARgenomeGenerate.toSortedList()

  output:
  file("${genomeDir}/*") into runSTARgenomeGenerateOutput
  file("${genomeDir}") into runSTARgenomeGenerate_to_runSTAR2pass

  script:
  genomeDir="genomeGenerate"

  """
  module list
  mkdir ${genomeDir}/
  STAR \
    --runMode genomeGenerate \
    --genomeDir ${genomeDir} \
    --genomeFastaFiles ${params.ref_fasta} \
    --sjdbFileChrStartEnd ${sjdb_files.join(' ')} \
    --sjdbGTFfile ${params.gene_gtf} \
    --sjdbOverhang ${params.read_length - 1} \
    --runThreadN \$NSLOTS \
    --limitSjdbInsertNsj ${params.STAR.limitSjdbInsertNsj}
  """
}

process runSTAR2pass {
  tag "Running STAR second pass on ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}"

  module params.modules.star
  module params.modules.samtools

  input:
  set indivID, sampleID, libraryID, rgID, platform_unit, platform,
    platform_model, run_date, center, R1, R2 \
    from readInput_to_runSTAR2pass
  file(genomeDir) from runSTARgenomeGenerate_to_runSTAR2pass

  output:
  // Output to MultiQC
  file(star_log_file) into runSTAR2pass_to_runMultiQCSample
  // Output to RSeQC
  set indivID, sampleID, file(outfile_bam) \
    into (
      runSTAR2pass_to_RSeQC_bam_stat,
      runSTAR2pass_to_RSeQC_geneBody_coverage,
      runSTAR2pass_to_RSeQC_junction_annotation,
      runSTAR2pass_to_RSeQC_junction_saturation,
      runSTAR2pass_to_RSeQC_tin,
      runSTAR2pass_to_RSeQC_inner_distance,
      runSTAR2pass_to_RSeQC_clipping_profile,
      runSTAR2pass_to_RSeQC_infer_experiment,
      runSTAR2pass_to_RSeQC_insertion_profile,
      runSTAR2pass_to_RSeQC_deletion_profile,
      runSTAR2pass_to_RSeQC_read_distribution,
      runSTAR2pass_to_RSeQC_QC,
      runSTAR2pass_to_RSeQC_read_duplication,
      runSTAR2pass_to_RSeQC_read_NVC,
      runSTAR2pass_to_RSeQC_read_quality
    )
  // Output to RSEM
  set indivID, sampleID, file(star_transcriptome_bam) \
    into runSTAR2pass_to_runRSEM
  // Files to be stored, but not passed on to other processes
  set file(outfile_bai), file(outfile_bambai) \
    into runSTAR2passOutput

  script:
  star_output_path = "Processing/Libraries/${libraryID}/${rgID}/STAR_2Pass"
  star_outFileNamePrefix = "${star_output_path}/${sampleID}."
  star_coordinate_bam =
    "${star_outFileNamePrefix}Aligned.sortedByCoord.out.bam"
  star_transcriptome_bam = 
    "${star_outFileNamePrefix}Aligned.toTranscriptome.out.bam"
  star_log_file = "${star_outFileNamePrefix}Log.final.out"
  outfile_bam = "${sampleID}.bam"
  outfile_bai = "${sampleID}.bai"
  outfile_bambai = "${sampleID}.bam.bai"

  """
  module list

  # Extract the 'mem_total' qsub parameter, which should be an integer followed
  # by the letter 'G' (e.g., '30G') or an empty string if not specified
  MEM_TOTAL="\$(qstat -j \$JOB_ID | grep -oP "mem_total=[^,]+" | cut -f2 -d'=')"

  # Create @RG (read group) header line for STAR output
  # Note: --outSAMattrRGline flag replaces call to Picard AddOrReplaceReadGroups
  # Read group fields are discussed at:
  # https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
  RGline="ID:${rgID} SM:${sampleID} LB:${libraryID}"
  RGline="\${RGline} PL:${platform} PU:${platform_unit}"

  # Perform second-pass STAR alignment
  mkdir -p ${star_output_path}/
  STAR \
    --genomeDir ${genomeDir} \
    --readFilesIn "${R1}" \
      \$([[ "${params.paired_end}" == "true" ]] && echo "${R2}") \
    --runThreadN \$NSLOTS \
    --outFileNamePrefix ${star_outFileNamePrefix} \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --outFilterMultimapNmax ${params.STAR.outFilterMultimapNmax} \
    --outFilterType BySJout \
    --outSAMunmapped Within \
    --readFilesCommand zcat \
    --outSAMattrRGline \${RGline} \
    \$([[ \$MEM_TOTAL != "" ]] && \
      echo "--limitBAMsortRAM \$((\${MEM_TOTAL/G/} * 1024**3))")

  # Rename coordinate-sorted BAM file
  mv -v "${star_coordinate_bam}" "${outfile_bam}"
  # Index coordinate-sorted BAM file
  samtools index "${outfile_bam}"
  # Copy the .bam.bai index to .bai index for convenience
  cp -av ${outfile_bambai} ${outfile_bai}
  """
}

process runRSEM {
  tag "Running RSEM on ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSEM"

  module params.modules.R
  module params.modules.rsem

  input:
  set indivID, sampleID, star_transcriptome_bam \
    from runSTAR2pass_to_runRSEM

  output:
  // RSEM output files used to create SummarizedExperiment objects
  file("${sampleID}.genes.results") into runRSEM_genes_to_createSE
  file("${sampleID}.isoforms.results") into runRSEM_isoforms_to_createSE
  // Files to be stored, but not passed on to other processes
  file(outfile_plot) into runRSEMOutput

  script:
  outfile_plot = "${sampleID}_RSEM.pdf"

  """
  module list

  # Extract the 'mem_total' qsub parameter, which should be an integer followed
  # by the letter 'G' (e.g., '30G') or an empty string if not specified
  MEM_TOTAL="\$(qstat -j \$JOB_ID | grep -oP "mem_total=[^,]+" | cut -f2 -d'=')"

  # Compute expression values using RSEM
  rsem-calculate-expression \
    --calc-ci --estimate-rspd --no-bam-output --bam \
    \$([[ "${params.paired_end}" == "true" ]] && echo "--paired-end") \
    --forward-prob \
      \$([[ "${params.stranded}" == "true" ]] && echo 0 || echo 0.5) \
    -p \$NSLOTS \
    ${star_transcriptome_bam} \
    ${params.RSEM.reference} \
    ${sampleID}

  rsem-plot-model ${sampleID} ${outfile_plot}
  """
}

// Run FastQC //////////////////////////////////////////////////////////////////

process runFastQC {
  tag "Running FastQC for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastQC"

  module params.modules.fastqc

  input:
  set indivID, sampleID, libraryID, rgID, platform_unit, platform,
    platform_model, run_date, center, R1, R2 \
    from readInput_to_runFastQC

  output:
  file("*.zip") into runFastQC_to_runMultiQCFastq

  script:
  """
  module list

  # Run FastQC
  fastqc -t 1 -o . \
    ${R1} \
    \$([[ "${params.paired_end}" == "true" ]] && echo "${R2}")
  """
}

// Run RSeQC ///////////////////////////////////////////////////////////////////
// Each Python script is run separately, because they do not need to be run
// seqcombining them into a single
// process is a bottleneck (assuming that the concurrent job limit is high
// enough that more than one script can run at a time): the scripts do not need
// to be run sequentially.

process runRSeQC_bam_stat {
  tag "Running RSeQC bam_stat for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_bam_stat
  output:
  file("${sampleID}.bam_stat.txt") into runRSeQC_bam_stat_to_runMultiQCSample
  script:
  """
  module list
  bam_stat.py -i ${bam} > ${sampleID}.bam_stat.txt
  """
}
process runRSeQC_clipping_profile {
  tag "Running RSeQC clipping_profile for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_clipping_profile
  output:
  file("${sampleID}.*") into runRSeQC_clipping_profile_Output
  script:
  """
  module list
  clipping_profile.py -i ${bam} -o ${sampleID} \
    -s \$([[ "${params.paired_end}" == "true" ]] && echo "PE" || echo "SE")
  """
}
process runRSeQC_deletion_profile {
  tag "Running RSeQC deletion_profile for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_deletion_profile
  output:
  file("${sampleID}.*") into runRSeQC_deletion_profile_Output
  script:
  """
  module list
  deletion_profile.py -i ${bam} -l ${params.read_length} -o ${sampleID}
  """
}
process runRSeQC_geneBody_coverage {
  tag "Running RSeQC geneBody_coverage for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_geneBody_coverage
  output:
  file("${sampleID}.geneBodyCoverage.txt") \
    into runRSeQC_geneBody_coverage_to_runMultiQCSample
  file("${sampleID}.*") into runRSeQC_geneBody_coverage_Output
  script:
  """
  module list
  geneBody_coverage.py -i ${bam} -r ${params.gene_bed} -o ${sampleID}
  """
}
process runRSeQC_infer_experiment {
  tag "Running RSeQC infer_experiment for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_infer_experiment
  output:
  file("${sampleID}.infer_experiment.txt") \
    into runRSeQC_infer_experiment_to_runMultiQCSample
  script:
  """
  module list
  infer_experiment.py -i ${bam} -r ${params.gene_bed} > \
    ${sampleID}.infer_experiment.txt
  """
}
process runRSeQC_inner_distance {
  tag "Running RSeQC inner_distance for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_inner_distance
  output:
  file("${sampleID}.inner_distance_freq.txt") \
    into runRSeQC_inner_distance_to_runMultiQCSample
  file("${sampleID}.*") into runRSeQC_inner_distance_Output
  script:
  """
  module list
  inner_distance.py -i ${bam} -r ${params.gene_bed} -o ${sampleID}
  """
}
process runRSeQC_insertion_profile {
  tag "Running RSeQC insertion_profile for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_insertion_profile
  output:
  file("${sampleID}.*") into runRSeQC_insertion_profile_Output
  script:
  """
  module list
  insertion_profile.py -i ${bam} -o ${sampleID} \
    -s \$([[ "${params.paired_end}" == "true" ]] && echo "PE" || echo "SE")
  """
}
process runRSeQC_junction_annotation {
  tag "Running RSeQC junction_annotation for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_junction_annotation
  output:
  file("${sampleID}.junction_annotation.txt") \
    into runRSeQC_junction_annotation_to_runMultiQCSample
  file("${sampleID}.*") into runRSeQC_junction_annotation_Output
  script:
  """
  module list
  junction_annotation.py -i ${bam} -r ${params.gene_bed} -o ${sampleID} 2> \
    ${sampleID}.junction_annotation.txt
  """
}
process runRSeQC_junction_saturation {
  tag "Running RSeQC junction_saturation for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_junction_saturation
  output:
  file("${sampleID}.junctionSaturation_plot.r") \
    into runRSeQC_junction_saturation_to_runMultiQCSample
  file("${sampleID}.*") into runRSeQC_junction_saturation_Output
  script:
  """
  module list
  junction_saturation.py -i ${bam} -r ${params.gene_bed} -o ${sampleID}
  """
}
process runRSeQC_read_distribution {
  tag "Running RSeQC read_distribution for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_read_distribution
  output:
  file("${sampleID}.read_distribution.txt") \
    into runRSeQC_read_distribution_to_runMultiQCSample
  script:
  """
  module list
  read_distribution.py -i ${bam} -r ${params.gene_bed} > \
    ${sampleID}.read_distribution.txt
  """
}
process runRSeQC_read_duplication {
  tag "Running RSeQC read_duplication for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_read_duplication
  output:
  file("${sampleID}.pos.DupRate.xls") \
    into runRSeQC_read_duplication_to_runMultiQCSample
  file("${sampleID}.*") into runRSeQC_read_duplication_Output
  script:
  """
  module list
  read_duplication.py -i ${bam} -o ${sampleID}
  """
}
process runRSeQC_read_GC {
  tag "Running RSeQC read_GC for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_QC
  output:
  file("${sampleID}.GC.xls") into runRSeQC_read_GC_to_runMultiQCSample
  file("${sampleID}.*") into runRSeQC_read_GC_Output
  script:
  """
  module list
  read_GC.py -i ${bam} -o ${sampleID}
  """
}
process runRSeQC_read_NVC {
  tag "Running RSeQC read_NVC for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_read_NVC
  output:
  file("${sampleID}.*") into runRSeQC_read_NVC_Output
  script:
  """
  module list
  read_NVC.py -i ${bam} -o ${sampleID}
  """
}
process runRSeQC_read_quality {
  tag "Running RSeQC read_quality for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_read_quality
  output:
  file("${sampleID}.*") into runRSeQC_read_quality_Output
  script:
  """
  module list
  read_quality.py -i ${bam} -o ${sampleID}
  """
}
process runRSeQC_tin {
  tag "Running RSeQC tin for ${sampleID}"
  publishDir "${params.output_dir}/${indivID}/${sampleID}/RSeQC"
  module "${params.modules.python3}:${params.modules.rseqc}"
  input:
  set indivID, sampleID, bam from runSTAR2pass_to_RSeQC_tin
  output:
  file("${sampleID}.summary.txt") into rseqc_tin_to_createSE
  file("${file(bam).baseName}.*") into rseqc_tin_Output
  script:
  """
  module list
  tin.py -i ${bam} -r ${params.gene_bed} > ${sampleID}.summary.txt
  """
}

process runMultiQCFastq {
  tag "Generating FASTQ-level summary and QC plots"
  publishDir "${params.output_dir}/Output/QC/Fastq"

  module params.modules.python2
  module params.modules.multiqc

  input:
  // Method toSortedList() used to ensure task will be cached, not resubmitted
  val fastqc_files from runFastQC_to_runMultiQCFastq.flatten().toSortedList()

  output:
  set file("fastq_multiqc.html"), file("fastq_multiqc_data/*") \
    into runMultiQCFastqOutput
  file("fastq_multiqc_data/multiqc_fastqc.txt") \
    into runMultiQCFastq_to_createSE

  script:
  """
  module list
  # Create temporary MultiQC input file
  multiqc_input_tempfile=\$(mktemp)
  echo -e "${fastqc_files.join('\n')}" > \$multiqc_input_tempfile
  # Run MultiQC
  multiqc -n fastq_multiqc --file-list \$multiqc_input_tempfile
  """
}

process runMultiQCSample {
  tag "Generating sample-level summary and QC plots"
  publishDir "${params.output_dir}/Output/QC/Sample"

  module params.modules.python2
  module params.modules.multiqc

  input:
  // Method toSortedList() used to ensure task will be cached, not resubmitted
  val rseqc_bam_stat_files \
    from runRSeQC_bam_stat_to_runMultiQCSample.toSortedList()
  val rseqc_geneBody_coverage_files \
    from runRSeQC_geneBody_coverage_to_runMultiQCSample.toSortedList()
  val rseqc_infer_experiment_files \
    from runRSeQC_infer_experiment_to_runMultiQCSample.toSortedList()
  val rseqc_inner_distance_files \
    from runRSeQC_inner_distance_to_runMultiQCSample.toSortedList()
  val rseqc_junction_annotation_files \
    from runRSeQC_junction_annotation_to_runMultiQCSample.toSortedList()
  val rseqc_junction_saturation_files \
    from runRSeQC_junction_saturation_to_runMultiQCSample.toSortedList()
  val rseqc_read_distribution_files \
    from runRSeQC_read_distribution_to_runMultiQCSample.toSortedList()
  val rseqc_read_duplication_files \
    from runRSeQC_read_duplication_to_runMultiQCSample.toSortedList()
  val rseqc_read_GC_files \
    from runRSeQC_read_GC_to_runMultiQCSample.toSortedList()
  val star_files \
    from runSTAR2pass_to_runMultiQCSample.toSortedList()

  output:
  set file("sample_multiqc.html"), file("sample_multiqc_data/*") \
    into runMultiQCSampleOutput
  file("sample_multiqc_data/multiqc_rseqc_*.txt") \
    into runMultiQCSample_rseqc_to_createSE
  file("sample_multiqc_data/multiqc_star.txt") \
    into runMultiQCSample_star_to_createSE

  script:
  """
  module list

  # Create temporary MultiQC input file
  multiqc_input_tempfile=\$(mktemp)
  for filename in \
    '${rseqc_bam_stat_files.join("' '")}' \
    '${rseqc_geneBody_coverage_files.join("' '")}' \
    '${rseqc_infer_experiment_files.join("' '")}' \
    '${rseqc_inner_distance_files.join("' '")}' \
    '${rseqc_junction_annotation_files.join("' '")}' \
    '${rseqc_junction_saturation_files.join("' '")}' \
    '${rseqc_read_distribution_files.join("' '")}' \
    '${rseqc_read_duplication_files.join("' '")}' \
    '${rseqc_read_GC_files.join("' '")}' \
    '${star_files.join("' '")}'
  do
    echo "\$filename" >> \$multiqc_input_tempfile
  done
  # Run MultiQC
  multiqc -n sample_multiqc --file-list \$multiqc_input_tempfile
  """
}

// Combine results into SummarizedExperiment objects ///////////////////////////

process createSE {
  tag "Combining results into SummarizedExperiment objects"
  publishDir "${params.output_dir}/Output/Expression"

  module params.modules.R

  input:
  // Method toSortedList() used to ensure task will be cached, not resubmitted
  val multiqcfastqc_file from runMultiQCFastq_to_createSE
  val rseqc_tin_files from rseqc_tin_to_createSE.toSortedList()
  val multiqc_rseqc_files \
    from runMultiQCSample_rseqc_to_createSE.flatten().toSortedList()
  val multiqc_star_file from runMultiQCSample_star_to_createSE
  val rsem_genes_files from runRSEM_genes_to_createSE.toSortedList()
  val rsem_isoforms_files from runRSEM_isoforms_to_createSE.toSortedList()

  output:
  file("*.rds") into runCreateSEOutput

  script:
  rds_files = [
    gene: "${params.prefix}_Gene_Expression.rds",
    isoform: "${params.prefix}_Isoform_Expression.rds"
  ]
  """
  #!/usr/bin/env Rscript

  # Load and list packages
  library(biomaRt)
  library(GenomicFeatures)
  library(SummarizedExperiment)
  sessionInfo()

  # Convenience function to read whitespace-separated values file with a header,
  # and leaving character values and column names as-is
  read.wsv <- function (file, ...) {
    read.table(
      file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE,
      ...
    )
  }

  output <- list()

  # Parse input file ###########################################################

  output[["inputfile"]] <- read.delim(
    "${params.infile}", stringsAsFactors=FALSE
  )
  n.samples <- nrow(output[["inputfile"]])
  # Aggregate by individual ID, sample ID, library ID
  output[["inputfile"]] <- aggregate(
    output[["inputfile"]],
    by=output[["inputfile"]][c("INDIVIDUAL_ID","SAMPLE_ID","LIBRARY_ID")],
    FUN=paste, collapse=","
  )
  # Remove aggregation columns and copy SAMPLE_ID column into row names
  output[["inputfile"]] <- output[["inputfile"]][-(1:3)]
  rownames(output[["inputfile"]]) <- output[["inputfile"]][["SAMPLE_ID"]]

  # Parse MultiQC FastQC output table ##########################################

  # FILE FORMAT
  # 'Filename' column: (base) FASTQ filename.
  # 'Sample' column: (base) FASTQ filename, without ".fastq.gz" extension.
  # The remaining columns have the following types:
  # 'pass', 'fail', or 'warn':
  #   adapter_content, sequence_duplication_levels, per_base_sequence_quality,
  #   sequence_length_distribution, basic_statistics, per_sequence_gc_content,
  #   per_base_n_content, per_base_sequence_content, overrepresented_sequences,
  #   per_tile_sequence_quality, per_sequence_quality_scores
  # Percentages [0..100]:
  #   %GC, total_deduplicated_percentage
  # Floating point values:
  #   avg_sequence_length
  # Integer values (represented as floating-point with trailing .0):
  #   Sequences flagged as poor quality, Total Sequences
  # Range of integers (mm-nn):
  #   Sequence length
  # Text:
  #   Encoding, File type

  # Read file, placing "Sample" column in row names
  output[["multiqcfastqc"]] <- read.wsv(
    "${multiqcfastqc_file}", sep="\\t", row.names=1
  )

  # Initialize a data frame with FASTQ filenames in the row names
  aggregate.by <- data.frame(row.names=output[["multiqcfastqc"]][["Filename"]])
  # Extract the read ID (R1 or R2) from FASTQ filenames, and use filenames
  # and read IDs to look up the sample IDs from input table
  aggregate.by[["sample"]] <- NA_character_
  aggregate.by[["read"]] <- sub(
    ".+_(R[12])_.+", "\\\\1", rownames(aggregate.by)
  )

  reads <- if ("${params.paired_end}" == "true") c("R1", "R2") else c("R1")
  for (read in reads) {
    i <- which(aggregate.by[["read"]] == read)
    aggregate.by[["sample"]][i] <- output[["inputfile"]][["SAMPLE_ID"]][
      match(rownames(aggregate.by)[i], basename(output[["inputfile"]][[read]]))
    ]
  }
  # Collapse the MultiQC results first by sample ID, then by read
  output[["multiqcfastqc"]] <- aggregate(
    output[["multiqcfastqc"]], by=aggregate.by, FUN=c, simplify=FALSE
  )

  # Simplify columns with multiple sets of FASTQ files
  # using an appropriate function
  columns.to.simplify <- list(
    sum = c("Sequences flagged as poor quality"),
    mean = c("avg_sequence_length", "%GC", "total_deduplicated_percentage")
  )
  for (f in names(columns.to.simplify)) {
    for (j in columns.to.simplify[[f]]) {
      output[["multiqcfastqc"]][[j]] <- sapply(
        output[["multiqcfastqc"]][[j]], f
      )
    }
  }
  # Aggregate data frame again by sample ID
  # (so R1 and R2 metrics are on the same row)
  output[["multiqcfastqc"]] <- aggregate(
    output[["multiqcfastqc"]][-(1:2)],
    by=output[["multiqcfastqc"]]["sample"], FUN=c, simplify=FALSE
  )
  # Move the sample IDs into the row names of the data frame
  rownames(output[["multiqcfastqc"]]) <- output[["multiqcfastqc"]][["sample"]]
  output[["multiqcfastqc"]][["sample"]] <- NULL
  # Add a prefix to the column names
  colnames(output[["multiqcfastqc"]]) <- paste(
    "FastQC", colnames(output[["multiqcfastqc"]])
  )

  # Parse tables of TIN metrics from RSeQC #####################################

  # FILE FORMAT
  # 'Bam_file' column: (base) BAM filename.
  # The remaining columns have the following types:
  # Floating point values:
  #   TIN(mean), TIN(median), TIN(stdev)

  tinfiles <- c('${rseqc_tin_files.join("','")}')
  output[["tinfiles"]] <- c(
    # Get tab-delimited column names from first line of first file
    readLines(tinfiles[1], n=1),
    # Read tab-delimited contents of second lines of all files
    sapply(tinfiles, scan, what="", n=1, sep="\\n", skip=1, quiet=TRUE)
  )
  # Parse the character vector into a data frame
  output[["tinfiles"]] <- read.wsv(
    textConnection(output[["tinfiles"]]), row.names=1
  )
  # Add a prefix to the column names
  colnames(output[["tinfiles"]]) <- paste(
    "RSeQC", colnames(output[["tinfiles"]])
  )

  # Parse remaining MultiQC output files #######################################

  # Read tables, placing first column ("Sample") in row names
  # and adding prefix to column names
  multiqc.rseqc.filenames <- c('${multiqc_rseqc_files.join("','")}')
  for (filename in multiqc.rseqc.filenames) {
    output[[basename(filename)]] <- read.wsv(filename, row.names=1)
    colnames(output[[basename(filename)]]) <- paste(
      "RSeQC", colnames(output[[basename(filename)]])
    )
  }
  output[["starfile"]] <- read.wsv("${multiqc_star_file}", row.names=1)
  colnames(output[["starfile"]]) <- paste(
    "STAR", colnames(output[["starfile"]])
  )

  # Reformat each data frame ###################################################

  for (i in which(names(output) != "inputfile")) {
    # Reorder rows of each data frame to match order of rows in input file
    # Note: a while loop is used because pmatch() returns NA where there are
    #       multiple matches, e.g.,
    #         pmatch("Sample_1", c("Sample_1.bam", "Sample_10.bam")) == NA
    #       Each iteration excludes any values that were previously matched.
    j <- rep(NA_integer_, n.samples)
    k <- seq(n.samples)
    while (any(is.na(j))) {
      j[is.na(j)] <- pmatch(
        rownames(output[["inputfile"]])[is.na(j)], rownames(output[[i]])[k]
      )
      k[j[!is.na(j)]] <- NA_integer_
    }
    output[[i]] <- output[[i]][j, , drop=FALSE]
    # 1. Replace "%" with "percent", e.g., "%GC" -> "percent GC"
    # 2. Replace whitespace/punctuation in column names with underscores
    # 3. Remove trailing underscores from column names
    colnames(output[[i]]) <- gsub("% ?", "percent ", colnames(output[[i]]))
    colnames(output[[i]]) <- gsub(
      "_\$", "", gsub("[' \\"\\\\(\\\\)\\\\-]", "_", colnames(output[[i]]))
    )
  }

  # Combine output into a single data frame ####################################

  # Initialize with sample names from input file, then add each table in turn
  SE.colData <- data.frame(row.names=rownames(output[["inputfile"]]))
  for (i in seq_along(output)) {
    SE.colData <- cbind(SE.colData, output[[i]])
  }

  # Open connection to Biomart #################################################

  # Determine from parameters which Biomart dataset and hostname should be used
  biomart <- list(
    database = "ENSEMBL_MART_ENSEMBL",
    dataset = paste(
      tolower(
        sub("^(.)[^ ]+ (.+)\$", "\\\\1\\\\2", "${params.genome.species}")
      ),
      "gene_ensembl", sep="_"
    ),
    host = with(
      listEnsemblArchives(), url[version == "${params.genome.ensembl}"]
    )
  )
  # Note: this connection is done in a repeat loop because it doesn't always
  #       work on the first try; a maximum number of attempts is included to
  #       keep it from endlessly looping if Biomart is down
  cat(
    "Opening connection to Biomart",
    with(
      biomart,
      sprintf(
        "(database %s, dataset %s, host %s, Ensembl version %s).\\n",
        database, dataset, host, "${params.genome.ensembl}"
      )
    )
  )
  max.attempts <- 10
  attempt <- 1
  repeat {
    cat(sprintf("Connection attempt #%d of %d...", attempt, max.attempts))
    mart <- try(
      with(biomart, useMart(biomart=database, dataset=dataset, host=host)),
      silent=TRUE
    )
    if (inherits(mart, "Mart")) {
      cat("successful.\\n")
      break
    }
    cat("\\n")
    if (attempt == max.attempts) {
      stop("Could not connect to Biomart after ", max.attempts, " attempts")
    }
    attempt <- attempt + 1
  }

  # Create SummarizedExperiment objects ########################################

  # Copy Nextflow variables into R variables
  rsem.filenames <- list(
    gene=c('${rsem_genes_files.join("','")}'),
    isoform=c('${rsem_isoforms_files.join("','")}')
  )
  biomart.attributes <- list(
    gene = c('${params.createSE.biomart_attributes.gene.join("','")}'),
    isoform = c('${params.createSE.biomart_attributes.isoform.join("','")}')
  )
  rds.files <- c(gene='${rds_files.gene}', isoform='${rds_files.isoform}')

  # RSEM output columns related to gene/isoform-level annotation
  rsem.annotation.columns <- c(
    "gene_id", "transcript_id(s)", "length", "effective_length"
  )
  # Names of RSEM output columns holding gene/isoform IDs
  rsem.id.columns <- c(gene="gene_id", isoform="transcript_id")

  # Create a TxDb object from the GTF file; this will be used to populate the
  # rowRanges of the SummarizedExperiment objects
  txdb <- makeTxDbFromGFF("${params.gene_gtf}")

  for (type in c("gene", "isoform")) {
    # Reorder RSEM output files to match order of sample annotation
    filename.suffix <- paste0("\\\\.", type, "s", "\\\\.results\$")
    rsem.filenames[[type]] <- rsem.filenames[[type]][
      match(
        rownames(SE.colData),
        sub(filename.suffix, "", basename(rsem.filenames[[type]]))
      )
    ]
    # Read the RSEM output files into a list of matrices
    SE.assays <- list()
    n <- length(rsem.filenames[[type]])
    for (i in seq(n)) {
      cat("Processing", type, "level RSEM output file", i, "of", n, "\\r")
      rsem.output <- read.wsv(
        rsem.filenames[[type]][i], row.names=rsem.id.columns[type]
      )
      # Copy each RSEM column (minus any annotation columns) into assay list
      for (j in setdiff(colnames(rsem.output), rsem.annotation.columns)) {
        SE.assays[[j]] <- cbind(SE.assays[[j]], rsem.output[[j]])
      }
    }
    cat("\\n")
    # Name the rows and columns of each matrix
    SE.assays <- lapply(SE.assays, `rownames<-`, rownames(rsem.output))
    SE.assays <- lapply(SE.assays, `colnames<-`, rownames(SE.colData))

    # Retrieve annotation from Biomart
    BM <- getBM(attributes=biomart.attributes[[type]], mart=mart)
    # Collapse by first column (either gene or transcript ID)
    BM <- aggregate(BM[-1], BM[1], unique)
    # Reorder to match order of features from RSEM output
    BM <- BM[match(rownames(SE.assays[[1]]), BM[[1]]), ]
    # Combine Biomart annotation with data from TxDb object
    SE.rowRanges <- switch(type, gene=genes(txdb), isoform=transcripts(txdb))
    mcols(SE.rowRanges) <- c(mcols(SE.rowRanges), BM[-1])

    # Assemble variables into a SummarizedExperiment and write to RDS file
    dataset <- SummarizedExperiment(
      assays=SE.assays, colData=DataFrame(SE.colData), rowRanges=SE.rowRanges
    )
    saveRDS(dataset, rds.files[[type]])
  }
  """
}

workflow.onComplete {
  log.info ""
  log.info "========================================="
  log.info "Duration:       $workflow.duration"
  log.info "========================================="
}

// FUNCTIONS ///////////////////////////////////////////////////////////////////

// Read input file and save it into list of lists //////////////////////////////
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for (s in p) {
    file << "${s.key}:\t${s.value}\n"
  }
}
