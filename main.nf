// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Print help message if requested
if (params.h) {
    log.info """
    =========================================
    MTB Pipeline (DSL2 Version)
    =========================================
    --SEQ       Sequencing technology (ILL or ONT) [default: ILL]
    --ref       Reference Genome to use [default: M._tuberculosis_H37Rv_2015-11-13]
    --join      Perform joint analysis [default: true]
    --sj        List of samples to use for joint analysis [optional]
    --proj      Name of the project for joint analysis [default: def]
    --pharma    Perform Drug resistance analysis at Custom frequencies [default: false]
    --tdrug     Custom frequencies % [default: 10]
    --WHO       Path to the WHO catalogue file formatted
    --extra     Perform extra analysis [default: true]
    """
    exit 0
}

// Import sub-workflows
include { PREPROCESS }          from './subworkflows/preprocess.nf'
include { CORE_ANALYSIS }       from './subworkflows/core_analysis.nf'
include { DOWNSTREAM_ANALYSIS } from './subworkflows/downstream_analysis.nf'

// Import standalone modules for the pharma-only run
include { PHARMA; MUT_GATHER; WHO; OUT_WHO } from './modules/all_modules.nf'

workflow {
    
    // --- Standalone Pharma Mode ---
    // Bypasses the rest of the pipeline if --pharma is provided
    if (params.pharma) {
        ch_tabs = Channel.fromPath('Called/*corrected.tab').collect()
        PHARMA(ch_tabs, params.tdrug, params.pgene)
        MUT_GATHER(ch_tabs)
        WHO(MUT_GATHER.out, params.dhead, params.WHO)
        OUT_WHO(WHO.out, params.headWHO)
        
        // Terminate pipeline execution early
        exit 0 
    }

    // --- Standard Pipeline Execution ---
    
    // Log starting parameters for traceability
    log.info """\
    =========================================
    Starting MTB Pipeline...
    Reads       : ${params.reads}
    Reference   : ${params.ref}
    Technology  : ${params.SEQ}
    Results Dir : ${params.results}
    =========================================
    """

    // Create the initial input channel based on the sequencing technology
    if (params.SEQ == "ILL") {
        ch_raw_reads = Channel.fromFilePairs(params.reads + '*_R{1,2}*.fastq.gz')
            .map { id, file -> tuple((id - ~/_.*/), file) }
    } else {
        ch_raw_reads = Channel.fromPath(params.reads + '/*fastq.gz')
            .map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }
    }

    // 1. Pre-processing Sub-workflow
    // Handles reading, standardized naming, and optional Kraken filtering
    PREPROCESS(
        ch_raw_reads, 
        params.SEQ, 
        params.kraken, 
        params.krakendb,
        params.minbqual,
        params.RP,
        params.minphred20
    )

if (params.historical_dir) {
        // Use groupTuple to ensure .bam and .bai for the same ID travel together
        ch_historical_bams = Channel.fromPath("${params.historical_dir}/*bam*")
            .map { file -> tuple((file.name.replaceAll(/\.bam.*$/, '') - ~/_.*/), file) }
            .groupTuple()

        ch_historical_ptables = Channel.fromPath("${params.historical_dir}/*table.tab")
            .map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }

        ch_historical_var_low = Channel.fromPath("${params.historical_dir}/*variants_cf1*001.tab")
            .map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }

        ch_historical_var_std = Channel.fromPath("${params.historical_dir}/*variants_cf4*")
            .map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }

        ch_historical_corrected = Channel.fromPath("${params.historical_dir}/*corrected.tab")
            .map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }
    } else {
        ch_historical_bams      = Channel.empty()
        ch_historical_ptables   = Channel.empty()
        ch_historical_var_low   = Channel.empty()
        ch_historical_var_std   = Channel.empty()
        ch_historical_corrected = Channel.empty()
    }

    // 2. Core Analysis Sub-workflow
    // Takes the clean reads from PREPROCESS and performs mapping, GATK refinement, and variant calling
    CORE_ANALYSIS(
        PREPROCESS.out.ready_reads, 
		ch_historical_bams,      // Inject old BAMs
        ch_historical_ptables,   // Inject old Position Tables
        ch_historical_var_std,   // Inject old Standard Variants
		ch_historical_var_low,   // Inject old Low-freq Variants
        params.SEQ, 
        params.ref,
        params.ascii,
        params.minbqual,
        params.mincovf,
        params.mincovr,
        params.minphred20,
        params.join,
        params.sj,
        params.proj
    )

    // 3. Downstream Analysis Sub-workflow
    // Takes BAMs and variants from CORE_ANALYSIS to perform depth calculation, deletions, and WHO reporting
    DOWNSTREAM_ANALYSIS(
        CORE_ANALYSIS.out.bam,
        CORE_ANALYSIS.out.var_low,
        CORE_ANALYSIS.out.map_strain,
		ch_historical_var_low,   // Inject old Low-freq Variants
        ch_historical_corrected, // Inject old Corrected Mutations
        params.SEQ,
        params.ref,
        params.bed,
        params.bedix,
        params.tgene,
        params.extra,
        params.pgene,
        params.tdrug,
        params.dhead,
        params.who_cat,
        params.head_who
    )

	if (params.historical_dir) {
        
        // Extract only the raw files (dropping the sample ID) and flatten them into a single stream
        ch_files_to_archive = CORE_ANALYSIS.out.new_bams.map { it.drop(1) }.flatten()
            .mix( CORE_ANALYSIS.out.new_ptables.map { it.drop(1) }.flatten() )
            .mix( CORE_ANALYSIS.out.new_var_std.map { it.drop(1) }.flatten() )
            .mix( CORE_ANALYSIS.out.new_var_low.map { it.drop(1) }.flatten() )
            .mix( DOWNSTREAM_ANALYSIS.out.new_corrected.map { it.drop(1) }.flatten() )

        // Send all newly generated files to the historical directory
        ARCHIVE_HISTORICAL(ch_files_to_archive)
		}
}

process ARCHIVE_HISTORICAL {
    tag "Archiving ${file_to_save.name}"
    
    // We remove the "archived_" prefix when publishing to the final folder
    publishDir "${params.historical_dir}", mode: 'link', overwrite: true, saveAs: { filename -> filename.replace('archived_', '') }
    
    input:
    path file_to_save
    
    output:
    path "archived_${file_to_save}"
    
    script:
    """
    ln ${file_to_save} archived_${file_to_save}
    """
}