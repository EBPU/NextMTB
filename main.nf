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
        Channel.fromPath("${params.historical_dir}/**")
            .map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }
            .branch {
                bams:      it[1].name.endsWith('.bam')
                ptables:   it[1].name.endsWith('table.tab')
                var_low:   it[1].name.contains('variants_cf1') && it[1].name.endsWith('001.tab')
                var_std:   it[1].name.contains('variants_cf4')
                corrected: it[1].name.endsWith('corrected.tab')
            }
            .set { ch_historical }
    } else {
        // Create empty channels if no historical directory is provided
        ch_historical = [
            bams:      Channel.empty(),
            ptables:   Channel.empty(),
            var_low:   Channel.empty(),
            var_std:   Channel.empty(),
            corrected: Channel.empty()
        ]
    }

    // 2. Core Analysis Sub-workflow
    // Takes the clean reads from PREPROCESS and performs mapping, GATK refinement, and variant calling
    CORE_ANALYSIS(
        PREPROCESS.out.ready_reads, 
		ch_historical.bams,      // Inject old BAMs
        ch_historical.ptables,   // Inject old Position Tables
        ch_historical.var_std,   // Inject old Standard Variants
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
		ch_historical.var_low,   // Inject old Low-freq Variants
        ch_historical.corrected, // Inject old Corrected Mutations
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
}