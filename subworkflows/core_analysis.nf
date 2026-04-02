// Include modules
include { MAPPING; MAPPING_ONT; REFINE; REFINE_ONT; PILE; PILE_ONT; LIST; VARIANTS_LOW; VARIANTS; STATS; STRAIN; MAP_STRAIN; JOIN } from '../modules/all_modules.nf'

workflow CORE_ANALYSIS {
    take:
        ch_reads
		ch_historical_bams
		ch_historical_ptables
		ch_historical_var_std
        SEQ
        ref
        ascii
        minbqual
        mincovf
        mincovr
        minphred20
        run_join
        sj
        proj

    main:
        // Branch input based on sequencing technology
        ch_reads.branch {
            illumina: SEQ == 'ILL'
            nanopore: SEQ == 'ONT'
        }.set { processing_branch }

        // --- Illumina Pipeline ---
        MAPPING(processing_branch.illumina, ref)
		MAPPING_ONT(processing_branch.nanopore, ref)

		ch_bams_to_refine = MAPPING.out.bam
            .mix(MAPPING_ONT.out.bam)
            .mix(ch_historical_bams)
            .unique { it[0] }

        REFINE(ch_bams_to_refine, ref)
        PILE(REFINE.out.gatk, ref)

        // --- ONT Pipeline ---
        REFINE_ONT(MAPPING_ONT.out.bam, ref, ascii)
        PILE_ONT(REFINE_ONT.out.gatk, ref, minbqual)

        // Merge pileup outputs
        ch_mpile = PILE.out.mpile.mix(PILE_ONT.out.mpile)

        // Generate position tables
        LIST(ch_mpile, minbqual, ref)

		ch_ptables = LIST.out.list
            .mix(ch_historical_ptables)
            .unique { it[0] }

        // Call variants
        VARIANTS_LOW(ch_ptables, ref)
        VARIANTS(ch_ptables, mincovf, mincovr, minphred20, ref)

		STATS(ch_bams_to_refine.join(ch_ptables, by: 0), mincovf, mincovr, minphred20)
        STRAIN(ch_ptables)

        // Map strain and statistics together
		ch_map_strain_input = STATS.out.stats.join(STRAIN.out.strain, by: 0)
            .map { id, file1, file2 -> tuple(file1, file2) }
            .collect()
        MAP_STRAIN(ch_map_strain_input)

		// --- Joint Analysis ---
        if (run_join) {
            // Mix new standard variants with historical ones for the joint calling
            ch_call_mixed = VARIANTS.out.var
                .mix(ch_historical_var_std)
                .unique { it[0] }
                .map { id, file -> file }
                .collect()
            
            ch_list_mixed = ch_ptables.map { id, file -> file }.collect()
                
            JOIN(ch_call_mixed, ch_list_mixed, Channel.fromPath(sj, checkIfExists: true).collect(), minbqual, minphred20, proj, ref)
        }

    emit:
        // Export channels needed for downstream analysis
        bam = ch_bams_to_refine
        var_low = VARIANTS_LOW.out.var_low
        map_strain = MAP_STRAIN.out
}	