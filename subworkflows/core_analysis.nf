// Include modules
include { MAPPING; MAPPING_ONT; REFINE; REFINE_ONT; PILE; PILE_ONT; LIST; VARIANTS_LOW; VARIANTS; STATS; STRAIN; MAP_STRAIN; JOIN } from '../modules/all_modules.nf'

workflow CORE_ANALYSIS {
    take:
        ch_reads
        SEQ
        ref
        ascii
        minbqual
        mincovf
        mincovr
        minphred20
        join
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
        REFINE(MAPPING.out.bam, ref)
        PILE(REFINE.out.gatk, ref)

        // --- ONT Pipeline ---
        MAPPING_ONT(processing_branch.nanopore, ref)
        REFINE_ONT(MAPPING_ONT.out.bam, ref, ascii)
        PILE_ONT(REFINE_ONT.out.gatk, ref, minbqual)

        // Merge mapping outputs handling both new bams and previously existing ones in the directory
        ch_new_mapped = MAPPING.out.bam.mix(MAPPING_ONT.out.bam)
        ch_old_mapped = Channel.fromPath('Bam/*bam*').map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }.groupTuple()
        ch_mapped_bam = ch_new_mapped.mix(ch_old_mapped).unique { it[0] }

        // Merge pileup outputs
        ch_mpile = PILE.out.mpile.mix(PILE_ONT.out.mpile)

        // Generate position tables
        LIST(ch_mpile, minbqual, ref)
        
        // Handle existing position tables
        ch_old_list = Channel.fromPath('Position_Tables/*table.tab').map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }.groupTuple()
        ch_ptables = LIST.out.list.mix(ch_old_list).unique { it[0] }

        // Call variants
        VARIANTS_LOW(LIST.out.list, ref)
        VARIANTS(LIST.out.list, mincovf, mincovr, minphred20, ref)

        // Generate statistics and classify strain
        STATS(ch_mapped_bam.join(ch_ptables, by: 0), mincovf, mincovr, minphred20)
        STRAIN(ch_ptables)

        // Map strain and statistics together
        ch_map_strain_input = STATS.out.stats.join(STRAIN.out.strain, by: 0).map { id, file1, file2 -> tuple(file1, file2) }.collect()
        MAP_STRAIN(ch_map_strain_input)

        // Perform joint analysis if requested
        if (join) {
            ch_call = VARIANTS.out.var.mix(Channel.fromPath('Called/*variants_cf4*').map { file -> tuple((file.getSimpleName() - ~/_.*/), file) })
                .unique { it[0] }.map { id, file -> file }.collect()
            
            ch_list = LIST.out.list.mix(Channel.fromPath('Position_Tables/*').map { file -> tuple((file.getSimpleName() - ~/_.*/), file) })
                .unique { it[0] }.map { id, file -> file }.collect()
                
            JOIN(ch_call, ch_list, Channel.fromPath(sj, checkIfExists: true).collect(), minbqual, minphred20, proj, ref)
        }

    emit:
        // Export channels needed for downstream analysis
        bam = ch_mapped_bam
        var_low = VARIANTS_LOW.out.var_low
        var_standard = VARIANTS.out.var
        map_strain = MAP_STRAIN.out
}