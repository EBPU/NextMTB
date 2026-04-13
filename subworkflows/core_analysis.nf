// Include modules
include { MAPPING; MAPPING_ONT; REFINE; REFINE_ONT; PILE; PILE_ONT; LIST; VARIANTS_LOW; VARIANTS; STATS; STRAIN; MAP_STRAIN; JOIN } from '../modules/all_modules.nf'

workflow CORE_ANALYSIS {
	take:
		ch_reads
		ch_hist_bams
		ch_hist_ptables
		ch_hist_var_std
		ch_hist_var_low
		seq_type
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

		ch_hist_ptable_ids  = ch_hist_ptables.map  { it[0] }.toList()
		ch_hist_var_std_ids = ch_hist_var_std.map  { it[0] }.toList()
		ch_hist_var_low_ids = ch_hist_var_low.map  { it[0] }.toList()


		// Branch input based on sequencing technology
		ch_reads.branch {
			illumina: seq_type == 'ILL'
			nanopore: seq_type == 'ONT'
		}.set { processing_branch }

		// --- Illumina Pipeline ---
		MAPPING(processing_branch.illumina, ref)
		MAPPING_ONT(processing_branch.nanopore, ref)

		ch_all_bams = MAPPING.out.bam
			.mix(MAPPING_ONT.out.bam)
			.mix(ch_hist_bams)
			.unique { it[0] }
		
		ch_bams_to_refine = ch_all_bams
			.join(ch_hist_ptables, by: 0, remainder: true) // Keep all BAMs, even those without historical data
			.filter { it.size() == 3 && it[1] != null && it[2] == null }
            .map { tuple(it[0], it[1]) }

		
		ch_bams_to_refine.branch {
            illumina: seq_type == 'ILL'
            nanopore: seq_type == 'ONT'
        }.set { bams_for_gatk }

		REFINE(bams_for_gatk.illumina, ref)
        PILE(REFINE.out.gatk, ref)

        REFINE_ONT(bams_for_gatk.nanopore, ref, ascii)
        PILE_ONT(REFINE_ONT.out.gatk, ref, minbqual)

		// Merge pileup outputs
		ch_mpile = PILE.out.mpile.mix(PILE_ONT.out.mpile)

		// Generate position tables
		LIST(ch_mpile, minbqual, ref)

		ch_all_ptables = LIST.out.list
					.mix(ch_hist_ptables)
					.unique { it[0] }

		ch_ptables_for_var_std = ch_all_ptables
            .join(ch_hist_var_std, by: 0, remainder: true)
            .filter { it.size() == 3 && it[1] != null && it[2] == null }
            .map { tuple(it[0], it[1]) }


		VARIANTS(ch_ptables_for_var_std, mincovf, mincovr, minphred20, ref)

        // ONLY call low freq variants if they DO NOT exist historically
        ch_ptables_for_var_low = ch_all_ptables
            .join(ch_hist_var_low, by: 0, remainder: true)
            .filter { it.size() == 3 && it[1] != null && it[2] == null }
            .map { tuple(it[0], it[1]) }

        VARIANTS_LOW(ch_ptables_for_var_low, ref)

        // Recombine all variants
        ch_all_var_std = VARIANTS.out.var.mix(ch_hist_var_std).unique { it[0] }
        ch_all_var_low = VARIANTS_LOW.out.var_low.mix(ch_hist_var_low).unique { it[0] }

        // --- Statistics and Strain Classification ---
        // Stats uses all BAMs and all Position Tables
        STATS(ch_all_bams.join(ch_all_ptables, by: 0), mincovf, mincovr, minphred20)
        STRAIN(ch_all_ptables)

        ch_map_strain_input = STATS.out.stats.join(STRAIN.out.strain, by: 0)
            .map { id, file1, file2 -> tuple(file1, file2) }
            .collect()
            
        MAP_STRAIN(ch_map_strain_input)

        // --- Joint Analysis ---
        if (run_join) {

			ch_intersected = ch_all_var_std.join(ch_all_ptables, by: 0)

			// Now you can split them back safely, knowing they perfectly match
			ch_call_mixed  = ch_intersected.map { id, var, ptable -> var }.collect()
			ch_list_mixed = ch_intersected.map { id, var, ptable -> ptable }.collect()
						
            //ch_call_mixed = ch_all_var_std.map { id, file -> file }.collect()
            //ch_list_mixed = ch_all_ptables.map { id, file -> file }.collect()


                
            JOIN(ch_call_mixed, ch_list_mixed, Channel.fromPath(sj, checkIfExists: true).collect(), minbqual, minphred20, proj, ref)
        }

    emit:
        bam = ch_all_bams
        var_low = ch_all_var_low
        map_strain = MAP_STRAIN.out

		new_bams = MAPPING.out.bam
        new_ptables = LIST.out.list
        new_var_std = VARIANTS.out.var
        new_var_low = VARIANTS_LOW.out.var_low
}