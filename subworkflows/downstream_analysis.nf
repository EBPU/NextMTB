// Include modules
include { DEL; DEL_ONT; OUT_DEL; DEPTH; OUT_DEPTH; MUT_CORRECTION_DEL; MUT_CORRECTION; MUT_GATHER; PHARMA; WHO; OUT_WHO; FINAL_OUT } from '../modules/all_modules.nf'

workflow DOWNSTREAM_ANALYSIS {
    take:
        ch_bam
        ch_var_low
        ch_map_strain
		ch_hist_var_low
        ch_hist_corrected
        seq_type
        ref
        bed
        bedix
        tgene
        run_extra
        pgene
        tdrug
        dhead
        who_cat
        head_who

    main:
        // --- Genome Breadth & Depth Analysis ---

		ch_hist_corrected_ids = ch_hist_corrected.map { it[0] }.toList()

        DEPTH(ch_bam, tgene)
        
        // Collect depth outputs
        ch_depth_collected = DEPTH.out.map { id, file -> file }.collect()
        OUT_DEPTH(ch_depth_collected)

        // Generate final output if reference is H37Rv
        if (ref == "M._tuberculosis_H37Rv_2015-11-13") {
            FINAL_OUT(OUT_DEPTH.out, ch_map_strain)
        }

		ch_var = ch_var_low
            .mix(ch_hist_var_low)
            .unique { it[0] }
		

        // --- Extra Analysis (Deletions) OR Standard Mutation Correction ---
        if (run_extra) {
            
            // Branch bams for deletion calling based on technology
            ch_bam.branch {
                illumina: seq_type == 'ILL'
                nanopore: seq_type == 'ONT'
            }.set { ch_bam_del }

            // Call deletions
            DEL(ch_bam_del.illumina, ref, bed, bedix)
            DEL_ONT(ch_bam_del.nanopore, ref, bed, bedix)
            
            // Mix deletion outputs
            ch_deletion = DEL.out.mix(DEL_ONT.out)
            
            // Summarize deletions
            OUT_DEL(ch_deletion.map { id, file -> file }.collect())

            // Correct mutations including deletions
            ch_var_del = ch_var.join(ch_deletion, by: 0)

			ch_var_del_to_correct = ch_var_del
				.join(ch_hist_corrected, by: 0, remainder: true)
				.filter { it.size() == 4 && it[1] != null && it[3] == null }
				.map { tuple(it[0], it[1], it[2]) }
			

			MUT_CORRECTION_DEL(ch_var_del_to_correct)
            ch_new_corrected = MUT_CORRECTION_DEL.out
            
        } else {
            // ONLY correct mutations if they DO NOT exist historically
            ch_var_to_correct = ch_var_low
				.join(ch_hist_corrected, by: 0, remainder: true)
                .filter { it.size() == 3 && it[1] != null && it[2] == null }
                .map { tuple(it[0], it[1]) }

            MUT_CORRECTION(ch_var_to_correct)
            ch_new_corrected = MUT_CORRECTION.out
        }

        // --- Pharmacoresistance and WHO Catalogue Analysis ---
        // Combine newly corrected mutations with historical ones
        ch_all_corrected = ch_new_corrected
            .mix(ch_hist_corrected)
            .unique { it[0] }
            .map { id, file -> file }
            .collect()

        MUT_GATHER(ch_all_corrected)
        PHARMA(ch_all_corrected, tdrug, pgene)
        WHO(MUT_GATHER.out, dhead, who_cat)
        OUT_WHO(WHO.out, head_who)

	emit:
        new_corrected = ch_new_corrected
}