// Include modules
include { DEL; DEL_ONT; OUT_DEL; DEPTH; OUT_DEPTH; MUT_CORRECTION_DEL; MUT_CORRECTION; MUT_GATHER; PHARMA; WHO; OUT_WHO; FINAL_OUT } from '../modules/all_modules.nf'

workflow DOWNSTREAM_ANALYSIS {
    take:
        ch_bam
        ch_var_low
        ch_map_strain
        SEQ
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
        DEPTH(ch_bam, tgene)
        
        // Collect depth outputs
        ch_depth_collected = DEPTH.out.map { id, file -> file }.collect()
        OUT_DEPTH(ch_depth_collected)

        // Generate final output if reference is H37Rv
        if (ref == "M._tuberculosis_H37Rv_2015-11-13") {
            FINAL_OUT(OUT_DEPTH.out, ch_map_strain)
        }

        // Incorporate old variants to ensure completeness
        ch_old_var = Channel.fromPath('Called/*variants_cf1*001.tab').map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }.groupTuple()
        ch_var = ch_var_low.mix(ch_old_var).unique { it[0] }

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
            MUT_CORRECTION_DEL(ch_var_del)
            ch_mut = MUT_CORRECTION_DEL.out
            
        } else {
            // Correct mutations without deletions
            MUT_CORRECTION(ch_var)
            ch_mut = MUT_CORRECTION.out
        }

        // Incorporate existing corrected mutations
        ch_old_mut = Channel.fromPath('Called/*corrected.tab').map { file -> tuple((file.getSimpleName() - ~/_.*/), file) }
        ch_mut_gathered = ch_mut.mix(ch_old_mut).unique { it[0] }.map { id, file -> file }.collect()

        // --- Pharmacoresistance and WHO Catalogue Analysis ---
        MUT_GATHER(ch_mut_gathered)
        PHARMA(ch_mut_gathered, tdrug, pgene)
        WHO(MUT_GATHER.out, dhead, who_cat)
        OUT_WHO(WHO.out, head_who)
}