// Include modules
include { COLLECT_READS; COLLECT_READS_ONT; KRAKEN; BRACKEN; BRACKNOUT; KRAKEN_FILTER; KRAKEN_STATS } from '../modules/all_modules.nf'

workflow PREPROCESS {
    take:
        ch_raw_reads
        SEQ
        kraken
        krakendb
        minbqual
        RP
        minphred20

    main:
        // Branch reads based on sequencing technology
        ch_raw_reads.branch {
            illumina: SEQ == 'ILL'
            nanopore: SEQ == 'ONT'
        }.set { branched_reads }

        // Standardize Illumina reads
        COLLECT_READS(branched_reads.illumina, seq_type, minbqual, RP, minphred20)
        
        // Standardize ONT reads
        COLLECT_READS_ONT(branched_reads.nanopore, seq_type, minbqual, RP, minphred20)

        // Merge standardized reads into a single channel
        ch_collected_reads = COLLECT_READS.out.mix(COLLECT_READS_ONT.out)

        // Run Kraken and Bracken if requested and sequencing is Illumina
        if (kraken && seq_type == 'ILL') {
            KRAKEN(ch_collected_reads, krakendb)
            BRACKEN(KRAKEN.out.kreport, krakendb)
            
            // Collect all bracken outputs to generate the summary
            BRACKNOUT(BRACKEN.out.bout.map { id, file -> file }.collect(sort:true))

            // Join reads with kraken output for filtering
            joined_kraken_ch = ch_collected_reads.join(KRAKEN.out.kraken)
            KRAKEN_FILTER(joined_kraken_ch, seq_type, minbqual, RP, minphred20)
            
            // Set the filtered reads as the final output
            final_reads = KRAKEN_FILTER.out.reads
            
            // Collect stats and generate summary
            kraken_stats = KRAKEN_FILTER.out.stats.map { id, file -> tuple(file) }.collect()
            KRAKEN_STATS(kraken_stats)
        } else {
            // Bypass Kraken entirely
            final_reads = ch_collected_reads
        }

    emit:
        // Expose the ready-to-use reads
        ready_reads = final_reads
}