/*
* Collect the reads
*/

process COLLECT_READS {
cpus 1
tag "$replicateId"
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
input:
	tuple val(replicateId), path(reads) 
	val SEQ
	val minbqual
	val r
	val minphred20
output:
	tuple val(replicateId), path('*150bp_R1.fastq.gz'), path('*150bp_R2.fastq.gz')
script:
"""

for i in *_R2*.fastq.gz; do basename=`ls \$i | cut -d "_" -f 1`;mv \$i \${basename}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R2.fastq.gz; done
for i in *_R1*.fastq.gz; do basename=`ls \$i | cut -d "_" -f 1`; mv \$i \${basename}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R1.fastq.gz; done

"""
}

process COLLECT_READS_ONT {
cpus 1
tag "$replicateId"
input:
	tuple val(replicateId), path(reads) 
	val SEQ
	val minbqual
	val r
	val minphred20
output:
	tuple val(replicateId), path('*nbp_ONT.fastq.gz')
script:
"""
for i in *.fastq.gz; do basename=`ls \$i | cut -d "_" -f 1`; mv \$i \${basename}_${SEQ}-Q${minbqual}-RP${r}-PH${minphred20}_nbp_ONT.fastq.gz; done
"""
}




process KRAKEN {
cpus 16
memory '150GB'
tag "$replicateId"
publishDir "kraken"

input:
	tuple val(replicateId), path(reads1), path(reads2)
	path(krakendb)
output:
	tuple val(replicateId), path("*.kreport"), emit: kreport
	tuple val(replicateId), path("*.kraken"), path("*.kreport"), emit: kraken

script:
"""
mkdir kraken
kraken2 --db $krakendb --threads ${task.cpus} --use-names --gzip-compressed --output kraken/${replicateId}.kraken --report kraken/${replicateId}.kreport --paired $reads1 $reads2
mv kraken/*.kreport .
mv kraken/*.kraken .
#touch ${replicateId}.kreport
"""
}


/*
* Kraken Filter
*/

process KRAKEN_FILTER {
conda "/idle/ric.cirillo/dimarco.federico/envs/tools"
cpus 8
tag "$replicateId"
input:
	tuple val(replicateId), path(R1), path(R2)
	tuple val(replicateId), path(kraken), path(kreport)
	val SEQ
	val minbqual
	val r
	val minphred20
output:
	//tuple val(replicateId), path('*150bp_R1.fastq.gz'), path('*150bp_R2.fastq.gz')
	tuple val(replicateId),path("${R1}"), path("${R2}")
script:
"""

#R1=\$(ls ${replicateId}_*R1*.fastq.gz)
#R2=\$(ls ${replicateId}_*R2*.fastq.gz)

FILE1=\$(basename ${R1} .gz)
FILE2=\$(basename ${R2} .gz)

mkdir samp
extract_kraken_reads.py -s1 ${R1} -s2 ${R2} -t 1762 -k "${kraken}" --include-children --include-parents -o "samp/\${FILE1}" -o2 "samp/\${FILE2}" -r "${kreport}" --fastq-output > /dev/null

#gzip "samp/${replicateId}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R1.fastq"
#gzip "samp/${replicateId}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R2.fastq"

#pigz samp/${replicateId}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R1.fastq
#pigz samp/${replicateId}_ILL-Q${minbqual}-RP${r}-PH${minphred20}_150bp_R2.fastq

pigz samp/\${FILE1} &
pigz samp/\${FILE2} &

wait

rm ${R1} ${R2}

mv samp/* .

"""
}




/*
* BRAKEN
*/



process BRACKEN {
cpus 16
tag "$replicateId"
publishDir "bracken", mode:"copy"

input:
    tuple val(replicateId), path(report)
	path(krakendb)
output:
    tuple val(replicateId), path("*.bout"), emit : bout
	tuple val(replicateId),path("*.report"), emit: breport
	val 'done', emit:done
script:
"""
mkdir bracken
bracken -d $krakendb -i $report -o bracken/${replicateId}.bout -w bracken/${replicateId}.report -r 150
mv bracken/* .

#touch ${replicateId}.report
#touch ${replicateId}.bout

"""
}


process BRACKNOUT {
publishDir "OUTPUT", mode:'copy', pattern: 'bracken_summary.csv'
conda '/idle/ric.cirillo/dimarco.federico/envs/prokka'
input:
	path(BRK)

output:
	path("bracken_summary.csv")
script:
"""
for i in *bout; do awk -F '\t' '{print gensub(".bout","","g",FILENAME),\$0}' OFS='\t' \${i} |tr '\t' ';' | sort -k 8 -n -t ';' | tail -n 1 | cut -f1,2,8; done > bracken_summary.csv
"""

}


/*
* Bam
*/

process MAPPING {
//errorStrategy 'ignore'
cpus 8
memory "20GB"
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Bam", mode:'copy', pattern: "*bam*"
input:
	tuple val(replicateId), path(reads1), path(reads2)
        val(ref)
output:
	tuple val(replicateId), path("Bam"), emit: BAM
	tuple val(replicateId), path("*bam*"), emit: bam
	val 'done', emit:done
script:
"""


USER=a perl /opt/conda/bin/MTBseq --step TBbwa --threads ${task.cpus} --ref ${ref}|| echo "processed \$?"
ln -s Bam/*bam* .



"""
}




process MAPPING_ONT {
//errorStrategy 'ignore'
cpus 8
memory "20GB"
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Bam", mode:'copy', pattern: "*bam*"
input:
	tuple val(replicateId), path(reads)
        val(ref)
output:
	tuple val(replicateId), path("Bam"), emit: BAM
	tuple val(replicateId), path("*bam*"), emit: bam
	val 'done', emit:done
script:
"""
ss=\$(ls -1 *fastq.gz | cut -f2 -d '_')
mkdir Bam
bwa mem -t ${task.cpus} /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta *.fastq.gz > Bam/${replicateId}.sam 2>> Bam/${replicateId}.bamlog
samtools view -@ ${task.cpus} -b -T /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta -o Bam/${replicateId}.bam Bam/${replicateId}.sam 2>> Bam/${replicateId}.bamlog
samtools sort -@ ${task.cpus} -T /tmp/${replicateId}.sorted -o Bam/${replicateId}.sorted.bam Bam/${replicateId}.bam 2>> Bam/${replicateId}.bamlog
samtools index -b Bam/${replicateId}.sorted.bam 2>> Bam/${replicateId}.bamlog
samtools rmdup Bam/${replicateId}.sorted.bam Bam/${replicateId}.nodup.bam 2>> Bam/${replicateId}.bamlog
samtools index -b Bam/${replicateId}.nodup.bam 2>> Bam/${replicateId}.bamlog

rm Bam/${replicateId}.sam Bam/${replicateId}.bam Bam/${replicateId}.sorted.bam Bam/${replicateId}.sorted.bam.bai

mv Bam/${replicateId}.nodup.bam Bam/${replicateId}_\${ss}_nBP.bam
mv Bam/${replicateId}.nodup.bam.bai Bam/${replicateId}_\${ss}_nBP.bam.bai
ln -s Bam/*bam* .
"""
}







/*
* GATK
*/


process REFINE {
cpus 8
memory "20GB"
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "GATK_Bam", mode:'copy', pattern: "*gatk*"
input:
        tuple val(replicateId), path(bam1)
        val(ref)
output:
        tuple val(replicateId), path("*gatk*"), emit: gatk
        tuple val(replicateId), path("GATK_Bam"), emit: GATK_Bam
	val 'done', emit:done
script:
"""
mkdir GATK_Bam
mkdir Bam
mv *bam* Bam/
USER=a perl /opt/conda/bin/MTBseq --step TBrefine --threads ${task.cpus} --ref ${ref} || echo "processed \$?"
ln -s GATK_Bam/* .

"""
}


process REFINE_ONT {
cpus 8
memory "20GB"
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "GATK_Bam", mode:'copy', pattern: "*gatk*"
input:
        tuple val(replicateId), path(bam1)
        val(ref)
        path(ascii)
output:
        tuple val(replicateId), path("*gatk*"), emit: gatk
        tuple val(replicateId), path("GATK_Bam"), emit: GATK_Bam
	val 'done', emit:done
script:
"""
mkdir GATK_Bam
mkdir Bam
mkdir temp_Bam
ss=\$(ls -1 *bam | cut -f2 -d '_' | sort -u)
mv *bam* Bam/
cat <(samtools view -H Bam/${replicateId}_\${ss}_nBP.bam) <(paste <(samtools view Bam/${replicateId}_\${ss}_nBP.bam | cut -f1-10 ) <(samtools view Bam/${replicateId}_\${ss}_nBP.bam | cut -f 11 | tr "\$(cat ${ascii})" "K")) | samtools view -b -o temp_Bam/${replicateId}_\${ss}_dump.bam -
picard AddOrReplaceReadGroups I=temp_Bam/${replicateId}_\${ss}_dump.bam O=temp_Bam/${replicateId}_\${ss}_final.bam RGPU=unit1 RGID=11 RGLB=LaneX RGSM=AnySampleName RGPL=illumina 2>> GATK_Bam/${replicateId}_\${ss}.gatk.bamlog || echo "processed \$?"
samtools index temp_Bam/${replicateId}_\${ss}_final.bam

gatk3 -Xmx30g --analysis_type RealignerTargetCreator --reference_sequence /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta --input_file temp_Bam/${replicateId}_\${ss}_final.bam --downsample_to_coverage 10000 --num_threads ${task.cpus} --out GATK_Bam/${replicateId}_\${ss}.gatk.intervals 2>> GATK_Bam/${replicateId}_\${ss}.gatk.bamlog
gatk3 -Xmx30g --analysis_type IndelRealigner --reference_sequence /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta --input_file temp_Bam/${replicateId}_\${ss}_final.bam --defaultBaseQualities 4 --targetIntervals GATK_Bam/${replicateId}_\${ss}.gatk.intervals --noOriginalAlignmentTags --out GATK_Bam/${replicateId}_\${ss}.realigned.bam 2>> GATK_Bam/${replicateId}_\${ss}.gatk.bamlog
gatk3 -Xmx30g --analysis_type BaseRecalibrator --reference_sequence /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta --input_file GATK_Bam/${replicateId}_\${ss}.realigned.bam --knownSites /opt/conda/share/mtbseq-1.0.4-2/var/res/MTB_Base_Calibration_List.vcf --maximum_cycle_value 400000  --num_cpu_threads_per_data_thread ${task.cpus} --out GATK_Bam/${replicateId}_\${ss}.gatk.grp 2>>GATK_Bam/${replicateId}_\${ss}.gatk.bamlog
gatk3 -Xmx30g -T --analysis_type PrintReads --reference_sequence /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta --input_file GATK_Bam/${replicateId}_\${ss}.realigned.bam --BQSR GATK_Bam/${replicateId}_\${ss}.gatk.grp --num_cpu_threads_per_data_thread ${task.cpus} --out GATK_Bam/${replicateId}_\${ss}_nBP.gatk.bam 2>> GATK_Bam/${replicateId}_\${ss}.gatk.bamlog
samtools index GATK_Bam/${replicateId}_\${ss}_nBP.gatk.bam
rm GATK_Bam/*.realigned.*
rm -r temp_Bam
ln -s GATK_Bam/* .
"""
}







/*
* Pileup
*/

process PILE {
memory "20GB"
cpus 8
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Mpileup", mode:'copy', pattern: "*mpileup*"
input:
        tuple val(replicateId), path(bam)
        val(ref)
output:
        tuple val(replicateId), path("*mpileup*"), emit: mpile
        tuple val(replicateId), path("Mpileup"), emit: MPILE
	val 'done', emit:done
script:
"""
mkdir Mpileup
mkdir GATK_Bam
mv *gatk* GATK_Bam
USER=a perl /opt/conda/bin/MTBseq --step TBpile --threads 8 --ref ${ref} || echo "processed \$?"
ln -s Mpileup/* .
"""
}



process PILE_ONT {
memory "20GB"
cpus 8
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Mpileup", mode:'copy', pattern: "*mpileup*"
input:
        tuple val(replicateId), path(bam)
        val(ref)
        val(minbqual)
output:
        tuple val(replicateId), path("*mpileup*"), emit: mpile
        tuple val(replicateId), path("Mpileup"), emit: MPILE
	val 'done', emit:done
script:
"""
mkdir Mpileup
mkdir GATK_Bam
ss=\$(ls -1 *bam | cut -f2 -d '_' | sort -u)
mv *gatk* GATK_Bam
samtools mpileup -B -A -x -Q ${minbqual} -f /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta -o Mpileup/${replicateId}_\${ss}_nBP.gatk.mpileup GATK_Bam/${replicateId}_\${ss}_nBP.gatk.bam
ln -s Mpileup/* .
"""
}

/*
* Position Table
*/

process LIST {
memory "20GB"
cpus 8
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Position_Tables", mode:'copy', pattern: "*position_table*"
input:
        tuple val(replicateId), path(pile)
	val(minbq)
        val(ref)
output:
        tuple val(replicateId), path("*position_table*"), emit: list
        tuple val(replicateId), path("Position_Tables"), emit: LIST
script:
"""
mkdir Mpileup
mv *mpileup* Mpileup/
mkdir Position_Tables
USER=a perl /opt/conda/bin/MTBseq --step TBlist --threads 8 --minbqual $minbq --ref ${ref}|| echo "processed \$?"
ln -s Position_Tables/* .
"""
}

/*
* Variants_low
*/

process VARIANTS_LOW {
memory "20GB"
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Called", mode:'copy', pattern: "*tab"
input:
        tuple val(replicateId), path(list)
        val(ref)
output:
        tuple val(replicateId), path("*variants*tab"), emit: var_low
        tuple val(replicateId), path("Called"), emit: VAR_LOW
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Called
USER=a perl /opt/conda/bin/MTBseq --step TBvariants --ref ${ref} --mincovf 1 --mincovr 1 --lowfreq_vars --minfreq 5 --minphred20 1 || echo "processed \$?"
ln -s Called/* .
"""
}

/*
* Variants
*/


process VARIANTS {
memory "20GB"
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Called", mode:'copy', pattern: "*tab"
input:
        tuple val(replicateId), path(list)
	val(mincovf)
	val(mincovr)
	val(minphred)
        val(ref)
output:
        tuple val(replicateId), path("*tab"), emit: var
        tuple val(replicateId), path("Called"), emit: VAR
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Called
USER=a perl /opt/conda/bin/MTBseq --step TBvariants --ref ${ref}  --mincovf $mincovf --mincovr $mincovr --minphred20 $minphred || echo "processed \$?"
ln -s Called/* .
"""
}



process STATS {
memory "20GB"
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Statistics", mode:'copy', pattern: "*tab"
input:
        tuple val(replicateId), path(bam), path(list)
		val(mincovf)
        val(mincovr)
        val(minphred)
output:
        tuple val(replicateId), path("*tab"), emit: stats
        tuple val(replicateId), path("Statistics"), emit: STATS
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Bam
mv *bam* Bam/
mkdir Statistics
USER=a perl /opt/conda/bin/MTBseq --step TBstats  --mincovf $mincovf --mincovr $mincovr --minphred20 $minphred || echo "processed \$?"
mv Statistics/Mapping_and_Variant_Statistics.tab Statistics/${replicateId}_Mapping_and_Variant_Statistics.tab
ln -s Statistics/* .
"""
}



process JOIN {
memory "60GB"
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
publishDir "Joint", mode: "copy", pattern: "Joint/*"
publishDir "Amend", mode: "copy", pattern: "Amend/*"
publishDir "Groups", mode: "copy", pattern: "Groups/*"
input:
        path(cc)
        path(ll)
        path(sample_joint)
	val(minbqual)
	val(minphred20)
        val(proj)
        val(ref)
output:
        path("Joint/*"), emit: joint
        path("Amend/*"), emit: amend
        path("Groups/*"), emit: groups
script:
def joint_select = sample_joint.name != 'placehold' ? "| awk 'NR==FNR{a[\$1];next}(\$1 in a)' ${sample_joint} -" : ''
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Called
mv *variants_cf4* Called/
mkdir Joint
mkdir Amend
mkdir Groups

ls -1 Called/*_variants_cf4* | cut -f2 -d'/' | cut -f1,2 -d '_' | tr '_' '\\t' | sort -r | sort -u -k1,1 $joint_select  > sample_joint



USER=a perl /opt/conda/bin/MTBseq --step TBjoin --continue --ref ${ref} --samples sample_joint --distance 5 --project ${proj} --minbqual ${minbqual} --minphred20 ${minphred20} || echo "processed \$?"
"""

}

process STRAIN {
memory "20GB"
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
tag "$replicateId"
publishDir "Classification", mode:'copy', pattern: "*tab"
input:
        tuple val(replicateId), path(list)
output:
        tuple val(replicateId), path("*.tab"), emit: strain
		tuple val(replicateId), path("Classification"), emit: STRAIN
script:
"""
mkdir Position_Tables
mv *position_table* Position_Tables
mkdir Classification
USER=a perl /opt/conda/bin/MTBseq --step TBstrains || echo "processed \$?"
mv Classification/Strain_Classification.tab Classification/${replicateId}_Strain_Classification.tab
ln -s Classification/* .
"""
}




process MAP_STRAIN{
memory '5GB'
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
publishDir "OUTPUT", mode: 'copy', pattern: "Mapping_Classification.tab"
input:
	path(stats)
output:
	path("Mapping_Classification.tab")
script:
"""
for i in \$(ls -1 *Mapping_and_Variant_Statistics* | cut -f1 -d '_' | sort -u)
do
paste \${i}_Mapp* \${i}_Strain* | cut --complement -f1,25  >> Mapping_Classification.tab
done

sort -u -r -k1 Mapping_Classification.tab | uniq >Mapping_Classification.tab.tmp
mv Mapping_Classification.tab.tmp Mapping_Classification.tab

chmod 666 Mapping_Classification.tab
"""
}


process DEL {
tag "$replicateId"
cpus 1
errorStrategy 'ignore'
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory '20GB'
input:
	tuple val(replicateId), path(bam)
	val(ref)
	path(bed)
	path(bedix)
output:
	tuple val(replicateId), path("*.dels")
script:
"""
delly call -o ${replicateId}.bcf -g /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta ${replicateId}*.bam

bcftools view -i '((SVTYPE="DEL" || SVTYPE="INS") && FILTER="PASS" && (END - POS + 1) <= 100000) || (SVTYPE="INS" && INFO/INSLEN <= 100000)' ${replicateId}.bcf | bcftools filter -e 'FORMAT/FT!="PASS"' |bcftools query -f '%CHROM\\t%POS\\t%END\\t%SVTYPE\\t%REF\\t[%DR]\\t[%DV]\\t[%RR]\\t[%RV]\\t%PRECISE\\t%INFO/INSLEN\\n' | awk 'BEGIN { OFS="\\t" } { len = (\$4 == "DEL") ? (\$3 - \$2 - 1) : \$11; print \$1, \$2, \$3, \$4, \$5, \$6 + \$8, \$7 + \$9, \$10, (\$7+\$9) / (\$6 + \$7 +\$8 +\$9), len }' > variants.bed

# Find overlaps with genes
bedtools intersect -a variants.bed -b ${bed} -wa -wb -loj > overlaps.txt

# Adjust start and end positions to match overlapping part of gene and include variant length
awk 'BEGIN { OFS=";" } { if (\$12 == -1) { print \$2 + 1 , \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10 } else { start = (\$2 > \$12) ? \$2 : \$12; end = (\$3 < \$13) ? \$3 : \$13; print start + 1, end, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$14 } }' overlaps.txt > ${replicateId}.dels
"""
}


process DEL_ONT {
tag "$replicateId"
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory '20GB'
input:
	tuple val(replicateId), path(bam)
	val(ref)
	path(bed)
	path(bedix)
output:
	tuple val(replicateId), path("*.dels")
script:
"""
delly lr -o ${replicateId}.bcf -g /opt/conda/share/mtbseq-1.0.4-2/var/ref/${ref}.fasta ${replicateId}*.bam

bcftools annotate -a ${bed} -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') ${replicateId}.bcf | grep -P "=DEL|=INS" | grep "PASS" | grep -i -v "lowQual" | cut -f2,4,5,8,10 |  awk '{if (\$4 ~ "SVLEN") {print \$0,gensub(".*END=([0-9]*);.*","\\\\1","g",\$4),gensub(".*SVLEN=([0-9]*);.*","\\\\1","g",\$4)} else {print \$0, gensub(".*END=([0-9]*);.*","\\\\1","g",\$4),0}}' | awk '{ if (\$4 ~ "GENE") {print \$1,\$6,gensub(".*SVTYPE=([A-Z]*);.*","\\\\1","g",\$4),\$2,gensub(".*:([^:]*):([^:]*\$)","\\\\1;\\\\2","g",gensub(":0:0\$","","g",\$5)),\$7,gensub("(^[PI]).*","\\\\1","g",\$4),gensub(".*GENE=","","g",\$4)} else {print \$1,\$6,gensub(".*SVTYPE=([A-Z]*);.*","\\\\1","g",\$4),\$2,gensub(".*:([^:]*):([^:]*\$)","\\\\1;\\\\2","g",gensub(":0:0\$","","g",\$5)),\$7,gensub("(^[PI]).*","\\\\1","g",\$4)}}' OFS=';' |  awk -F ';' '{if (\$7+\$2-\$1 < 100000) {print \$1,\$2, \$7+\$2-\$1,\$3,\$4,\$5,\$6,\$6/(\$6+\$5),\$8,\$9}}' OFS=';' > ${replicateId}.dels
"""     
}


process OUT_DEL {
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory '20GB'
publishDir 'OUTPUT', mode: "copy", pattern: "DELETIONS.tab"
input:
	path(del)
output:
	path("DELETIONS.tab")
script:
"""

awk '{print gensub("\\.dels","","g",FILENAME),\$0}' OFS='\\t' *.dels | tr ';' '\\t' > DELETIONS.tab


sort -k1,1 -k2,2n -u -o DELETIONS.tab DELETIONS.tab

chmod 666 DELETIONS.tab

"""
}



process DEPTH {
cpus 8
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory "40GB"
tag "$replicateId"
input:
	tuple val(replicateId), path(bam)
	path(tgene)
output:
	tuple val(replicateId), path("gbcov_*")
script:
"""
mosdepth -b target_genes.bed -n -T ${task.cpus} ${replicateId} ${replicateId}*.bam

zcat ${replicateId}.thresholds.bed.gz | tail -n +2 | awk -v value="GB" '\$4 == value' > ${replicateId}.t.bed
zcat ${replicateId}.thresholds.bed.gz | tail -n +2 | awk -v value="GB" '\$4 != value' | sort --ignore-case -k4,4 >> ${replicateId}.t.bed

zcat ${replicateId}.regions.bed.gz | tail -n +2 | awk -v value="GB" '\$4 == value' > ${replicateId}.r.bed 
zcat ${replicateId}.regions.bed.gz | tail -n +2 | awk -v value="GB" '\$4 != value' | sort --ignore-case -k4,4 >> ${replicateId}.r.bed

cat <(cat ${replicateId}.t.bed | awk '{print \$1,\$2,\$3,\$4,\$5/(\$3-\$2)}' OFS='\\t' | sed 's/GB/GB_perc/g') <(cat ${replicateId}.r.bed) | cut -f4,5 | awk -F '\t' '{for (i=1; i<=NF; i++) {a[NR,i] = \$i}}; NF>p { p = NF }; END {for (j=1; j<=p; j++) {str=a[1,j]; for (i=2; i<=NR; i++){str=str";"a[i,j];};print str}}' > gbcov_${replicateId}

rm ${replicateId}*

"""

}



process OUT_DEPTH {
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory '20GB'
publishDir 'OUTPUT', mode: "copy", pattern: "GB_cov.csv"
input:
	path(cov)
output:
	path("GB_cov.csv")
script:
"""
for i in gbcov_*;do awk '{if (NR == 1) {print "000000",\$0} else {print gensub(".*_","","g",FILENAME),\$0} }' OFS=';' \${i}; done | sort -u -k1 -t ';' | sed 's/000000/Samples/g' > GB_cov.csv

sort -u -o GB_cov.csv GB_cov.csv
sort -k2 -r -t ';' -o GB_cov.csv GB_cov.csv

chmod 666 GB_cov.csv
"""
}


process FINAL_OUT {
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory '10GB'
publishDir 'OUTPUT', mode: "copy", pattern: "FINAL_OUT.csv"
input:
	path(map_strain)
	path(depth)
output:
	path("FINAL_OUT.csv")
script:
"""
#!/usr/bin/env Rscript

library(tidyverse)

m=read_delim('Mapping_Classification.tab')
g=read_delim('GB_cov.csv') %>% 
  select(Samples,Genome_Breadth=GB_perc,Genome_Depth=GB)



m %>% 
  mutate(across(everything(),function(x){str_remove_all(x,"'")})) %>% 
  select(Samples=SampleID...1,Map_reads=`Mapped Reads`,Map_reads_perc=`% Mapped Reads`,`Homolka species`:`Beijing quality (easy)`)->m1


g %>% left_join(m1) %>% write_delim('FINAL_OUT.csv',delim=';')

Sys.chmod('FINAL_OUT.csv', mode = "0777")

"""
}

process MUT_CORRECTION_DEL {
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory "20GB"
tag "$replicateId"
publishDir "Called", mode:'copy', pattern: "*corrected.tab"
input:
	tuple val(replicateId), path(var), path(del)
output:
	tuple val(replicateId), path('*corrected.tab')

script:
"""
#!/usr/bin/env Rscript

library(tidyverse)
library(Biostrings)



l=list.files()



l<-l[str_detect(l,'variants_cf1_cr1_fr5_ph1_outmode001.tab\$')]

for (i in l){
 print(i)
  s=str_extract(i,'[^_]*_[^_]*')


  read_delim(i,delim='\\t',show_col_types = FALSE) %>%
        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character)) %>%
        filter(Subst!=' ',Type=='SNP') %>%
    extract(Subst,c('refA','pA','varA','codR','cV1','cV2','cV3'),'([A-z]*)([0-9]*)([A-z]*) \\\\((.*)/(.)(.)(.)\\\\)') %>%
    group_by(Gene,pA) %>%
    filter(n()>1,sd(Freq)<15) %>%
    mutate(across(starts_with('cV'),function(x){ifelse(sum(LETTERS%in%x)>0,x[x%in%LETTERS],x)})) %>%
    unite('codV',cV1:cV3,sep='') %>%
    {if(dim(.)[1]>0) mutate(.,varA=plyr::mapvalues(as.character(translate(DNAString(codV[1]))),from = names(AMINO_ACID_CODE),to=AMINO_ACID_CODE,warn_missing = F)) else
            mutate(.,varA='A')}%>%
    arrange(`#Pos`) %>%
    mutate(Ref=ifelse(n()==2 &max(diff(`#Pos`))==2,codR,(str_c(Ref,collapse=''))),
           Allel=ifelse(n()==2 &max(diff(`#Pos`))==2,codV,(str_c(Allel,collapse='')))#,
           #`#Pos`=str_c(`#Pos`,collapse=',')
    ) %>%
    #mutate(across(c('Ref','Allel'),function(x){ifelse(str_detect(Gene,'c'),as.character(reverseComplement(DNAString(x[1]))),x[1])})) %>%
    unite('Subst',refA,pA,varA,sep='') %>% mutate(Subst=paste(Subst,' (',codR,'/',codV,')',sep='')) %>% select(-codR,-codV) %>%
    mutate(across(where(is.logical),function(x){x=as.character(x)})) %>%
    rbind(read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                       across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
        filter(Subst!=' ',Type=='SNP') %>%
            extract(Subst,c('refA','pA','varA','codR','cV1','cV2','cV3'),'([A-z]*)([0-9]*)([A-z]*) \\\\((.*)/(.)(.)(.)\\\\)') %>%
            group_by(Gene,pA) %>%
            filter(n()<2|sd(Freq)>=15) %>%  unite('codV',cV1:cV3,sep='') %>%
            unite('Subst',refA,pA,varA,sep='') %>% mutate(Subst=paste(Subst,' (',codR,'/',codV,')',sep='')) %>%
            #mutate(`#Pos`=as.character(`#Pos`)) %>%
            select(-codR,-codV),
          read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
            #mutate(`#Pos`=as.character(`#Pos`)) %>%
            filter(Subst==' ',Type=='SNP'),
          read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
        filter(Type!='SNP') %>% arrange(Type,`#Pos`) %>%
            group_by(Type) %>%
            mutate(p=cumsum(c(1,diff(`#Pos`)>1))) %>%
            group_by(Type,Subst,Gene,GeneName,Product,ResistanceSNP,PhyloSNP,InterestingRegion,p) %>%
            summarise(`#Pos`=min(`#Pos`),Insindex=0,Ref='_',Allel=str_to_upper(Type),across(CovFor:Cov,mean)) %>% arrange(Type,p) %>% ungroup() %>%
            select(-p)
    ) ->a
    #%>% write_delim(paste(i,'corrected.tab',sep='_'),delim='\\t')



  if (file.exists(paste("${replicateId}",'.dels',sep=''))) {

    a %>% bind_rows(read_delim(paste("${replicateId}",'.dels',sep=''),show_col_types = FALSE,delim=';',
                         col_names = c('Start','End','Type','Ref','RefR','VarR','Precision','Freq','Length','Gene'))%>% 
                         filter(as.character(Precision)=="1")%>%
                         separate(Length,c('Length','Gene'),sep=';') %>%
                {if(dim(.)[1]>0) mutate(.,Insindex=0,Ref='_',Allel=Type,Type=str_to_title(Type),Subst=" ",GeneName='-', Product=" ",Freq=Freq*100,Qual20=RefR+VarR) %>%
                select(`#Pos`=Start,Insindex,Ref,Type,Allel,Subst,Gene,GeneName,Product,Freq,Qual20)})->a}
        a%>%bind_rows(filter(unique(a%>%filter(Type != 'SNP') %>% arrange(Type, '#Pos') %>% add_count(across(everything()))) %>% ungroup() %>% mutate(Allel= case_when(n %% 3!=0 ~ 'LOF', TRUE ~as.character(Type))), Allel=='LOF'))->a
  a%>%mutate(n=ifelse(is.na(CovFor),'long',as.character(n)))->a

  a%>%filter(str_detect(Subst,'^[^_]+[1-9]+_'))%>%
  mutate(Ref=as.character(ifelse(str_detect(Subst,'^[^_]+[1-9]+_'),'_',Ref)),
        Allel=as.character(ifelse(str_detect(Subst,'^[^_]+[1-9]+_'),'STOP',Allel)))->b

   a%>%filter(str_detect(Subst,'^Met1[A-z]+|^Val1[A-z]+')) %>% 
    filter(!str_detect(Subst,'^(Val1|Met1)(Val|Met)')) %>% 
    mutate(Ref='START',Allel='STARTLOSS')->c
  a %>%bind_rows(b)%>%
  bind_rows(c)%>%
write_delim(paste(i,'corrected.tab',sep='_'),delim='\\t')
Sys.chmod(paste(i,'corrected.tab',sep='_'), mode = "0777")
}
"""
}


process MUT_CORRECTION {
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory "20GB"
tag "$replicateId"
publishDir "Called", mode:'copy', pattern: "*corrected.tab"
input:
	tuple val(replicateId), path(var)
output:
	tuple val(replicateId), path('*corrected.tab')

script:
"""
#!/usr/bin/env Rscript

library(tidyverse)
library(Biostrings)



l=list.files()



l<-l[str_detect(l,'variants_cf1_cr1_fr5_ph1_outmode001.tab\$')]

for (i in l){
 print(i)
  s=str_extract(i,'[^_]*_[^_]*')


  read_delim(i,delim='\\t',show_col_types = FALSE) %>%
        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character)) %>%
        filter(Subst!=' ',Type=='SNP') %>%
    extract(Subst,c('refA','pA','varA','codR','cV1','cV2','cV3'),'([A-z]*)([0-9]*)([A-z]*) \\\\((.*)/(.)(.)(.)\\\\)') %>%
    group_by(Gene,pA) %>%
    filter(n()>1,sd(Freq)<15) %>%
    mutate(across(starts_with('cV'),function(x){ifelse(sum(LETTERS%in%x)>0,x[x%in%LETTERS],x)})) %>%
    unite('codV',cV1:cV3,sep='') %>%
    {if(dim(.)[1]>0) mutate(.,varA=plyr::mapvalues(as.character(translate(DNAString(codV[1]))),from = names(AMINO_ACID_CODE),to=AMINO_ACID_CODE,warn_missing = F)) else
            mutate(.,varA='A')}%>%
    arrange(`#Pos`) %>%
    mutate(Ref=ifelse(n()==2 &max(diff(`#Pos`))==2,codR,(str_c(Ref,collapse=''))),
           Allel=ifelse(n()==2 &max(diff(`#Pos`))==2,codV,(str_c(Allel,collapse='')))#,
           #`#Pos`=str_c(`#Pos`,collapse=',')
    ) %>%
    #mutate(across(c('Ref','Allel'),function(x){ifelse(str_detect(Gene,'c'),as.character(reverseComplement(DNAString(x[1]))),x[1])})) %>%
    unite('Subst',refA,pA,varA,sep='') %>% mutate(Subst=paste(Subst,' (',codR,'/',codV,')',sep='')) %>% select(-codR,-codV) %>%
    mutate(across(where(is.logical),function(x){x=as.character(x)})) %>%
    rbind(read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                       across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
        filter(Subst!=' ',Type=='SNP') %>%
            extract(Subst,c('refA','pA','varA','codR','cV1','cV2','cV3'),'([A-z]*)([0-9]*)([A-z]*) \\\\((.*)/(.)(.)(.)\\\\)') %>%
            group_by(Gene,pA) %>%
            filter(n()<2|sd(Freq)>=15) %>%  unite('codV',cV1:cV3,sep='') %>%
            unite('Subst',refA,pA,varA,sep='') %>% mutate(Subst=paste(Subst,' (',codR,'/',codV,')',sep='')) %>%
            #mutate(`#Pos`=as.character(`#Pos`)) %>%
            select(-codR,-codV),
          read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
            #mutate(`#Pos`=as.character(`#Pos`)) %>%
            filter(Subst==' ',Type=='SNP'),
          read_delim(i,delim='\\t',show_col_types = FALSE) %>%
                        mutate(across(c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.numeric),
                               across(-c('#Pos','Insindex','CovFor','CovRev','Qual20','Freq','Cov'),as.character))%>%
        filter(Type!='SNP') %>% arrange(Type,`#Pos`) %>%
            group_by(Type) %>%
            mutate(p=cumsum(c(1,diff(`#Pos`)>1))) %>%
            group_by(Type,Subst,Gene,GeneName,Product,ResistanceSNP,PhyloSNP,InterestingRegion,p) %>%
            summarise(`#Pos`=min(`#Pos`),Insindex=0,Ref='_',Allel=str_to_upper(Type),across(CovFor:Cov,mean)) %>% arrange(Type,p) %>% ungroup() %>%
            select(-p)
    ) ->a
    #%>% write_delim(paste(i,'corrected.tab',sep='_'),delim='\\t')



  if (file.exists(paste("${replicateId}",'.dels',sep=''))) {

    a %>% bind_rows(read_delim(paste("${replicateId}",'.dels',sep=''),show_col_types = FALSE,delim=';',
                         col_names = c('Start','End','Type','Ref','RefR','VarR','Precision','Freq','Length','Gene'))%>% separate(Length,c('Length','Gene'),sep=';') %>%
                {if(dim(.)[1]>0) mutate(.,Insindex=0,Ref='_',Allel=Type,Type=str_to_title(Type),Subst=" ",GeneName='-', Product=" ",Freq=Freq*100,Qual20=RefR+VarR) %>%
                select(`#Pos`=Start,Insindex,Ref,Type,Allel,Subst,Gene,GeneName,Product,Freq,Qual20)})->a}
        a%>%bind_rows(filter(unique(a%>%filter(Type != 'SNP') %>% arrange(Type, '#Pos') %>% add_count(across(everything()))) %>% ungroup() %>% mutate(Allel= case_when(n %% 3!=0 ~ 'LOF', TRUE ~as.character(Type))), Allel=='LOF'))->a
  
  a%>%filter(str_detect(Subst,'^[^_]+[1-9]+_'))%>%
  mutate(Ref=as.character(ifelse(str_detect(Subst,'^[^_]+[1-9]+_'),'_',Ref)),
        Allel=as.character(ifelse(str_detect(Subst,'^[^_]+[1-9]+_'),'STOP',Allel)))->b

   a%>%filter(str_detect(Subst,'^Met1[A-z]+|^Val1[A-z]+')) %>% 
    filter(!str_detect(Subst,'^(Val1|Met1)(Val|Met)')) %>% 
    mutate(Ref='START',Allel='STARTLOSS')->c
  a %>%bind_rows(b)%>%
  bind_rows(c)%>%
write_delim(paste(i,'corrected.tab',sep='_'),delim='\\t')
Sys.chmod(paste(i,'corrected.tab',sep='_'), mode = "0777")
}
"""
}

process MUT_GATHER{
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory "10GB"
input:
        path(tabs)
output:
        path('all_cf1N.tab')
script:
"""
#!/usr/bin/env Rscript
library(tidyverse)
files <- list.files(path = ".", pattern = "corrected.tab")
lapply(files, function(x){
    read_delim(x,delim='\\t',
    col_types='nncccnnnnncccccccc') %>%
    mutate(File=str_remove_all(x,'_.*'))%>%
    relocate(File)
})->l
bind_rows(l)%>%
write_delim('all_cf1N.tab',delim='\\t',col_names = F,na='0')
Sys.chmod('all_cf1N.tab', mode = "0777")
"""
}



process PHARMA {
cpus 1
container "library://allen13x/mtbseq/nf_mtbseq:1.0.1"
memory "10GB"
publishDir "Called", mode: 'copy', pattern: 'pharma_gene*'
publishDir "OUTPUT", mode: 'copy', pattern: '*format*'
input:
        path (tabs)
        val (thresholds)
        path (genes)

output:
        path("pharma_gene*"), emit:pgene
        path("*format*"), emit:format
script:
"""
#!/usr/bin/env Rscript
library(tidyverse)

t=${thresholds}

# Read in the data
#read all the file that ends with "corrected.tab"
files <- list.files(path = ".", pattern = "corrected.tab")

lapply(files, function(x){
    read_delim(x,delim='\\t',
    col_types='nncccnnnnncccccccc') %>%
    filter(Freq>=t) %>%
    mutate(ID=str_remove_all(x,'_.*'))
})->l

genes<-read_delim('gene_drug.csv',col_names=c('Gene_name','Start','Stop','Gene','S'))

#from genes get a vector containing all the numbers in the intervavals between Start and Stop

lapply(genes\$Start,function(x){
    seq(x,genes\$Stop[genes\$Start==x][1])
})->l2

unlist(l2)->gene_intervals



bind_rows(l) %>%
    filter(`#Pos` %in% gene_intervals) %>%
    relocate(ID)->p


genes %>% 
  group_by(Gene_name) %>% 
  mutate(i=str_c(seq(Start,Stop),collapse=',')) %>% 
  separate_rows(i,sep=',') %>% 
  transmute(Name=Gene_name,
            Range=ifelse(Start>Stop,
                         paste('r[',Stop,'-',Start,']',sep=''),
                         paste('r[',Start,'-',Stop,']',sep='')),
            Locus=Gene,`#Pos`=as.numeric(i),
            Start1=ifelse(S=='rev',max(Start,Stop),min(Start,Stop)),
            Start=Start,Stop=Stop,S=S) %>% 
  right_join(p) %>% 
  select(Name,Range,Locus,`#Pos`,ID,Start1,Freq,Ref,Allel,Type,Subst,Start,Stop,S) %>%
  mutate(
    Type=ifelse(Type=='SNP',Type,paste('__',Type,sep='')),
    Mut=case_when((str_detect(Name,'_ups')& S=='rev')~paste0(Ref,'-',abs(`#Pos`-Start)+1,Allel,'-',`#Pos`),
                  (str_detect(Name,'_ups'))~paste0(Ref,'-',abs(`#Pos`-Stop)+1,Allel,'-',`#Pos`),
                       (Subst!=' ')~str_remove(Subst,' .*'),
                  (Type=='SNP')~paste0(`#Pos`,'_',Ref,'>',Allel,'_',Type),
                  (Type!='SNP')~paste0(`#Pos`,'_',str_to_upper(Type))),
    Mut=ifelse(Freq<75,paste0(Mut,':',Freq),Mut)) %>% 
  ungroup() %>% 
  select(Name,Range,LocusName=Locus,ID,Mut) %>% 
  group_by(Name,Range,LocusName,ID) %>% 
  summarise(Mut=str_c(Mut,collapse=' ')) %>%ungroup() %>%  
  select(Name,ID,Mut) %>% 
  spread(Name,Mut) %>% rename(Name=ID)->tail


genes %>% 
  transmute(Name=Gene_name,
            Range=ifelse(Start>Stop,
                         paste('r[',Stop,'-',Start,']',sep=''),
                         paste('r[',Start,'-',Stop,']',sep='')),
            LocusName=Gene)->head


as.data.frame(t(head[,-1]))->h
colnames(h)<-head\$Name
h %>% rownames_to_column('Name') %>% as_tibble()->head



bind_rows(list(head,tail)) %>% 
  write_delim(paste0('pharma_gene_format_',t,'.tab'),delim='\\t',na = '0')
Sys.chmod(paste0('pharma_gene_format_',t,'.tab'), mode = "0777")


p%>%
write_delim('pharma_gene.tab',delim='\\t',col_names=F)
Sys.chmod("pharma_gene.tab", mode = "0777")

"""
}



process  WHO {
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory "10GB"
input:
        path(tab)
        path(head)
        path(WHO)
output:
        path('WHO_raw.csv')
script:
"""
#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
from functools import reduce


#Check wheter a subs is Syn
def syn(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res = temp.match(s).groups()
    if res[0] == res[2]:
        print(s)

#Check wheter a subs is NoSyn
def nsyn(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res = temp.match(s).groups()
    if res[0] != res[2]:
        print(s)
        
def nsyn2(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res1 = temp.match(s).group(1)
    res2 = temp.match(s).group(2)
    res3 = temp.match(s).group(3)
    if res1 != res3:
        print(s)
        
# Check whether there is STOP codon       
def nsyn3(s):
    import re
    temp = re.compile("([a-zA-Z_]+)([0-9]+)([a-zA-Z_]+)")
    #check whether s IS NOT white space
    if s and not s.isspace():
        res1 = temp.match(s).group(1)
        res2 = temp.match(s).group(2)
        res3 = temp.match(s).group(3)
        if res1 != res3:
            return(True)
        else:
            return(False)
    else:
        return(False)

# remove content between () in s
def cleanvar(s):
    import re
    return(re.sub(r"\\([^()]*\\)", "", s))  

# convert tab to space
def tab2space(s):
    import re
    return(s.replace('\\t', ' '))

def removespace(s):
    import re
    return(s.replace(' ', ''))

def remove_(s):
    S=s.split('_')
    return(S[0])

def remove1_(s):
    S=s.split('_')
    return(S[1])

def first(s):
    return(s[0])

def last(s):
    return(s[-1])

def one2tree(s):
    from Bio.PDB.Polypeptide import one_to_three
    s=s.rstrip()
    if s.isupper() and not s.isnumeric() :
        try:
            return(one_to_three(s))
        except:
            return(s)
    else:
        return(s)
    
def codon(s):
    import re
    temp = re.findall(r'\\d+', s)
    return(temp[0])


myfile = open("head", "r")
head1 = myfile.read()
head_list1=head1.split('\\t')
head_list1[16]=head_list1[16].strip()
head_list1.insert(0,"File")
all=pd.read_csv('all_cf1N.tab',sep='\\t',header=None,names=head_list1)


#clean substitution and add a new column SubstClean
all['SubstClean']=all.Subst.apply(cleanvar)
all['SubstClean2'] = all.SubstClean.apply(removespace)
#all['ID']=all.File.apply(remove_)
all['ID']=all.File
#all['NSYN'] = all.SubstClean2.apply(nsyn3)
all.rename({'#Pos': 'Genome position'}, axis=1,inplace=True)
all['genome_index']=all['Genome position'].astype(int)
all.genome_index.astype(int)
whoG=pd.read_csv('${WHO}',sep='\\t')

whoG.genome_index=whoG.genome_index.astype(str).str.split(',')
whoG = whoG.explode('genome_index').reset_index(drop=True)
whoG.genome_index=whoG.genome_index.astype(int)
whoG['Allel']=whoG['final_annotation.AlternativeNucleotide'].str.upper()
whoG['Ref']=whoG['final_annotation.ReferenceNucleotide'].str.upper()
whoG['common_name']=whoG['final_annotation.Gene'] + '_' + whoG['final_annotation.TentativeHGVSNucleotidicAnnotation']



# 1) match esatti (no wildcard)
whoG_no_wc = whoG[(whoG['Allel'] != 'ANY') & (whoG['Ref'] != 'ANY')]
exact = pd.merge(all, whoG_no_wc, on=['genome_index','Allel','Ref'], how='inner')

# 2) anti-join per ottenere le righe di `all` che NON hanno match esatto
tmp = all.merge(whoG_no_wc, on=['genome_index','Allel','Ref'], how='left', indicator=True)
non_matched = tmp[tmp['_merge'] == 'left_only'][all.columns]   # prendo solo le colonne originali di `all`

# 3) mantieni solo quelle nel range da riempire con wildcard
non_matched_range = non_matched[non_matched['genome_index'].between(761083, 761160)]

# 4) match di quelle con i wildcard (solo su genome_index)
whoG_wc = whoG[(whoG['Allel'] == 'ANY') & (whoG['Ref'] == 'ANY')]
wildcard_matches = pd.merge(non_matched_range, whoG_wc, on='genome_index', how='inner')
# 5) Rimuovi sinonime
wildcard_matches = wildcard_matches[wildcard_matches['SubstClean2'].apply(nsyn3)]

all_whoG = pd.concat([exact, wildcard_matches], ignore_index=True)

all_whoG['File'] = all_whoG['File'].astype(str)
all_whoG.sort_values('File',inplace=True)


all_whoG.filter(regex=r'(File|_Conf_Grade)')
drugs_conf=all_whoG.filter(regex=r'(_Conf_Grade)').columns
# create array of drugs
drugs=drugs_conf.str.split('_').str[0]

cutoff=10
A = {}
Asubst = {}
Asubst1 = {}
Asubst5 = {}
Ali=[]
Alisubst=[]
Alisubst5=[]
Alisubst1=[]
for i in drugs_conf:
    t = i.split('_')[0]
    A[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File',i]]
    Asubst[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Freq','Type','LOF','Genome position','Subst',i]]
    Asubst1[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim') | (all_whoG[i]=='3) Uncertain significance')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Freq','Type','LOF','Genome position','Subst',i]]
    Asubst5[t] = all_whoG[(all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Type','Freq','LOF','Genome position','Subst',i]]
    Ali.append(A[t])
    Alisubst.append(Asubst[t])
    Alisubst5.append(Asubst5[t])
    Alisubst1.append(Asubst1[t])


A5=pd.concat(Alisubst5)

A5.fillna(0,inplace=True)


A5a=A5.drop_duplicates()
A5a['variant']=np.where(A5a['variant'].str.contains("lof"),A5a['variant'] + '('+ A5a['Genome position'].astype(str)+':'+ A5a['Type'] + "-" + A5a['LOF'].astype(str) + ')',A5a['variant'])
A5a['variant'] = np.where(
    A5a['variant'].str.contains("rrdr"),
    A5a['variant'] + "_" + A5a['Subst'].str.replace(' .*', '', regex=True),
    A5a['variant']
)
A5a['variant']=np.where(A5a['Freq'] > 75, A5a['variant'],A5a['variant'] + ':' + A5a['Freq'].astype(str))

A5a=A5a.drop(['LOF', 'Type','Genome position','Subst'], axis=1)
A5b=A5a.groupby(['File','RIF_Conf_Grade','INH_Conf_Grade', 'EMB_Conf_Grade', 'PZA_Conf_Grade', 'LEV_Conf_Grade',       'MXF_Conf_Grade', 'BDQ_Conf_Grade', 'LZD_Conf_Grade', 'CFZ_Conf_Grade',       'DLM_Conf_Grade', 'AMI_Conf_Grade', 'STM_Conf_Grade', 'ETH_Conf_Grade',       'KAN_Conf_Grade', 'CAP_Conf_Grade'])['variant'].apply(', '.join).reset_index()


#A5b=A5.groupby(['File','RIF_Conf_Grade','INH_Conf_Grade', 'EMB_Conf_Grade', 'PZA_Conf_Grade', 'LEV_Conf_Grade',       'MXF_Conf_Grade', 'BDQ_Conf_Grade', 'LZD_Conf_Grade', 'CFZ_Conf_Grade',       'DLM_Conf_Grade', 'AMI_Conf_Grade', 'STM_Conf_Grade', 'ETH_Conf_Grade',       'KAN_Conf_Grade', 'CAP_Conf_Grade'])['variant'].apply(', '.join).reset_index()


A5b.insert(1, "ID", A5b.File.str.split('_').str[0], True)
#T5b['ID']=T5b.File.str.split('_').str[0]
#A5b['new'] = A5b['ID'].str.extract('(\\d+)').astype(int)
#https://stackoverflow.com/questions/66134896/python-pandas-sort-an-alphanumeric-dataframe
#A5b = A5b.sort_values(by=['new'], ascending=True).drop('new', axis=1)
A5b = A5b.sort_values(by=['ID'], ascending=True)#.drop('new', axis=1)

A5b.to_csv('WHO_raw.csv',index=False, sep=';')
os.chmod('WHO_raw.csv', 0o777)
"""
}

process OUT_WHO {
cpus 1
container 'library://allen13x/mtbseq/nf_mtbseq:1.0.1'
memory "20G"
publishDir "OUTPUT", mode: 'copy', pattern: 'res_WHO.csv'
publishDir "OUTPUT", mode: 'copy', pattern: 'res_who_long.csv'
input:
        path(WHO)
		path(head)
output:
        path('res_WHO.csv')
        path('res_who_long.csv')
script:
"""
#!/usr/bin/env Rscript
library(tidyverse)

read_delim('${WHO}',delim=';') %>%
  select(ID) %>% distinct() %>%
  left_join(
    (read_delim('${WHO}',delim=';') %>%
       select(ID, ends_with('Grade'),variant) %>%
       gather(drug,grade,ends_with('Grade')) %>%
       filter(grade!='0') %>%
       mutate(drug=str_remove_all(drug,'_.*')) %>%
       group_by(ID,drug,grade) %>%
       summarise(v=str_c(variant,collapse=', ')))) %>%
  mutate(grade=str_replace_all(grade, ' ','_')) %>%
  group_by(ID) %>% 
  complete(drug,grade) %>% 
  ungroup() %>% 
  unite('drug',drug,grade,sep='_') %>%
  spread(drug,v) %>%
  mutate(Interpretation_126=case_when((!if_all(matches('^RIF_1|^RIF_2|^RIF_6'),is.na)&
                                       !if_all(matches('^INH_1|^INH_2|^INH_6'),is.na)&
                                       !if_all(matches('^LEV_1|^LEV_2|^LEV_6|^MXF_1|^MXF_2|^MXF_6'),is.na)&
                                       !if_all(matches('^BDQ_1|^BDQ_2|^BDQ_6'),is.na))~'XDR',
                                    (!if_all(matches('^RIF_1|^RIF_2|^RIF_6'),is.na)&
                                       !if_all(matches('^INH_1|^INH_2|^INH_6'),is.na)&
                                       !if_all(matches('^LEV_1|^LEV_2|^LEV_6|^MXF_1|^MXF_2|^MXF_6'),is.na))~'Pre-XDR',
                                    (!if_all(matches('^RIF_1|^RIF_2|^RIF_6'),is.na)&
                                       !if_all(matches('^INH_1|^INH_2|^INH_6'),is.na))~'MDR',
                                    (!if_all(matches('^RIF_1|^RIF_2|^RIF_6'),is.na))~'RR',
                                    (!if_all(matches('^INH_1|^INH_2|^INH_6'),is.na))~'HR')) %>% 
  select(ID,starts_with('Interpretation'),matches('\\\\)')) ->temp
  bind_rows(read_delim('${head}',delim=';'),temp) %>%
  write_delim('res_WHO.csv',delim=';',na='')
  Sys.chmod("res_WHO.csv", mode = "0777")

read_delim('res_WHO.csv') %>% 
  gather(WHO_grade,Mut,-ID,-Interpretation_126) %>% 
  separate_rows(Mut,sep=', ') %>% 
  filter(!is.na(Mut)) %>% 
  write_excel_csv('res_who_long.csv',delim=';')
  Sys.chmod("res_who_long.csv", mode = "0664")
"""

}
