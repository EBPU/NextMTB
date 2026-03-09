nextflow.enable.dsl = 2
autoMounts = true
/*
 * Define the default parameters
 */ 
	params.reads	= "$baseDir/samp"
	params.results	= "OUTPUT"
	params.SEQ	= "ILL"
	params.minbqual	= "13"
	params.RP	= "0"
	params.minphred20	= "4"
	params.mincovf	= "4"
	params.mincovr	= "4"
	params.ref="M._tuberculosis_H37Rv_2015-11-13"
	params.reff = "${baseDir}/REF/${params.ref}.fasta"
	params.bed = "$baseDir/REF/h37rv_ups_ordered.bed.gz"
	params.bedix= "$baseDir/REF/h37rv_ups_ordered.bed.gz.tbi"
	params.kraken = true
	params.tgene="$baseDir/REF/target_genes.bed"
	params.pharma=false
	params.pgene="$baseDir/REF/gene_drug.csv"
	params.tdrug=10
	params.dhead="/$baseDir/REF/head"
	params.WHO="/$baseDir/REF/WHO_custom.csv"
	params.extra=true
	params.join=true
	params.sj="$baseDir/REF/placehold"
	params.proj='def'
	params.ascii="$baseDir/REF/ascii_string"
	params.h=false
	params.headWHO="$baseDir/REF/Header_WHO.csv"


/* 
 * Import modules 
 */

include{COLLECT_READS;
	COLLECT_READS_ONT;
	KRAKEN;
	BRACKEN;
	BRACKNOUT;
	KRAKEN_FILTER;
	MAPPING;
	MAPPING_ONT;
	REFINE;
	REFINE_ONT;
	PILE;
	PILE_ONT;
	LIST;
	VARIANTS_LOW;
	VARIANTS;
	STATS;
	STRAIN;
	MAP_STRAIN;
	JOIN;
	DEL;
	DEL_ONT;
	OUT_DEL;
	DEPTH;
	OUT_DEPTH;
	MUT_CORRECTION_DEL;
	MUT_CORRECTION;
	MUT_GATHER;
	PHARMA;
	WHO;
	OUT_WHO;
	FINAL_OUT} from "$baseDir/module.nf"
/* 
 * main pipeline logic
 */


workflow {
if (params.h){
log.info """
  --SEQ		    Sequencing technology:
				          ILL: Illumina [default]
				          ONT: Oxford Nanopore technology  
  --ref        Reference Genome to use:
                  M._abscessus_CIP-104536T_2014-02-03
                  M._chimaera_DSM44623_2016-01-28
                  M._fortuitum_CT6_2016-01-08
                  M._tuberculosis_H37Rv_2015-11-13 [default]
  --join       Perform joint analysis [default: true]
  --sj         List of samples to use for joint analysis [optional]
  --proj       Name of the project for joint analysis [default: def]
  --pharma     Perform Drug resistance analysis at Custom frequencies (to use after first round of analysis) [default: false]
  --tdrug      Custom frequencies % [default: 10]
  --WHO        Path to the WHO catalogue file formatted (Default = Catalogue v1)
  --extra      Perform extra analysis [default: true]:
                 Deletion/Insertion detection with Delly2
                 Sequencing depth with Mosdepth
                 Pharma analysis at 10% and Comparison with WHO catalogue              
"""
}else{
if (params.pharma){

PHARMA(Channel.fromPath('Called/*corrected.tab').collect(),params.tdrug,params.pgene)
MUT_GATHER(Channel.fromPath('Called/*corrected.tab').collect())
WHO(MUT_GATHER.out,params.dhead,params.WHO)
OUT_WHO(WHO.out,params.headWHO)

}
else{
log.info """\
================================
reads   	: $params.reads + '*_R{1,2}*.fastq.gz'
reference	: $params.ref
bed			: $params.bed
SEQ		: $params.SEQ
minbq		: $params.minbqual
REP REG		: $params.RP
minphred20	: $params.minphred20
mincovF		: $params.mincovf
mincovR		: $params.mincovr
results		: $params.results
Interesing genes: $params.bed
Target genes: $params.tgene
Drug genes: $params.pgene
Mutation Threshold: $params.tdrug
WHO Catalogue: $params.WHO
SampleList: $params.sj
"""

if(params.SEQ == "ILL"){

reads_ch=channel.fromFilePairs(params.reads + '*_R{1,2}*.fastq.gz').map{id,file ->tuple((id - ~/_.*/),file)}
//reads_ch.view()
COLLECT_READS(reads_ch,params.SEQ,params.minbqual,params.RP,params.minphred20)
collected=COLLECT_READS.out


KRAKEN(collected)
BRACKEN(KRAKEN.out.kreport)
brackenOUT=BRACKEN.out.breport
brackenOUT=brackenOUT.concat(channel.fromPath("bracken/*.report").map{file->tuple(file.getSimpleName(),file)}).unique{it[0]}
brackenOUTB=BRACKEN.out.bout
brackenOUTB=brackenOUTB.concat(channel.fromPath("bracken/*.bout").map{file->tuple(file.getSimpleName(),file)}).unique{it[0]}
MULTIQC(brackenOUT.map{id,file->file}.collect(sort:true),fastqcOUT.map{id,file->file}.collect(sort:true))
BRACKNOUT(brackenOUTB.map{id,file->file}.collect(sort:true))

if (params.kraken){
KRAKEN_FILTER(collected,KRAKEN.out.kraken)
collected=KRAKEN_FILTER.out
}

MAPPING(COLLECT_READS.out,params.ref)
mapped=MAPPING.out
REFINE(MAPPING.out.bam,params.ref)
refined=REFINE.out
PILE(REFINE.out.gatk,params.ref)
piled=PILE.out
}
else{
reads_ch=channel.fromPath(params.reads + '/*fastq.gz').map{file ->tuple((file.getSimpleName() - ~/_.*/),file)}
//reads_ch.view()
COLLECT_READS_ONT(reads_ch,params.SEQ,params.minbqual,params.RP,params.minphred20)
collected=COLLECT_READS_ONT.out
MAPPING_ONT(COLLECT_READS_ONT.out,params.ref)
mapped=MAPPING_ONT.out
REFINE_ONT(MAPPING_ONT.out.bam,params.ref,params.ascii)
refined=REFINE_ONT.out
PILE_ONT(REFINE_ONT.out.gatk,params.ref,params.minbqual)
piled=PILE_ONT.out}
old_mapped=channel.fromPath('Bam/*bam*').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}.groupTuple()
new_mapped=mapped.bam
mapped_bam=new_mapped.concat(old_mapped).unique{it[0]}
LIST(piled.mpile,params.minbqual,params.ref)
old_list=channel.fromPath('Position_Tables/*table.tab').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}.groupTuple()
new_list=LIST.out.list
ptables=new_list.concat(old_list).unique{it[0]}
VARIANTS_LOW(LIST.out.list,params.ref)
VARIANTS(LIST.out.list,params.mincovf,params.mincovr,params.minphred20,params.ref)
STATS(mapped_bam.join(ptables,by: 0),params.mincovf,params.mincovr,params.minphred20)
STRAIN(ptables)
map_strain=STATS.out.stats.join(STRAIN.out.strain,by:0).map{id,file1,file2 -> tuple(file1,file2)}.collect()
//old_map=channel.fromPath('OUTPUT/Mapping_Classification.tab')
map_strain=map_strain.collect()
MAP_STRAIN(map_strain)
DEPTH(mapped_bam,params.tgene)
depth=DEPTH.out.map{id,file->file}
//old_cov=channel.fromPath('OUTPUT/GB_cov.*')
depth=depth.collect()
OUT_DEPTH(depth)
var=VARIANTS_LOW.out.var_low
old_var=Channel.fromPath('Called/*variants_cf1*001.tab').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}.groupTuple()
var=var.concat(old_var).unique{it[0]}

if (params.ref == "M._tuberculosis_H37Rv_2015-11-13"){

FINAL_OUT(OUT_DEPTH.out,MAP_STRAIN.out)

}

if (params.extra){
if (params.SEQ == "ILL"){	
DEL(mapped_bam,params.ref,params.bed,params.bedix)
deletion=DEL.out}
else{
DEL_ONT(mapped_bam,params.ref,params.bed,params.bedix)
deletion=DEL_ONT.out}
//var.view()
delly=deletion.map{id,file -> file}
//old_del=channel.fromPath('OUTPUT/DELETIONS.*')
delly=delly.collect()
OUT_DEL(delly)
var_del=var.join(deletion,by:0)
//var_del.view()
MUT_CORRECTION_DEL(var_del)
mut=MUT_CORRECTION_DEL.out
old_mut=Channel.fromPath('Called/*corrected.tab').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}
mut=mut.concat(old_mut).unique{it[0]}.map{id,file->file}.collect()
MUT_GATHER(mut)
PHARMA(mut,"10",params.pgene)
WHO(MUT_GATHER.out,params.dhead,params.WHO)
OUT_WHO(WHO.out,params.headWHO)
}
else{
MUT_CORRECTION(var)
mut=MUT_CORRECTION.out
old_mut=Channel.fromPath('Called/*corrected.tab').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}
mut=mut.concat(old_mut).unique{it[0]}.map{id,file->file}.collect()
//mut.view()
MUT_GATHER(mut)
PHARMA(mut,"10",params.pgene)
WHO(MUT_GATHER.out,params.dhead,params.WHO)
OUT_WHO(WHO.out,params.headWHO)
}

if (params.join){
call=VARIANTS.out.var
old_call=Channel.fromPath('Called/*variants_cf4*').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}
call=call.concat(old_call).unique{it[0]}.map{id,file->file}.collect()
list=LIST.out.list
old_list=Channel.fromPath('Position_Tables/*').map{file -> tuple ((file.getSimpleName())- ~/_.*/,file)}
list=list.concat(old_list).unique{it[0]}.map{id,file->file}.collect()
JOIN(call,list,channel.fromPath(params.sj,checkIfExists:true).collect(),params.minbqual,params.minphred20,params.proj,params.ref)
}


}
}}
