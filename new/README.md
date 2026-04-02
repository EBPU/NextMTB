# NF_TBSEQ

Presenting an enhanced version of the MTBseq pipeline, now implemented in Nextflow. Additional features include the integration of deletion/insertion detection using delly2, breadth coverage analysis for target genes with Mosdepth, and the incorporation of drug resistance pattern detection and interpretation through the utilization of the WHO catalogue as a point of reference (available at: https://www.who.int/publications/i/item/9789240028173).

This pipeline exploits the following software components:

MTBseq: https://github.com/ngs-fzb/MTBseq_source

Mosdepth: https://github.com/brentp/mosdepth

Delly2: https://github.com/dellytools/delly

## Requirements
  - Singularity
  - Nextflow

## Quick Usage
Put the .fastqs in a folder named "samp" and then use:
```
nextflow run https://github.com/Allen13x/NF_TBSEQ -latest -r main -with-singularity library://allen13x/mtbseq/nf_mtbseq:1.0.0 -resume --join false --extra false --reads $(pwd)/samp/
```
This will perform the basic single-sample-level characterizion from MTBseq pipeline

## Parameters

```
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
```
## Customization
To use on other species, it is necessary to rebuild the singularity image, adding the new reference genome to the /opt/conda/share/mtbseq-1.0.4-2/var/ref/ folder and rebuild the indexes with the commands:
```
samtools faidx /opt/conda/share/mtbseq-1.0.4-2/var/ref/NEWGENOME.fasta
bwa index /opt/conda/share/mtbseq-1.0.4-2/var/ref/NEWGENOME.fasta
picard CreateSequenceDictionary -R /opt/conda/share/mtbseq-1.0.4-2/var/ref/NEWGENOME.fasta
```

## Future Plans

-	Implentation for Hybrid reads
-	Graphic tree



MIT License

Copyright (c) [2023] [Federico Di Marco]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
