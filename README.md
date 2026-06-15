# NextMTB

A Nextflow implementation of the MTBseq pipeline for whole-genome sequencing analysis of *Mycobacterium tuberculosis* and related species. The pipeline supports both Illumina paired-end and Oxford Nanopore Technology (ONT) reads, and extends the core MTBseq workflow with:

<img width="8000" height="4500" alt="MainFigure" src="https://github.com/user-attachments/assets/bd89945e-2f36-45ba-b873-bdd0e2c0e252" />

- **Taxonomic classification and read filtering** using Kraken2 and Bracken (optional)
- **Structural variant (deletion/insertion) detection** using Delly2
- **Target gene coverage analysis** using Mosdepth
- **Drug resistance prediction** with pharmacological analysis and WHO catalogue interpretation

**Tools used:**

| Tool | Purpose |
|------|---------|
| [MTBseq v1.0.4](https://github.com/ngs-fzb/MTBseq_source) | Core mapping, variant calling, strain typing, and joint analysis |
| [Kraken2](https://github.com/DerrickWood/kraken2) + [Bracken](https://github.com/jenniferlu717/Bracken) | Taxonomic classification and mycobacterial read filtering |
| [BWA](https://github.com/lh3/bwa) + [Samtools](https://github.com/samtools/samtools) | Alignment for ONT reads |
| [GATK / Picard](https://gatk.broadinstitute.org/) | BAM refinement |
| [Delly2](https://github.com/dellytools/delly) | Structural variant detection |
| [Mosdepth](https://github.com/brentp/mosdepth) | Depth of coverage for target genes |

## Requirements

- [Singularity](https://sylabs.io/singularity/)
- [Nextflow](https://www.nextflow.io/)

## Pipeline Overview

```
FASTQ input (Illumina or ONT)
        │
        ▼
┌───────────────────────────────────────────────────┐
│  1. Read preparation (COLLECT_READS)               │
│     Rename reads to MTBseq naming convention       │
└───────────────────────┬───────────────────────────┘
                        │
                        ▼  [if --kraken true (default)]
┌───────────────────────────────────────────────────┐
│  2. Taxonomic filtering (Kraken2 + Bracken)        │
│     Classify reads; retain Mycobacteria (taxid     │
│     1762) and report per-sample statistics         │
└───────────────────────┬───────────────────────────┘
                        │
                        ▼
┌───────────────────────────────────────────────────┐
│  3. Mapping (MTBseq TBbwa / BWA for ONT)           │
│     Align reads to the reference genome            │
└───────────────────────┬───────────────────────────┘
                        │
                        ▼
┌───────────────────────────────────────────────────┐
│  4. BAM refinement (MTBseq TBrefine / GATK)        │
│     Indel realignment and base quality             │
│     score recalibration                            │
└───────────────────────┬───────────────────────────┘
                        │
                        ▼
┌───────────────────────────────────────────────────┐
│  5. Pileup and variant calling (MTBseq)            │
│     TBpile → TBlist → TBvariants                   │
│     Position tables + variant call files           │
└───────────────────────┬───────────────────────────┘
                        │
                        ▼
┌───────────────────────────────────────────────────┐
│  6. Statistics & strain typing                     │
│     Mapping/variant statistics (MTBseq TBstats)    │
│     Strain classification (MTBseq TBstrains)       │
│     Combined Mapping_Classification.tab output     │
└───────────────────────┬───────────────────────────┘
                        │
          ┌─────────────┴──────────────┐
          ▼                            ▼
  [if --extra true]           [if --join true]
┌──────────────────────┐   ┌──────────────────────────┐
│  7. Extra analysis   │   │  8. Joint analysis        │
│                      │   │     (MTBseq TBjoin)       │
│  - Delly2: INDELs    │   │     Multi-sample SNP      │
│  - Mosdepth: target  │   │     matrix and phylo-     │
│    gene coverage     │   │     genetic grouping      │
│  - Drug resistance   │   └──────────────────────────┘
│    (PHARMA + WHO     │
│     catalogue)       │
└──────────────────────┘
```

### Output directories

| Directory | Contents |
|-----------|---------|
| `Bam/` | Per-sample BAM files |
| `GATK_Bam/` | Refined BAM files |
| `Position_Tables/` | Per-sample position tables |
| `Called/` | Per-sample variant call files |
| `Classification/` | Per-sample strain classification |
| `Joint/` | Joint-analysis SNP matrices and alignments |
| `Amend/` | Amended position tables from joint analysis |
| `Groups/` | Sample groupings from joint analysis |
| `OUTPUT/` | Summary files: `Mapping_Classification.tab`, coverage (`GB_cov.*`), deletions (`DELETIONS.*`), drug resistance, WHO comparison |
| `Kraken_Stats/` | Per-sample Mycobacteria read counts (`*_MycoReads.csv`, `Kraken_reads_summary.csv`) |
| `bracken/` | Bracken abundance reports |

## Quick Usage

### Basic single-sample analysis (Illumina)

Place paired-end FASTQ files (matching `*_R1*.fastq.gz` / `*_R2*.fastq.gz`) in a folder named `samp`, then run:

```bash
nextflow run https://github.com/EBPU/NextMTB \
  -latest -r main \
  -with-singularity library://allen13x/mtbseq/nf_mtbseq:1.0.2 \
  -resume \
  --join false \
  --extra false \
  --reads $(pwd)/samp/
```

### Full analysis with extra modules and joint calling

```bash
nextflow run https://github.com/EBPU/NextMTB \
  -latest -r main \
  -with-singularity library://allen13x/mtbseq/nf_mtbseq:1.0.2 \
  -resume \
  --reads $(pwd)/samp/ \
  --extra true \
  --join true \
  --proj my_project
```

### ONT reads

```bash
nextflow run https://github.com/EBPU/NextMTB \
  -latest -r main \
  -with-singularity library://allen13x/mtbseq/nf_mtbseq:1.0.2 \
  -resume \
  --SEQ ONT \
  --reads $(pwd)/samp/ \
  --extra true \
  --join false
```

### Re-run drug resistance analysis at a custom allele frequency

After a completed run, use the corrected mutation tables already present in `Called/`:

```bash
nextflow run https://github.com/EBPU/NextMTB \
  -latest -r main \
  -with-singularity library://allen13x/mtbseq/nf_mtbseq:1.0.0 \
  -resume \
  --pharma true \
  --tdrug 5
```

## Parameters

```
  --reads      Path to the directory containing input FASTQ files
                 [default: <baseDir>/samp]
  --results    Name of the top-level output directory
                 [default: OUTPUT]
  --SEQ        Sequencing technology:
                 ILL: Illumina paired-end [default]
                 ONT: Oxford Nanopore Technology
  --ref        Reference genome to use:
                 M._tuberculosis_H37Rv_2015-11-13 [default]
                 M._abscessus_CIP-104536T_2014-02-03
                 M._chimaera_DSM44623_2016-01-28
                 M._fortuitum_CT6_2016-01-08
  --kraken     Run Kraken2/Bracken taxonomic filtering [default: true]
  --krakendb   Path to the Kraken2 database (required when --kraken true)
  --join       Perform joint multi-sample analysis [default: true]
  --sj         Path to a file listing samples to include in joint analysis
                 (one sample per line; omit to use all samples) [optional]
  --proj       Project name used by MTBseq TBjoin [default: def]
  --extra      Perform extra analyses [default: true]:
                 - Deletion/insertion detection with Delly2
                 - Target gene coverage with Mosdepth
                 - Drug resistance analysis at 10% frequency
                 - Comparison with WHO catalogue
  --pharma     Re-run drug resistance analysis at a custom allele frequency
                 using previously called mutations in Called/ [default: false]
  --tdrug      Allele frequency threshold (%) for --pharma mode [default: 10]
  --WHO        Path to a formatted WHO catalogue CSV file
                 [default: <baseDir>/REF/WHO_custom.csv]
  --tgene      Path to a BED file of target genes for Mosdepth coverage
                 [default: <baseDir>/REF/target_genes.bed]
  --pgene      Path to a gene-to-drug mapping CSV file
                 [default: <baseDir>/REF/gene_drug.csv]
  --minbqual   Minimum base quality score for MTBseq [default: 13]
  --minphred20 Minimum Phred 20 base count per position [default: 4]
  --mincovf    Minimum forward strand coverage for variant calling [default: 4]
  --mincovr    Minimum reverse strand coverage for variant calling [default: 4]
  --RP         Repetitive region masking flag for MTBseq [default: 0]
  -h / --h     Print this help message and exit
```

## Resuming and Incremental Runs

The pipeline is designed to be resumed or extended with new samples. On each run, Nextflow automatically discovers previously generated files in `Bam/`, `Position_Tables/`, `Called/`, and similar directories and merges them with new results before running downstream steps. Use the `-resume` flag to skip already-completed tasks.

## Customization

### Adding a new reference genome

To analyse a species not included in the pre-built Singularity image, rebuild the image after adding the reference to the MTBseq reference directory and regenerating the required indexes:

```bash
samtools faidx /opt/conda/share/mtbseq-1.0.4-2/var/ref/NEWGENOME.fasta
bwa index /opt/conda/share/mtbseq-1.0.4-2/var/ref/NEWGENOME.fasta
picard CreateSequenceDictionary \
  -R /opt/conda/share/mtbseq-1.0.4-2/var/ref/NEWGENOME.fasta \
  -O /opt/conda/share/mtbseq-1.0.4-2/var/ref/NEWGENOME.dict
```

Then pass `--ref NEWGENOME` (without the `.fasta` extension) when running the pipeline.

## Future Plans

- Support for hybrid (Illumina + ONT) reads
- Graphical phylogenetic tree output

---

MIT License

Copyright (c) 2023 Federico Di Marco

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
