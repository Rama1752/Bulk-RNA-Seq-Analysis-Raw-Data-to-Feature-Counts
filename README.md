# Bulk RNA-Seq Analysis: Raw Data to Feature Counts

🎯 Project Overview
This project demonstrates the initial stages of a Bulk RNA-Seq workflow, focusing on processing raw sequencing data through to feature count generation. We analyzed the effect of LCOR overexpression and how it affects gene transcription in two cell lines that differ in nuclear receptor status:
- **MCF7 (NR-positive)** — nuclear receptor–positive
- **MDA-MB-231 (NR-negative)** — nuclear receptor–negative  

LCOR acts both as a transcriptional corepressor and activator, and this dataset enables comparison of its regulatory activity in nuclear receptor–dependent versus independent contexts.

🧪 Sample Information
| Sample ID | Cell Line  | Condition | SRA Accession |
|-----------|------------|-----------|---------------|
| Sample 1  | MDA-MB-231 | LCOR OE   | SRR32858437   |
| Sample 2  | MDA-MB-231 | Wild Type | SRR32858438   |
| Sample 3  | MCF7       | LCOR OE   | SRR32858439   |
| Sample 4  | MCF7       | Wild Type | SRR32858440   |

**Note:** This dataset has no biological replicates. For robust differential expression analysis, replicates are essential.

📂 Dataset Information:
- **GEO Accession:** [GSE292767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292767)
- **Experiment Type:** Expression profiling by high-throughput sequencing
- **Description:** RNA-Seq experiment in breast cancer cell lines to check the effect of LCOR overexpression in MCF7 and MDA-MB-231.

## Quick Navigation

- [1. Downloading necessary tools and make directories](#1-downloading-necessary-tools-and-make-directories)
- [2. Fetching SRA Files](#2-fetching-sra-files)
- [3. FastQC - Quality Control](#3-fastqc---quality-control)
- [4. MultiQC](#4-multiqc)
- [5. Trimming Reads (optional)](#5-trimming-reads-optional)
- [6. Post-trimming Quality Control](#6-post-trimming-quality-control)
- [7. Reference Genome Preparation (Indexing)](#7-reference-genome-preparation-indexing)
- [8. Alignment/Mapping](#8-alignmentmapping)
- [9. BAM index file](#9-bam-index-file)
- [10. Assessing Alignment Quality](#10-assessing-alignment-quality)
- [11. Convert GTF to BED](#11-convert-gtf-to-bed)
- [12. Determine Library Strandedness](#12-determine-library-strandedness)
- [13. Feature Counting (Read Quantification)](#13-feature-counting-read-quantification)

---

## Workflow

---

## 1. Downloading necessary tools and make directories
- Download all the necessary tools and make directories

```bash

conda install -y -c bioconda -c conda-forge fastqc multiqc sra-tools hisat2 samtools trimmomatic subread qualimap rseqc bedops

mkdir -p SRA_files FASTQ_files FASTQC_reports Multiqc_reports reference aligned_reads quants rnaseq_qc_results

```
---

## 2. Fetching SRA Files
- Download SRA files using the SRA toolkit.
- Convert the SRA files into FASTQ files using fastq-dump.
- Save the file in .gz format i.e. gzip file so that it would take less space.
- **Run fastq-dump for each sample.**
```bash
#Download SRA files
prefetch SRR32858437 SRR32858438 SRR32858439 SRR32858440 --progress

#Convert SRA files to FASTQ files
fastq-dump --outdir FASTQ_files --gzip --skip-technical \
--readids --read-filter pass --dumpbase --split-3 --clip \
SRR32858437/SRR32858437.sra
```
---

## 3. FastQC - Quality Control
- Check the quality of raw sequencing reads using FastQC.

```bash

#Run FASTQC
fastqc FASTQ_files/*.fastq.gz -o FASTQC_reports/ --threads 8

```
---

## 4. MultiQC 
- Merge all the fastqc reports to get a summarised report of it.

```bash

multiqc FASTQC_reports/ -o Multiqc_reports

```

## 5. Trimming Reads (optional)
- Remove adapter contamination and poor-quality bases if present.

```bash
trimmomatic SE -threads 8 -phred33 \
  FASTQ_files/SRR32858437.fastq.gz \
  FASTQ_files/SRR32858437_trimmed.fastq.gz \
  TRAILING:10
```  
---

## 6. Post-trimming Quality Control
- Re-run FastQC after trimming to ensure cleaning steps were effective.
- Same as step-3

---

## 7. Reference Genome Preparation (Indexing)
- Download HISAT2 prebuilt GRCh38 genome index and ensemble gtf annotation.

```bash

wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz -C reference/

```
- Download Ensembl GTF annotation:

```bash

wget -P reference/ https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip reference/Homo_sapiens.GRCh38.115.gtf.gz 
```
---

## 8. Alignment/Mapping
 - Rename the files for better understanding.
 - Align reads with the Human Genome and convert the SAM file to BAM file.
 - **Repeat alignment step for each sample, changing the file names accordingly.**
   
```bash

#Renaming the files:
mv FASTQ_files/SRR32858437.fastq.gz FASTQ_files/MDA_MB_231_LCOR_OE.fastq.gz
mv FASTQ_files/SRR32858438.fastq.gz FASTQ_files/MDA_MB_231_WT.fastq.gz
mv FASTQ_files/SRR32858439.fastq.gz FASTQ_files/MCF7_LCOR_OE.fastq.gz
mv FASTQ_files/SRR32858440.fastq.gz FASTQ_files/MCF7_WT.fastq.gz

#Aligning the fastq files with genome
hisat2 -q -x reference/grch38/genome -U FASTQ_files/MDA_MB_231_LCOR_OE.fastq.gz | \
  samtools sort -o aligned_reads/MDA_MB_231_LCOR_OE.bam

```
---

## 9. BAM index file
- Creates .bai index file for fast random access.
- **Run for each sample.**
```bash

samtools index aligned_reads/MDA_MB_231_LCOR_OE.bam

```
---

## 10. Assessing Alignment Quality
- Generate reports on mapping performance and quality.
- **Run for each sample.**
```bash

#QC Check
qualimap rnaseq -bam aligned_reads/MDA_MB_231_LCOR_OE.bam -gtf reference/Homo_sapiens.GRCh38.115.gtf \
 -outdir rnaseq_qc_results/MDA_MB_231_LCOR_OE --java-mem-size=10G
 
```

---

## 11. Convert GTF to BED
- Create a .bed file from Human .gtf file which is required to check strandedness of RNA-Seq data using RSeQC.

```bash

gtf2bed < reference/Homo_sapiens.GRCh38.115.gtf > reference/Homo_sapiens.GRCh38.115.bed

```
---

## 12. Determine Library Strandedness
- RSeQC is used to check the strandedness of the RNA-Seq data.
- Checking strandedness is necessary cause it assigns the reads to the correct gene, when genes overlap on opposite strands.
  
```bash

infer_experiment.py -i aligned_reads/MDA_MB_231_LCOR_OE.bam \
  -r reference/Homo_sapiens.GRCh38.115.bed

```
### Interpreting Strandedness Results:
- **Fraction of reads explained by "++,--":** Matches forward strand
- **Fraction of reads explained by "+-,-+":** Matches reverse strand

Common interpretations:
- ~0.5/0.5: **Unstranded** library (use `-s 0` in featureCounts)
- ~0.01/0.99 or ~0.99/0.01: **Stranded** library (use `-s 1` or `-s 2`)

---

## 13. Feature Counting (Read Quantification)
- Count the reads mapping to genes/features.  

```bash

featureCounts -S 2 -a reference/Homo_sapiens.GRCh38.115.gtf \
  -o quants/featurecounts.txt -T 8 aligned_reads/*.bam

```

---

> *Click on any step in [Quick Navigation](#quick-navigation) to jump directly to that section!*
