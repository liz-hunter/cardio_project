# Bacterial Expression Pipeline
**Didn't end up using this analysis - coverage wasn't high enough to provide meaningful output, the most highly expressed genes were rRNAs with both STAR and HISAT2**

## STAR

#### Creating the Genome Index

**STAR/2.5.3a, gcc/6.1.0**
```
STAR --runThreadN 16 \
#uses about 6G mem for every thread
--runMode genomeGenerate \
--genomeDir star_index \
#must make this directory before running
--sjdbGTFfeatureExon CDS \
#if your gtf file doesn't use "exon" you have to specify
--genomeFastaFiles Genome/genome.fa \
--sjdbGTFfile Genome/genome.gtf \
--genomeSAindexNbases 9 \
#MUST use this for small genomes or the alignment will fail, remove for normal sized stuff (for 1MB = 9, 100kB = 7, formula in manual)
--sjdbOverhang 100
#100 is default should be fine for almost any assembly but is readlength-1
```

- can rename contigs by editing the chrName.txt file (preserve order)
- other than that do not mess with these files
- used 9 for alphas and beta, 8 for bac (genomeSAindexNbases)
- reordered the both sets of alpha contigs (now named numerically by size - largest to smallest)

#### Mapping

```
STAR --runThreadN 16 \
--genomeDir star_index \
#genome index directory
--readFilesIn reads1F,reads2F reads1R,reads2R \
--outFileNamePrefix results/prefix \
#must make the results directory
--outSAMtype BAM Unsorted SortedByCoordinate\
#creates a sorted and unsorted BAM
--quantMode TranscriptomeSAM GeneCounts
#generates a gene count table as well as a BAM aligned to transcriptome
--twopassMode Basic
#updates the index as it goes through the first round of mapping then does it again, slow but more sensitive
```
- this is kinda slow for big datasets, especially with the twopassMode enabled (multiple days)
- add --alignIntronMax 1 for bacteria (no splices)
- want column 2 from the ReadsPerGene.out.tab file, columns 3 and 4 are + and - strand specific counts

#### Convert to Usable Format with R

##### Merge the KASS and Prokka annotations with the ordered counts

```
KASS <- read.table("KASS.txt",
                               header = TRUE)
counts <- read.table("cardio_collapsed.txt",
                               header = TRUE)
prokka <- read.table("cardio_prokka.txt",
                               header = TRUE)
merged <- merge(counts,KASS,by="GENE_ID")
merged
merged_all <- merge(merged_beta,prokka,by="GENE_ID")
write.table(merged_all, "merged_all.txt", sep='\t')
```

##### Convert to DeSeq2 format

```
ff <- list.files( path = "./counts", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
counts <- as.data.frame( sapply( counts.files, function(x) x[ , number ] ) )
ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/counts/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1
```

## HISAT2

#### Index the genomes

**HISAT2/2.0.4**
```
hisat2-build genome.fasta index/index.fna
```

##### Mapping
```
hisat2 -p 16 --un-conc results/aligned_hisat -x index.fna \
-1 forward.fastq \
-2 reverse.fastq \
-S results/hisat.sam
```

##### Sort and stuff
```
#!/bin/bash
#SBATCH -J samtools
#SBATCH -t 05:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=END
#SBATCH --mail-user=lizhunter@uri.edu

module load SAMtools/1.9-foss-2018b

samtools view -S -b hisat.sam > hisat.bam | samtools sort -n hisat.bam -o hisat_sorted.bam
```

#### featureCount

**Subread/1.6.3**
```
featureCounts -p -T [threads] -t transcript -g gene_id -a prokka.gtf -o bac_count.txt bac.bam
```
- can take bam or sam sorted or unsorted
- make sure the gtf file contig names match the bam/sam file
- works with gffread converter (gtf format)

### Normalize
- length / 1000 (kb)n = length_kb
- raw counts/length_kb = RPK
- sum(RPK) / 1,000,000 = scaling factor
- RPK / scaling factor = TPM (transcripts per million)
