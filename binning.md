CAT/Metabat run on the URI-Bluewaves server, all mapping on Brown Oscar (EPSCoR Condo)

# **CAT** (used for Eukaryotic Binning)

- Download CAT_pack in directory
- Check ``` CAT_pack/CAT --help ```
- Download the taxonomy database ``` wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20190719.tar.gz ```
  - (takes hours)
  - unpack ``` tar -xvzf CAT_prepare_20190719.tar.gz ```
  - get the most recent one from tbb.bio.uu.nl/bastiaan/CAT_prepare/

### Run CAT contigs

**Python-3.5.2, prodigal/2.6.2, DIAMOND/0.9.23**
```
CAT-master/CAT_pack/CAT contigs -c <input fasta> -d CAT-master/CAT_prepare_20190719/2019-07-19_CAT_database -t CAT-master/CAT_prepare_20190719/2019-07-19_taxonomy
```
~8-10hr runtime, needs a lot of memory (~100G)

### Run CAT add_names
```
CAT-master/CAT_pack/CAT add_names -i out.CAT.ORF2LCA.txt -o tax_named.txt -t CAT-master/CAT_prepare_20190719/2019-07-19_taxonomy
```
(super fast) 

##### Sorted output with textwrangler
##### Pulled headers with fasomerecords

``` fasomerecords <in.fasta> <list.txt> <out.fasta> ```

# RNA-Coverage Binning (also used for Eukaryotic Binning)

### Make bowtie2 databases for the transcriptome (binned) and genome (unbinned) assemblies 

**bowtie2/2.3.5.1**
```
bowtie2-build --threads 16 -f assembly.fasta database.out
```

### Map trimmed RNA reads to the transcriptome
```
bowtie2 -p 16 -q --very-sensitive \
--al-conc RNA_mapped_reads \
-x transcriptome_db \
-1 RNA_1F.fq,RNA_2F.fq,RNA_3F.fq \
-2 RNA_1R.fq,RNA_2R.fq,RNA_3R.fq
```
### Map the read output from this (RNA_mapped_reads) to the full, unbinned genome assembly 

### Create a genome coverage file with bedtools

**bedtools/2.26.0**
```
bedtools genomecov -bg -ibam DNA_mapped_sorted.bam -g assembly.fasta > out.bed
```

### Sort by contigs with high coverage from the binned transcriptome reads in R

**R version 1.2.5019**
```
library('dplyr')

out <- read.delim("out_2.bed", stringsAsFactors=FALSE)
#read in bedfile from genomecov -bg with added headers (SEQID, START, END, COV)
sorted <- out %>% group_by(SEQID) %>% filter(COV == max(COV))
filter <- select(sorted,-c(START,END))
filter_byCOV <- filter[order(filter$COV),]
outlier_removed <- filter_byCOV[-3034,] #3034 is the column number to remove
#pulls out the highest coverage per contig
#optionally create a file without the unnecessary columns and sorts it by coverage

mean(sorted$COV)
median.default(sorted$COV)
#averages from the sorted dataset

write.table(filter_byCOV, "cov.txt", sep='\t', quote = FALSE)
#cumulative fraction chart
plot(ecdf(outlier$COV))

cardio_CAT <- read.csv("genome_cov/cardio_CAT_v2_head.txt", sep="", stringsAsFactors=FALSE)
View(cardio_CAT)
super <- merge(sorted, cardio_CAT, by=c('SEQID') )
#read in and compare to headers from the CAT sorted bin

sorted_cutoff <- sorted[sorted$COV>=50,]
#pull contigs by a cutoff value (in this case 50)
write.table(sorted_cutoff, "cov.txt", sep='\t', quote = FALSE)
```


# **Metabat** (used for Prokaryotic Binning)

### Make Bowtie database from Assembly

**bowtie2/2.3.0**
``` 
bowtie2-build --threads 16 -f all_contigs.fasta all_contigs_db ```
```

### Map Trimmed Reads to the Assembly Database to get a sam file

```
bowtie2 -p 16 -q --very-sensitive \
--al-conc aligned_reads_out.fasta \
-x full_assembly.fasta \
--no-unal -S out.sam \
-1 trimmed_1F.fastq,trimmed_2F.fastq,trimmed_3F.fastq \
-2 trimmed_1R.fastq,trimmed_2R.fastq,trimmed_3R.fastq \
```

### Sam to Bam + Bam Sort

**samtools/1.9**
```
samtools view -S -b out.sam > out.bam && samtools sort out.bam -o out_sorted.bam
```

### Metabat

**metabat/2.12.1, checkm/1.0.5**
```
runMetaBat.sh -m 1500 contigs.fasta Bam/*.bam && \
#makes the depth file

metabat2 -i contigs.fasta -m 1500 -a contigs.fasta.depth.txt -o /Bins/bin1
#makes bins

checkm lineage_wf -t 20 -x fa /Bins/ Bins/Checkm/
#gives lineage guess for the bins
```

#### Can use Checkm functionality to check for SSUs
``` checkm ssu_finder -x .fa bin1.22.fa . ssu/ ```
