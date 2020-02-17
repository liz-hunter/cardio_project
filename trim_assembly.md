# QC, Trimming, and Assembly

All run on the Brown University Oscar Server - EPSCoR Condo

### Initial QC - FastQC

**fastqc 0.11.8**
```
fastqc *.fastq
```

### TrimmomaticPE

**trimmomatic/0.36**
```
TrimmomaticPE -threads 16 F1.fastq R1.fastq -baseout out1.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100 HEADCROP:5 
```

**Repeat FastQC to ensure quality and removal of adapters**

### Spades

Used larger Kmers - more optimal for 150x150 chemistry Illumina runs

**spades/3.13.0**
```
spades.py --memory 250 --threads 16 -k 55,77,99,127 \
--pe1-1 trimmed_F1.1.fq \
--pe1-2 trimmed_R1.fq \
--pe2-1 trimmed_F2.fq \
--pe2-2 trimmed_R2.fq \
-o /outdir/path/
```
Needs a lot of memory (250G) and about 10-12 days 
Have to use export OMP_NUM_THREADS=16 to use multiple threads 
And also specify the path to spades.py (/gpfs/runtime/opt/spades/3.13.0/bin/ on Bluewaves)

Can pick up spades with the continue function if it runs out of memory or fails..
```
spades.py --continue -o /outdir/path
```

### Check assembly stats with BBmap

**bbmap 38.23**
```
stats.sh in=contigs.fasta
```
