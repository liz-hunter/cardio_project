# Assembly and QC 

#### Initial QC - FastQC

```
#!/bin/bash
#SBATCH -J fastQC
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 6:00:00
##SBATCH --account=epscor-condo

module load fastqc
fastqc *.fastq
```

#### TrimmomaticPE

``` #!/bin/bash
#SBATCH -J fastQC
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 6:00:00
#SBATCH --account=epscor-condo

module load trimmomatic/0.36

TrimmomaticPE -threads 16 F1.fastq R1.fastq -baseout out1.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100 HEADCROP:5 
```

**Repeat FastQC to ensure quality and removal of adapters**

#### Spades
Used larger Kmers - more optimal for 150x150 chemistry Illumina runs

```
#!/bin/bash
#SBATCH -J spades
#SBATCH --account=epscor-condo
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH --mem=250G
#SBATCH -n 16

module load spades/3.13.0
export OMP_NUM_THREADS=16

/gpfs/runtime/opt/spades/3.13.0/bin/spades.py --memory 250 --threads 16 -k 55,77,99,127 \
--pe1-1 trimmed_F1.1.fq \
--pe1-2 trimmed_R1.fq \
--pe2-1 trimmed_F2.fq \
--pe2-2 trimmed_R2.fq \
-o /outdir/path/
```

Can pick up spades with the continue function if it runs out of memory or fails..
```
#!/bin/bash
#SBATCH -J spades
#SBATCH --account=epscor-condo
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH --mem=250G
#SBATCH -n 16

module load spades/3.13.0

export OMP_NUM_THREADS=16

/gpfs/runtime/opt/spades/3.13.0/bin/spades.py --continue \
-o /outdir/path
```

#### Check assembly stats with BBmap
```
#!/bin/bash
#SBATCH -J bbstat
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:00:00
##SBATCH --account=epscor-condo

module load bbmap/38.23
stats.sh in=contigs.fasta
```
