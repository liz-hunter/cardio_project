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

TrimmomaticPE -threads 16 /users/ehunter6/data/ehunter6/CardioDNA/DNA/XRICL_20190308_K00134_IL100122607_ATGTCA-BP4-27_L008_R1.fastq /users/ehunter6/data/ehunter6/CardioDNA/DNA/XRICL_20190308_K00134_IL100122607_ATGTCA-BP4-27_L008_R2.fastq -baseout BP4_trim.fastq ILLUMINACLIP:/users/ehunter6/data/ehunter6/scripts/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100 HEADCROP:5 ```

**Repeat FastQC**

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
#SBATCH --mail-type=END
#SBATCH --mail-user=liz.hunter1122@gmail.com

module load spades/3.13.0
export OMP_NUM_THREADS=16

/gpfs/runtime/opt/spades/3.13.0/bin/spades.py --memory 250 --threads 16 -k 55,77,99,127 \
--pe1-1 /users/ehunter6/data/ehunter6/CardioDNA/Bowtie2/Bin/BP4bin_reads.1.fq \
--pe1-2 /users/ehunter6/data/ehunter6/CardioDNA/Bowtie2/Bin/BP4bin_reads.2.fq \
--pe2-1 /users/ehunter6/data/ehunter6/CardioDNA/Bowtie2/Bin/BP6bin_reads.1.fq \
--pe2-2 /users/ehunter6/data/ehunter6/CardioDNA/Bowtie2/Bin/BP6bin_reads.2.fq \
--pe3-1 /users/ehunter6/data/ehunter6/CardioDNA/Bowtie2/Bin/BP9bin_reads.1.fq \
--pe3-2 /users/ehunter6/data/ehunter6/CardioDNA/Bowtie2/Bin/BP9bin_reads.2.fq \
--pe4-1 /users/ehunter6/data/ehunter6/CardioDNA/Bowtie2/Bin/BP11bin_reads.1.fq \
--pe4-2 /users/ehunter6/data/ehunter6/CardioDNA/Bowtie2/Bin/BP11bin_reads.2.fq \
-o /users/ehunter6/data/ehunter6/CardioDNA/AllAssembly2
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
#SBATCH --mail-type=END
#SBATCH --mail-user=liz.hunter1122@gmail.com

module load spades/3.13.0

export OMP_NUM_THREADS=16

/gpfs/runtime/opt/spades/3.13.0/bin/spades.py --continue \
-o /users/ehunter6/data/ehunter6/CardioDNA/AllAssembly2
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
