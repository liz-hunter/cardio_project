# CAT

- Download CAT_pack in directory
- Check ``` CAT_pack/CAT --help ```
- Download the taxonomy database ``` wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20181212.tar.gz ```
  - (takes hours)
  - unpack ``` tar -xvzf CAT_prepare_20181212.tar.gz ```
  - get the most recent one from tbb.bio.uu.nl/bastiaan/CAT_prepare/

### Run CAT contigs

```
#!/bin/bash
#SBATCH -J catbat
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -n 16

module load python-igraph/0.7.0-foss-2016b-Python-3.5.2
module load prodigal/2.6.2
module load DIAMOND/0.9.23-foss-2016b

CAT-master/CAT_pack/CAT contigs -c <input fasta> -d CAT-master/CAT_prepare_20190108/2019-01-08_CAT_database -t CAT-master/CAT_prepare_20190108/2019-01-08_taxonomy
```

~8hr runtime, needs a lot of memory

### Run CAT add_names
```
#!/bin/bash
#SBATCH -J catbat
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -n 16

module load python-igraph/0.7.0-foss-2016b-Python-3.5.2
module load prodigal/2.6.2
module load DIAMOND/0.9.23-foss-2016b

CAT-master/CAT_pack/CAT add_names -i out.CAT.ORF2LCA.txt -o tax_named.txt -t CAT-master/CAT_prepare_20190108/2019-01-08_taxonomy
```

##### Sorted output with textwrangler
##### Pulled headers with fasomerecords

``` fasomerecords <in.fasta> <list.txt> <out.fasta> ```

# **Run Metabat**
##### Metagenomic Binning

### Make Bowtie database from Assembly

``` #!/bin/bash
#SBATCH -J Bow2Build
#SBATCH --account=epscor-condo
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --constraint=intel
#SBATCH --mail-type=END
#SBATCH --mail-user=liz.hunter1122@gmail.com

module load bowtie2/2.3.0

bowtie2-build --threads 16 -f /users/ehunter6/data/ehunter6/CardioDNA/UnbinAssembly/all_contigs.fasta all_contigs.fna ```

### Map Trimmed Reads to the Assembly Database
(the Sam files are the necessary output here)

```
#!/bin/bash
#SBATCH -J Bow_BP4
#SBATCH --account=epscor-condo
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --constraint=intel
#SBATCH --mail-type=END
#SBATCH --mail-user=liz.hunter1122@gmail.com

module load bowtie2/2.3.0

bowtie2 -p 16 -q --very-sensitive \
--al-conc /users/ehunter6/data/ehunter6/CardioDNA/UnbinAssembly/Meta/BP4_reads.fa \
-x /users/ehunter6/data/ehunter6/CardioDNA/UnbinAssembly/Meta/all_contigs.fna \
--no-unal -S /users/ehunter6/data/ehunter6/CardioDNA/UnbinAssembly/Meta/BP4.sam \
-1 /users/ehunter6/data/ehunter6/CardioDNA/Trim/BP4/BP4_trim_1P.fastq \
-2 /users/ehunter6/data/ehunter6/CardioDNA/Trim/BP4/BP4_trim_2P.fastq
```

### Sam to Bam

```
#!/bin/bash
#SBATCH -J samtools
#SBATCH --account=epscor-condo
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=END
#SBATCH --mail-user=liz.hunter1122@gmail.com

module load samtools/1.9

samtools view -S -b BP4.sam > BamFiles/BP4_paired.bam
```
### Bam Sort

```
#!/bin/bash
#SBATCH -J samtools
#SBATCH --account=epscor-condo
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=END
#SBATCH --mail-user=liz.hunter1122@gmail.com

module load samtools

samtools sort BamFiles/BP4_paired.bam -o Bam/BP4_paired_sorted.bam
```

### Metabat
(have to switch to Bluewaves for this)
```
#!/bin/bash
#SBATCH --time=100:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20

module load metabat/2.12.1-foss-2018a
module load checkm/1.0.5

echo "START"
date

runMetaBat.sh -m 1500 -o /data3/lanelab/liz/CardioDNA/Metabat/Bin1/bin1 /data3/lanelab/liz/CardioDNA/Metabat/contigs.fasta /data3/lanelab/liz/CardioDNA/Metabat/Bam/*.bam
checkm lineage_wf -t 20 -x fa /data3/lanelab/liz/CardioDNA/Metabat/Bin1 /data3/lanelab/liz/CardioDNA/Metabat/Bin1/Checkm/

echo "DONE"
date
```

```
#!/bin/bash
#SBATCH --time=100:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20

echo "START"
date

module load metabat/2.12.1-foss-2018a
module load checkm/1.0.5

metabat2 -i /data3/lanelab/liz/CardioDNA/Metabat/contigs.fasta -m 1500 -a /data3/lanelab/liz/CardioDNA/Metabat/contigs.fasta.depth.txt -o /data3/lanelab/liz/CardioDNA/Metabat/Bin1/bin1
checkm lineage_wf -t 20 -x fa /data3/lanelab/liz/CardioDNA/Metabat/Bin1/ /data3/lanelab/liz/CardioDNA/Metabat/Bin1/Checkm/

echo "DONE"
date
```

#### Can use Checkm functionality to check for SSUs
``` checkm ssu_finder -x .fa bin1.22.fa . ssu/ ```
