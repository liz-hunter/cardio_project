# CAT/BAT

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
#SBATCH --mail-type=END
#SBATCH --mail-user=liz.hunter1122@gmail.com

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
