# Running Blast on a Server (without a taxonomy database)

### Split file into pieces and rename numerically
```awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%[chunks of how many sequences]==0){file=sprintf("%d_split.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; close(file)}' < [input.fasta]
```
```
n=1; for x in *.fa; do mv $x  $n.fa; n=$(($n+1)); done
```

### Run blast array on Bluewaves/Oscar
**blast/2.2.31+**
```
blastp -query "$SLURM_ARRAY_TASK_ID".fa -db /data3/shared/ncbi-nr/nr \
-outfmt '6 qseqid sseqid sacc evalue pident' \
-max_target_seqs 10 \
-evalue 1e-2 \
-num_threads 16
```

### Pull down outfile and assign taxonomy in R with taxonomizr 
```
#first time do all this shit 

#install.packages("taxonomizr")
#install.packages("readr")

#prepareDatabase('accessionTaxa.sql')

#getAccession2taxid(outDir = ".",
baseUrl = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/",
types = c("prot"))
#read.accession2taxid(taxaFiles = "prot.accession2taxid.gz", sqlFile = "prot_accessionTaxa.sql", vocal = TRUE,
indexTaxa = FALSE, overwrite = FALSE)
#need this second block for protein accessions, but still need to use the nucleotide files to convert taxaID to names
###########################

library(taxonomizr)
library(readr)

#load in blast accession results (single column with header ACCESSION)
neph_accession <- read_csv("apico_v3/accession.txt")
print(neph_accession)
#create vector
neph_vector <- neph_accession$ACCESSION
print(neph_vector)
#assign taxaId
taxaId<-accessionToTaxa(neph_vector,"prot_accessionTaxa.sql",version='base')
print(taxaId)

#assign taxonomy
phylum<-getTaxonomy(taxaId, sqlFile = "accessionTaxa.sql",
                      desiredTaxa = c("phylum"))
species<-getTaxonomy(taxaId, sqlFile = "accessionTaxa.sql",
                      desiredTaxa = c("species"))

#merge vectors
neph_accession$PHYLUM <- phylum
neph_accession$SPECIES <- species
neph_accession$TAX_ID <- taxaId
print(neph_accession)
#read seqids
seq_ids <- read_csv("apico_v3/seqid.txt")
seq_vector <- seq_ids$SEQID
print(seq_vector)
#merge vectors
neph_accession$SEQID <- seq_vector
print(neph_accession)

#print
write.table(neph_accession, "blast_tax.txt", sep='\t', quote = FALSE)
```
