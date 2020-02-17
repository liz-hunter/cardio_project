# Bacterial Gene Trees

### Install all the things
```
conda create -c conda-forge -c bioconda -c astrobiomike gtotree
conda create -n trees
conda install -c bioconda ncbi-genome-download
conda install -c etetoolkit ete3 ete_toolchain
```
(Gotta put these in separate environments there are package conflicts)

### Make a list of the taxids you want

- Use this to pull the taxids for all descendent taxa
- Put them into a list
- Use: gimme_taxa.py (from user kblin on github)
```python gimme_taxa.py -o alpha.txt 28211```
  - will grab all taxa id that include alphaproteobacteria (28211) in their taxonomy chain

### Pull the complete genomes with that taxonomy from NCBI

```
ncbi-genome-download --dry-run --taxid alpha_taxid.txt --format fasta --assembly-level complete > alpha_test.txt
```
- use the --dry-run flag and wc -l to double check what you are actually getting, and that output will give you the accessions to use with GToTree

### Use GToTree

```
GToTree -a all_accession.txt -f fasta_files.txt -H Universal -t -L Species,Strain -m genome_to_id_map.tsv -j 4 -o GToTree_out
```
- can use the H flag to specify the hmm model to use (gtt-hmms to view all) or use the Universal
- fasta_files.txt is a list of the relative paths to the genomes you are providing (provide full assembly - not proteins or predicted CDS, of your genomes of interest and an outgroup - used E. coli here)
- accessions.txt is a list of genbank accession IDs (ex. GCF_000006965.1) you want to include, and it will automatically pull what it needs from NCBI
- genome_to_id_map.tsv is to provide new names for the sequences you are adding, just simple 2 column format
-  -t and -L will automatically and label everything it pulls from NCBI accession_Genus_species,strain in the finished tree
- -j is how many jobs to run in parallel


# Eukaryotic Apicoplast Trees

#### Runs the blast, pulls the results column, deduplicates it, and grabs the sequences from the master file (locally)
```
for f in $(cat list.txt); do blastp -query apicoplastgenes/"$f".faa -db 1258_db/1258_PROKKA_01052020.faa -outfmt 6 -evalue .1 -max_target_seqs 3 > "$f"_1258_blast.txt && cut -f2 "$f"_1258_blast.txt | sort -u > "$f"_1258_head.txt && fasomerecords 1258_PROKKA_01052020.faa "$f"_1258_head.txt done_1258/"$f"_1258_final.txt ; done
```

#### Concatenates the files of interest
```
for f in $(cat list.txt); do cat apicoplastgenes/"$f".faa done_1212/"$f"_1212_final.txt done_1258/"$f"_1258_final.txt done_1263/"$f"_1263_final.txt > concat/"$f"_cat_apico.fasta; done
```

#### Remove duplicate sequences (with different headers) if necessary and rename any duplicate headers (adjusted the e-value for blast and didn't end up needing to do this)
```
for f in $(cat list.txt); do sed 's/*//g' concat_final/"$f"_cat_apico.fasta > concat_final/"$f"_cat_apico_destar.fasta ; done
```
(remove the * in the fasta first because this doesn't play nice with cdhit apparently)
```
for f in $(cat list.txt); do cd-hit -i concat_final/"$f"_cat_apico_destar.fasta -o concat_final/cdhit/"$f"_cat_apico_destar_100.fasta -c 1.0 -t 1 -n 5 && \
perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' concat_final/cdhit/"$f"_cat_apico_destar_100.fasta > concat_final/dedup/"$f"_cat_apico_destar_100_dedup.fasta ; done
```

#### Move the final concatenated fastas up to a server (used Bluewaves)

### Run the alignments, trimming, and RAxML trees in parallel (on a server)

**mafft/7.215, trimAl/1.4.1-GCC-8.2.0-2.31.1, raxml/8.2.3**
```
FILE=$(head -n $SLURM_ARRAY_TASK_ID list.txt | tail -n 1)

mafft --auto concat_dedup/${FILE}_cat_apico_dedup.fasta > aligned/${FILE}_aligned.fasta && \
trimal -in aligned/${FILE}_aligned.fasta -out aligned_trimmed/${FILE}_aligned_trimmed.fasta -gt .7 -cons 60 && \
raxmlHPC-PTHREADS -m PROTGAMMALG -T 8 -f a -x 1034 -p 1034 -N 100 -s aligned_trimmed/${FILE}_aligned_trimmed.fasta -n ${FILE}.tre
```
- use an array to run this in parallel for every file
- will take a couple hours (longer if more trees/bigger proteins - 35 here, relatively small job)
- need to provide a list.txt file with all the variable names

### Parse these manually for paralogs
- remove paralogs (should have one copy of every gene in the final fastas)
- make sure all the headers are identical for each species (for every protein across all files)
- realign and trim new parsed fastas 
- use a less stringent trimal script (-gt .01)

### Concatenate alignments with catfasta2phyml program 
``` perl catfasta2phyml.pl -c -f -v *_individual_alignments.fasta > concat_alignment.fasta ```

### Run final tree (RAxML)

```raxmlHPC-PTHREADS -m PROTGAMMALG -T 8 -f a -x 1034 -p 1034 -N 100 -s concat_alignment2.fasta -n concat_alignment2.tre```
