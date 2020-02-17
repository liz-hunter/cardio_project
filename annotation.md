# Annotation of C. cionae

## Run Augustus

## MAKER

#### Install (local)

```
conda create --name maker
conda activate maker
conda install -c bioconda maker
conda install -c bioconda blast-legacy
maker -h #check installation
maker -CTL #generate option files
export REPEATMASKER_LIB_DIR=/Users/AlgaeLab/miniconda2/pkgs/repeatmasker-4.0.9_p2-pl526_0/share/RepeatMasker/Libraries

# needed to install legacy blast to get blastall
# skipped genemark and WUBLAST - optional
# exported the repeatmasker library to overcome ERROR: Could not determine if RepBase is installed, not sure if this solved the issue as intended but RepBase now costs $1300 and it did mask without it
```
 #### Setting up the option files

 - maker -CTL will generate the three option files
  - maker_exe.log should have all the paths to the executable programs
  - maker_opts.log
    - genome = genome.fasta
    - est = CDS (nucleotide) from transcriptome
    - protein = protein homology evidence from a library of references ex. plasmodiuma.fasta,plasmodiumb.fasta,plasmodiumc.fasta
    - augustus_species = species parameter file from augustus (must put in the species folder of the augustus installation - /Users/AlgaeLab/miniconda2/envs/maker/config/species)
  - maker_bopts.log contains blast options - left defaults

#### Run MAKER
```
conda activate maker
export REPEATMASKER_LIB_DIR=/Users/AlgaeLab/miniconda2/pkgs/repeatmasker-4.0.9_p2-pl526_0/share/RepeatMasker/Libraries
maker --base cardio_maker_v1 #base name for all the output files
#completed in 48hrs for 1700 contigs on AlgaeLab (6.6Mb)
```

#### Output
```
grep -c "FAILED" master_datastore_index.log
#confirm all contigs ran successfully

fasta_merge -d master_datastore_index.log
#merge all the many output files for every contig
```

#### Reannotation with SNAP
```
mkdir SNAP
gff3_merge -d master_datastore_index.log
maker2zff cardio_maker_v1.all.gff
# merge gff3s to a single file and convert to zff
# add params -x 0.25 -l 50 to maker2zff pull only "confident" gene models (AED <.25 and length >50)

fathom genome.ann genome.dna -gene-stats
fathom genome.ann genome.dna -validate > validate.txt
# remove genuine errors from these files (textwrangler)

fathom genome.ann genome.dna -categorize 1000
# break into pieces
fathom uni.ann uni.dna -export 1000 -plus
# convert uni to plus stranded

mkdir params
cd params
forge ../export.ann ../export.dna
cd ..
hmm-assembler.pl my-genome params > my-genome.hmm
# forge makes a lot of files, put params in their own directory
# don't forget to move up a directory before running the hmm-assembler
# generates a nice hmm model for your genome to retrain MAKER with

```
- this process directly from the SNAP user manual

#### Generate AED score cumulative fraction values
- use AED_cdf_generator.pl
- the b value sets the interval and .025 is the recommended default

```
perl AED_cdf_generator.pl -b 0.025 maker.gff
```
- The lower the better, AED values represent how well a prediction matches the model (with 0 being a perfect match)
- Below .5 is ideal for annotation done with a transcriptome reference
- can also filter for AEDs below a certain threshold with quality_filter.pl

```
perl quality_filter.pl -options transcripts.gff
```
- options for quality_filter:
  - -d Prints transcripts with an AED <1 (MAKER default)
  - -s Prints transcripts with an AED <1 and/or Pfam domain if in gff3 (MAKER Standard)
  - -a <number between="" 0="" and="" 1=""> Prints transcripts with an AED < the given value

#### Rerun MAKER
- can reuse the CTL files from the the first run but disable
- edit hmm into the maker_opts file (snaphmm=my-genome.hmm)
- also add in the maker_gff file from the first pass
- relaunch with a new --base name
- redo AED assessment and trim genes that are over .5

#### Pull rRNAs
- kmer binning with bbduk using a database of eukaryotic ribosomes (Oscar)
```
module load bbmap/38.23
bbduk.sh in=sequences.fasta \
outm=ribo.fasta \
outu=nonribo.fasta \
k=31 ref=ribokmers.fa.gz
```
- break the output into smaller files (needs to be <1MB per file) and rename
```
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%[chunks of how many sequences]==0){file=sprintf("%d_split.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; close(file)}' < [input.fasta]
```
```
n=1; for x in *.fa; do mv $x  $n.fa; n=$(($n+1)); done
```
- run through the RNammer web server (this entire program is an antique and impossible to install locally without legacy everything and perl from '96)

#### Pull tRNAs (if necessary to do separately)
```
conda activate maker
tRNAscan-SE -o out_tRNA.txt -m summary.txt sequences.fasta
```

#### maker_opts.log
```
#-----Genome (these are always required)
genome=RNA_CAT_combo_v2_1kb.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff=/Users/AlgaeLab/Desktop/LIZ_sucks/Maker/Cardio/RNA_CAT_combo/cardio_combo_maker_out_run2/cardio_combo_maker_out_run2.all.gff #MAKER derived GFF3 file
est_pass=1 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=1 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=1 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=1 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=1 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=1 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=1 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=/Users/AlgaeLab/Desktop/LIZ_sucks/Maker/Cardio/data/Cardio_clctrim_nulc_renamed.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=../RNA_CAT_combo/cardio_combo_maker_out_run1/SNAP/cardio_combo.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=Cardiosporidium_cionae_combo #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=4 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```
