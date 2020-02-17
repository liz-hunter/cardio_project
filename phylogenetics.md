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
- Use: python gimme_taxa.py -o alpha.txt 28211
  - will grab all taxa id that include alphaproteobacteria in their taxonomy chain

```
#!/usr/bin/env python
from __future__ import print_function

__author__ = "Joe R. J. Healey; Nick Youngblut"
__version__ = "1.1"
__title__ = "gimme_taxa"
__license__ = "Apache2.0"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"


# Check script compatibilities and module requirements
import sys
import argparse
try:
    from ete3 import NCBITaxa
except ImportError:
    msg = """
The ete3 import failed, the module doesn't appear to be installed
(at least in the PYTHONPATH for this python binary").
Try running:
 $ python -m pip install ete3
or
 $ conda install -c etetoolkit ete3 ete_toolchain
"""
    print(msg)
    sys.exit(1)

def get_args():
    """Parse command line arguments
    """
    desc = 'Perform queries against the NCBI Taxa database'
    epi = """DESCRIPTION:
    This script lets you find out what TaxIDs to pass to ngd, and will write
    a simple one-item-per-line file to pass in to it. It utilises the ete3
    toolkit, so refer to their site to install the dependency if it's not
    already satisfied.
    You can query the database using a particular TaxID, or a scientific name.
    The primary function of the script is to return all the child taxa of the
    specified parent taxa. If specified with -v verbose flags however, the
    script will also print out some information about the lineages etc.

    WARNING: This script is still somewhat experimental
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('taxid', metavar='taxid', type=str,
			help='A comma-separated list of TaxIDs and/or taxon names. (e.g. 561,2172)')
    parser.add_argument('-v', '--verbose', action='count', default=0,
			help='Verbose behaviour. Supports 3 levels at present: Off = 0, Info = 1, Verbose = 2. (default: %(default)s)')
    parser.add_argument('-d', '--database', type=str, default=None,
			help='NCBI taxonomy database file path. If "None", it will be downloaded (default: %(default)s)')
    parser.add_argument('-u', '--update', action='store_true', default=False,
                        help='Update the local taxon database before querying. Recommended if not used for a while. (default: %(default)s)')
    parser.add_argument('-j', '--just-taxids', action='store_true', default=False,
                        help='Just write out a list of taxids an no other information (default: %(default)s)')
    parser.add_argument('-i', '--taxon-info', action='store_true', default=False,
                        help='Just write out rank & lineage info on the provided taxids (default: %(default)s)')
    parser.add_argument('-o', '--outfile', action='store',
			help='Output file to store the descendent TaxIDs for the query.')
    return parser.parse_args()

def pretty(d, indent=0):
    """A prettier way to print nested dicts
    """
    for key, value in d.items():
        print('  ' * indent + str(key))
        if isinstance(value, dict):
            pretty(value, indent+1)
        else:
            sys.stderr.write('  ' * (indent+1) + str(value) + '\n')

def desc_taxa(taxid, ncbi, outFH, just_taxids=False):
    """Write descendent taxa for taxid
    """
    # Main feature of the script is to get all taxa within a given group.
    descendent_taxa = ncbi.get_descendant_taxa(taxid)
    descendent_taxa_names = ncbi.translate_to_names(descendent_taxa)

    if just_taxids:
        for taxid in descendent_taxa:
            outFH.write(str(taxid) + '\n')
    else:
        for dtn, dt in zip(descendent_taxa_names, descendent_taxa):
            x = [str(x) for x in [taxid, dt, dtn]]
            outFH.write('\t'.join(x) + '\n')

def taxon_info(taxid, ncbi, outFH):
    """Write info on taxid
    """
    taxid = int(taxid)
    tax_name = ncbi.get_taxid_translator([taxid])[taxid]
    rank = list(ncbi.get_rank([taxid]).values())[0]
    lineage = ncbi.get_taxid_translator(ncbi.get_lineage(taxid))
    lineage = ['{}:{}'.format(k,v) for k,v in lineage.items()]
    lineage = ';'.join(lineage)
    x = [str(x) for x in [tax_name, taxid, rank, lineage]]
    outFH.write('\t'.join(x) + '\n')

def name2taxid(taxids, ncbi):
    """Converting taxon names to taxids
    """
    new_taxids = []
    for taxid in taxids:
        try:
            new_taxids.append(ncbi.get_name_translator([taxid])[taxid][0])
        except KeyError:
            try:
                new_taxids.append(int(taxid))
            except ValueError:
                msg = 'Error: cannot convert to taxid: {}'
                raise ValueError(msg.format(taxid))

    return new_taxids

def main():
    """Make queries against NCBI Taxa databases
    """
    # Get commandline args
    args = get_args()

    # Instantiate the ete NCBI taxa object
    ncbi = NCBITaxa(dbfile=args.database)
    ## dbfile location
    if args.verbose > 1:
        sys.stderr.write('Taxa database is stored at {}\n'.format(ncbi.dbfile))

    # Update the database if required.
    if args.update is True:
        if args.verbose > 1:
            msg = 'Updating the taxonomy database. This may take several minutes...\n'
            sys.stderr.write(msg)
        ncbi.update_taxonomy_database()

    # If names were provided in taxid list, convert to taxids
    args.taxid = args.taxid.replace('"', '').replace("'", '').split(',')
    args.taxid = name2taxid(args.taxid, ncbi)

    # Output
    if args.outfile is None:
        outFH = sys.stdout
    else:
        outFH = open(args.outfile, 'w')
    ## header
    if args.taxon_info:
        outFH.write('\t'.join(['name', 'taxid', 'rank', 'lineage']) + '\n')
    elif not args.just_taxids:
        outFH.write('\t'.join(['parent_taxid',
                               'descendent_taxid',
                               'descendent_name']) + '\n')
    ## body
    for taxid in args.taxid:
        if args.taxon_info:
            taxon_info(taxid, ncbi, outFH)
        else:
            desc_taxa(taxid, ncbi,  outFH, args.just_taxids)

    outFH.close()


if __name__ == "__main__":
	main()

```

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


# Gene Trees of Any Variety

#### Runs the blast, pulls the results column, deduplicates it, and grabs the sequences from the master file (locally)
```
for f in $(cat list.txt); do blastp -query apicoplastgenes/"$f".faa -db 1258_db/1258_PROKKA_01052020.faa -outfmt 6 -evalue .1 -max_target_seqs 3 > "$f"_1258_blast.txt && cut -f2 "$f"_1258_blast.txt | sort -u > "$f"_1258_head.txt && fasomerecords 1258_PROKKA_01052020.faa "$f"_1258_head.txt done_1258/"$f"_1258_final.txt ; done
```

#### Concatenates the files of interest (locally)
```
for f in $(cat list.txt); do cat apicoplastgenes/"$f".faa done_1212/"$f"_1212_final.txt done_1258/"$f"_1258_final.txt done_1263/"$f"_1263_final.txt > concat/"$f"_cat_apico.fasta; done
```

##### Remove duplicate sequences (with different headers) if necessary and rename any duplicate headers so RAxML doesn't shit a brick
```
for f in $(cat list.txt); do sed 's/*//g' concat_final/"$f"_cat_apico.fasta > concat_final/"$f"_cat_apico_destar.fasta ; done
```
(remove the * in the fasta first because this doesn't play nice with cdhit apparently)
```
for f in $(cat list.txt); do cd-hit -i concat_final/"$f"_cat_apico_destar.fasta -o concat_final/cdhit/"$f"_cat_apico_destar_100.fasta -c 1.0 -t 1 -n 5 && \
perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' concat_final/cdhit/"$f"_cat_apico_destar_100.fasta > concat_final/dedup/"$f"_cat_apico_destar_100_dedup.fasta ; done
```

##### Move the final concatenated fastas up to a server

#### Run the alignments, trimming, and RAxML trees in parallel (on a server)
```
#!/bin/bash
#SBATCH -J mafft_raxml
#SBATCH -t 3-00:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --array=1-35%6

echo "START"
date

module load mafft/7.215
module load trimAl/1.4.1-GCC-8.2.0-2.31.1
module load raxml/8.2.3

FILE=$(head -n $SLURM_ARRAY_TASK_ID list.txt | tail -n 1)

mafft --auto concat_dedup/${FILE}_cat_apico_dedup.fasta > aligned/${FILE}_aligned.fasta && \
trimal -in aligned/${FILE}_aligned.fasta -out aligned_trimmed/${FILE}_aligned_trimmed.fasta -gt .7 -cons 60 && \
raxmlHPC-PTHREADS -m PROTGAMMALG -T 8 -f a -x 1034 -p 1034 -N 100 -s aligned_trimmed/${FILE}_aligned_trimmed.fasta -n ${FILE}.tre

echo "DONE"
date
```

#### Parse these manually for paralogs
- remove paralogs (should have one copy of every gene in the final fastas)
- realign those that were edited to remove paralogs
- make sure all the headers are identical for each species (for every protein)
- concatenate with catfasta2phyml program
- run tree (IQtree = fast but RAxML is better)
