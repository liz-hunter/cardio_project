### Splits fastas into chunks to run arrays
```awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%50==0){file=sprintf("%d_split.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < input.fasta```

### Renames split files sequentially 
```n=1; for x in *.fa; do mv $x  $n.fa; n=$(($n+1)); done```

### Replaces fasta headers with names from a tab separated file

```awk -f replace.awk list.txt file.fasta > output.txt```

**replace.awk:**
```
NR == FNR {
  rep[$1] = $2
  next
} 

{
    for (key in rep) {
      gsub(key, rep[key])
    }
    print
}
```

### Append a unique number to duplicate headers
```
perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' file.fa > outputfile.txt
```

### Make fastas single-line (preserves headers)
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < input.fa > output.fa
```
### Show physical memory use by directory 
``` du -s * ```

### Count fastq reads
``` awk '{s++}END{print s/4}' file.fastq```

### Cool scripts:
- fasomerecords (to pull fasta records - from santiagosnchez)
- interleave-fasta.py (create an interleaved fasta file with forward and reverse reads - from sebhtml)
- taxnameconvert.pl (rename trees with a tab separated file - from cibiv)
- catfasta2phyml (concatenates multi-gene alignments by header - from nylander)
- shannon.py (computes shannon entropy from multi-sequence alignments)
```
# This script will calculate Shannon entropy from a MSA.

# Dependencies:

# Biopython, Matplotlib [optionally], Math

"""
Shannon's entropy equation (latex format):
    H=-\sum_{i=1}^{M} P_i\,log_2\,P_i
    Entropy is a measure of the uncertainty of a probability distribution (p1, ..... , pM)
    https://stepic.org/lesson/Scoring-Motifs-157/step/7?course=Bioinformatics-Algorithms&unit=436
    Where, Pi is the fraction of nuleotide bases of nuleotide base type i,
    and M is the number of nuleotide base types (A, T, G or C)
    H ranges from 0 (only one base/residue in present at that position) to 4.322 (all 20 residues are equally
    represented in that position).
    Typically, positions with H >2.0 are considerered variable, whereas those with H < 2 are consider conserved.
    Highly conserved positions are those with H <1.0 (Litwin and Jores, 1992).
    A minimum number of sequences is however required (~100) for H to describe the diversity of a protein family.
"""
import os
import sys
import warnings
import traceback

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "ShannonMSA"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Compute per base/residue Shannon entropy of a Multiple Sequence Alignment.')

        parser.add_argument('-a',
                            '--alignment',
                            action='store',
                            required=True,
                            help='The multiple sequence alignment (MSA) in any of the formats supported by Biopython\'s AlignIO.')
        parser.add_argument('-f',
                            '--alnformat',
                            action='store',
                            default='fasta',
                            help='Specify the format of the input MSA to be passed in to AlignIO.')
        parser.add_argument('-v',
                            '--verbose',
                            action='count',
                            default=0,
                            help='Verbose behaviour, printing parameters of the script.')
        parser.add_argument('-m',
                            '--runningmean',
                            action='store',
                            type=int,
                            default=0,
                            help='Return the running mean (a.k.a moving average) of the MSAs Shannon Entropy. Makes for slightly smoother plots. Providing the number of points to average over switches this on.')
        parser.add_argument('--makeplot',
                            action='store_true',
                            help='Plot the results via Matplotlib.')
    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()



def parseMSA(msa, alnformat, verbose):
    """Parse in the MSA file using Biopython's AlignIO"""

    from Bio import AlignIO
    alignment = AlignIO.read(msa, alnformat)

    # Do a little sanity checking:
    seq_lengths_list = []
    for record in alignment:
       seq_lengths_list.append(len(record))

    seq_lengths = set(seq_lengths_list)

    if verbose > 0: print("Alignment length is:" + str(list(seq_lengths)))

    if len(seq_lengths) != 1:
        sys.stderr.write("Your alignment lengths aren't equal. Check your alignment file.")
        sys.exit(1)

    index = range(1, list(seq_lengths)[0]+1)

    return alignment, list(seq_lengths), index

##################################################################
# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i (http://imed.med.ucm.es/Tools/svs_help.html)
# Gaps and N's are included in the calculation
##################################################################

def shannon_entropy(list_input):
    """Calculate Shannon's Entropy per column of the alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)"""

    import math
    unique_base = set(list_input)
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base) # Number of residues of type i
        P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)

    sh_entropy = -(sum(entropy_list))

    return sh_entropy


def shannon_entropy_list_msa(alignment):
    """Calculate Shannon Entropy across the whole MSA"""

    shannon_entropy_list = []
    for col_no in xrange(len(list(alignment[0]))):
        list_input = list(alignment[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))

    return shannon_entropy_list


def plot(index, sel, verbose):
    """"Create a quick plot via matplotlib to visualise"""
    import matplotlib.pyplot as plt

    if verbose > 0: print("Plotting data...")

    plt.plot(index, sel)
    plt.xlabel('MSA Position Index', fontsize=16)
    plt.ylabel('Shannon Entropy', fontsize=16)

    plt.show()


def running_mean(l, N):
    sum = 0
    result = list(0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result

def main():
    """Compute Shannon Entropy from a provided MSA."""

    # Parse arguments
    args = parseArgs()

    # Convert object elements to standard variables for functions
    msa = args.alignment
    alnformat = args.alnformat
    verbose = args.verbose
    makeplot = args.makeplot
    runningmean = args.runningmean

# Start calling functions to do the heavy lifting

    alignment, seq_lengths, index = parseMSA(msa, alnformat, verbose)
    sel = shannon_entropy_list_msa(alignment)

    if runningmean > 0:
        sel = running_mean(sel, runningmean)

    if makeplot is True:
        plot(index, sel, verbose)

    if verbose > 0: print("Index" + '\t' + "Entropy")
    for c1, c2 in zip(index, sel):

        print(str(c1) + '\t' + str(c2))



if __name__ == '__main__':
    main() 
```    
    

- AEDgenerator.pl (generates AED values from maker annotation output)
``` 
#!/usr/bin/perl -w 
use strict;
use FileHandle;
use Getopt::Std;
use vars qw($opt_b);
getopts('b:');

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\nAED_cdf_generator.pl: returns values for a cumulative fraction plot based on AED.
USAGE: AED_cdf_generator.pl -b <bin size> <maker1.gff ,maker2.gff ...>\n
Options -b <number between 0 and 1> sets the y axis intervals. 0.025 works well.\n\n";



die($usage) unless ($ARGV[0] && $opt_b);
my $dividend = 1/$opt_b;
my ($test_value) = ($dividend  =~/\d+\.(\d+)/);

die('make sure your bin size is a dividend of 1') if $test_value;

my $BIN_UNIT = $opt_b;
my %HASH;
my $BIN_S = 0;

my $BL = set_bins();
foreach my $F (@ARGV){
    parse($F);
}
report();
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{
    my $places = 0 x $BL;
    my $max_bin = '1.' . $places;
#    print "max bin $max_bin\n";
    my @files = sort @ARGV;
    print "AED\t";
    print join("\t", @files);

    foreach my $x (sort keys (%HASH)){
	print "\n$x";
	foreach my $f (sort keys (%{$HASH{$x}})){
	    my $total = $HASH{$max_bin}{$f};
	    my $frac_below = $HASH{$x}{$f}/$total;
	    # more rounding to deal with floating point stuff
	    my $rfrac_below = sprintf("%.3f", $frac_below);
	    print "\t$rfrac_below";
	}
    }
    print "\n";
}
#-----------------------------------------------------------------------------
sub set_bins{
    my ($y) = $BIN_UNIT =~ /\d?\.(\d+)/;
    my $bl = length($y);

    while ($BIN_S <= 1){
        #I used sprintf to deal with floating point issue in perl
       
	my $x = sprintf("%.${bl}f", $BIN_S);
	foreach my $file (@ARGV){
	    $HASH{$x}{$file} = 0;
	}
	$BIN_S = $x +$BIN_UNIT;
    }
	return($bl);
}
#-----------------------------------------------------------------------------
sub parse {
    my $file = shift;	
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @cols = split(/\t/, $line);
	if ($cols[2] eq 'mRNA' && $line =~ /_AED=(\d\.\d+?)\;/){
	    load($1, $file);
	}
	else {
	    next;
	}
    }
    
    $fh->close();
}
#-----------------------------------------------------------------------------
sub load{
    my $value = shift;
    my $file  = shift;

    #go through each fo the values in the bin hash and 
    #add 1 if the value passed into the routine in less 
    #than or equal to the bin value	
    foreach my $x (keys %HASH){
	if ($value <= $x || $value eq $x){
	    $HASH{$x}{$file}++;
	}
    }
}
#-----------------------------------------------------------------------------
```

