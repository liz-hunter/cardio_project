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
