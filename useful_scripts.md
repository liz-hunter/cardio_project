### Rename fasta headers numerically
```awk '/^>/{print ">" ++i; next}{print}' < file.fasta ```

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

### Make fastas single-line (preserves headers)
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < [input.fa] > [output.fa]```
