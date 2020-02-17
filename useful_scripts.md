### Rename fasta headers numerically
```awk '/^>/{print ">" ++i; next}{print}' < file.fasta ```

### Splits fastas into chunks to run arrays
```awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%50==0){file=sprintf("%d_split.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < input.fasta```

### Renames split files sequentially 
```n=1; for x in *.fa; do mv $x  $n.fa; n=$(($n+1)); done```
