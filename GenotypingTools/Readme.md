## Genotyping Tools

### Installation of miRNA_Search

```
gcc miRNA_Search.c -o MIRNA_SEARCH
```

### Running the genotyping

```
./MIRNA_SEARCH <INPUT> <SIGFILES> > <OUTPUT>
```

**<INPUT>**
The input is a text file of reads with one read per line. 
The quality and sequence identifiers need to me removed.

**<SIGFILES>**
This is a fasta file with unique signature motifs specific to 
each genotype. Example:

```
>RL5A-G1
ATTTTATCATTAAACCCCAATATT
>RL5A-G2
ATGTCCACGTTGAATTCAACTTAG
>RL5A-G3
GCTTCCAAGTTGACGTAAGATCGT
```

MIRNA_SEARCH will determine how many reads have each of these motifs.


