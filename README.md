#DNA Sequence Alignment

`DNAalign.py` takes 2 DNA sequences, in the form of text files, and calculates a similarity score that measures how similar the 2 DNA sequences are and an alignment of the same 2 DNA sequences.

` DNAalign.py` uses a dynamic programming algorithm in order to efficiently compute the similarity score for 2 DNA sequences.

To run the program:

```
Python DNAalign.py DNASeq1.txt DNASeq2.txt out.txt
cat out.txt
```  

DNA Similarity score results will look like:

```
Similarity score: 6.85
Sequence alignment1: -G-TG-A-ACG-C-TGGCGG-CGTG-CTAAAA
Sequence alignment2: AGCT-AATACCCCAT-----ACGT-TC-----
```

Run on `Python version 2.7.9`
