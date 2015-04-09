# CDSgffSNP2fa
A perl script to generate fasta sequences of coding regions specified in gff file for a population sample.

Inputs
======

1. gff3 file containing ONLY CDS annotations.
2. single [POPBAM SNP](https://github.com/dgarriga/POPBAM) file, compressed with [BGZIP](http://samtools.sourceforge.net/tabix.shtml) and indexed with [TABIX](http://samtools.sourceforge.net/tabix.shtml).
3. text file with list of sample names, one name per line.
4. reference genome sequence in fasta format.


Output
======

Writes single file (transcriptID.fa) for each transcript. Files contain fasta formated sequences for each sample.
