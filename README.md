## This script combines sequences and metadata from a local database with GISAID data to produce input files for Nextstrain pipeline
To speed up querying of sequences in GISAID fasta file, sequences are stored in form SeqIO.index_db database.

## Input files:
+ sequences.fq.gz from GISAID bulk download
+ metadata.tsv from GISAID bulk download
+ local SQL database with sequences and metadata

## Usage:
+ python3 filter_GISAID.py > output.fasta

## Cleaning sequences file
The sequences file from GISAID has to be cleanded from duplicated sequences and whitespaces in sequence IDs:
+ sed -e 's/[|].*$//' -e 's/ /_/g' sequences.fasta | seqkit rmdup -n | bgzip -c > sequences_clean.fasta.bgz 
 
