#!/bin/bash
# Uses grep to pull out identifier lines in fasta file and write them to a csv file.

echo 'accession_number,strain,year,month,day,country' > ../data/genbank_data/genbank_B_HA_seqs.csv
grep '>' ../data/genbank_data/genbank_B_HA_seqs.fasta >> ../data/genbank_data/genbank_B_HA_seqs.csv