#!/bin/bash
#SBATCH --job-name=B_seq_blast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o B_seq_blast.out        
#SBATCH -e B_seq_blast.err        
#SBATCH --time=500:00:00

# Blasts all sequences in genbank_B_HA_seqs.fasta against reference Vic and Yam strains.
module load blastplus

# First, replace all commas in fasta identifiers with _ (to avoid messing up blast)
# Also, replace all /'s (From B strain names)

cat ../data/genbank_B_HA_seqs.fasta | tr ',' '_' | tr '/' '-' | tr ' ' '|' > temp_query.fasta

cat ../data/reference_strains.fasta | tr ',' '_' | tr '/' '-' > temp_subject.fasta

blastn -task blastn -query temp_query.fasta -out ../results/B_seq_blast/B_seq_blast_temp.csv -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -subject temp_subject.fasta 

rm temp_query.fasta
rm temp_subject.fasta



# Add column names, undo slash and space replacements
echo 'query_id,subject_id,percent_identity,alignment_length,mismatches,gap_opens,query_start,query_end,subject_start,subject_end,evalue,bit_score' > ../results/B_seq_blast/B_seq_blast.csv
cat ../results/B_seq_blast/B_seq_blast_temp.csv | tr '-' '/' | tr '|' ' ' >> ../results/B_seq_blast/B_seq_blast.csv


rm ../results/B_seq_blast/B_seq_blast_temp.csv