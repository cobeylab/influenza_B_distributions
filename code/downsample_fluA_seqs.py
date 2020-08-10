"""
Samples up to 100 H3N2 and H1N1 sequences for each year for analysis of stalk divergence
"""

from Bio import SeqIO
from sequence_divergence import extract_year
from random import sample
import sys

def add_subtype_to_id(record_list, subtype):
    for record in record_list:
        record_id = record.id
        record_id = record_id.split('|')
        record_id[2] = subtype
        record.id = '|'.join(record_id)
        record.description = ''
    return record_list
    


def downsample(seq_file, max_seqs_by_year):
    
    seqs_by_year = {}
    
    downsampled_seqs = []
    
    with open(seq_file, 'rU') as h:
        for record in SeqIO.parse(h, "fasta"):
            year = extract_year(record.id)

            if year not in seqs_by_year.keys():
                seqs_by_year[year] = []
            seqs_by_year[year].append(record)
            
    for y in [y for y in seqs_by_year.keys() if y != 'NA']:
        year_seqs = seqs_by_year[y]
        n_seqs = len(year_seqs)
        
        if n_seqs > max_seqs_by_year:
            retained_seqs = sample(year_seqs, max_seqs_by_year)
            downsampled_seqs = downsampled_seqs + retained_seqs
        else:
            downsampled_seqs = downsampled_seqs + year_seqs

    return(downsampled_seqs)
    
    
def main(argv):
    seq_file_H3N2_HA = '../data/sequence_divergence/gisaid_H3N2_HA.fasta'
    seq_file_H1N1_HA = '../data/sequence_divergence/gisaid_H1N1_HA.fasta'
    seq_file_H3N2_NA =  '../data/sequence_divergence/gisaid_H3N2_NA.fasta'
    seq_file_H1N1_NA = '../data/sequence_divergence/gisaid_H1N1_NA.fasta'
#    seq_file_Vic = '../data/sequence_divergence/gisaid_Vic_HA.fasta'
#    seq_file_Yam = '../data/sequence_divergence/gisaid_Yam_HA.fasta'
    
    
    max_seqs_by_year = 100
    
    downsampled_H3N2_HA_seqs = downsample(seq_file_H3N2_HA, max_seqs_by_year)
    downsampled_H3N2_HA_seqs = add_subtype_to_id(downsampled_H3N2_HA_seqs, 'H3N2')
    
    with open("../data/sequence_divergence/H3N2_downsampled_HAs.fasta", "w") as output_handle:
        SeqIO.write(downsampled_H3N2_HA_seqs, output_handle, "fasta")
    
    downsampled_H1N1_HA_seqs = downsample(seq_file_H1N1_HA, max_seqs_by_year)
    downsampled_H1N1_HA_seqs = add_subtype_to_id(downsampled_H1N1_HA_seqs, 'H1N1')
    
    with open("../data/sequence_divergence/H1N1_downsampled_HAs.fasta", "w") as output_handle:
        SeqIO.write(downsampled_H1N1_HA_seqs, output_handle, "fasta")
        
        
    downsampled_H3N2_NA_seqs = downsample(seq_file_H3N2_NA, max_seqs_by_year)
    downsampled_H3N2_NA_seqs = add_subtype_to_id(downsampled_H3N2_NA_seqs, 'H3N2')
    
    with open("../data/sequence_divergence/H3N2_downsampled_NAs.fasta", "w") as output_handle:
        SeqIO.write(downsampled_H3N2_NA_seqs, output_handle, "fasta")
        
    
    downsampled_H1N1_NA_seqs = downsample(seq_file_H1N1_NA, max_seqs_by_year)
    downsampled_H1N1_NA_seqs = add_subtype_to_id(downsampled_H1N1_NA_seqs, 'H1N1')
    
    with open("../data/sequence_divergence/H1N1_downsampled_NAs.fasta", "w") as output_handle:
        SeqIO.write(downsampled_H1N1_NA_seqs, output_handle, "fasta")

    
#    downsampled_Vic_seqs = downsample(seq_file_Vic, max_seqs_by_year)
#    downsampled_Vic_seqs = add_subtype_to_id(downsampled_Vic_seqs, 'Victoria')
#    
#    with open("../data/sequence_divergence/Vic_downsampled_HAs.fasta", "w") as output_handle:
#        SeqIO.write(downsampled_Vic_seqs, output_handle, "fasta")
#        
#    
#    downsampled_Yam_seqs = downsample(seq_file_Yam, max_seqs_by_year)
#    downsampled_Yam_seqs = add_subtype_to_id(downsampled_Yam_seqs, 'Yamagata')
#    
#    with open("../data/sequence_divergence/Yam_downsampled_HAs.fasta", "w") as output_handle:
#        SeqIO.write(downsampled_Yam_seqs, output_handle, "fasta")

    
    
    
    
if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)


    