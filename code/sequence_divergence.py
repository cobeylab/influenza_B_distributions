#!/usr/bin/python

from Bio import SeqIO
from random import sample

import sys
sys.path.insert(0, './')
from copy import deepcopy
from bs4 import BeautifulSoup
from parse_numbering_annotation import pull_seqs_from_HTML, extract_year



# Site numbering for flu A based on the reference column of numbering annotation tables
# https://www.fludb.org/brc/haNumbering.spg?method=ShowCleanInputPage&decorator=influenza
# Influenza A: "H3" reference sequence (A/AICHI/2/68)
# For B, based on Kirkpatrick et al. 2018 (Scientific Reports), Fig. S2 compared to our alignment
# Numbering starts at 1
reference_site_partitions = {'influenza_A':{'stalk' : range(1,52) + range(278,515),
                                            'head' : range(52,278)},
                             'influenza_B':{'stalk': range(16,57) + range(311,551),
                                            'head' : range(57,311),
                                            '120_loop' : [90,92,131,133,137,144,152],
                                            '150_loop' : [156] + range(159,167),
                                            '160_loop' : range(177,188),
                                            '190_helix' : range(214,221),
                                            # Full concatenated seq. including head and stalk
                                            'head_plus_stalk' : range(16,551)}
}
# Sequence from initial methionine up to beginning of stalk                             
B_HA_lead_seq_sites = range(1,16)
    


def has_frameshift(seq):
    frameshift = False
    for codon_position in range(len(seq) / 3):
        codon = seq[codon_position * 3: codon_position * 3 + 3]
        n_gaps = len([base for base in codon if base == '-'])

        # Sequence has a frameshift mutation if any codon has 1 or 2 gaps (but not 3)
        if n_gaps == 1 or n_gaps == 2:
            frameshift = True

    return frameshift

# Reads alignment into dictionary of years > lineages > sequences
def read_alignment(alignment_file, reference_site_partitions, flu_type):
    alignment = {}
    
    # If no region partitions given, assume it's NA.
    if reference_site_partitions is None:
        regions = 'NA'
    else:
        regions = reference_site_partitions['influenza_B'].keys()    

    if flu_type == 'B':
        year_dic = {'Victoria':{},'Yamagata':{}}
    else:
        assert flu_type == 'A'
        year_dic = {'H1N1':{},'H3N2':{}}
    
    with open(alignment_file, 'rU') as h:
        for record in SeqIO.parse(h, "fasta"):
            year = extract_year(record.id)
            isolate_name = record.id.split('|')[1]       
            

            excluded_seq = isolate_name in ['B/Ann_Arbor/1994',
                                            'B/Kanagawa/73',
                                            'B/Catalonia/NSVH100750997/2018',
                                            'B/Catalonia/NSVH100773835/2018'] and alignment_file.find('NA_') > -1
            
            # Remove two non-sensical NA sequences
            if not excluded_seq:
                # Correcting two isolates misidentified as B/Victoria
                if isolate_name in ['B/Hong_Kong/548/2000','B/Victoria/504/2000']:
                    lineage = 'Yamagata'
                else:
                    lineage = record.id.split('|')[2]
                if year not in alignment.keys():
                    alignment[year] = deepcopy(year_dic)
    
                if isolate_name not in alignment[year][lineage]:
                    sequence = str(record.seq)
                    if regions == 'NA':
                        alignment[year][lineage][isolate_name] = {regions: sequence}
                    else:
                        alignment[year][lineage][isolate_name] = {}
                        for region in regions:
                            region_sites = reference_site_partitions['influenza_B'][region]
                            # Adjust region sites for python numbering
                            region_sites = [s - 1 for s in region_sites]
                            region_seq = [sequence[i] for i in range(len(sequence)) if i in region_sites]
                            region_seq = ''.join(region_seq)
                            
                            alignment[year][lineage][isolate_name][region] = region_seq
                        
    return alignment
                
     

def get_reference_seq(lineage, alignment, region):
    if lineage == 'Yamagata':
        return alignment['1988']['Yamagata']['B/Yamagata/16/88'][region]
    else:
        assert lineage == 'Victoria'
        return alignment['1987']['Victoria']['B/Victoria/2/87'][region]
    
def pairwise_distance(seq1, seq2, normalize = True):
    """
    :param seq1: string with amino acid sequence
    :param seq2: string with amino acid sequence
    :param sites: list of sites (in python indexing) to be considered
    :return: amino acid divergence between seqs. (n diffs. / number of sites excluding gaps and ambiguities)
    """
    
    assert len(seq1) == len(seq2)
    
    sites = range(len(seq1))
    
    # Number of sites with gaps or ambiguities in one sequence or the other
    invalid_sites = [i for i in sites if seq1[i] in {'-','X'} or seq2[i] in {'-','X'}]

    # Number of sites that are different between sequence and reference sequence
    diff_sites = [i for i in sites if seq1[i] != seq2[i]]

    # Exclude invalid sites (with gaps or ambiguities)
    diff_sites = [site for site in diff_sites if site not in invalid_sites]

    if (len(sites) - len(invalid_sites)) == 0:
        return 'NA'
    else:
        if normalize is True:
            # Normalize difference by length of sequence excluding invalid sites
            distance = float(len(diff_sites)) / (len(sites) - len(invalid_sites))
        else:
            distance = len(diff_sites)

    return distance

assert(pairwise_distance('MAANRSCNASCNPSC','MXANPSCDAS-NPSA') == float(3)/13)


def year_pairwise_diffs(alignment, year):
    
    # Names of lineages or subgroups present in alignment object
    lineages = alignment[year].keys()
    
    # Regions in the alignment ('head' & 'stalk' for HA or just 'NA' for NA)
    if len(alignment[year][lineages[0]]) > 0:
        regions = alignment[year][lineages[0]].itervalues().next().keys()
    else:
        regions = alignment[year][lineages[1]].itervalues().next().keys()
    
    empty_pwdist_dic = {}
    for region in regions:
        empty_pwdist_dic[region] = {}
    
    
    pwdists = {lineages[0]: deepcopy(empty_pwdist_dic),
               lineages[1]: deepcopy(empty_pwdist_dic),
               lineages[0] + '_vs_' + lineages[1]: deepcopy(empty_pwdist_dic)}
    
    
    lin1_strains = alignment[year][lineages[0]].keys()
    lin2_strains = alignment[year][lineages[1]].keys()
    
    if len(lin1_strains) > 100:
        lin1_strains = sample(lin1_strains, 100)
    if len(lin2_strains) > 100:
        lin2_strains = sample(lin2_strains, 100)
        
    for lineage in lineages:
        for region in regions:
            if region in ['120_loop','150_loop','160_loop','190_helix']:
                normalize = False
            else:
                normalize = True
                
            if lineage == lineages[0]:
                strains_list = lin1_strains
            else:
                strains_list = lin2_strains
            
            if len(strains_list) > 1: 
                for focal_strain in strains_list:
                    other_strains = [s for s in strains_list if s != focal_strain]
                    if len(other_strains) > 0:
                        focal_seq =  alignment[year][lineage][focal_strain][region]
                        for strain in other_strains:
                            pair = focal_strain + ';' + strain
                            # Check if pair already exists, if so skip
                            existing_pairs = [set(p.split(';')) for p in pwdists[lineage][region].keys()]
                            
                            pair_exists = set([focal_strain, strain]) in existing_pairs
                            
                            if pair_exists is False:
                                other_seq = alignment[year][lineage][strain][region]
                                dist = pairwise_distance(focal_seq, other_seq,
                                                         normalize)
                                pwdists[lineage][region][pair] = dist
                                
    # Between lineage differences
    if len(lin1_strains) > 0 and len(lin2_strains) > 0:
        for region in regions:
            if region in ['120_loop','150_loop','160_loop','190_helix']:
                normalize = False
            else:
                normalize = True
                
            for lin1_strain in lin1_strains:
                for lin2_strain in lin2_strains:
                    # Excludes conflicting lineage assignments
                    if lin1_strain != lin2_strain: 
                        pair = lineages[0] + ':' + lin1_strain + ';' + lineages[1] + ':' + lin2_strain
                        lin1_seq = alignment[year][lineages[0]][lin1_strain][region]
                        lin2_seq = alignment[year][lineages[1]][lin2_strain][region]
                        pwdists[lineages[0] + '_vs_' + lineages[1]][region][pair] = pairwise_distance(lin1_seq,
                               lin2_seq, normalize)
    return pwdists

def divergence_from_reference(alignment, year):
    
     # Regions in the alignment ('head' & 'stalk' for HA or just 'NA' for NA)
    if len(alignment[year]['Victoria']) > 0:
        regions = alignment[year]['Victoria'].itervalues().next().keys()
    else:
        regions = alignment[year]['Yamagata'].itervalues().next().keys()

    
    empty_pwdist_dic = {}
    for region in regions:
        empty_pwdist_dic[region] = {}
 
    
    divref = {'Victoria': deepcopy(empty_pwdist_dic),
              'Yamagata': deepcopy(empty_pwdist_dic)}
    
    vic_strains = alignment[year]['Victoria'].keys()
    yam_strains = alignment[year]['Yamagata'].keys()
    
    #if len(vic_strains) > 100:
    #    vic_strains = sample(vic_strains, 100)
    #if len(yam_strains) > 100:
    #    yam_strains = sample(yam_strains, 100)
        
    for lineage in ['Victoria','Yamagata']:
        for region in regions:
            if region in ['120_loop','150_loop','160_loop','190_helix']:
                normalize = False
            else:
                normalize = True
            if lineage == 'Victoria':
                strains_list = vic_strains
                ref_seq = get_reference_seq('Victoria', alignment, region)
            else:
                strains_list = yam_strains
                ref_seq = get_reference_seq('Yamagata', alignment, region)
            if len(strains_list) > 0: 
                for strain in strains_list:
                    strain_seq = alignment[year][lineage][strain][region]
                    divergence = pairwise_distance(strain_seq, ref_seq,
                                                   normalize)
                    divref[lineage][region][strain] = divergence
    return divref
         
def find_PNGS_sites(seq):
    """
    :param seq: an amino acid sequence
    :return: number of glycosylation motifs, N-X-[ST], where X is not proline
    """
    PNGS_sites = []

    for i in range(len(seq) - 3):
        motif = seq[i:i+4]
        if motif[0] == 'N':
            if motif[1] not in  ['P','X','-'] and motif[2] in ['S','T']:
                PNGS_sites.append(i)

    return PNGS_sites

   
def main(argv):
    
    B_HA_alignment_file = '../results/sequence_divergence/B_HA_alignment.fasta'
    B_NA_alignment_file = '../results/sequence_divergence/B_NA_alignment.fasta'
    A_NA_alignment_file = '../results/sequence_divergence/A_NA_alignment.fasta'
    
    
    B_HA_alignment = read_alignment(B_HA_alignment_file, reference_site_partitions, 'B')
    B_NA_alignment = read_alignment(B_NA_alignment_file, None, 'B')
    A_NA_alignment = read_alignment(A_NA_alignment_file, None, 'A')
    
    
    # Parse influenza A annotation HTMLs from IRC (Burke & Smith numbering)
    annotation_file_H3N2 = '../results/sequence_divergence/H3N2_H3_numbering.htm'
    with open(annotation_file_H3N2, 'rU') as fp:
        parsed_H3N2_file = BeautifulSoup(fp, "html.parser")
        
    annotation_file_H1N1 = '../results/sequence_divergence/H1N1_H3_numbering.htm'
    with open(annotation_file_H1N1, 'rU') as fp:
        parsed_H1N1_file = BeautifulSoup(fp, "html.parser")
        
    H3N2_seqs = pull_seqs_from_HTML(parsed_H3N2_file, 'influenza_A','H3N2', 
                                    reference_site_partitions)
    
    H1N1_seqs = pull_seqs_from_HTML(parsed_H1N1_file, 'influenza_A', 'H1N1',
                                    reference_site_partitions)
    
    # Merge influenza A sequences into a single object 
    # Same structure as the processed influenza B alignment from read_alignment
    A_HA_alignment = {}
    for year in list(set(H3N2_seqs.keys() + H1N1_seqs.keys())):
        A_HA_alignment[year] = {}
        if year not in H3N2_seqs.keys():
            H3N2_dict = {'H3N2':[]}
        else:
            H3N2_dict = H3N2_seqs[year]
        if year not in H1N1_seqs:
            H1N1_dict = {'H1N1':[]}
        else:
            H1N1_dict = H1N1_seqs[year]
        A_HA_alignment[year].update(H3N2_dict)
        A_HA_alignment[year].update(H1N1_dict)
        

    with open('../results/sequence_divergence/pairwise_divergence.csv','w') as results:
        results.write('segment,year,pair,pair_type,pairwise_divergence\n')
        for alignment in [A_HA_alignment, A_NA_alignment, B_HA_alignment, B_NA_alignment]:
            for year in sorted(alignment.keys()):
                print year
                divergence = year_pairwise_diffs(alignment, year)
                    
                for pair_type in divergence.keys():
                    for region in divergence[pair_type].keys():
                        if len(divergence[pair_type][region]) > 0:
                            for pair in divergence[pair_type][region].keys():
                                value = str(divergence[pair_type][region][pair])
                                results.write(','.join([region,year,pair,
                                                        pair_type, value]))
                                results.write('\n')
                
    with open('../results/sequence_divergence/divergence_from_reference.csv','w') as results:
        results.write('segment,year,lineage,strain,divergence_from_reference\n')
        for alignment in [B_HA_alignment, B_NA_alignment]:
            for year in sorted(alignment.keys()):
                divref = divergence_from_reference(alignment, year)
                for lineage in divref.keys():
                    for region in divref[lineage].keys():
                        for strain in divref[lineage][region].keys():
                            row = ','.join([region,year,lineage, strain,
                                        str(divref[lineage][region][strain])]) + '\n'
                            results.write(row)
                    
    with open('../results/sequence_divergence/glycosylation_over_time.csv','w') as glyc_results:
        glyc_results.write('segment,year,lineage,strain,glycosylation_site\n')
        for alignment in [B_HA_alignment, B_NA_alignment]:
            if alignment is B_HA_alignment:
                region = 'head_plus_stalk'
                label = 'HA'
            else:
                region = 'NA'
                label = 'NA'
            for year in sorted(alignment.keys()):
                for lineage in alignment[year].keys():
                    for strain in alignment[year][lineage].keys():
                        full_seq = alignment[year][lineage][strain][region]
                        PNGS_sites = find_PNGS_sites(full_seq)
                        # Convert site numbers to alignment ones by adding length of lead signal sequence
                        PNGS_sites = [s + len(B_HA_lead_seq_sites) for 
                                      s in PNGS_sites]
                        # Convert to numbering starting at 1
                        PNGS_sites = [s + 1 for s in PNGS_sites]
                        
                        row = [','.join([label,year,lineage,strain,str(s)]) for 
                               s in PNGS_sites]
                        row = '\n'.join(row) + '\n'
                        glyc_results.write(row)
                        
    
if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)

