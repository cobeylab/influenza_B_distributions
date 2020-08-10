#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Parses output of the IRD annotation based on the Burke and Smith numbering scheme
"""

import sys
sys.path.insert(0, './')
import re

def extract_year(seq_id):
    if len(seq_id.split('|')) >= 3:
        date = seq_id.split('|')[3]
        year = date.split('-')[0]
        year = year.split('_')[0]
        
        if len(year) == 0:
            year = 'NA'
    else:
        year = 'NA'
    return year


def find_seq_ids_from_html(parsed_HTML_file):
    isolates = []

    div_tags = parsed_HTML_file.find_all('div')
    
    for div in div_tags:
        if len(div.find_all('input')) > 0:        
            isolate_name = [x for x in div if x.find('EPI_ISL') > - 1]
            if len(isolate_name) > 0:
                assert len(isolate_name) == 1
                isolates.append(str(isolate_name[0]))
    return isolates
                                          

def get_seq_from_numbering_table(table, flu_type, reference_site_partitions):
    head_seq = ''
    stalk_seq = ''
    
    # Find all table headers
    table_headers = table.find_all('th')
    table_headers = [re.search(r'>[^<]*', str(x)).group() for 
                     x in table_headers]
    
    # Table has either 4 or 6 columns
    if len(table_headers) == 4:
        # 4 columns when best match is the chosen reference
        query_state_var = 2 # column with amino acid in query sequence
        ref_site_var = 1 # column with position in reference sequence
    else:
        # 6 columns when best match is not chosen reference
        assert len(table_headers) == 6
        query_state_var = 3
        ref_site_var = 2
        

    table_rows = table.find_all('tr')

    # For each site in query sequence...
    for row in table_rows[1:len(table_rows)]:
        row_values = [re.search(r'>[^<]*',str(x)).group().replace('>','') 
        for x in row.find_all('td')]
        
        reference_site = row_values[ref_site_var]
        if reference_site != '-':
            reference_site = int(reference_site)
        
            query_state = row_values[query_state_var]
        
            # If corresponding site of the reference seq. is in the head...
            if reference_site in reference_site_partitions[flu_type]['head']:
                head_seq = head_seq + query_state
            # Else if it's in the stalk, append query state to stalk sequence
            elif reference_site in reference_site_partitions[flu_type]['stalk']:
                stalk_seq = stalk_seq + query_state
    return {'head': head_seq, 'stalk': stalk_seq}
 
       
def pull_seqs_from_HTML(parsed_HTML_file, flu_type, subtype,
                        reference_site_partitions):
    isolates = find_seq_ids_from_html(parsed_HTML_file)
    
    tables = parsed_HTML_file.find_all('table')
    
    sequences = {}
    
    for i in range(len(tables)):
        isolate_seqs = get_seq_from_numbering_table(tables[i], flu_type,
                                           reference_site_partitions)
        isolate_name = isolates[i]
        year = extract_year(isolate_name)

        if year not in sequences.keys():
            sequences[year] = {subtype:{}}
        
        sequences[year][subtype][isolate_name] = isolate_seqs

    return sequences