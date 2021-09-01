# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 12:59:35 2021

@author: jfgui
"""

import collections

#%%

'''
Counting DNA Nucleotides:
Create Dictionary of Nucleotide:count key-pairs
'''

def count_nt(sequence_unstripped):
    sequence = sequence_unstripped.strip('\n')            
    nt_counts = {'A':0,
              'C':0,
              'G':0,
              'T':0}
    for nt in sequence:
        nt_counts[nt] += 1
    return nt_counts


#%%

'''
Transcribing DNA into RNA:
Creates RNA strand complementary to non-coding DNA strand
'''

def transcribe(sequence_unstripped):
    sequence = sequence_unstripped[0].strip('\n')
    RNAString = ''
    for nt in sequence:
        if nt == 'T':
            RNAString += 'U'
        else:
            RNAString += nt
    return RNAString

#%%

'''
(Not a Rosalind Problem)
Creates reverse complement strand of DNA
'''

def rev_comp(seq):
    # sequence = seq_unstripped[0].strip('\n')
    rev_comp = ''
    for nt in seq:
        if nt == 'A':
            rev_comp = 'T' + rev_comp
        elif nt == 'T':
            rev_comp = 'A' + rev_comp
        elif nt == 'C':
            rev_comp = 'G' + rev_comp
        else:
            rev_comp = 'C' + rev_comp
    return rev_comp


#%%

'''
(Not a Rosaling Problem)
Get GC content of specific sequence
'''

def get_gc_cont(seq):
    # seq = seq_raw[0].strip('\n')
    gc_count = 0
    for nt in seq:
        if nt == 'C' or nt == 'G':
            gc_count += 1
    total_nt = len(seq)
    gc_cont = gc_count / total_nt
    return gc_cont


#%%

'''
(Not a Rosalind Problem)
Define dictionary of Fasta:Sequence key-pairs
'''

def create_fasta_dictionary(fasta_seq):
    fasta_dict = {}
    for i in fasta_seq:
        if i[0] == '>':
            fasta_name = i[1:-1]
            fasta_dict[fasta_name] = ''
            seq_start = fasta_seq.index(i) + 1
            for j in fasta_seq[seq_start:]:
                if j[0] != '>':
                    fasta_dict[fasta_name] += j.strip('\n')
                else:
                    break
    return fasta_dict


#%%

'''
Computing GC Content:
Finds strand in FASTA file with highest GC content
'''

def high_gc_cont(fasta_seq):
    highest_gc_content = 0
    seq_dict = create_fasta_dictionary(fasta_seq)
    for i in seq_dict:
        current_gc_content = get_gc_cont(seq_dict[i])
        if current_gc_content > highest_gc_content:
            highest_gc_content_name = i
            highest_gc_content = round(current_gc_content, 6)
            highest_gc = (highest_gc_content_name + '\n' + str(highest_gc_content*100))
    return highest_gc


#%%

'''
Complementing a Strand of DNA
'''

def get_comp(seq):
    # sequence = seq_unstripped[0].strip('\n')
    rev_comp = ''
    for nt in seq:
        if nt == 'A':
            rev_comp += 'T'
        elif nt == 'T':
            rev_comp += 'A'
        elif nt == 'C':
            rev_comp += 'G'
        else:
            rev_comp += 'C'
    return rev_comp


#%%

'''
Counting Point Mutations:
Find Hamming Distance between two strands
'''

with open('C:/Users/jfgui/Dropbox/Programming/rosalind_data/rosalind_hamm.txt',
          'r') as rosalind_hamming:
    hamming_strands = rosalind_hamming.readlines()
    hamstrings_unsplit = hamming_strands[0].split('\n')
    hamstrings = [hamming_strands[0], hamming_strands[1]]


def find_hamming_distance(strands):
    hamming_distance = 0
    hamstring1 = strands[0].strip('\n')
    hamstring2 = strands[1].strip('\n')
    for i in range(len(hamstring1)):
        if hamstring1[i] != hamstring2[i]:
            hamming_distance += 1
    return hamming_distance


#%%

'''
Translatin RNA into Protein
'''

with open('C:/Users/jfgui/Dropbox/Programming/rosalind_data/rosalind_prot.txt',
          'r') as rosalind_protein:
    rna_sequence = rosalind_protein.readlines()
    
codon_map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


def translate_rna(seq_unstrip):
    seq = seq_unstrip[0].strip('\n')
    protein_seq = ''
    codon_seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
    for codon in codon_seq:
        protein_seq += codon_map[codon]
    protein_seq = protein_seq[:-4]
    return protein_seq


#%%

'''
Finding a Motif in DNA
'''

with open('C:/Users/jfgui/Dropbox/Programming/rosalind_data/rosalind_subs.txt',
          'r') as rosalind_motif:
    motif_full = rosalind_motif.readlines()

def find_motif(motif_full):
    motif_index = ''
    seq = motif_full[0].strip('\n')
    motif = motif_full[1].strip('\n')
    for i in range(len(seq)):
        if seq[i:i+len(motif)] == motif:
            motif_index += str(i+1) + ' '
    return motif_index


#%%

'''
Rabbits and Recurrence Relations:
Recurrence Relation with immortal rabbit population
'''

# Simple Fibonacci calculator
def fibonaccis_seq(iters):
    new, old = 1, 1
    for i in range(iters):
        temp = new
        new += old
        old = temp
    return new

# Finds population of rabbits after n months
# When rabbits have n pairs of offspring
def rabbit_population(months, offspring):
    mat_pop, juv_pop = 1, 1
    for i in range(months - 1):
        juv_pop, mat_pop = mat_pop, mat_pop + (juv_pop* offspring)
    return juv_pop


#%%

'''
Calculating Protein Mass:
Calculating the atomic mass of protein from sequence of amino acids
'''

with open('C:/Users/jfgui/Dropbox/Programming/rosalind_data/rosalind_prtm.txt', 'r') as prob:
    aa_seq = prob.readlines()

mass_h2o = 18.01056
aa_mass_dict = {'A':71.03711,
'C':103.00919,
'D':115.02694,
'E':129.04259,
'F':147.06841,
'G':57.02146,
'H':137.05891,
'I':113.08406,
'K':128.09496,
'L':113.08406,
'M':131.04049,
'N':114.04293,
'P':97.05276,
'Q':128.05858,
'R':156.10111,
'S':87.03203,
'T':101.04768,
'V':99.06841,
'W':186.07931,
'Y':163.06333
}

def find_mass_protein(aa_seq_raw):
    aa_seq = aa_seq_raw[0].strip('\n')
    mass_protein = 0
    for aa in range(len(aa_seq)):
        mass_protein += aa_mass_dict[aa_seq[aa]]
    return mass_protein


#%%

'''
Mendel's First Law:
Mendelian Probability of populations
'''

with open('C:/Users/jfgui/Dropbox/Programming/rosalind_data/rosalind_iprb.txt', 'r') as prob:
    prob_vars = prob.readlines()

prob_vars = prob_vars[0].split(' ')
dom, het, rec = int(prob_vars[0]), int(prob_vars[1]), int(prob_vars[2])

def prob_dominant(dom, het, rec):
    total_al = dom + het + rec
    h_h = (het/total_al)*((het-1)/(total_al-1))
    r_r = (rec/total_al)*((rec-1)/(total_al-1))
    h_r = ((rec/total_al)*(het/(total_al-1))) + ((het/total_al)*(rec/(total_al-1)))
    p = 1 - (r_r + (h_h * 0.25) + (h_r * 0.5))
    return p


#%%

'''
Consensus and Profile:
Finding consensus sequence
'''

with open('C:/Users/jfgui/Dropbox/Programming/rosalind_data/rosalind_cons.txt', 'r') as prob:
    fasta_cons = prob.readlines()

def consensus_and_profile(fasta_cons):  
    fasta_cons_dict = collections.OrderedDict(create_fasta_dictionary(fasta_cons))
    fasta_cons_dict_items = list(fasta_cons_dict.values())
    
    sequence_matrix = []
    
    for i in range(len(fasta_cons_dict)):   #Create list of sublists containing single nucleotides as strings
        sequence_matrix.append([])   # Add empty sublist
        for j in fasta_cons_dict_items[i]:
            sequence_matrix[i].append(j)
            
    profile_matrix_A = [0]*len(sequence_matrix[0])
    profile_matrix_C = [0]*len(sequence_matrix[0])
    profile_matrix_G = [0]*len(sequence_matrix[0])
    profile_matrix_T = [0]*len(sequence_matrix[0])
    
    
    for seq in range(len(sequence_matrix)):
        for nt in range(len(sequence_matrix[seq])):
            if sequence_matrix[seq][nt] == 'A':
                profile_matrix_A[nt] += 1
            elif sequence_matrix[seq][nt] == 'C':
                profile_matrix_C[nt] += 1
            elif sequence_matrix[seq][nt] == 'G':
                profile_matrix_G[nt] += 1
            elif sequence_matrix[seq][nt] == 'T':
                profile_matrix_T[nt] += 1
    
    nts = ['A', 'C', 'G', 'T']
    profile_matrix = [profile_matrix_A, profile_matrix_C, profile_matrix_G, profile_matrix_T]
    ans_prof_string = ''
    for i in range(len(profile_matrix)):
        ans_prof_string+=nts[i]+': '
        for j in profile_matrix[i]:
            ans_prof_string += str(j) + ' '
        ans_prof_string = ans_prof_string.rstrip()
        ans_prof_string += '\n'
    
    locus_nt_count = [[]]
    for i in range(len(sequence_matrix[0])):
        locus_nt_count.append([])
        for j in range(len(profile_matrix)):
            locus_nt_count[i].append(profile_matrix[j][i])
    locus_nt_count = locus_nt_count[0:-1]
    consensus_string = ''
    for i in range(len(locus_nt_count)):
        max_nt_count = max(locus_nt_count[i])
        for j in range(len(profile_matrix)):
            if profile_matrix[j][i] == max_nt_count:
                consensus_string += nts[j]
                break
    return consensus_string, ans_prof_string


#%%

'''
Mortal Fibonacci Rabbits:
Recurrence Relation with mortal rabbit population
'''

def rabbit_pop_mortal(months, lifespan):
    rabbits = [1,1]
    month = 2
    while month < months:
        if month < lifespan:
            rabbits.append((rabbits[-2]+rabbits[-1]))
            month += 1
        elif month == lifespan:
            rabbits.append(rabbits[-1]+rabbits[-2]-rabbits[-lifespan])
            month += 1
        else:
            rabbits.append((rabbits[-1]+rabbits[-2]-rabbits[-1-lifespan]))
            month += 1
    return rabbits[-1]


#%%





























    
    