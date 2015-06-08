#!/usr/bin/python
from __future__ import division
from __future__ import with_statement


aa_colors = {'start'    : 'c',
             'basic'    : 'b',
             'acidic'   : 'r',
             'polar'    : 'g',
             'nonpolar' : 'y',
             'stop'     : 'k'}

nucleotides = ['a', 'c', 'g', 'u']
start_codon = 'aug'
stop_codons = ['uaa', 'uag', 'uga']

class Codon:
    def __init__(self, name, letter, category):
        self.name = name
        self.letter = letter
        self.category = category
        self.color = aa_colors[self.category]

codons = {}

alanine       = Codon('alanine',       'A', 'nonpolar')
arginine      = Codon('arginine',      'R', 'basic')
asparagine    = Codon('asparagine',    'N', 'polar')
aspartic      = Codon('aspartic',      'D', 'acidic')
cysteine      = Codon('cysteine',      'C', 'polar')
glutamic      = Codon('glutamic',      'E', 'acidic')
glutamine     = Codon('glutamine',     'Q', 'polar')
glycine       = Codon('glycine',       'G', 'polar')
histidine     = Codon('histidine',     'H', 'basic')
isoleucine    = Codon('isoleucine',    'I', 'nonpolar')
leucine       = Codon('leucine',       'L', 'nonpolar')
lysine        = Codon('lysine',        'K', 'basic')
methionine    = Codon('methionine',    'M', 'start')
phenylalanine = Codon('phenylalanine', 'F', 'nonpolar')
proline       = Codon('proline',       'P', 'nonpolar')
serine        = Codon('serine',        'S', 'polar')
threonine     = Codon('threonine',     'T', 'polar')
tryptophan    = Codon('tryptophan',    'W', 'nonpolar')
tyrosine      = Codon('tyrosine',      'Y', 'polar')
valine        = Codon('valine',        'V', 'nonpolar')

codons[start_codon] = methionine
for nt in nucleotides:
    codons['ac'+nt] = threonine
codons['aac'] = asparagine
codons['aau'] = asparagine
codons['aag'] = lysine
codons['aaa'] = lysine
codons['agc'] = serine
codons['agu'] = serine
for nt in nucleotides:
    codons['uc'+nt] = serine
codons['aga'] = arginine
codons['agg'] = arginine
for nt in nucleotides:
    codons['cg'+nt] = arginine
for nt in nucleotides:
    codons['gu'+nt] = valine
for nt in nucleotides:
    codons['gc'+nt] = alanine
codons['gau'] = aspartic
codons['gac'] = aspartic
codons['gaa'] = glutamic
codons['gag'] = glutamic
for nt in nucleotides:
    codons['gg'+nt] = glycine
codons['uuu'] = phenylalanine
codons['uuc'] = phenylalanine
codons['uua'] = leucine
codons['uug'] = leucine
for nt in nucleotides:
    codons['cu'+nt] = leucine
codons['uau'] = tyrosine
codons['uac'] = tyrosine
codons['ugu'] = cysteine
codons['ugc'] = cysteine
codons['ugg'] = tryptophan
for nt in nucleotides:
    codons['cc'+nt] = proline
codons['cau'] = histidine
codons['cac'] = histidine
codons['caa'] = glutamine
codons['cag'] = glutamine
codons['auu'] = isoleucine
codons['auc'] = isoleucine
codons['aua'] = isoleucine

stop = Codon('stop', '!', 'stop')
for c in stop_codons:
    codons[c] = stop
