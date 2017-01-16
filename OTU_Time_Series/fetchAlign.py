#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 18:29:36 2016

@author: login
"""

import random, os, sys
random.seed(137)
from cogent import LoadSeqs, DNA, LoadTree
from cogent.db.ncbi import EFetch
from cogent.app.muscle import align_unaligned_seqs
from cogent.app.fasttree import build_tree_from_alignment

def fetch_ncbi_data(ofile,s):
    # get the seqs from Genbank
    input = [e.split() for e in s.strip().split('\n')]
    id_list = [t[0] for t in input]
    names = [t[1] for t in input]
    ef = EFetch(id=','.join(id_list), rettype='fasta')
    data = ef.read().strip()
    
    # title lines are too long, replace by genus_species
    rL = list()
    for i,e in enumerate(data.split('\n\n')):
        old_title, seq = e.strip().split('\n',1)
        new_title = '>' + names[i]
        seq = seq[:500]
        rL.append('\n'.join([new_title,seq]))
    FH = open(ofile,'w')
    FH.write('\n\n'.join(rL))
    FH.close()

def mutagenize(seq, mrate=1):
    L = list(seq)
    D = { 'A':'CGT', 'C':'AGT', 'G':'ACT', 'T':'ACG' }
    N = int(mrate / 100.0 * len(seq))
    X = len(seq)
    for i in range(N):
        j = random.choice(range(X))
        nt = L[j]
        if not nt in 'ACGT':  continue
        L[j] = random.choice(D[nt])
    return ''.join(L)

def distribute_seqs(ifile,ofile):
    # set up our samples
    FH = open(ifile,'r')
    data = FH.read().strip().split('\n\n')
    FH.close()
    seqs = list()
    for e in data:
        title,seq = e.split('\n',1)
        seqs.append(''.join(seq.split()))
    
    outgroup = '>Thermotoga\n' + seqs.pop()
    
    A = {0:5,1:5,2:0,3:1,4:0,5:1,6:1,7:1}  # A has lots of Firmicutes
    B = {0:0,1:1,2:5,3:5,4:1,5:0,6:1,7:1}  # B has Bacteroidetes
    C = {0:1,1:0,2:1,3:0,4:5,5:5,6:1,7:1}  # C has enterics
    dL = [A,B,C]
    L = list()
    
    for distr, sample in zip(dL,list('ABC')):
        counter = 1
        for k in distr:
            seq = seqs[k]
            n = distr[k]
            for i in range(n):
                if n == 1:  mrate = 5
                else:       mrate = random.choice((1,2,3))
                copy = mutagenize(seq[:],mrate)
                name = sample + str(counter)
                L.append(DNA.makeSequence(copy,name))
                counter += 1
    FH = open(ofile,'w')
    L = [seqs_.toFasta() for seqs_ in L]
    L.insert(0,outgroup)
    FH.write('\n\n'.join(L))
    FH.close()

def align_seqs(ifile,ofile):
    seqs = LoadSeqs(ifile, moltype=DNA, aligned=False)
    aln = align_unaligned_seqs(seqs, DNA)
    aln.writeToFile(ofile)
    return aln

def get_tree(ifile):
    aln = LoadSeqs(ifile, moltype=DNA, aligned=True)
    tr = build_tree_from_alignment(aln,moltype=DNA)
    return tr

#===============================================
s = '''
AY005045.1     Streptococcus_mitis_bv2
D83363.1       Staphylococcus_epidermidis_14990
L14639.1       Capnocytophaga_gingivalis
AB053940.1     Tannerella_forsythensis_HA3
EU009197.1     Shigella_sonnei_FBD023
AB435616.1     Serratia_marcescens_JCM24201
AB302401.1     Pseudomonas_cinnamophila
AF411020.1     Achromobacter_xylosoxidans_AU1011
AJ401017.1     Thermotoga_maritima_SL7
'''

fn1 = 'rRNA_gb.fasta'
fn2 = 'samples.fasta'
fn3 = 'samples.aln.fasta'
fn4 = 'samples.tree'

if not os.path.exists(fn1) or False:  
    fetch_ncbi_data(fn1,s)
if not os.path.exists(fn2) or False:  
    distribute_seqs(fn1,fn2)
    
if not os.path.exists(fn3) or True:  
    aln = align_seqs(fn2,fn3)
    tr = get_tree(fn3)
    # re-root manually
    print tr.asciiArt()
    n = tr.getNodeMatchingName('Thermotoga')
    for a in n.ancestors():  print a.Name
    
    tr2 = tr.rootedAt(n.ancestors()[0].Name)
    tree_str = tr2.getNewick(with_distances=True)
    FH = open(fn4,'w')
    FH.write(tree_str + '\n')
    FH.close()

tr = LoadTree(fn4)
print tr.asciiArt()