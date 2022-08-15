#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 16:04:22 2021

@author: ccck67
"""
import os, re, sys
from os.path import exists
from Bio.PDB import *


def make_backbone(test_name):
    fw = open(test_name+'CA.pdb','w+')
    fr = open(test_name+os.path.basename(os.path.normpath(test_name))+'.pdb', 'r')
    for record in fr:
        if(re.search(r'^ATOM\s+\d+\s+CA\s+', record)):
            fw.write(record)
    fw.close()
    fr.close()

def get_coords(backbone_flname):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(backbone_flname[:-4],backbone_flname)
    coords=[]
    for model in structure:
        chains = model.get_chains()
        ids = []
        for chain in chains:
            ids.append(chain.get_id())
        chain = model[ids[0]]
        for residue in chain:
            for atom in residue:
                coords.append(atom.get_coord())
    return coords

def read_ss(fp_flname):
    fp = []
    with open(fp_flname,'r') as fin:
        for line in fin:
            fp.append(line)
    return fp[-1]

def count_ss(fp_flname):
    ss = read_ss(fp_flname)
    ss_count = []
    count = 1
    i = 0
    while i<len(ss)-1:
        if ss[i+1] == ss[i]:
            count += 1
            i += 1
        else:
            ss_count.append([ss[i], count])
            count = 1
            i += 1
    ss_count.append([ss[-1], count])
    return ss_count

def split_coords(backbone_flname, fp_flname):
    coords = get_coords(backbone_flname)
    ss_count = count_ss(fp_flname)
    split_coords = []
    idx=0
    for i in ss_count:
        split_coords.append(coords[idx:idx+i[1]])
        idx+=i[1]
    return split_coords

def write_coords(test_name):
    coords = split_coords(test_name+'CA.pdb', test_name+'fingerPrint.dat')
    with open(test_name+'coordinates.dat','w+') as fout:
        for subsection in coords:
            for point in subsection:
                fout.write(' '.join(map(str,point)))
                fout.write('\n')
            fout.write('\n')
        fout.write('End Chain 1')
    
if __name__ == '__main__':
    make_backbone(sys.argv[1])
    fp_file_exists = exists(sys.argv[1]+'fingerPrint.dat')
    if fp_file_exists:        
        write_coords(sys.argv[1])
    else:
        print('Missing fingerPrint.dat, run fingerPrint.py first')
    













































