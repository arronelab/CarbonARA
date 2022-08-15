#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:22:41 2022

@author: ccck67
"""
import sys
import os


simplify_dict = {'H': 'H',
                 'B': 'S',
                 'E': 'S',
                 'G': 'H',
                 'I': 'H',
                 'T': '-',
                 'S': '-',
                 '-': '-',
                 ' ': '-'
                 }

def make_string(lst):
    return ''.join(lst)

def read_fasta_file(fasta_flname):
    lines = []
    with open(fasta_flname) as input_data:
        for line in input_data:
            lines.append(line.strip())
    return lines[1]
        

def read_dssp_file(dssp_flname):
    lines=[]
    with open(dssp_flname) as input_data:
        # Skips text before the beginning of the interesting block:
        for line in input_data:
            if line.strip() == '#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA':  # Or whatever test is needed
                break
        # Reads text until the end of the block:
        for line in input_data:  # This keeps reading the file
            lines.append(simplify_dict[line[16]])
    return make_string(lines)


def make_fingerprint(test_name):
    dssp = read_dssp_file(test_name+os.path.basename(os.path.normpath(test_name))+'_dssp.txt')
    fasta = read_fasta_file(test_name+os.path.basename(os.path.normpath(test_name))+'.fasta')
    assert len(dssp) == len(fasta)
    with open(test_name+'fingerPrint.dat','w+') as fout:
        fout.write('1\n\n')
        fout.write(fasta+'\n\n')
        fout.write(dssp)

if __name__ == '__main__':
    make_fingerprint(sys.argv[1])