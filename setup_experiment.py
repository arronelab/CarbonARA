#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 15:18:34 2022

@author: ccck67
"""
import sys
import os
from os.path import exists

fp_file_exists = exists(sys.argv[1]+'fingerPrint1.dat')
if not fp_file_exists:
    print('Missing fingerPrint.dat, run fingerPrint.py first')
coord_file_exists = exists(sys.argv[1]+'coordinates1.dat')
if not coord_file_exists:
    print('Missing coordinates.dat, run coordinates.py first')
fix_all = input('Are you varying all subsections of your structure? (Y/N)\n')
is_mix = input('Is your SAXS data for one singular structure? (Y/N)\n')

def write_sh_file(test_name,project_name):
    output_filepath = os.path.basename(os.path.normpath(test_name))+'_'+project_name+'.sh'
    print(output_filepath)
    with open(output_filepath, 'w+') as fout:
        fout.write('#!/bin/bash')
        fout.write('\nScatterFile={t_n}{bt_n}Saxs.dat'.format(t_n=test_name,bt_n=os.path.basename(os.path.normpath(test_name))))
        fout.write('\nfileLocs={t_n}'.format(t_n=test_name))
        fout.write('\ninitialCoordsFile={}coordinates.dat'.format(test_name))
        if is_mix=='Y':
            fout.write('\nnoStructures=1')
        else:
            return('This functionality is still in progress')
        fout.write('\npairedPredictions=none')
        if fix_all == 'N':
            fout.write('\nfixedsections={}/varyingSectionSecondary.dat'.format(test_name))
        else:
            fout.write('\nfixedsections=none')
        # no options for crystal symmetry or or hydro cover in prototype
        fout.write('\ncrystalSymmetry=none')
        fout.write('\nwithinMonomerHydroCover=none')
        fout.write('\nbetweenMonomerHydroCover=none')
        fout.write('\nkmin=0.01;')
        fout.write('\nkmax=0.1;')
        fout.write('\nmaxNoFitSteps=10')
        fout.write('\n\nif [ ! -d {t_n}{p_n} ]; then'.format(t_n=test_name,p_n=project_name))
        fout.write('\n  mkdir -p {t_n}{p_n};'.format(t_n=test_name,p_n=project_name))
        fout.write('\nfi')
        fout.write('\n\nfor i in {1..1}')
        fout.write('\ndo')
        fout.write('\n    predictStructure $ScatterFile $fileLocs')
        fout.write(' $initialCoordsFile $pairedPredictions $fixedsections $noStructures')
        fout.write(' $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax')
        fout.write(' $maxNoFitSteps {t_n}{p_n}/mol$i'.format(t_n=test_name,p_n=project_name))
        fout.write(' {t_n}{p_n}/scatter$i.dat'.format(t_n=test_name,p_n=project_name))
        fout.write(' {t_n}mixtureFile.dat'.format(t_n=test_name))
        fout.write('\ndone')

if __name__ == "__main__":
    write_sh_file(sys.argv[1],sys.argv[2])
