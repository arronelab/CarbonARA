#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:16:28 2022

@author: ccck67
"""
import sys

def ask_varying_sectionSecondary(test_name):
    fix_all = input('Would you like to vary all subsections of your structure? (Y/N)\n')
    if fix_all == 'Y':
        pass
    else:
        varying_sections_string = input('List the subsections you would like to vary. E.g. 1 2 3 5 8\n')
        varying_sections_list = varying_sections_string.split()
        with open(test_name+'varyingSectionSecondary.dat','w+') as f:
            for i in varying_sections_list[:-1]:
                f.write(i+'\n')
            f.write(varying_sections_list[-1])
            
if __name__ == '__main__':
    ask_varying_sectionSecondary(sys.argv[1])
    