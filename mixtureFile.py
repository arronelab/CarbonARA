#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 15:28:57 2022

@author: ccck67
"""
import numpy as np
import sys

def ask_is_mix(test_name):
    is_mix = input('Is your SAXS data for one singular structure? (Y/N)\n')
    if is_mix == 'Y':
        with open(test_name+'mixtureFile.dat','w+') as f:
            f.write('1.0')
    else:
        with open(test_name+'mixtureFile.dat','w+') as f:
            for i in np.arange(0.0,1.0,0.1):
                f.write(str(np.round(1-i,2))+ ' ' + str(np.round(i,2))+'\n')
            f.write(str(0.0) + ' ' + str(1.0))
            
if __name__ == '__main__':
    ask_is_mix(sys.argv[1]) 
