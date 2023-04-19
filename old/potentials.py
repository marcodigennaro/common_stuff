# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:56:54 2021

@author: MDI0316
"""

import numpy as np

def morse(r, eps=1, alp=1, r_eq=1):
    return -eps*( 1- (1-np.exp( -alp*(r-r_eq) ) )**2 )

def lennard_jones(r, eps=1, sig=1):
    return 4*eps*( (sig/r)**12 - (sig/r)**6 )

