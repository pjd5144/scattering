# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:35:22 2019

@author: pjdudenas
"""
#import matplotlib.pyplot as plt
from scattering import reduction as rd


image1 = rd('PAAGNP1_A0p160_2m.edf')
image1.load()
image1.raw_plot((12,10))
image1.geometry(1990,622,1320)
image1.plot((12,10))
