# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 11:34:10 2018

@author: Lucas
"""

class airfoil:
    def __init__(self,r,rmax,i, alfa):
        self.i = i
        self.r = r
        self.rmax = rmax
        self.n = len(self.rmax)
        number= []
        for a in  range (len(rmax)):
            number.append(a)    
        foil = 'profile'
        ext = '.out'
        if r <= rmax[i]:
            self.Airfoil = foil+str(number[i])+ext
            self.alpha = alfa[i]
        else:
            self.i = i+1
            self.Airfoil = foil+str(number[i])+ext
            self.alpha = alfa[i]
        #print('Airfoil',i)
        #print('Angulo de Ataque: ',alfa[i])
