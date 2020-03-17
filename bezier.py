# -*- coding: utf-8 -*-
"""
Created on Wed May 30 23:48:57 2018

@author: Lucas
"""
import math as m
import numpy as np
import os
import matplotlib.pyplot as plt# importa a biblioteca matplotlib
def fat(n):
     if n == 0:
         return 1
     else:
         return n * fat(n-1)


def bezier_profiles(cpx,cpy,order,npts):

#
#  Some input data for the stretching function...
    nocp      = order + 1
    #npts      = 60
    stre_coef = 95.0
    dx        = 1.0 / npts
    etau      = list()
    x         = list()
#

#  Generating the profile bunching  - ATANH  distribution
#
    for i in range(npts+1):
        x.append(float(0.0 + i * dx))
        stre  = 1.0 + m.tanh(stre_coef*(x[i]-1.0)*m.pi/180.0) / m.tanh(stre_coef*m.pi/180.0)
        etau.append(float(stre+0.0))

#
# Defining the coefficients of the PASCAL TRIANGLE.....
    lcni = list()
    for i in range(0, order+1, 1):
        aux = fat(order)/(fat(i)*fat(order-i))
        lcni.append(int(aux))
    duaux = [[0. for i in range(npts+1)] for k in range(order+1)]
    for i in range(0,npts+1,1):
        for k in range(0, order+1,1):
            duaux[int(k)][int(i)] = lcni[k]*(m.pow(etau[i],k))*(m.pow((1-etau[i]),(order-k)))
#
#   Here I am putting the dul array into other arrays using the numpy lib.
    du            = np.matrix(duaux).T
    cpx2          = np.reshape(cpx, (nocp,1))
    cpy2          = np.reshape(cpy, (nocp,1))
    x      = du*cpx2#mudei x_bezier
    y      = du*cpy2#y_bezier

    
    #y = list()
    '''
   f1 = open('profile.out','w')
    for i in reversed(range(0,npts+1)):                                # Upper Surface
        f1.write('%9.6f %9.6f\n' % (x_bezier[i],y_bezier[i]))
        x.append(float(x_bezier[i]))
        y.append(float(y_bezier[i]))

    for i in range(1,npts+1):                                   # Lower Surface
        f1.write('%9.6f %9.6f\n' % (x_bezier[i],-y_bezier[i]))
        x.append(float(x_bezier[i]))
        y.append(float(-y_bezier[i]))
        
    f1.close()
    '''
    return x,y
'''
cpx=[0.1,.4,.9,1.25,1.5]

cpy = [0.1,0.2,.25,0.18,0]
order = 4
npts = 20
a = bezier_profiles(cpx,cpy,order,npts)
plt.plot(a[0],a[1])
print(a[0])
print(len(a[0]),
          a[1])
plt.plot(cpx,cpy,'ro')
r = a[0]
c= a[1]
from xlwt import Workbook
wb = Workbook()
Sheet1 = wb.add_sheet('Sheet1')
for i in range(0,len(r)): 
       Sheet1.write(i+1,0,r[i])
       Sheet1.write(i+1,1,c[i])
wb.save('koch4.xls')
os.startfile('koch4.xls')
'''