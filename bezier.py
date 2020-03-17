
import math as m
import numpy as np
import os


import matplotlib.pyplot as plt
def fat(n):
     if n == 0:
         return 1
     else:
         return n * fat(n-1)


def bezier_profiles(cpx,cpy,order,npts):

    nocp      = order + 1
    stre_coef = 95.0
    dx        = 1.0 / npts
    etau      = list()
    x         = list()


    for i in range(npts+1):
        x.append(float(0.0 + i * dx))
        stre  = 1.0 + m.tanh(stre_coef*(x[i]-1.0)*m.pi/180.0) / m.tanh(stre_coef*m.pi/180.0)
        etau.append(float(stre+0.0))


    lcni = list()
    for i in range(0, order+1, 1):
        aux = fat(order)/(fat(i)*fat(order-i))
        lcni.append(int(aux))
    duaux = [[0. for i in range(npts+1)] for k in range(order+1)]
    for i in range(0,npts+1,1):
        for k in range(0, order+1,1):
            duaux[int(k)][int(i)] = lcni[k]*(m.pow(etau[i],k))*(m.pow((1-etau[i]),(order-k)))
#
    du            = np.matrix(duaux).T
    cpx2          = np.reshape(cpx, (nocp,1))
    cpy2          = np.reshape(cpy, (nocp,1))
    x      = du*cpx2
    y      = du*cpy2


    return x,y
