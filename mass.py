# -*- coding: utf-8 -*-
"""
Created on Mon May 28 09:28:17 2018

@author: Lucas
"""
def inte(x,y):
            n = len(x)-1#ver
            A = 0
            for i in range(0,n):
                multi = x[i+1]-x[i]
                a = abs(multi*(y[i]+y[i+1])/2)
                A = A+a
            return A
def S2(x,y):
    n = len(x)-1
    A=0
    for i in range(0,n):
       multi = (x[i+1]-x[i])/3
       if i ==0:
           A = A + abs(multi*y[i])
       elif i % 2 == 0:
           A = A +abs(2*multi*y[i])
       elif i == n:
           A = A + abs(y[n])
       elif i % 2 != 0:
           A = A + abs(4*multi*y[i])
    return A

def xt(airfoil,c,esp,r):
        xi = []
        yi = []
        xo = []
        yo = []
        
        es = esp/100
        prof_name = airfoil#+'.dat'
        f=open(prof_name,'r')
        #verificar a ordem do ponto lendo de x
        for line in f:
            X = float(line.strip().split()[0])
            Y = float(line.strip().split()[1])
            Xo = X*c
            Yo = Y*c
            Xi = X*c*(1-es)
            Yi = Y*c*(1-es)
            xi.append(Xi)#xrot)
            yi.append(Yi)#yrot)
            Ai = inte(xi,yi)
            #Ai = S2(xi,yi)
            xo.append(Xo)
            yo.append(Yo)               
            Ao =inte(xo,yo)
            #Ao = S2(xo,yo)
            A = Ao-Ai
            v = A*r
        return v #xi,yi,xo,yo

class mass:
    def __init__ (self,airfoil,c,esp,r):
        self.airfoil = airfoil
        self.c = c
        self.r = r
        '''
        '''
        Vol = xt( airfoil, c, esp, r)      
        self.V = Vol
            
       
        
