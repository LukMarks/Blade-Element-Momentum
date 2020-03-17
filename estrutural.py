# Author: Lucas Marques Moreno
# Date: 9/2018 
# Description: Class designed
# to performe stress calculations

import numpy as np
import os 
import matplotlib.pyplot as plt

def xt(airfoil,c):
        x = []
        y = []
        yu = []
        xu = []
        yl =[]
        xl = []
        prof_name = airfoil
        f=open(prof_name,'r')

        Tks = 0
        for line in f:
            X = float(line.strip().split()[0])*c
            Y = float(line.strip().split()[1])*c
            x.append(X)
            y.append(Y)



        for i in range(1,len(x)):
            if x[i]<x[i-1]:
                cp = i+1

        for i in range(0,cp):
            xu.append(x[i])
            yu.append(y[i])

        y.reverse()
        for i in range(cp, len(x)):
            xl.append(x[i])
            yl.append(y[i])

        for i in range(0,len(yl)):
                for j in range(0,len(yu)):
                    tks = abs(yu[j]-yl[i])
                    if Tks<tks:
                        Tks=tks
                        ref1 = y[i]
                        ref2 = y[j]
                        x_t = xl[i]
                        #print(ref1,ref2)
        '''
        for g in range (0,len(x)):
            if y[g] == ref1:
                x_t = xl[g]
            elif y[g] == ref2:
                x_t = x[g]
        '''
        
        return Tks,x_t,x,y
def cone (Di,Df,esp,L,rl):
    
    b = Di/2
    a = (Df-Di)/(2*L)
    
    r = a*rl+b
    
    
    de = 2*r
    di = de-2*esp
    return de, di 
    

class struct:
    def __init__(self, torque, thrust,r,c,m_s,omega,profile,Di,mode):
        self.c = c 
        self.m_s = m_s
        self.profile = profile
        self.omega = omega
        self.torque = torque
        self.thrust = thrust
        self.mode = mode
        self.Di = Di
        self.rs=[]
        #r.reverse()
        if self.mode == 'tubo':
            d_list = []
            for j in range(0,len(c)-1):
                D = xt(profile,c[j])
                d_list.append(D[0])
            De = (10)*1e-3#min(d_list)*0.85
            print('Maximum external diameter: ',round(De*1000,2),'mm')
        self.normal = []
        self.mflet = []
        r2 = r
        #r2.reverse()
        #m_s.reverse()
        #thrust.reverse()
        self.d = []
        for i in range(0,len(r)):
    
            
            ace_rad = r2[i]*omega**2
            weight = m_s[i]*ace_rad            
            if self.mode == 'tubo':
                De,Di = cone(10e-3,10e-3,1e-3,1.42,r[i])
                self.d.append(De*1e3)
                A_sec = (np.pi/4)*((De**2)-(Di**2))
                Normal = weight/A_sec
                self.normal.append(Normal)
                I = (np.pi/64)*((De**4)-(Di**4))
                M_thrust = thrust[j]*r2[j]
                M_flet = (Normal+((torque[j]/I)*De/2)+((M_thrust/I)*De/2))
                self.mflet.append(round((M_flet/1e6),2))
        #r.reverse()
               
            
            