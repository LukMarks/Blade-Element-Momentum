# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 16:34:56 2018

@author: 14.02236-2
"""
#from xlwt import Workbook
import numpy as np
import os
#import xlrd
import matplotlib.pyplot as plt

def x_t(airfoil):
        x = []
        y = []
        yu = []
        xu = []
        yl =[]
        xl = []
        prof_name = airfoil
        f=open(prof_name,'r')
        #verificar a ordem do ponto lendo de x
        Tks = 0
        for line in f:
            X = float(line.strip().split()[0])
            Y = float(line.strip().split()[1])
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
        
        #print('thickness: ',round(100*Tks,1),'%')
        #print(x_t)
        return tks, x_t

def xt(airfoil,c,phi,x_t,tks,i):
        x = []
        y = []
        yu = []
        xu = []
        yl =[]
        xl = []           
        prof_name = airfoil
        f=open(prof_name,'r')
        if i == 0:
            centerx = c*.25
            Z_local = centerx
        else:
            centerx = 0
            Z_local = 0.25*c-centerx
        for line in f:
            X = (float(line.strip().split()[0]))*c-Z_local
            Y = float(line.strip().split()[1])*c#-Z_local*.1
         
            xrot = X*np.cos(phi)+Y*np.sin(phi)#-Z_local*np.cos(phi) #x_t*.25*c*np.cos(phi)
            yrot = Y*np.cos(phi)-X*np.sin(phi)#-Z_local*np.sin(phi)+tks*.25*c#*np.sin(phi)
            x.append(xrot)
            y.append(yrot)
        return x,y
class Solid:
    def __init__(self,c,perfil,phi,r,i):
        
        self.c = c
        self.perfil = perfil
        self.phi = phi
        self.r=r
        self.i = i
        
        tks,spc= x_t(self.perfil)
        G = xt(self.perfil,self.c,self.phi,spc,tks,self.i)
        self.x = G[0]
        self.y = G[1]
        n = self.r*1000
        nome = 'R'+str(int(n))+'.xls'
        '''
        wb = Workbook()
        Sheet1 = wb.add_sheet('Sheet1')
        #Sheet1.write(0,0,'X')  
        #Sheet1.write(0,1,'Y') 
        #Sheet1.write(0,2,'r')
        for i in range(0,len(self.x)):
           Sheet1.write(i,0,self.x[i])
           Sheet1.write(i,1,self.y[i])
           Sheet1.write(i,2,self.r)
        wb.save(nome)
        '''
        nome2 = 'R'+str(int(n))+'.txt'
        
        self.title = nome2
        f=open(nome2,'w')
        for i in range(0,len(self.x)):
         f.write('%s \t' % (str(self.x[i])))
         f.write('%s \t' % (str(self.y[i])))
         f.write('%s \n' % (str(self.r)))
        f.close()
