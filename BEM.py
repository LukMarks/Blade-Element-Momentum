# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 16:21:08 2018

@author: Lucas
"""
# BEM
import numpy as np# importa a biblioteca numpy
import matplotlib.pyplot as plt# importa a biblioteca matplotlib
import math as m
import os
import xlrd
import time
from corda import corda
from airfoil import airfoil
from ve import ve
from xlwt import Workbook
from xfoil import xfoil
from mass import mass, S2
from solid import Solid
from estrutural import struct
from bezier import bezier_profiles
from pot import Pot
#=====================path Xfoil=================
start_time = time.time()
global path
path  ='.\\xfoil_output'
#============== Resultados em Excel ==============
Excel =False
################ Comando da Correção #################
Correcao = True # Ativa a correção determinada
############# Metodos de correção #######
Hansen = True
#Liebeck = False
############ Chave para ativação do Xfoil ##############
Xfoil_use = True  
############### Configuração do Xfoil###############
npp = 220 # numero de paineis
inter =200# numero de interações do Xfoil
#=============== Modos de Uso ==============
MultiFoil = True
#N_perfil = 0
pos = 0
#R_max = [.4,.42,.85,1.05,1.3]#[m] raios de troca de perfil
R_max =[.4,.42,.85,.95,1.5] ##[.4,.42,.85,1.05,1.6]#[m] raios de troca de perfil
print(len(R_max))
#Ang_at = [5,5.5,7,6.5,6]#[°] angulo de ataque do perfil
Ang_at =[3.,3.,8.5,8.6,9.]#[0.,3.316,6.471,6.,6.]#[6.5,7,6.5,6.5,6.5]#[°] angulo de ataque do perfil
#============ Corda ========================
Curva_Bezier = False
plot_bezier = False

#tri - Triangular
Mode = 'trapz'# modo de distribuição de corda
c_pa = 0.25#[m] corda constante da pá
pcorda = [0.2,1.25]#[0.4, 1.25]#[m] trecho de corda constante 
esp = 100 #[%] Espessura da parede
den = 420 #[kg/m3] |Nylon (poliamida)| massa especifica do material
#=========== Plot =========================
plot = True
ss = True
#escreve pontos das secções
#============== Estudo de Eficiencia ===============
EE = False # habilita o calcuclo para varias velocidades
V_inicial =7.#[m/s] Velocidade inicial
V_final = 18.#[m/s] velocidade final
passo_v = 1.# passo para Variação de velocidades
#================ Estudo Estrutural ====================
struc = True
#================ Entrada do Excel ======================
cpx = [0.1,.5,.75,1.2,1.6]#[0.1,.4,.75,.9,1.6]#coordenadas dos pontos Em X
order = len(cpx)-1
npts = 20
cpy = [.15,0.24999986,0.34668816,0.34996463,.1]#[.1,.25,.35,.25,.01]#[0.1,.25,.31,.25,.01]
prop = bezier_profiles(cpx,cpy,order,npts)
r = []#[m]lista do raio da hélice
c=[]#[m] lista com os comprimentos da corda
for i in range(0,len(prop[0])):
    r.append(float(prop[0][i]))
    c.append(float(prop[1][i]))
#r = prop[0]
#c = prop[1]


#################### Dados Entrada ########################
#-------------------- Dados Natuais -----------------------
g = 9.81#[m/s²] aceleração da gravidade
p = 1.08#[kg/m³] massa especifica do ar(1.225)
u = 1.4607*10**(-5)#[N.s/m²] viscosidade dinamica do a 18.20*10**(-6) #
#-------------------- Dados da Hélice -------------------
AeroFolio='naca2412.dat'
n_pontos = 20
passo = 0.01#[m] Distancia de divisão de interação
Vo =12.#[m/s]
rpm = 94.73845432#112.#[rpm] rotação da hélice
c_inicial =0.15#[m]corda inicial
cmin = 0.#[m]corda inicial corda minima
Ac = 0.1#[]Constante de inicio correção
D = 2*cpx[len(cpx)-1]#[m] Diametro da hélice
print(D)
D_hub = 0.2#[m] Diametro do Hub
alfaD = 5#[º]Angulo de ataque
B = 2.#numero de pás
Re_min = 1.5e3#[] Minimo numero de Reynolds para convergencia

if Xfoil_use == False:
    Cl=1.#[] Coeficiente de sustentação
    Cd=0.1#[] coeficiente de de arrasto


############### Calculo do Fator de indução ############
theta=[]#[º] Lista do angulo de torção da pá
efi = []#[] lista das eficiencias
pt=[]#[N] lista para força tangenciais
pn=[]#[N] lista para força Normais
vs=[]#[m/s] velocidade na componente seno de phi
vc=[]#[m/s] velocidade na componente cosseno de phi
v_abs=[]#[m/s] velocidade corrigida
v_rel = []#[m/s] velocidade não corrigida
reynolds = []
mach=[]
CL = []#[] lista com os coeficientes de sustentação
CD = []#[] lista com os coeficientes de arrasto
m1 = []
t =[]
ratio =[]# lista para relação CL/Cd
ratio32 = []# lista para relação CL^3/2/CD
a = [0]
a_l = [0]
w = rpm*2*np.pi/60#[rad/s] velocidade angular
rps = rpm/60#[rot/s] Frequencia angular
Phi = []

Jo = Vo/(D*rps)
cy = 0
print('Coeficiente de Avanço:',round(Jo,2),'[m/rot]')
torque = 0#[N] Valor inicial da força radial
Thrust = 0#[N] Valor inicial da força axial
Ri = D_hub/2
R_hub = D_hub/2
Rf = D/2
r_foil = []
if Curva_Bezier == True:
    order = len(cpx)-1
    npts = int((Rf-Ri)/passo)
    print(npts)
    Corda = bezier_profiles(cpx,cpy,order,npts)
    cy_min = int((R_min-Ri)/passo)
PnT = 0
PtT = 0
V_som = 343#[m/s] velocidade do som
reuso=0
R_max.append(Rf)
print('')
print('INICIANDO CALCULO DO FATOR DE INDUÇÃO')
if MultiFoil == False:
    alphaD = alfaD
s=[]
Vol = 0
prof = []
if Correcao == True:
    ac  = Ac
    #while Ri<=Rf:
    m_s = []
    j = 0
    for j in range(0,len(r)):
        R = r[j] #Ri
        if MultiFoil == True:
                    perfil = airfoil(R,R_max,pos,Ang_at)
                    pos = perfil.i
                    AeroFolio = perfil.Airfoil
                    r_foil.append(AeroFolio)
                    alphaD = perfil.alpha

        alpha = alphaD*np.pi/180
        #A = a.pop()
        A = a[len(a)-1]
        #a.append(A)
        #A_l = a_l.pop()
        A_l = a_l[len(a_l)-1]#a_l.pop()
        #a_l.append(A_l)
        Vt = w*R
        Tan_phi = ((1+A)*Vo)/((1-A_l)*Vt)
        phi = np.arctan(Tan_phi)
        Phi.append(phi)
        Theta = phi-alpha
        
        theta.append(Theta*180/np.pi)
        Vrel = (Vo**2+(Vt)**2)**(1/2)
        v_rel.append(Vrel)
        if Curva_Bezier == True:
                if R<R_min:
                    C = cmin
                   
                    cy = cy_min
                else:
                    C = Corda[cy]
                    print(C)
                    cy = cy+1
                    

        #if Curva_Bezier == False:
            #Corda=corda(R,0,pcorda[0],pcorda[1],Rf,c_pa,Mode)
        C = c[j]
        VrelS = Vo*(1+A)
        VrelC = Vt*(1-A_l)
        vs.append(VrelS)
        vc.append(VrelC)
        V_abs = np.sqrt(VrelS**2+VrelC**2)
        v_abs.append(V_abs)
        Re = p*V_abs*C/u
        reynolds.append(Re)
        Mach = V_abs/V_som#[] numero de Mach
        mach.append(Mach)
        #S = S2(r,c)# C*np.cos(alpha)*R
        #s.append(S)
        J = Vo/Vt
        f = (B/2)*(1/J)*(1+J**2)**(1/2)*(1-(R/Rf))
        #F = (2/np.pi)*(np.arccos(np.exp(-f)))**(-1)
        F = (2/np.pi)*(np.arctan((np.exp(2*f)-1)**(1/2)))
        #=================
        volu = mass(AeroFolio, C, esp, passo)
        vol = volu.V
        Vol = Vol+vol
        m_s.append(Vol)
        #===================

        if Xfoil_use == True:
            if Re>=Re_min:

                Foil=xfoil(AeroFolio,alphaD,round(Re,0),round(Mach,0),inter,npp,path,CL,CD)
                prof.append(AeroFolio)                 
                Cl=Foil[0]
                CL.append(Cl)
                Cd=Foil[1]
                CD.append(Cd)
                reuso = reuso + Foil[2]
                print('reuso:',reuso)
                Ratio = Cl/Cd
                ratio.append(Ratio)
                Ratio32 = Cl**(3/2)/Cd
                ratio32.append(Ratio32)
            else:
                Cl=0.0001#CL[len(CL)-1]
                CL.append(Cl)
                Cd=0.00001#CD[len(CD)-1]
                CD.append(Cd)
                Ratio = Cl/Cd
                ratio.append(Ratio)
                Ratio32 = Cl**(3/2)/Cd
                ratio32.append(Ratio32)

        L = (1/2)*p*C*Cl*Vrel**2
        Dr = (1/2)*p*C*Cd*Vrel**2
        Cn = Cl*np.cos(phi)+Cd*np.sin(phi)
        Ct = Cl*np.sin(phi)-Cd*np.cos(phi)
        sigma = C*B/(2*np.pi*R)
        I1 = 4*np.sin(phi)**2
        I2 = sigma*Cn
        A = 1/((I1/I2)-1)
        a.append(A)
        I3 = 4*np.sin(phi)*np.cos(phi)
        I4 = sigma*Ct
        A_l = 1/((I3/I4)+1)
        a_l.append(A_l)
        if Hansen == True:

            if A>ac:
                K_h = 4*F*np.sin(phi)**2/(sigma*Cn)
                A = (1/2)*(2+K_h*(1-2*ac)-((K_h*(1-2*ac)+2)**2+4*(K_h*ac**2-1))**(1/2))
                a.pop()
                a.append(A)
            
            if A<=(1/3):
                Ct = 4*A*(1-A)*F
            else:
                Ct = 4*A*(1-(1/4)*(5-3*A)*A)*F
            
        '''        
        elif Liebeck == True:
            
             if A>ac:
                   k_l = Cn/(4*np.sin(phi)**2)
                   A = sigma*k_l/(F-sigma*k_l)
                   a.pop()
                   a.append(A)
                   kl_l = Ct/(4*np.cos(phi)*np.sin(phi))
                   A_l = sigma*kl_l/(F+sigma*kl_l)
                   a_l.pop()
                   a_l.append(A_l)


        '''
        #Pn = (1/2)*p*S*Cn*Vrel**2
        #Pt = (1/2)*p*S*Ct*Vrel**2
        
        Pn = (1/2)*p*C*Cn*Vrel**2 #mudar isso
        Pt = (1/2)*p*C*Ct*Vrel**2
        
        pn.append(Pn)
        pt.append(Pt)
        if ss == True:
            SS = Solid(C,AeroFolio,phi,R,j)
            plt.figure(1)
            plt.plot(SS.x,SS.y)
#           plt.title(SS.title)
        progresso = 100*j/len(r)
        print('Progresso: ',round(progresso,2),'%')
    

if Correcao == False:
    while Ri<=Rf:
        R = Ri
        r.append(R)
        A = a.pop()
        a.append(A)
        A_l = a_l.pop()
        a_l.append(A_l)
        Vt = w*R
        Tan_phi = ((1-A)*Vo)/((1-A_l)*Vt)
        Vrel = (Vo**2+(Vt)**2)**(1/2)
        v_rel.append(Vrel)
        Corda=corda(Kc,Kq,Kl,c_inicial,R,cmin,R_min)
        C = Corda.c
        c.append(C)
        S = R*C*np.cos(alpha)
        J = Vo/Vt
        f = (B/2)*(1/J)*(1+J**2)**(1/2)*(1-(R/Rf))
        F = (2/np.pi)*(np.cos(np.exp(-f)))**(-1)
        L = (1/2)*p*S*Cl*Vrel**2
        D = (1/2)*p*S*Cd*Vrel**2
        phi = np.arctan(Tan_phi)
        Cn = Cl*np.cos(phi)+Cd*np.sin(phi)
        Ct = Cl*np.sin(phi)-Cd*np.cos(phi)
        sigma = C*B/(2*np.pi*R)
        I1 = 4*np.sin(phi)**2
        I2 = sigma*Cn
        A = 1/((I1/I2)-1)
        a.append(A)
        I3 = 4*np.sin(phi)*np.cos(phi)
        I4 = sigma*Ct
        A_l = 1/((I3/I4)+1)
        a_l.append(A_l)

        Pn = (1/2)*p*S*Cn*Vrel**2
        Pt = (1/2)*p*S*Ct*Vrel**2
        pn.append(Pn)
        pt.append(Pt)
        VrelS = Vo*(1+A)
        VrelC = Vt*(1-A_l)
        vs.append(VrelS)
        vc.append(VrelC)
        V_abs = (VrelS**2+VrelC**2)**(1/2)
        v_abs.append(V_abs)
        Theta = phi-alpha
        theta.append(Theta*180/np.pi)
        Rey = p*Vrel*C/u
        reynolds.append(Rey)
        progresso = (Ri/Rf)*100
        print('Progresso: ',round(progresso,2),'%')
        Ri = Ri+passo

print('')
print('CALCULO DO FATOR DE INDUÇÃO COMPLETO')
#==========Desenho da Hélise==========
hub_x=[]
hub_y=[]
fi = len(r)-1
l_esq=[0,c[0]]
r_esq=[r[0],r[0]]
l_dir=[0,c[len(c)-1]]
r_dir=[r[fi],r[fi]]


#####
EF = 0
Ct = 0
Cm = 0
i=0
n = len(r)-1
##############Calculo das Forças#####################
rn = []
Pot=0
print('')
print('INICIANDO CALCULO DAS FORÇAS')
print('pni',len(pn))
while i <n:
    Rn=r[i]
    rn.append(Rn)
    Pn = pn[i]
    Pni = pn[i+1]
    Pt = pt[i]
    Pti = pt[i+1]

    R = r[i]
    Ri = r[i+1]
    ######CALCULO DA FORÇA AXIAL#####
    Yn = (Pni-Pn)/(Ri-R)
    Sn = (Pn*Ri-Pni*R)/(Ri-R)
    PN = Yn*R + Sn

    T = (1/2)* Yn*(Ri**2-R**2)+Sn*(Ri-R)
    t.append(T)
    Thrust = (Thrust+T)
    

    ######CALCULO DA FORÇA RADIAL#####

    Yt = (Pti-Pt)/(Ri-R)
    St = (Pt*Ri-Pti*R)/(Ri-R)
    PT = Yt*R + St
    M = (1/3)* Yt*(Ri**3-R**3)+(1/2)*St*(Ri**2-R**2)
    m1.append(M)
    torque = (torque+M)
    #print(torque)
    progresso = (i/(1+n))*100
    print('Progresso: ',round(progresso,2),'%')
    #i=i+1



    ######CALCULO DA EFICIENCIA#####
    '''
    if i >= 1:
        POT = torque*w
        ct = Thrust/(p*v_abs[i]**2*(r[i])**4)
        cm = POT/(p*v_abs[i]**3*(r[i])**5)
        dV = -Vo+((Vo**2 + 4*0.5*Thrust*(1/(p*S)))**(1/2))
        E = 1 +(dV/(2*Vo))
        EFI =1/E*100
        efi.append(EFI)
    else:
        efi.append(0)
   '''
    i=i+1

F_axial = B*Thrust
Momento = B*torque
Pot_voo = F_axial*Vo

Area = S2(r,c)


POT = Momento*w
dV = -Vo+((Vo**2 + 4*0.5*F_axial*(1/(p*B*Area)))**(1/2))
E = 1 +(dV/(2*Vo))
EFI =1/E*100





massa = den*Vol*B

a.reverse()
a.pop()
a.reverse()

a_l.reverse()
a_l.pop()
a_l.reverse()


#CD.reverse()


print('')
print('CALCULO DAS FORÇAS COMPLETO')


print('')
print('Estudo de Eficiencia')
if EE == True:
    vv=[]
    te=[]
    ee=[]
    me=[]
    pe=[]
    print(r,Phi,c)
    while V_inicial<=V_final:
        V = V_inicial
        vv.append(V)
        PP = Pot(r,Phi,c,V_inicial)
        pe.appen(PP)
        #Est = ve(a,a_l,Phi,r,V,c,Area,Xfoil_use,r_foil)
        #te.append(Est.Tr)
        #ee.append(Est.efi)
        #me.append(Est.m)
        #pe.append(Est.p)
        V_inicial = V_inicial + passo_v 
       
        plt.plot(vv,pe)

###########Resultados########################
#print(len(r),len(Phi))
end_time = time.time()
tempo = end_time-start_time
minutes, seconds= divmod(tempo, 60)
print('tempo:',round(minutes,0),':', round(seconds,0))
print('\nRESULTADOS:')



print('\nForça Axial = ',round(F_axial,2),'[N]')

print('\nMomento = ',round(Momento,2),'[Nm]')


print('\nMassa Total',round(massa,3),'[kg]')

print('\nPotencia requerida',round(POT),'[W]')

test = 1-(Pot_voo/POT)**-1
#print(test)
print('\nPotencia de Vôo',round(Pot_voo),round(Pot_voo/test,2),'[W]')

print('\nEficiencia = ',round(test*100,2),'[%]')

print('\nNumero de vezes de reuso',reuso,'\n')


if plot == True:

    plt.figure(1)
    plt.plot(r,c)
    plt.plot(r_esq,l_esq)
    plt.plot(r_dir,l_dir)
    plt.plot(hub_x,hub_y)
    if plot_bezier == True:
        plt.plot(cpx,cpy,'ro')

    plt.ylabel('Corda [m]')
    plt.xlabel('Raio da pá [m]')
    plt.title('Distribuição De Cordas')
    plt.plot(cpx,cpy,'ro')
    plt.grid()

    plt.figure(2)
    plt.plot(r,theta)
    plt.ylabel('Angulo de torção theta [º]')
    plt.xlabel('Raio da pá[m]')
    plt.title('ThetaXR')
    plt.grid()

    plt.figure(3)
    plt.plot(r,reynolds)
    plt.ylabel('Numero de Reynolds')
    plt.xlabel('Raio da pá[m]')
    plt.title('Reynolds X R')
    plt.grid()

    plt.figure(4)
    plt.plot(r,a)
    plt.ylabel('Fator de Indução Axial')
    plt.xlabel('Raio da pá[m]')
    plt.title('Fator de Indução Axial X Raio da Pá')
    plt.grid()

    plt.figure(5)
    plt.plot(r,a_l)
    plt.ylabel('Fator de Indução Radial')
    plt.xlabel('Raio da pá[m]')
    plt.title('Fator de Indução Radial X Raio da Pá')
    plt.grid()

    if Xfoil_use == True:

        plt.figure(6)
        plt.plot(r,CL)
        plt.ylabel('Coeficiente de sustentação')
        plt.xlabel('Raio da Pá [m]')
        plt.title('Cl X Raio da Pá')
        plt.grid()

        plt.figure(7)
        plt.plot(r,CD)
        plt.ylabel('Coeficiente de Arrasto')
        plt.xlabel('Raio da Pá [m]')
        plt.title('Cd X Raio da Pá')
        plt.grid()


        plt.figure(9)
        plt.plot(r,ratio)
        plt.ylabel('Cl/Cd')
        plt.xlabel('Raio da pá[m]')
        plt.title('Cl/Cd X R')
        plt.grid()

        plt.figure(10)
        plt.plot(r,ratio32)
        plt.ylabel('Cl^(3/2)/Cd')
        plt.xlabel('Raio da pá[m]')
        plt.title('Cl^(3/2)/Cd X Raio da Pá')
        plt.grid()

    plt.figure(8)
    plt.plot(r,mach)
    plt.ylabel('numero de mach')
    plt.xlabel('Raio da pá[m]')
    plt.title('Mach X Raio da Pá')
    plt.grid()



    plt.figure(11)
    plt.plot(rn,t)
    plt.ylabel('Força Axial [N]')
    plt.xlabel('Raio da pá[m]')
    plt.title('Força Axial X R')
    plt.grid()

    plt.figure(12)
    plt.plot(rn,m1)
    plt.ylabel('Força Radial [N]')
    plt.xlabel('Raio da pá[m]')
    plt.title('Força Axial X R')
    plt.grid()

    plt.figure(13)
    plt.plot(r,v_abs)
    plt.ylabel('Velocidade[m/s]')
    plt.xlabel('Raio da pá[m]')
    plt.title('Velocidade corrigida X R')
    plt.grid()

    plt.figure(14)
    plt.plot(r,v_rel)
    plt.ylabel('Velocidade[m/s]')
    plt.xlabel('Raio da pá[m]')
    plt.title('Velocidade não corrigida X R')
    plt.grid()

    plt.figure(15)
    plt.plot(v_abs,v_rel)
    plt.ylabel('Velocidade não corrigida[m/s]')
    plt.xlabel('Velocidade corrigida[m/s]')
    plt.title('Velocidade não corrigida X Velocidade corrigida')
    plt.grid()
    
    if struc == True:
        
        Di = 9e-3
        mode_S = 'tubo'     
        print(r)
        Estru = struct(m1,t,r,c,m_s,w,AeroFolio,Di,mode_S)
 
        
        plt.figure(90)
        #plt.plot(Estru.rs,Estru.mflet)
        #r.pop()
        plt.plot(r,Estru.mflet)
        xlim = max(Estru.mflet)*1.1
        plt.autoscale(enable=True, axis='both', tight=None)
        plt.grid()
        plt.ylabel('Tensao [MPa]')
        plt.xlabel('Raio da pá [m]')
        plt.title('Esforços da Longarina')
        plt.figure(91)
        plt.plot(r, Estru.d)
        plt.grid()
        plt.ylabel('Diametro Externo [mm]')
        plt.xlabel('Raio da pá [m]')
        plt.title('Da longarina da Longarina')
        
        
        Max =0
        r_max =0
        for i in range(0,len(r)):
            if Estru.mflet[i]>=Max:
                Max = Estru.mflet[i]                
                r_max = r[i]
        
        print('\nEsforço maximo: ', round(Max,2),'MPa')
        print('\nPosicao: ',round(r_max,2),'m')
        print('2',r)
    plt.show()
    '''
    plt.figure(16)
    plt.plot(rn,efi)
    plt.ylabel('Eficiencia da Hélice [%]')
    plt.xlabel('Raio da pá[m]')
    plt.title('Eficiencia da Hélice X Raio da Pá')
    plt.grid()
    
    if EE == True:
        plt.figure(17)
        plt.plot(vv,ee)
        plt.ylabel('Eficiencia da Hélice [%]')
        plt.xlabel('Velocidade de Voo [m/s]')
        plt.title('Eficiencia da Hélice X Velocidade')
        plt.grid()
        
        plt.figure(18)
        plt.plot(vv,te)
        plt.ylabel('tração [N]')
        plt.xlabel('Velocidade de Voo [m/s]')
        plt.title('tracao X Velocidade')
        plt.grid()
    '''
    
#================ Resultados Em Excel =================
if Excel == True:
   print('Exportanto Resultados para Excel')
   #writer = pd.ExcelWriter('output.xlsx')
   Fx = round(F_axial,2)
   Fr = round(Momento,2)
   Pot_ex = round(Pot_voo,2)
   EF_ex = round(EFI,2)
   Thrust_ex = round(Thrust,2)
   torque_ex = round(torque,2)
   massa_ex = round(massa,3)
   wb = Workbook()
   Sheet1 = wb.add_sheet('Sheet1')
   Sheet1.write(0,0,'Aerofólio utilizado:')
   Sheet1.write(0,1,AeroFolio)
   Sheet1.write(1,0,'Força Axial [N]')
   Sheet1.write(1,1,Fx)
   Sheet1.write(2,0,'Momento [Nm]')
   Sheet1.write(2,1,Fr)
   Sheet1.write(3,0,'Potencia de Vôo [W]')
   Sheet1.write(3,1,Pot_ex)
   Sheet1.write(4,0,'Eficiencia [%]')
   Sheet1.write(4,1,EF_ex)
#   Sheet1.write(5,0,'Massa [Kg]')
#   Sheet1.write(5,1,massa_ex)
   Sheet1.write(5,0,'Numero de vezes do reuso')
   Sheet1.write(5,1,reuso)
   Sheet1.write(6,0,'Traçao Final Por Pa [N]')
   Sheet1.write(6,1,Thrust_ex)
   Sheet1.write(7,0,'Torque Final por Pa [Nm]')
   Sheet1.write(7,1,torque_ex)
   Sheet1.col(0).width = 5500
   Sheet1.col(1).width = 5000
   Sheet1.write(0,4,'Raio[m]')
   Sheet1.write(0,5,'c [m]')
   Sheet1.write(0,6,'beta [º]')
   Sheet1.write(0,7,'Reynolds[]')
   Sheet1.write(0,8,'Velocidade[m/s]')
   for i in range(0,len(r)): 
       Sheet1.write(i+1,4,r[i])
       Sheet1.write(i+1,5,c[i])
       Sheet1.write(i+1,6,theta[i])
       Sheet1.write(i+1,7,reynolds[i])
       Sheet1.write(i+1,8,v_abs[i])
       Sheet1.write(i+1,9,prof[i])
       
  
   wb.save('Resultados.xls')
   os.startfile('Resultados.xls')
