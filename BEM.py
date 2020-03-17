# Author: Lucas Marques Moreno
# Date: 3/2018 
# Description: Implementation
# of B.E.M method

import numpy as np
import matplotlib.pyplot as plt
import math as m
import os
import xlrd
import time
from corda import corda
from airfoil import airfoil
from xlwt import Workbook
from mass import mass, S2
from solid import Solid
from estrutural import struct
from bezier import bezier_profiles


#----------------------------------
# import the right version of Xfoil 
# script.
#----------------------------------
import platform
sys = platform.system() # get the O.S
if sys =='Linux':
	from xfoil import xfoil
else:
	from xfoil_windows.py import xfoil

#=====================path Xfoil=================
start_time = time.time()
global path
path  ='.\\xfoil_output'

#============== Export result to a Excel spreadsheet ==============
Excel =False

################ Flow correction #################
Correcao = True 

Hansen = True # apply the Hansen correction 


############ Active Xfoil ##############
Xfoil_use = True

############### Config Xfoil options ###############
npp = 220 # number of painels
inter =200 # number of iterations

#=============== Use Modes ==============
MultiFoil = True # set the program to switch between profiles
pos = 0 # position where it should start
R_max =[.4,.42,.85,.95,1.5] #[m] radius distance at the change of airfoil
Ang_at =[3.,3.,8.5,8.6,9.] #[°] angle of attack of each profile

#============ Corda ========================
Curva_Bezier = False
plot_bezier = False

Mode = 'trapz' # patter of chord distribution, to see all modes check corda.py
c_pa = 0.25 #[m] constant chord of the blade (it might not be used)
pcorda = [0.2,1.25] #[m] section of constant chord 
esp = 100 #[%] airfoil wall tickness, used for mass calculation
den = 420 #[kg/m³] material's specific mass 

#=========== Enable Plot =========================
plot = True
ss = True

#================ Enable Structure Stress calculation ====================
struc = True

#================ Blade parameters ======================
cpx = [0.1,.5,.75,1.2,1.6] # [m] x components for bezier curve's points
order = len(cpx)-1
npts = 20 # number of sections
cpy = [.15,0.24999986,0.34668816,0.34996463,.1] #[m] y components for bezier curve's points

prop = bezier_profiles(cpx,cpy,order,npts) # define the blade's shape with a bezier curve

r = [] # [m]list for radius value of each section
c = [] # [m]list for chord value of each section

for i in range(0,len(prop[0])):
    r.append(float(prop[0][i]))
    c.append(float(prop[1][i]))
#end for

#################### Input Data ########################
#-------------------- Natural Properties -----------------------
g = 9.81 #[m/s²] gravity
p = 1.08 #[kg/m³] air specifc mass
u = 1.4607*10**(-5) #[N.s/m²] air dynamic viscosity 

#-------------------- Dados da Hélice -------------------
AeroFolio='naca2412.dat' # airfoil that will be used

passo = 0.01 #[m] iteration distance
Vo =12. #[m/s] flgiht speed
rpm = 94.73845432 #[rpm] Propeller rotantion
c_inicial = 0.15 #[m] initial chord (if bezier curves isn't been used)
cmin = 0. # [m] minimum chord 
Ac = 0.1 #[] Correction initial factor
D = 2*cpx[len(cpx)-1] #[m] Propeller diameter
D_hub = 0.2 #[m] Propeller's Hub diameter
alfaD = 5 #[º] angle of attack (if multifoil isn't been used)
B = 2. #[] number of blades
Re_min = 1.5e3 #[] minimum Reynolds number to converge

##if Xfoil isn't been used##
if Xfoil_use == False:
    Cl=1.#[] lift coefficient
    Cd=0.1#[] drag coefficient


############### Calculation of the induction factor ############

theta=[] #[º] twist angle of the blade
efi = [] #[] efficiency of each section
pt=[] #[N] Tangential force
pn=[] #[N] Normal force
vs=[] #[m/s] air speed*sin(phi)
vc=[] #[m/s] air speed*cos(phi)
v_abs=[] #[m/s] real speed
v_rel = [] #[m/s] relative speed
reynolds = []
mach=[]
CL = []#[] lift coefficient
CD = []#[] drag coefficient
m1 = []
t =[]
ratio =[]# Cl/Cd ratio
ratio32 = []# Cl^3/2/Cd ratio
a = [0]
a_l = [0]
w = rpm*2*np.pi/60# [rad/s] angular velocity
rps = rpm/60# [rot/s] angular frequency
Phi = []

Jo = Vo/(D*rps) #[] advance coefficient
cy = 0
print('Advance Coefficient:',round(Jo,2),'[m/rot]')

torque = 0#[N] Initial Axial Force value
Thrust = 0#[N] Initial Thrust value
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
V_som = 343 #[m/s] speed of the sound

reuso=0

R_max.append(Rf)
print('')
print('Starting Induce Factor Calculations')
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
        Mach = V_abs/V_som#[] Mach number
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
                print('Section(s) Failed to Converge:',reuso)
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
            
        
        Pn = (1/2)*p*C*Cn*Vrel**2 
        Pt = (1/2)*p*C*Ct*Vrel**2
        
        pn.append(Pn)
        pt.append(Pt)
        if ss == True:
            SS = Solid(C,AeroFolio,phi,R,j)
            plt.figure(80)
            plt.plot(SS.x,SS.y,label=('Section '+str(j)))
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.title('Sections Of The Blade')
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
        print('Progress: ',round(progresso,2),'%')
        Ri = Ri+passo

print('')
print('Completed Induce Factor Calculations')

#========== Propeller's Shape ==========
# populates some lists to illustrate the
# final shape of the propeller

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

############## Force Calculation #####################
rn = []
Pot=0
print('')
print('Starting Force Calculations')
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
    ###### Axial Force Calculation #####
    Yn = (Pni-Pn)/(Ri-R)
    Sn = (Pn*Ri-Pni*R)/(Ri-R)
    PN = Yn*R + Sn

    T = (1/2)* Yn*(Ri**2-R**2)+Sn*(Ri-R)
    t.append(T)
    Thrust = (Thrust+T)
    

    ###### Radial Force Calculation L#####

    Yt = (Pti-Pt)/(Ri-R)
    St = (Pt*Ri-Pti*R)/(Ri-R)
    PT = Yt*R + St
    M = (1/3)* Yt*(Ri**3-R**3)+(1/2)*St*(Ri**2-R**2)
    m1.append(M)
    torque = (torque+M)
    #print(torque)
    progresso = (i/(1+n))*100
    print('Progress: ',round(progresso,2),'%')

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

print('')
print('Completed Force Calculations')

########### Results #######################

end_time = time.time()
tempo = end_time-start_time
minutes, seconds= divmod(tempo, 60)
print('Time elapsed (m:s): ',round(minutes,0),':', round(seconds,0))
print('\nResults: ')



print('\nThurst: ',round(F_axial,2),'[N]')

print('\nMomentum: ',round(Momento,2),'[Nm]')


print('\nTotal Mass: ',round(massa,3),'[kg]')

efi = 1-(Pot_voo/POT)**-1

print('\nRequired Power: ',round(Pot_voo/efi,2),'[W]')

print('\nFlight Power: ',round(Pot_voo),'[W]')

print('\nEfficiency: ',round(efi*100,2),'[%]')

print('\nSection(s) failed to converge',reuso,'\n')


if plot == True:

    plt.figure(1)
    plt.plot(r,c)
    plt.plot(r_esq,l_esq)
    plt.plot(r_dir,l_dir)
    plt.plot(hub_x,hub_y)
    if plot_bezier == True:
        plt.plot(cpx,cpy,'ro')

    plt.ylabel('Chord [m]')
    plt.xlabel('Blade Radius [m]')
    plt.title('Blade Shape')
    plt.plot(cpx,cpy,'ro')
    plt.grid()

    plt.figure(2)
    plt.plot(r,theta)
    plt.ylabel('Twist Angle [º]')
    plt.xlabel('Blade Radius [m]')
    plt.title('Twist Angle')
    plt.grid()

    plt.figure(3)
    plt.plot(r,reynolds)
    plt.ylabel("Reynolds's Number")
    plt.xlabel('Blade Radius [m]')
    plt.title("Reynolds's Number")
    plt.grid()

    plt.figure(4)
    plt.plot(r,a)
    plt.ylabel('Axial Induction Factor')
    plt.xlabel('Blade Radius [m]')
    plt.title('Axial Induction Factor')
    plt.grid()

    plt.figure(5)
    plt.plot(r,a_l)
    plt.ylabel('Radial Induction Factor')
    plt.xlabel('Blade Radius [m]')
    plt.title('Radial Induction Factor')
    plt.grid()

    if Xfoil_use == True:

        plt.figure(6)
        plt.plot(r,CL)
        plt.ylabel('Lift Coefficient')
        plt.xlabel('Blade Radius [m]')
        plt.title('Lift Coefficient (Cl)')
        plt.grid()

        plt.figure(7)
        plt.plot(r,CD)
        plt.ylabel('Drag Coefficient')
        plt.xlabel('Blade Radius [m]')
        plt.title('Drag Coefficient (Cd)')
        plt.grid()


        plt.figure(9)
        plt.plot(r,ratio)
        plt.ylabel('Cl/Cd')
        plt.xlabel('Blade Radius [m]')
        plt.title('(Lift Drag Ratio (Cl/Cd)')
        plt.grid()

        plt.figure(10)
        plt.plot(r,ratio32)
        plt.ylabel('Cl^(3/2)/Cd')
        plt.xlabel('Blade Radius [m]')
        plt.title('Powered Lift Drag Ratio (Cl^(3/2)/Cd)')
        plt.grid()

    plt.figure(8)
    plt.plot(r,mach)
    plt.ylabel("Mach's number")
    plt.xlabel('Blade Radius [m]')
    plt.title("Mach's number")
    plt.grid()



    plt.figure(11)
    plt.plot(rn,t)
    plt.ylabel('Axial Force [N]')
    plt.xlabel('Blade Radius [m]')
    plt.title('Axial Force')
    plt.grid()

    plt.figure(12)
    plt.plot(rn,m1)
    plt.ylabel('Radial Force [N]')
    plt.xlabel('Blade Radius [m]')
    plt.title('Radial Force')
    plt.grid()

    plt.figure(13)
    plt.plot(r,v_abs)
    plt.ylabel('Speed[m/s]')
    plt.xlabel('Blade Radius [m]')
    plt.title('Absolute Speed')
    plt.grid()

    plt.figure(14)
    plt.plot(r,v_rel)
    plt.ylabel('Relative Speed[m/s]')
    plt.xlabel('Blade Radius [m]')
    plt.title('Relative Speed')
    plt.grid()

    plt.figure(15)
    plt.plot(v_abs,v_rel)
    plt.ylabel('Relative Speed [m/s]')
    plt.xlabel('Speed [m/s]')
    plt.title('Absolute and Relative Speed Deviation')
    plt.grid()
    
    if struc == True:
        
        Di = 9e-3 # [m] internal diameter
        mode_S = 'tubo'
        Estru = struct(m1,t,r,c,m_s,w,AeroFolio,Di,mode_S)
 
        
        plt.figure(90)
        plt.plot(r,Estru.mflet)
        xlim = max(Estru.mflet)*1.1
        plt.autoscale(enable=True, axis='both', tight=None)
        plt.grid()
        plt.ylabel('Resulting Stress [MPa]')
        plt.xlabel('Blade Radius [m]')
        plt.title('Beam Stress')
        plt.figure(91)
        plt.plot(r, Estru.d)
        plt.grid()
        plt.ylabel('External Diameter [mm]')
        plt.xlabel('Blade Radius [m]')
        plt.title('Beam External Diameter')
        
        
        Max =0
        r_max =0
        for i in range(0,len(r)):
            if Estru.mflet[i]>=Max:
                Max = Estru.mflet[i]                
                r_max = r[i]
        
        print('\nMaximum Stress: ', round(Max,2),'MPa')
    plt.show()

    
#================ Export results =================
if Excel == True:
   print('Exporting Results')
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
   Sheet1.write(0,0,'Airfoil:')
   Sheet1.write(0,1,AeroFolio)
   Sheet1.write(1,0,'Thrust [N]')
   Sheet1.write(1,1,Fx)
   Sheet1.write(2,0,'Momentum [Nm]')
   Sheet1.write(2,1,Fr)
   Sheet1.write(3,0,'Flight Power [W]')
   Sheet1.write(3,1,Pot_ex)
   Sheet1.write(4,0,'Eficiency [%]')
   Sheet1.write(4,1,EF_ex)
#   Sheet1.write(5,0,'Massa [Kg]')
#   Sheet1.write(5,1,massa_ex)
   Sheet1.write(5,0,'Section(s) Failed to converge')
   Sheet1.write(5,1,reuso)
   Sheet1.write(6,0,'Thrust per Blade [N]')
   Sheet1.write(6,1,Thrust_ex)
   Sheet1.write(7,0,'Momentum per Blade  [Nm]')
   Sheet1.write(7,1,torque_ex)
   Sheet1.col(0).width = 5500
   Sheet1.col(1).width = 5000
   Sheet1.write(0,4,'Radius[m]')
   Sheet1.write(0,5,'c [m]')
   Sheet1.write(0,6,'beta [º]')
   Sheet1.write(0,7,'Reynolds[]')
   Sheet1.write(0,8,'speed[m/s]')
   for i in range(0,len(r)): 
       Sheet1.write(i+1,4,r[i])
       Sheet1.write(i+1,5,c[i])
       Sheet1.write(i+1,6,theta[i])
       Sheet1.write(i+1,7,reynolds[i])
       Sheet1.write(i+1,8,v_abs[i])
       Sheet1.write(i+1,9,prof[i])
       
  
   wb.save('results.xls')
   #os.startfile('Resultados.xls')# auto opening the spreadsheet