import numpy as np
import matplotlib.pyplot as plt
import os
import platform
from scipy.integrate import quad
import subprocess as sp
from mpl_toolkits.mplot3d import Axes3D
import math

class blade:
    def __init__(self, v, rpm, B, d, r, c,alfa,airfoil, changes_section, g, p, u):
        
        self.flight_speed = v # flgiht speed
        self.rpm = rpm # rotation per minutes
        self.number_blades = B # quantity of blades
        self.diameter = d # propeller diameter 
        self.radius = r # radius's points for the calculations
        self.chord = c # chord's points for the calculations
        self.angle_of_attack = alfa # angle of attack for each airfoil
        self.airfoil = airfoil # airfoils that will be used
        self.g = g # gravity acceleration
        self.p = p #specific mass
        self.u = u #dynamic viscosity

        self.current_section = 0 #initial section
        self.changes_section = changes_section # point's that occurs a airfoil/angle of attack changes

        self.a = [] # induced factor
        self.a_l = [] # radial induced factor
        self.a_critic = 0.4 # critial value for induced factor

        self.w = self.rpm*2*np.pi/60 #angular speed

        self.phi = [] #[rad] twist angle for the velocity 
        self.theta = [] #[degree] twist angle of the blade

        self.reynolds = []
        self.mach = []
        self.v_sound = 343 #[m/s] speed of the sound

        self.pn = []
        self.pt = []

        return

    def config(self, correction = True, speed_test = False, export_sections = False, show_coefficient = False, theta_reference = None, max_ite = 100):

        # This function configures some secundaries features 
        
        self.correction = correction # Enabel Hansen's corrections calculations
        self.speed_test = speed_test # Compares and change the current angle of attack besaed in a theta angle reference
        self.export_sections = export_sections # Enabel an exportations of each sections in .dat files
        self.show_coefficient = show_coefficient # shows Cl and Cd values of each section in the terminal
        self.theta_ref = theta_reference # Assign a reference for the theta angle (speed test must be enable to work properly)
        self.max_ite = int(max_ite) # induced factors maximum number of iterations
        return

    def twist_angle_reference(self,theta_ref):

        self.theta_ref = theta_ref
        
        return

        
    def advance_ratio(self):
        
        self.advance_ratio = self.flight_speed/(self.diameter*(self.rpm/60))
       
        return self.advance_ratio

    def reynolds_number(self,v,c):

        self.re = self.p*v*c/self.u

        return self.re

    def mach_number(self,v,v_sound):
        self.ma = v/v_sound
        return self.ma

    def select_section(self,current_radius):
        
        # This functions manage the values for the current seciton
        # based in the values of changes_section list

        if current_radius <= self.changes_section[self.current_section]:
            self.current_airfoil = self.airfoil[self.current_section]
            self.current_alfa = self.angle_of_attack[self.current_section]
        else:
            self.current_section += 1
            self.current_airfoil = self.airfoil[self.current_section]
            self.current_alfa = self.angle_of_attack[self.current_section]     
        
        return

    def export(self,var, name = "value"):

        #this function export any internal variables

        file_name = name+".dat"
        f = open(file_name,'w')
        for value in var:
            f.write('%s \n' % (value))
        return
    
    def xfoil(self,inter = 50 , np = 220):

        # This function uses some Xfoils calculations to retrieve the values for Cl and Cd
        Coef = list()
        if os.path.isfile("xfoil_output.txt"):
            os.remove("xfoil_output.txt")
        itens  = [self.current_airfoil,self.current_alfa, self.current_alfa, '0', round(self.ma,2),self.re,inter,np]

        if platform.system() =='Linux':
               
                name = '\n'+itens[0]
                rey = '\n'+str(itens[5])
                alfa = '\n'+str(itens[1])
                inter = '\n'+str(itens[6])
                mach = '\n'+str(itens[4])
                np = '\n'+str(itens[7])
                path = 'load'+name+'\npane\nppar\nn'+np+'\n\n\noper\nvisc'+rey+'\nmach'+mach+'\niter'+inter+'\npacc\nxfoil_output.txt\n\naseq'+alfa+alfa+'\n0\n\nquit'
                n = path
                p = sp.Popen(['xfoil'],
		        stdin=sp.PIPE,
		        stdout=sp.PIPE,
		        stderr=sp.STDOUT)
                grep_stdout = p.communicate(input=n.encode())[0]
        else:
            os.system('del dump')
            os.system('del xfoil_input')
            os.system('del screen')
            f=open('xfoil_input.txt','w')
            f.write('%s \n' % ('load'))
            f.write('%s \n' % (itens[0]))
            f.write('%s \n' % ('pane'))
            f.write('%s \n' % ('ppar'))
            f.write('%s \n' % ('n'))
            f.write('%s \n' % (itens[7]))
            f.write('%s \n' % ('   '))
            f.write('%s \n' % ('   '))
            f.write('%s \n' % ('oper'))
            f.write('%s \n' % ('visc'))
            f.write('%s \n' % (itens[5]))
            f.write('%s \n' % ('mach'))
            f.write('%s \n' % (itens[4]))
            f.write('%s \n' % ('iter'))
            f.write('%s \n' % (itens[6]))
            f.write('%s \n' % ('pacc'))
            f.write('%s \n' % ('xfoil_output'))
            f.write('%s \n' % ('dump'))
            f.write('%s \n' % ('aseq'))
            f.write('%s \n' % (itens[1]))
            f.write('%s \n' % (itens[2]))
            f.write('%s \n' % (itens[3]))
            f.write('%s \n' % ('  '))
            f.write('%s \n' % ('quit'))
            f.close()
            os.system('xfoil.exe < xfoil_input.txt > screen')

        i    = 0
        aoa  = list()
        cl   = list()
        cd   = list()
        cdp  = list()
        cm   = list()
        xutr = list()
        xltr = list()

        if os.path.isfile("xfoil_output.txt"):
            f      = open('xfoil_output.txt','r')

            for line in f:
                if (i > 11):

                    aoa.append(line.strip().split()[0])
                    cl.append(line.strip().split()[1])
                    cd.append(line.strip().split()[2])
                    cdp.append(line.strip().split()[3])
                    cm.append(line.strip().split()[4])
                    xutr.append(line.strip().split()[5])
                    xltr.append(line.strip().split()[6])
                    #print((line.strip().split()[6]))
                    res = float(line.strip().split()[2])
                else:
                    res = 1.0
                i += 1

            f.close()

            #colocar maximo
        j = 0


        if len(cl)==0:
            self.Cl =0.0000001
            self.Cd =0.0000001
            

            j +=1

        else:
            self.Cl =float(cl[-1])
            self.Cd =float(cd[-1])

            j=j


        return

    def velocity_curve(self, Theta, theta_ref):

        # This functions makes a correction for
        # the 'new' angle of attack if speed test 
        # is enable

        gamma= (theta_ref - Theta)*180/np.pi
        print("Gamma angle: ",round(gamma,1))
        self.current_alfa += gamma
        print("New Alfa angle: ",round(self.current_alfa ,1))
        #phi = (theta_ref+ self.current_alfa)*np.pi/180  
        
        return 

    
    def induced_factor(self):

        # This functions starts the BEM Method 
        # with the induced factors of each section
        
        for i in range(len(self.radius)):
            converge = False
            count = 0
            self.a.append(.1)
            self.a_l.append(.01)

            while(not converge):
                if self.radius[i] != self.radius[-1]:
                    self.select_section(self.radius[i])
                
                print("progress: ",round(i/len(self.radius)*100,2), " %")
                #print(self.a)
                Vt = self.w * self.radius[i]
                Tan_phi = ((1+self.a[i])*self.flight_speed)/((1-self.a_l[i])*Vt)
                
                
                v_rel = (self.flight_speed**2+(Vt)**2)**(1/2)
                V_relS = self.flight_speed*(1+self.a[i])
                V_relC = Vt*(1-self.a_l[i])
                v_abs = (V_relS**2+V_relC**2)**(1/2)

                if Tan_phi < 1e-10 :
                    Tan_phi = 1e-6
                Phi = np.arctan(Tan_phi)
                self.phi.append(Phi)
                self.theta.append((Phi-self.current_alfa*(np.pi/180))*180/np.pi)


                if self.speed_test:
                    self.velocity_curve(self.phi[-1],self.theta_ref[i])
                

                re = self.reynolds_number(v_abs,self.chord[i])
                self.reynolds.append(re)

                ma = self.mach_number(v_abs,self.v_sound)
                self.mach.append(ma)

                Lambda = self.flight_speed/Vt

                f = (self.number_blades/2)*(1/Lambda)*(1+Lambda**2)**(1/2)*(1-(self.radius[i]/(self.diameter/2)))

                F = (2/np.pi)*(np.arctan((np.exp(2*f)-1)**(1/2)))

                self.xfoil()
                if self.show_coefficient:
                    print("Cl: ",self.Cl,"       " ,"Cd: ",self.Cd)
                #self.Cl = 1.
                #self.Cd = 1e-4
                L = (1/2)*self.p*self.chord[i]*self.Cl*v_rel**2
                Dr = (1/2)*self.p*self.chord[i]*self.Cd*v_rel**2 
                Cn = self.Cl*np.cos(Phi)-self.Cd*np.sin(Phi)
                Ct = self.Cl*np.sin(Phi)+self.Cd*np.cos(Phi)*self.radius[i]
                
                sigma = self.chord[i]*self.number_blades/(2*np.pi*self.radius[i])
                I1 = 4*np.sin(Phi)**2
                I2 = sigma*Cn
                a = (1/((I1/I2)-1))
                
                I3 = 4*np.sin(Phi)*np.cos(Phi)
                I4 = sigma*Ct
                a_l = (1/((I3/I4)+1))

                if self.a[i] >self.a_critic :
                    K_h = 4*F*np.sin(Phi)**2/(sigma*Cn)
                    a=((1/2)*(2+K_h*(1-2*self.a_critic)-((K_h*(1-2*self.a_critic)+2)**2+4*(K_h*self.a_critic**2-1))**(1/2)))


                self.a[i] = a
                self.a_l[i]

                if self.a[i] <= (1/3):
                    Ct = 4*self.a[i]*(1-self.a[i])*F
                else:
                    Ct = 4*self.a[i]*(1-(1/4)*(5-3*self.a[i])*self.a[i])*F

                Pn = ((1/2)*self.p*self.chord[i]*Cn*v_rel**2)
                Pt = ((1/2)*self.p*self.chord[i]*Ct*v_rel**2)


                if self.correction :

                    anew = Pn/(4*np.pi*self.radius[i]*self.p*v_rel**2+(1+a))
                    alnew = Pt/(4*np.pi*self.radius[i]**3*self.p*v_rel*(1+a)*self.rpm)

                    amiddle = (anew+a)/2
                    almiddle = (alnew +a_l)/2
                    print("\n=======================================")
                    print("a: ",a,"     ","a_l: ",a_l)
                    print("anew: ",anew,"     ","a_lnew: ",alnew)
                    print("amiddle: ",amiddle,"     ","almiddle: ",almiddle)
                    print("Da: ",abs(amiddle-a),"     ","Da_l: ",(almiddle-a_l))
                    print("=======================================")
                    if abs(amiddle-a) < 1e-2 and abs(almiddle-a_l) <1e-2:
                        converge = True
                    count +=1
                    if count ==  self.max_ite: 
                        converge = True
                    '''
                    if self.a[i] >self.a_critic :
                        self.a.pop()
                        K_h = 4*F*np.sin(Phi)**2/(sigma*Cn)
                        self.a.append((1/2)*(2+K_h*(1-2*self.a_critic)-((K_h*(1-2*self.a_critic)+2)**2+4*(K_h*self.a_critic**2-1))**(1/2)))

                    if self.a[i] <= (1/3):
                        Ct = 4*self.a[i]*(1-self.a[i])*F
                    else:
                        Ct = 4*self.a[i]*(1-(1/4)*(5-3*self.a[i])*self.a[i])*F'''

            self.pn.append(Pn)
            self.pt.append(Pt)

            #self.a.append(a)
            #self.a_l.append(a_l)
            
            #self.build_geometry(i)

            #if self.export_sections:
            #    self.export_section(i)
                

        return

    def forces(self):


        # This functions ends the BEM Method 
        # integranting the total thrust and
        # momentum

        self.thrust = 0
        self.momentum = 0
        self.T_tes = 0
        for i in range(len(self.radius)-1):
            #Axial forces
                An = (self.pn[i+1]-self.pn[i])/(self.radius[i+1]-self.radius[i])
                Bn = (self.pn[i+1]*self.radius[i+1]-self.pn[i]*self.radius[i])/(self.radius[i+1]-self.radius[i])
                Pn = An*self.radius[i] + Bn 

                T = (1/2)* An*(self.radius[i+1]**2-self.radius[i]**2)+Bn*(self.radius[i+1]-self.radius[i])
                self.thrust = (self.thrust+T)
                dt = 4*np.pi*self.radius[-1] * self.p * self.flight_speed*self.a[i]*(1+self.a[i])*(self.radius[i]-self.radius[i-1])
                self.T_tes = self.T_tes+Pn*(self.radius[i+1]-self.radius[i])

            #Radial forces

                At = (self.pt[i+1]-self.pt[i])/(self.radius[i+1]-self.radius[i])
                Bt = (self.pt[i+1]*self.radius[i+1]-self.pt[i]*self.radius[i])/(self.radius[i+1]-self.radius[i])
                PT = At*self.radius[i] + Bt
                M = (1/3)* At*(self.radius[i+1]**3-self.radius[i]**3)+(1/2)*Bt*(self.radius[i+1]**2-self.radius[i]**2)
                self.momentum = (self.momentum+M)

        self.thrust = self.thrust*self.number_blades
        self.momentum = self.momentum*self.number_blades        
        self.power_flight = self.thrust*self.flight_speed
        return

    def force_test(self):

        self.T = 0
        self.Q = 0

        for i in range(len(self.radius)):
            self.T = self.T + self.pn[i] *0.05
            self.Q = self.Q + self.pt[i] *0.05

        return


    def I1_prime(self, r, zeta, cl, cd):
        ksi = r/self.radius[-1]
        omega = self.w
        Lambda = self.flight_speed/(r*omega)
        index = int(len(self.radius)*ksi)
        tan_phiT = Lambda*(1+zeta/2)
        phiT = np.arctan(tan_phiT)
        phi = (np.arctan(tan_phiT/ksi))
        tan_phi = np.tan(phi)
        f = (self.number_blades/2)*(1-ksi)/np.sin(phiT)
        F = (np.arctan(np.exp(2*f)-1**.5))
        G = F*np.sin(phi)*np.cos(phi)*(r*omega)/self.flight_speed
        re = self.reynolds_number(self.flight_speed,self.chord[index])
        self.reynolds.append(re)
        ma = self.mach_number(self.flight_speed,self.v_sound)
        self.mach.append(ma)
        #self.xfoil()
        epsilon = cd[index]/cl[index]
        i1 = 4*ksi*G*(1-epsilon*tan_phi)
        return i1

    def I2_prime(self, r, zeta, cl, cd):
        ksi = r/self.radius[-1]
        omega = self.w
        Lambda = self.flight_speed/(r*omega)
        index = int(len(self.radius)*ksi)
        tan_phiT = Lambda*(1+zeta/2)
        phiT = np.arctan(tan_phiT)
        phi = (np.arctan(tan_phiT/ksi))
        tan_phi = np.tan(phi)
        f = (self.number_blades/2)*(1-ksi)/np.sin(phiT)
        F = (np.arctan(np.exp(2*f)-1**.5))
        G = F*np.sin(phi)*np.cos(phi)*(r*omega)/self.flight_speed
        re = self.reynolds_number(self.flight_speed,self.chord[index])
        self.reynolds.append(re)
        ma = self.mach_number(self.flight_speed,self.v_sound)
        self.mach.append(ma)
        #self.xfoil()
        epsilon = cd[index]/cl[index]
        i2 = Lambda*(self.I1_prime(r, zeta, cl, cd)/ksi)*(1+epsilon/tan_phi)*np.sin(phi)*np.cos(phi)
        return i2


    def J1_prime(self, r, zeta, cl, cd):
        ksi = r/self.radius[-1]
        omega = self.w
        Lambda = self.flight_speed/(r*omega)
        index = int(len(self.radius)*ksi)
        tan_phiT = Lambda*(1+zeta/2)
        phiT = np.arctan(tan_phiT)
        phi = (np.arctan(tan_phiT/ksi))
        tan_phi = np.tan(phi)
        f = (self.number_blades/2)*(1-ksi)/np.sin(phiT)
        F = (np.arctan(np.exp(2*f)-1**.5))
        G = F*np.sin(phi)*np.cos(phi)*(r*omega)/self.flight_speed
        re = self.reynolds_number(self.flight_speed,self.chord[index])
        self.reynolds.append(re)
        ma = self.mach_number(self.flight_speed,self.v_sound)
        self.mach.append(ma)
        #self.xfoil()
        epsilon = cd[index]/cl[index]
        j1 = 4*ksi*G*(1+epsilon/tan_phi)
        return j1


    def J2_prime(self, r, zeta, cl, cd):
        ksi = r/self.radius[-1]
        omega = self.w
        Lambda = self.flight_speed/(r*omega)
        index = int(len(self.radius)*ksi)
        tan_phiT = Lambda*(1+zeta/2)
        phiT = np.arctan(tan_phiT)
        phi = (np.arctan(tan_phiT/ksi))
        tan_phi = np.tan(phi)
        f = (self.number_blades/2)*(1-ksi)/np.sin(phiT)
        F = (np.arctan(np.exp(2*f)-1**.5))
        G = F*np.sin(phi)*np.cos(phi)*(r*omega)/self.flight_speed
        re = self.reynolds_number(self.flight_speed,self.chord[index])
        self.reynolds.append(re)
        ma = self.mach_number(self.flight_speed,self.v_sound)
        self.mach.append(ma)
        #self.xfoil()
        epsilon = cd[index]/cl[index]
        j2 = (self.J1_prime(r, zeta, cl, cd))*(1-epsilon*tan_phi)*np.cos(phi)**2
        return j2


    def I1(self, zeta, cl, cd):
        i1 = quad(self.I1_prime,self.radius[0],self.radius[-1], args=(zeta, cl, cd))      
        return i1[0]
    def I2(self, zeta, cl, cd):
        i2 = quad(self.I2_prime,self.radius[0],self.radius[-1], args=(zeta, cl, cd))      
        return i2[0]
    def J1(self, zeta, cl, cd):
        j1 = quad(self.J1_prime,self.radius[0],self.radius[-1], args=(zeta, cl, cd))      
        return j1[0]
    def J2(self, zeta, cl, cd):
        j2 = quad(self.J2_prime,self.radius[0],self.radius[-1], args=(zeta, cl, cd))      
        return j2[0]


    def lieback_optimun_design(self):
        # step 1
        zeta = 1
        F=[None]*len(self.radius)
        phi =[None]*len(self.radius)
        W = [None]*len(self.radius)
        Wc = [None]*len(self.radius)
        cl =[None]*len(self.radius)
        cd = [None]*len(self.radius)
        Chord = [None]*len(self.radius)
        Twist = [None]*len(self.radius)
        omega = self.w
        converge = False
        ite=1
        while not converge:
            print("======================\n NEW ITERATION\n ======================")
            for i in range(len(self.radius)):
                
                # step 2
                if self.radius[i] != self.radius[-1]:
                    self.select_section(self.radius[i])
                Lambda = self.flight_speed/(self.radius[i]*omega)
                x = self.w*self.radius[i]/self.flight_speed
                ksi = self.radius[i]/self.radius[-1]
                tan_phiT = Lambda*(1+zeta/2)
                phiT = np.arctan(tan_phiT)
                phi[i] = (np.arctan(tan_phiT/ksi))
                f = (self.number_blades/2)*(1-ksi)/np.sin(phiT)
                F[i] = (np.arctan(np.exp(2*f)-1**.5))
                
                # step 3

                re = self.reynolds_number(self.flight_speed,self.chord[i])
                self.reynolds.append(re)
                ma = self.mach_number(self.flight_speed,self.v_sound)
                self.mach.append(ma)
                self.xfoil()
                epsilon = self.Cd/self.Cl
                cd[i] = self.Cd
                cl[i] = self.Cl
                G = F[i]*np.sin(phi[i])*np.cos(phi[i])*(self.radius[i]*omega)/self.flight_speed
                wc = 4*np.pi*Lambda*G*self.flight_speed*self.radius[-1]*ksi/(self.Cl*self.number_blades)
                Wc[i] = (wc)
                
                # step 6

                a = (zeta/2)*(np.cos(phi[i])**2)*(1-epsilon*np.tan(phi[i]))
                al = (zeta/(2*x))*np.cos(phi[i])*np.sin(phi[i])*(1+epsilon/np.tan(phi[i]))
                w = self.flight_speed*(1+a)/np.sin(phi[i])
                wt = self.flight_speed*ksi*np.sin(phi[i])*np.cos(phi[i])
                circulation = 2*np.pi*self.radius[i]*F[i]*wt/self.number_blades
                W[i] = (w)
                
                # step 7

                chord = 2*circulation/(w*self.Cl)
                twist = phi[i]*180/np.pi + self.current_alfa
                Chord[i]=(chord)
                Twist[i]=(twist)
                # step 8
                
            print(ite)
            self.power_constraint=300
            Pc = 2*self.power_constraint/(self.p*(self.flight_speed)**2*np.pi*(self.radius[-1]**2))
            new_zeta = -(self.J1(zeta, cl, cd)/(self.J2(zeta, cl, cd))*2)+((self.J1(zeta, cl, cd)/(self.J2(zeta, cl, cd))*2)**2 + Pc/self.J2(zeta, cl, cd))**.5
            print(new_zeta)
            print(abs((new_zeta)/new_zeta))
            if(abs((new_zeta-zeta)/new_zeta) <1e-4):
                converge = True
            else:
                zeta = new_zeta
                ite +=1
        Tc = self.I1(zeta, cl, cd)*zeta-self.I2(zeta, cl, cd)*zeta**2
        T = self.p*(self.flight_speed**2)*np.pi*(self.radius[-1]**2)*Tc/2

        print(Chord,"\n")
        print(Twist,"\n")

        print(T)
        print(Tc/Pc)
        
        return




    def Ct(self):
        self.ct = self.thrust /(self.p*((self.rpm/60)**2)*(self.diameter**4))
        return self.ct

    def Cp(self):
        
        self.cp = (self.power_flight)/(self.p*((self.rpm/60)**3)*(self.diameter**5))
        return self.cp

    def Efficiency(self):
        self.power_momentum = self.momentum * self.rpm*2*np.pi/60
        self.efficiency = self.advance_ratio()*self.Ct()/self.Cp() #1-(self.power_flight /self.power_momentum)**-1
        self.test_efi = self.power_flight /self.power_momentum
        return  self.efficiency

    def plot_blade(self):
        hub_x=[]
        hub_y=[]
        l_esq=[0,self.chord[0]]
        r_esq=[self.radius[0],self.radius[0]]
        l_dir=[0,self.chord[-1]]
        r_dir=[self.radius[-1],self.radius[-1]]

        plt.figure(1)
        plt.plot(self.radius,self.chord)
        plt.plot(r_esq,l_esq)
        plt.plot(r_dir,l_dir)
        plt.plot(hub_x,hub_y)
        plt.show()
        return


    def airfoil_max_thickness(self,airfoil):

        # This functions read a airfoil file
        # and calculates the maximum thickness
        # and the chord's point where that occours

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
        return tks, x_t


    def rotational_matrix(self,airfoil,c,phi,x_t,tks,i):

        # This function read the original airfoils and
        # and change them for the blade's 
        # corresponding sections coordenates

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
            Y = float(line.strip().split()[1])*c
         
            xrot = X*np.cos(phi)+Y*np.sin(phi)
            yrot = Y*np.cos(phi)-X*np.sin(phi)
            x.append(xrot)
            y.append(yrot)
        return x,y


    def build_geometry(self,i):
        tks,pos= self.airfoil_max_thickness(self.current_airfoil)
        G = self.rotational_matrix(self.current_airfoil, self.chord[i], self.phi[i], pos, tks, i)
        self.x = G[0]
        self.y = G[1]
        self.z=[]
        for j in range(len(self.x)):
            self.z.append(self.radius)

        return
    
    def export_section(self,i):

        # This functions export the sections coordenates
        # in a .dat file
        
        n = self.radius[i]*1000
        name_file = 'R'+str(int(n))+'.txt'

        f=open(name_file,'w')
        for i in range(0,len(self.x)):
         f.write('%s \t' % (str(self.x[i])))
         f.write('%s \t' % (str(self.y[i])))
         f.write('%s \n' % (str(n/1000)))
        f.close()
        
        return

    def export_cad(self):
        pass
        return

    def cad_priveiw(self):
        fig = plt.figure(20)
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.x, self.y, self.z)
        plt.show()
        return