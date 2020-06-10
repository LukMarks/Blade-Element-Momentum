import numpy as np
import matplotlib.pyplot as plt
import os
import platform
import subprocess as sp
from mpl_toolkits.mplot3d import Axes3D



class blade:
    def __init__(self, v, rpm, B, d, r, c,alfa,airfoil, changes_section, g, p, u):
        
        self.flight_speed = v
        self.rpm = rpm
        self.number_blades = B
        self.diameter = d
        self.radius = r
        self.chord = c
        self.angle_of_attack = alfa
        self.airfoil = airfoil
        self.g = g # gravity acceleration
        self.p = p #specific mass
        self.u = u #dynamic viscosity

        self.current_section = 0 #initial section
        self.changes_section = changes_section

        self.a = [0] # induced factor
        self.a_l = [0] # radial imduced factor
        self.a_critic = 0.4

        self.w = self.rpm*2*np.pi/60 #angular speed

        self.phi = []
        self.theta = [] # twist angle

        self.reynolds = []
        self.mach = []
        self.v_sound = 343 #[m/s]

        self.pn = []
        self.pt = []

        return

    def config(self, correction = True, speed_test = False, export_sections = False, show_coefficient = False):
        
        self.correction = correction
        self.speed_test = speed_test 
        self.export_sections = export_sections
        self.show_coefficient = show_coefficient
        return

    def twist_angle_reference(self,theta_ref):

        self.theta_ref = theta_ref
        
        return

    def advance_ratio(self):
        
        self.Jo = self.flight_speed/(self.diameter*(self.rpm/60))
        
        return

    def reynolds_number(self,v,c):

        self.re = self.p*v*c/self.u

        return self.re

    def mach_number(self,v,v_sound):
        self.ma = v/v_sound
        return self.ma

    def select_section(self,current_radius):
        if current_radius <= self.changes_section[self.current_section]:
            self.current_airfoil = self.airfoil[self.current_section]
            self.current_alfa = self.angle_of_attack[self.current_section]
        else:
            self.current_section += 1
            self.current_airfoil = self.airfoil[self.current_section]
            self.current_alfa = self.angle_of_attack[self.current_section]     
        
        return
    
    def xfoil(self,inter = 200 , np = 220):

        Coef = list()
        os.remove("xfoil_output.txt")
        itens  = [self.current_airfoil,self.current_alfa, self.current_alfa, '0', round(self.ma,2),int(self.re),inter,np]

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
            os.system('del xfoil_output')
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

        gamma= theta_ref*np.pi/180 - Theta
        print("Gamma angle: ",round(gamma*180/np.pi,1))
        self.current_alfa += (gamma*180/np.pi)
        print("New Alfa angle: ",round(self.current_alfa ,1))
        phi = (theta_ref+ self.current_alfa)*np.pi/180  
        
        return phi

    def induced_factor(self):


        for i in range(len(self.radius)):
            if self.radius[i] != self.radius[-1]:
                self.select_section(self.radius[i])
            
            
            Vt = self.w * self.radius[i]
            Tan_phi = Tan_phi = ((1+self.a[i])*self.flight_speed)/((1-self.a_l[i])*Vt)
            Phi = np.arctan(Tan_phi)
            self.phi.append(Phi)
            self.theta.append((Phi-self.current_alfa*(np.pi/180))*180/np.pi)

            #todo add here function to speed variations test

            if self.speed_test:
                Phi = self.velocity_curve(self.theta[-1],self.theta_ref[i])

            self.phi.append(Phi)
            self.theta.append((Phi-self.current_alfa*(np.pi/180))*180/np.pi)

            v_rel = (self.flight_speed**2+(Vt)**2)**(1/2)
            V_relS = self.flight_speed*(1+self.a[i])
            V_relC = Vt*(1-self.a_l[i])

            v_abs = (V_relS**2+V_relC**2)**(1/2)

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
            L = (1/2)*self.p*self.chord[i]*self.Cl*v_rel**2
            Dr = (1/2)*self.p*self.chord[i]*self.Cd*v_rel**2 
            Cn = self.Cl*np.cos(Phi)+self.Cd*np.sin(Phi)
            Ct = self.Cl*np.sin(Phi)-self.Cd*np.cos(Phi)
            
            sigma = self.chord[i]*self.number_blades/(2*np.pi*self.radius[i])
            I1 = 4*np.sin(Phi)**2
            I2 = sigma*Cn
            self.a.append(1/((I1/I2)-1))
            
            I3 = 4*np.sin(Phi)*np.cos(Phi)
            I4 = sigma*Ct
            self.a_l.append(1/((I3/I4)+1))

            if self.correction :

                if self.a[i] >self.a_critic :
                    a.pop()
                    K_h = 4*F*np.sin(Phi)**2/(sigma*Cn)
                    a.append((1/2)*(2+K_h*(1-2*ac)-((K_h*(1-2*ac)+2)**2+4*(K_h*ac**2-1))**(1/2)))

                if self.a[i] <= (1/3):
                    Ct = 4*self.a[i]*(1-self.a[i])*F
                else:
                    Ct = 4*self.a[i]*(1-(1/4)*(5-3*self.a[i])*self.a[i])*F

            self.pn.append((1/2)*self.p*self.chord[i]*Cn*v_rel**2)
            self.pt.append((1/2)*self.p*self.chord[i]*Ct*v_rel**2)
            
            self.build_geometry(i)

            if self.export_sections:
                self.export_section(i)
                

        return

    def forces(self):
        self.thrust = 0
        self.momentum = 0

        for i in range(len(self.radius)-1):
            #Axial forces
                Yn = (self.pn[i+1]-self.pn[i])/(self.radius[i+1]-self.radius[i])
                Sn = (self.pn[i]*self.radius[i+1]-self.pn[i+1]*self.radius[i])/(self.radius[i+1]-self.radius[i])
                PN = Yn*self.radius[i] + Sn

                T = (1/2)* Yn*(self.radius[i+1]**2-self.radius[i]**2)+Sn*(self.radius[i+1]-self.radius[i])
                self.thrust = (self.thrust+T)

            #Radial forces

                Yt = (self.pt[i+1]-self.pt[i])/(self.radius[i+1]-self.radius[i])
                St = (self.pt[i]*self.radius[i+1]-self.pt[i+1]*self.radius[i])/(self.radius[i+1]-self.radius[i])
                PT = Yt*self.radius[i] + St
                M = (1/3)* Yt*(self.radius[i+1]**3-self.radius[i]**3)+(1/2)*St*(self.radius[i+1]**2-self.radius[i]**2)
                self.momentum = (self.momentum+M)

        self.thrust = self.thrust*self.number_blades
        self.momentum = self.momentum*self.number_blades        
        self.power_flight = self.thrust*self.flight_speed
        return

    def Ct(self):
        self.ct = self.thrust /(self.p*((self.rpm/60)**2)*(self.diameter**4))
        return self.ct

    def Cp(self):
        
        self.cp = (self.power_flight)/(self.p*((self.rpm/60)**3)*(self.diameter**5))
        return self.cp

    def Efficiency(self):
        self.power_momentum = self.momentum * self.rpm*2*np.pi/60
        self.efficiency = 1-(self.power_flight /self.power_momentum)**-1
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