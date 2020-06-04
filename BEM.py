import numpy as np
import matplotlib.pyplot as plt


class blade:
    def __init__(self, v, rpm, B, d, r, c, g, p, u):
        
        self.flgiht_speed = v
        self.rpm = rpm
        self.number_blades = B
        self.diameter = d
        self.radius = r
        self.chord = c
        self.g = g # gravity acceleration
        self.p = p #specific mass
        self.u = u #dynamic viscosity

  
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

    def config(self):
        pass

        return

    def advance_ratio(self):
        
        self.Jo = self.flgiht_speed*(self.diameter/(self.rpm/60))
        
        return

    def reynolds_number(self,v,c)

        re = self.p*v*c/self.u

        return re

    def mach_number(self,v,v_sound)
        ma = v/v_sound
        return ma
    
    def start(self):


        for i in range(self.radius):

            alpha = alphaD*np.pi/180
            
            #A = a[i]
            #a_l = a_l[i]#a[-1]
            
            vt = self.w * self.radius[i]
            Tan_phi = Tan_phi = ((1+a[i])*self.flgiht_speed)/((1-a_l[i])*Vt)
            Phi = np.arctan(Tan_phi)
            self.phi.append(Phi)
            self.theta.append((Phi-alpha)*180/np.pi)

            #todo add here function to speed variations test

            v_rel = (self.flgiht_speed**2+(Vt)**2)**(1/2)
            V_relS = Vo*(1+a[i])
            V_relC = Vt*(1-a_l[i])

            v_abs = (V_relS**2+V_relC**2)**(1/2)

            re = reynolds_number(v_abs,self.chord[i])
            self.reynolds.append(re)

            ma = mach_number(v_abs,self.v_sound)
            self.mach.append(ma)

            Lambda = self.flgiht_speed/vt

            f = (self.number_blades/2)*(1/Lambda)*(1+Lambda**2)**(1/2)*(1-(self.radius[]i/(self.diameter/2)))

            F = (2/np.pi)*(np.arctan((np.exp(2*f)-1)**(1/2)))

            #=============================================
            #todo add xfoil
            #=============================================

            L = (1/2)*self.p*self.chord *Cl*v_rel**2
            Dr = (1/2)*self.p*self.chord*Cd*v_rel**2 
            Cn = Cl*np.cos(Phi)+Cd*np.sin(Phi)
            Ct = Cl*np.sin(Phi)-Cd*np.cos(Phi)
            
            sigma = self.chord*self.number_blades/(2*np.pi*self.radius[i])
            I1 = 4*np.sin(Phi)**2
            I2 = sigma*Cn
            self.a.append(1/((I1/I2)-1))
            
            I3 = 4*np.sin(Phi)*np.cos(Phi)
            I4 = sigma*Ct
            self.a_l.append(1/((I3/I4)+1))

            #hansen correction

            if self.a[i] >self.a_critic :
                a.pop()
                K_h = 4*F*np.sin(Phi)**2/(sigma*Cn)
                a.append((1/2)*(2+K_h*(1-2*ac)-((K_h*(1-2*ac)+2)**2+4*(K_h*ac**2-1))**(1/2)))

            if a[i] <= (1/3):
                Ct = 4*a[i]*(1-[ai])*F
            else:
                Ct = 4*a[i]*(1-(1/4)*(5-3*a[i])*a[i])*F

            self.pn.append((1/2)*self.p*self.chord*Cn*v_rel**2)
            self.pt.append((1/2)*self.p*self.chord*Ct*v_rel**2)

        return

    def force(self):
        return

    def plot_blade(self):
        hub_x=[]
        hub_y=[]
        l_esq=[0,self.chord[0]]
        r_esq=[self.radius[0],self.radius[0]]
        l_dir=[0,self.chord[-1]]
        r_dir=[self.radius[-1],r_dir=[self.radius[-1]]

        plt.figure(1)
        plt.plot(r,c)
        plt.plot(r_esq,l_esq)
        plt.plot(r_dir,l_dir)
        plt.plot(hub_x,hub_y)

        return