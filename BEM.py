import numpy as np



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

        self.w = self.rpm*2*np.pi/60 #angular speed

        self.phi = []
        self.theta = [] # twist angle

        self.reynolds = []
        self.mach = []
        self.v_sound = 343 #[m/s]
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

            #v_rel = (self.flgiht_speed**2+(Vt)**2)**(1/2)
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

                    

        return