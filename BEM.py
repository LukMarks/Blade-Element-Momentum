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

        self.Phi = []
        self.theta =[] # twist angle

        return

    def config(self):
        pass

        return

    def advance_ratio(self):
        
        self.Jo = self.flgiht_speed*(self.diameter/(self.rpm/60))
        
        return
    
    def start(self):


        for i in range(self.radius):
            r = self.radius[i]
            alpha = alphaD*np.pi/180
            A = a[-1]
            a_l = a[-1]
            vt = self.w * r
            Tan_phi = Tan_phi = ((1+A)*self.flgiht_speed)/((1-A_l)*Vt)
            phi = np.arctan(Tan_phi)
            self.Phi.append(phi)


        return