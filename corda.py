from bezier import bezier_profiles
class corda:
    def __init__(self,xa,xi1,xf1,xi2,xf2,ci,Mode,Rf):
        
        self.xa = xa
        self.xi1 = xi1
        self.xf1 = xf1
        self.xi2 = xi2
        self.xf2 = xf2
        self.ci = ci
        self.Mode = Mode
        self.Rf = Rf
        
        def frange(start, stop, step):
         i = start
         while i < stop:
             yield i
             i += step
        
        if Mode == 'cte' :
            if xa <= xf2:
                self.c = ci
            else:
                R = xi2-xf2
                r = xa-xf2
                self.c = ci*r/R  
            
        
        if Mode == 'trapz' :
            
            if xa <= xf1:
                r = xa-xi1
                R = xf1-xi1
                self.c = ci*r/R
                
            
            elif xa >= xi2:
                #dx =-ci/(xf-xi)
                R = xi2-xf2
                r = xa-xf2
                self.c = ci*r/R               
                
            else:
                self.c=ci
        
        if Mode == 'quad':
            rz = xa/Rf
            if 0<=rz<.40:
                self.c = -1.4444*xa**2 + 1.5033*xa + 0.0125
            elif .40<=rz<.8:
                self.c = -0.6667*xa**2 + 1.06*xa -0.031
            elif rz>=.8:
                
                self.c = -0.6667*xa + 1.1
                
        if Mode == 'koch':
            self.c = (0.1767*xa**3 - 0.9575*xa**2 + 1.0554*xa + 0.0813)*0.7
            
            #self.c = (a1*xa**3 - a2*xa**2 + a3*xa + a4)*0.7
        
        if Mode == 'grf':
            i = xa/Rf

            if i <=.25:
                c_aux = (.08/.25)*i
                self.c = c_aux*Rf*1.5
            else:
                 c_aux =((-.05/.75)*(i-0.25)+.08)*1.5
                 self.c = c_aux*Rf
                 
        if Mode == 'grf2':
             i = xa/Rf
             if i <=.2:
                c_aux = ((4.7/61.5)/.2)*i
                self.c = c_aux*Rf*(1.25)
             elif .2 < i <=.93:
                 c_aux =(((2.2-4.7)/61.5)/(.93-.2))*(i-.2)+(4.7/61.5)
                 self.c = c_aux*Rf*(1.25)
             else:
                 c_aux =(((.4-2.2)/61.5)/(1-.93))*(i-.93)+(2.2/61.5)
                 self.c = c_aux*Rf*(1.25)
            
