#Author: Lucas Marques Moreno
#Date: 06/2020
import BEM as bem
import numpy as np
from bezier import bezier_profiles


# To start using BEM method we must define some values

# Atmosphere Constantes
g = 9.81 #[m/s²] gravity acceleration
p = 1.08 #[kg/m³] air's specific mass
u = 1.4607*10**(-5) #[N.s/m²]  air's dynamic viscosity

# Geometric Specifications
D = 2*cpx[len(cpx)-1] #[m]propeller diameter
B = 2. #[] blade quantity

R_max =[.4,.42,.85,.95,1.5] #[m] radius points to change airfoil
Ang_at =[3.,3.,8.5,8.6,9.] #[degree] angle of attack of each airfoil

cpx = [0.1,.5,.75,1.2,1.6] #[m] x componentes for bezier coordenates
order = len(cpx)-1
npts = 20 #[] number of points in bezier's curve
cpy = [.15,0.24999986,0.34668816,0.34996463,.1] #[m] y componentes for bezier coordenates
prop = bezier_profiles(cpx,cpy,order,npts)
r = [] # list for radius values
c=[] # list for chord values
for i in range(0,len(prop[0])):
    r.append(float(prop[0][i]))
    c.append(float(prop[1][i]))

airfoils = ["profile0.out",
            "profile1.out",
            "profile2.out",
            "profile3.out",
            "profile4.out"]

# Operation Settings
rpm = 94.73845432 # [rpm] rotation / minute
vo = 12 #[m/s] flight speed

propeller = bem.blade(vo, rpm, B, D, r,c,Ang_at,airfoils,R_max,g,p,u)
propeller.config()
propeller.induced_factor()
propeller.forces()
propeller.Efficiency()
print(propeller.thrust) 
print(propeller.efficiency)
print(propeller.Ct())
print(propeller.Cp())
print(propeller.advance_ratio())
