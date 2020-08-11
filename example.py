#Author: Lucas Marques Moreno
#Date: 06/2020
import BEM as bem
import numpy as np


# To start using BEM method we must define some values

# Atmosphere Constantes
g = 9.81 #[m/s²] gravity acceleration
p = 1.08 #[kg/m³] air's specific mass
u = 1.4607*10**(-5) #[N.s/m²]  air's dynamic viscosity

#============================================================
# Hansen Method
#============================================================

R_max =[.4,.42,.85,.95,1.5] #[m] radius points to change airfoil
Ang_at =[3.,3.,8.5,8.6,9.] #[degree] angle of attack of each airfoil

r = [0.075, 0.25, 1.5] # list for radius values
c=[0.152, 0.200, 0.1] # list for chord values

# Geometric Specifications
D = 2*r[-1]#[m]propeller diameter
B = 2. #[] blade quantity


airfoils = ["clarky.dat",
            "e193.dat",
            "naca2412.dat"]

# Operation Settings
rpm = 94.73845432 # [rpm] rounds / minute
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

#============================================================
# Adkins-Liebeck and Attickins Method optimun design 
#============================================================

R_max =[.4,.42,.85,.95,1.5] #[m] radius points to change airfoil
Ang_at =[3.,3.,8.5,8.6,9.] #[degree] angle of attack of each airfoil

r = [0.075, 0.25, 1.5] # list for radius values
c=[0.152, 0.200, 0.1] # list for chord values

# Geometric Specifications
D = 2*r[-1]#[m]propeller diameter
B = 2. #[] blade quantity


airfoils = ["clarky.dat",
            "e193.dat",
            "naca2412.dat"]

# Operation Settings
rpm = 94.73845432 # [rpm] rounds / minute
vo = 12 #[m/s] flight speed


propeller = bem.blade(vo, rpm, B, D, r,c,Ang_at,airfoils,R_max,g,p,u)
propeller.config(power_constraint = 300)
propeller.liebeck_optimun_design()
print(propeller.Thrust)


#============================================================
# Adkins-Liebeck and Attickins Method design - Analysis 
#============================================================

R_max =[.4,.42,.85,.95,1.5] #[m] radius points to change airfoil
Ang_at =[3.,3.,8.5,8.6,9.] #[degree] angle of attack of each airfoil

r = [0.075, 0.25, 1.5] # list for radius values
c=[0.152, 0.200, 0.1] # list for chord values

# Geometric Specifications
D = 2*r[-1]#[m]propeller diameter
B = 2. #[] blade quantity


airfoils = ["clarky.dat",
            "e193.dat",
            "naca2412.dat"]

# Operation Settings
rpm = 94.73845432 # [rpm] rounds / minute
vo = 12 #[m/s] flight speed

propeller = bem.blade(vo, rpm, B, D, r,c,Ang_at,airfoils,R_max,g,p,u)
propeller.config(max_ite = 10)
propeller.liebeck_analysis_design()
print(propeller.Thrust)
print(propeller.efficiency)