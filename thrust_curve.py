#Author: Lucas Marques Moreno
#Date: 06/2020
import BEM as bem
import numpy as np
from bezier import bezier_profiles
import matplotlib.pyplot as plt

# rpm 90 - 190 # 95% de eficiencia trasmissão faixa = 0 - 17 m/s opt = 10m/s
R_max = [.4,.42,.85,.95,1.5]
Ang_at = [3.,3.,8.5,8.6,9.]
cpx =  [0.1,.5,.75,1.2,1.6]#[0.1,.4,.75,.9,1.6]#coordenadas dos pontos Em X
order = len(cpx)-1
npts = 20
cpy = [.15,0.24999986,0.34668816,0.34996463,.1]#[.1,.25,.35,.25,.01]#[0.1,.25,.31,.25,.01]
prop = bezier_profiles(cpx,cpy,order,npts)
r = []#[m]lista do raio da hélice
c=[]#[m] lista com os comprimentos da corda
for i in range(0,len(prop[0])):
    r.append(float(prop[0][i]))
    c.append(float(prop[1][i]))
g = 9.81#[m/s²] aceleração da gravidade
p = 1.08#[kg/m³] massa especifica do ar(1.225)
u = 1.4607*10**(-5)#[N.s/m²] viscosidade dinamica do a 18.20*10**(-6) #
D = 2*cpx[len(cpx)-1]
B = 2.#numero de pás
rpm = 94.73845432

vo = 12

airfoils = ["profile0.out",
            "profile1.out",
            "profile2.out",
            "profile3.out",
            "profile4.out"]

propeller = bem.blade(vo, rpm, B, D, r,c,Ang_at,airfoils,R_max,g,p,u)
propeller.config()
propeller.induced_factor()
propeller.forces()
propeller.Efficiency()
print(propeller.thrust) # rpm 90 - 190 # 95% de eficiencia trasmissão faixa = 0 - 17 m/s opt = 10m/s
print(propeller.efficiency)
print(propeller.Ct())
print(propeller.Cp())
print(propeller.advance_ratio())
propeller.export(propeller.theta, name = "theta_ref")

theta_ref = []

file = open("theta_ref.dat",'r')
for line in file:
    theta_ref.append(float(line.strip().split()[0]))
print(theta_ref)

thrust = []
v = []

for vo in range(7,21):
    print(vo)
    curve_propeller = bem.blade(vo, rpm, B, D, r,c,Ang_at,airfoils,R_max,g,p,u)
    curve_propeller.config(speed_test = True, theta_reference = theta_ref)
    curve_propeller.induced_factor()
    curve_propeller.forces()
    thrust.append(curve_propeller.thrust)
    v.append(vo)

plt.figure(1)
plt.plot(v,thrust)
plt.grid()
plt.show()