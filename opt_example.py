#Author: Lucas Marques Moreno
#Date: 06/2020
import BEM as bem
import numpy as np
from bezier import bezier_profiles
from scipy.optimize import minimize
from time import time

# max and min values to optmization parameters
cpy = [(.15,0.25), # [m] y coordenates for 2nd point
       (.20,.35), # [m] y coordenates for 3rd point
       (.30,.35), # [m] y coordenates for 4th point
       (0,10), # [degree] angle of attack for the 1st airfoil
       (0,10), # [degree] angle of attack for the 2nd airfoil
       (0,10), # [degree] angle of attack for the 3rd airfoil
       (0,10), # [degree] angle of attack for the 4th airfoil
       (0,10), # [degree] angle of attack for the 5th airfoil
       (90.,190.)] # [rpm] propeller's rotational speed

# initial point's
x0 = np.array([.2,.3,.30,3.,3.,8.5,8.6,9.,90.]) 
x1 = np.array([.2,.3,.30,4.,4.,6.,7.5,8.,90.])
x2 = np.array([.15,.20,.32,1.,1.,2.,2.,2.,130.])

def objective(x):
    
    g = 9.81#[m/s²] aceleração da gravidade
    p = 1.08#[kg/m³] massa especifica do ar(1.225)
    u = 1.4607*10**(-5)#[N.s/m²] viscosidade dinamica do a 18.20*10**(-6) #

    transmission_efficiency = 0.95
    R_max =[.4, .42, .85, .95, 1.5]
    Ang_at =[x[3], x[4], x[5], x[6], x[7]]
    cpx = [0.1,.5,.75,1.2,1.6]#[0.1,.4,.75,.9,1.6]#coordenadas dos pontos Em X
    order = len(cpx)-1
    npts = 20
    cpy = [.15, x[0], x[1], x[2], .1]#[.1,.25,.35,.25,.01]#[0.1,.25,.31,.25,.01]
    prop = bezier_profiles(cpx,cpy,order,npts)
    r = []#[m]lista do raio da hélice
    c=[]#[m] lista com os comprimentos da corda
    for i in range(0,len(prop[0])):
        r.append(float(prop[0][i]))
        c.append(float(prop[1][i]))
    D = 2*cpx[len(cpx)-1]
    B = 2.#numero de pás

    vo = 10

    airfoils = ["profile0.out",
                "profile1.out",
                "profile2.out",
                "profile3.out",
                "profile4.out"]

    propeller = bem.blade(vo, x[-1], B, D, r,c,Ang_at,airfoils,R_max,g,p,u)
    propeller.config()
    propeller.induced_factor()
    propeller.forces()
    blade_efficiency = propeller.Efficiency()
    total_power = propeller.power_flight/(blade_efficiency*transmission_efficiency)
    return total_power


def force_minimum(x):
    g = 9.81#[m/s²] aceleração da gravidade
    p = 1.08#[kg/m³] massa especifica do ar(1.225)
    u = 1.4607*10**(-5)#[N.s/m²] viscosidade dinamica do a 18.20*10**(-6) #

    R_max =[.4, .42, .85, .95, 1.5]
    Ang_at =[x[3], x[4], x[5], x[6], x[7]]
    cpx = [0.1,.5,.75,1.2,1.6]#[0.1,.4,.75,.9,1.6]#coordenadas dos pontos Em X
    order = len(cpx)-1
    npts = 20
    cpy = [.15, x[0], x[1], x[2], .1]#[.1,.25,.35,.25,.01]#[0.1,.25,.31,.25,.01]
    prop = bezier_profiles(cpx,cpy,order,npts)
    r = []#[m]lista do raio da hélice
    c=[]#[m] lista com os comprimentos da corda
    for i in range(0,len(prop[0])):
        r.append(float(prop[0][i]))
        c.append(float(prop[1][i]))
    D = 2*cpx[len(cpx)-1]
    B = 2.#numero de pás

    vo = 10

    airfoils = ["profile0.out",
                "profile1.out",
                "profile2.out",
                "profile3.out",
                "profile4.out"]

    propeller = bem.blade(vo, x[-1], B, D, r,c,Ang_at,airfoils,R_max,g,p,u)
    propeller.config()
    propeller.induced_factor()
    propeller.forces()
    result = -60+ propeller.thrust
    return 
force_constraint = ({'type':'ineq', 'fun':force_minimum})

start_time = time()



print("1st Guess")

solution = minimize(objective, x0, method = 'SLSQP', bounds=cpy,
                    constraints = [force_constraint], options={'disp':True,'maxiter':200})


xopt = solution.x
Potf = solution.fun
P = objective(solution.x)
thrust_solution = force_minimum(solution.x)

print('parameter: ',xopt)
print('final Pot: ',Potf,'W')
print('final thrust: ', thrust_solution,'N')

time_elapsed = time()-start_time
minutes, seconds = divmod(time_elapsed, 60)
print('time elapsed m:s ',round(minutes,0),':', round(seconds,0))

print("2nd Guess")

solution = minimize(objective, x1, method = 'SLSQP', bounds=cpy,
                    constraints = [force_constraint], options={'disp':True,'maxiter':200})


xopt = solution.x
Potf = solution.fun
P = objective(solution.x)
thrust_solution = force_minimum(solution.x)

print('parameter: ',xopt)
print('final Pot: ',Potf,'W')
print('final thrust: ', thrust_solution,'N')

time_elapsed = time()-start_time
minutes, seconds = divmod(time_elapsed, 60)
print('time elapsed m:s ',round(minutes,0),':', round(seconds,0))

print("3rd Guess")

solution = minimize(objective, x2, method = 'SLSQP', bounds=cpy,
                    constraints = [force_constraint], options={'disp':True,'maxiter':200})


xopt = solution.x
Potf = solution.fun
P = objective(solution.x)
thrust_solution = force_minimum(solution.x)

print('parameter: ',xopt)
print('final Pot: ',Potf,'W')
print('final thrust: ', thrust_solution,'N')

time_elapsed = time()-start_time
minutes, seconds = divmod(time_elapsed, 60)
print('time elapsed m:s ',round(minutes,0),':', round(seconds,0))