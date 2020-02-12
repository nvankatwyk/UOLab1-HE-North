#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

#%%
# Functions

def dTlm_F(Tsat,To,Ti):
    num = (Tsat-To) - (Tsat-Ti)
    dom = np.log((Tsat-To) / (Tsat-Ti))
    return num/dom
    
    
def mu_w(T):
    A = -52.843
    B = 3703.6
    C = 5.866
    D = -5.879e-29
    E = 10
    mu = np.exp(A + B/T + C*np.log(T) + D*T**E) #P*s
    return mu

def cond_w(T):
    A = -.432
    B = .0057255
    C = -8.078e-6
    D = 1.861e-9
    E = 0
    cond = A + B*T + C*T**2 + D*T**3 + E*T**4 #W/m K
    return cond

def find_f(f, esp, Di, ReD):
    func = -2*np.log10((eps/Di)/3.7 + 2.51/ReD/np.sqrt(f)) - 1/np.sqrt(f)
    return func

def find_ho(var, Tsat, rho_w, mu_w, cond_w, h_fg, g, L_pipe, Ao, Tavg, Pr, q):
    ho = var[0]
    Ts = var[1]
    P = cond_w(Tavg) * L_pipe * (Tsat - Ts) / mu_w(Tavg) / h_fg / ((mu_w(Tavg)/rho_w(Tavg))**2/g)**(1/3)

    if P<=2530 and P>= 15.8:
        Nu = 1/P * (0.68*P+0.89)**0.82
        print("P is between")
    elif P<15.8:
        Nu = 0.943*P**(-1/4)
        print("P is less than 15.8")
    elif P>2530 and Pr>=1:
        Nu = 1/P * ((0.024*P-53)*np.sqrt(Pr)+89)**(4/3)
        print("P is greater than 2530")
    else:
        print("Conditions not met in find_ho")
        
    eq1 = Nu*cond_w(Tavg)/((mu_w(Tavg)/rho_w(Tavg))**2/g)**(1/3) - ho
    eq2 = (Tsat - Ts)/(1/(ho*Ao))-q
    
    return [eq1, eq2]

def find_hi(Tavg, NuD, L_tot, N_tubes, cond_w):
    return cond_w(Tavg)*NuD/L_tot*N_tubes

def find_fric(Rfoul, UA, hi, Ai, Do, Di, k_pipe, L_pipe, ho, Ao):
    eq = 1/hi/Ai + Rfoul/Ai + np.log(Do/Di)/2/np.pi/k_pipe/L_tot/N_tubes + 1/ho/Ao - 1/UA
    return eq

def rho_w(T):
    Tc = 647.096
    T = 1-T/Tc
    A = 17.863
    B = 58.606
    C = -95.396
    D = 213.89
    E = -141.26
    rho = A + B*T**.35 + C*T**(2/3) + D*T + E*T**(4/3)
    return rho*18.01528
    

#def find_To(To, Tsat, Ti, m_dot, Cp_w, Do, Di, hi, ):
#    
#    UA = (1/hi/Ai + Rfoul/Ai + np.log(Do/Di)/2/np.pi/k_pipe/L_tot/N_tubes + 1/ho/Ao)**-1
#    q = UA * F * dTlm_F(Tsat,To,Ti)
#    
#    eq1 = m_dot*Cp_w*(To - Ti) - q
#    
#    return eq1
    
#%%
#Constants
Do = .25 * .0254 #m               altered
N_tubes = 56                     #altered
#kr = .015
L_pipe = 4 * .0254 #m             altered
L_tot = 12 * .0254 #m             
Ao = np.pi*Do*N_tubes*L_tot      #altered
Ao_4 = np.pi*Do*N_tubes*L_pipe   #altered
w_tube = .022 * .0254 #m
Di = .25*.0254 - 2*w_tube #m     #altered
Ai = np.pi*N_tubes*L_tot*Di      #altered with ntubes and Di
k_pipe = 16.3  #W/m*K 
Cp_w = 75.364  #J/mol*K
#rho_w = 977.16 #kg/m^3
rho_w_gal = 3.699 #kg/gal
Tsat25 = 403.57 #K
Tsat40 = 414.65 #K
eps = .000045 #m
MW_w = .01801528 #kg/mol
h_fg25 = 39159/MW_w #J/mol       altered
h_fg40 = 38560/MW_w #J/mol       altered
Ac = np.pi*(Di/2)**2             #altered with Di
F = 1

g = 9.80665  #m/s^2
#%%
#Run Data
file = "data-files/round1/50GPM-40PSIG.csv"

P_run = 40

if P_run == 25:
    Tsat = Tsat25
    h_fg = h_fg25
elif P_run == 40:
    Tsat = Tsat40
    h_fg = h_fg40
    
  
# %%
#Import Data
data = pd.read_csv(file)
Tin_data = data.loc[:,"Inlet Water Temperature (C)"]
Tout_data = data.loc[:,"Outlet Water Temperature (C)"]
Vdot_data = data.loc[:,"Water Flowrate (GPM)"]

# %%
#Averages
Tin = Tin_data.mean() + 273.15 #K
Tout = Tout_data.mean() + 273.15 #K
Tavg = (Tin+Tout)/2
Vdot = Vdot_data.mean()


#%%
mdot = Vdot * rho_w_gal/60 #kg/s
ndot = mdot / .018015 #mol/s
q = ndot * Cp_w * (Tout - Tin) # W
dTlm = dTlm_F(Tsat,Tout,Tin)

UA = q/dTlm  #J/K

ReD = mdot*Di/mu_w(Tavg)/Ac/N_tubes

Pr = Cp_w*mu_w(Tavg)/cond_w(Tavg)/MW_w

fguess = .01
f = fsolve(find_f, fguess, args=(eps, Di, ReD))[0]              # POSSIBLE SOURCE OF ERROR

NuD = (f/8*(ReD-1000)*Pr)/(1+12.7*(f/8)**.5*(Pr**(2/3)-1))

hi = find_hi(Tavg, NuD, L_tot, N_tubes, cond_w) #W/m^2 K

ho_guess = 8000  #W/m^2 K
Ts_guess = 380  #K
guesses = [ho_guess,Ts_guess]
ho, Ts = fsolve(find_ho,guesses, args=(Tsat, rho_w, mu_w, cond_w, h_fg, g, L_pipe, Ao, Tavg, Pr, q))
ho_tot = ho*N_tubes/L_pipe*L_tot

foul_guess = 0.00007
fouling_factor = fsolve(find_fric,foul_guess, args=(UA, hi, Ai, Do, Di, k_pipe, L_pipe, ho_tot, Ao))[0]

print("\nFouling Factor for ",file," =", fouling_factor)

#%%

Tsat = 489.7 #K
Ti = 25 + 273.15 #K




















