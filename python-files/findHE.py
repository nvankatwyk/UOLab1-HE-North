#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
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
    

#%%
#Constants
Do = 3/8 * .0254 #m                altered
N_tubes = 210                     #altered
L_pipe = 6 * .0254 #m              altered
L_tot = (48-2) * .0254 #m          altered    
Ao = np.pi*Do*N_tubes*L_tot       
w_tube = .022 * .0254 #m
Di = Do - 2*w_tube #m
Ai = np.pi*N_tubes*L_tot*Di      
k_pipe = 16.3  #W/m*K 
Cp_w = 75.244  #J/mol*K
#rho_w = 977.16 #kg/m^3
rho_w_gal = 3.699 #kg/gal
Tsat25 = 403.57 #K
Tsat40 = 414.65 #K
eps = .000045 #m
MW_w = .01801528 #kg/mol    
Ac = np.pi*(Di/2)**2             
F = 1

g = 9.80665  #m/s^2
#%%
#Run Data
#file = "data-files/round1/50GPM-40PSIG.csv"

#P_run = 40
#
#if P_run == 25:
#    Tsat = Tsat25
#    h_fg = h_fg25
#elif P_run == 40:
#    Tsat = Tsat40
#    h_fg = h_fg40
    
  
# %%
#Import Data
#data = pd.read_csv(file)
#Tin_data = data.loc[:,"Inlet Water Temperature (C)"]
#Tout_data = data.loc[:,"Outlet Water Temperature (C)"]
#Vdot_data = data.loc[:,"Water Flowrate (GPM)"]

# %%
#Averages
Tin = 25 + 273.15 #K
Tout = 75 + 273.15 #K  GUESS
Vdot = 200
Tsat = 489.7 #K
h_fg = 38560/MW_w #J/mol
Rfoul = 1.6e-4


#%%
def find_To(To):
    To = To[0]
    Tavg = (Tin+To)/2
    m_dot = Vdot * rho_w_gal/60 #kg/s
    n_dot = m_dot / .018015 #mol/s
    q = n_dot * Cp_w * (To - Tin) # W
    
    ReD = m_dot*Di/mu_w(Tavg)/Ac/N_tubes
    
    Pr = Cp_w*mu_w(Tavg)/cond_w(Tavg)/MW_w
    
    fguess = .01
    f = fsolve(find_f, fguess, args=(eps, Di, ReD))[0]              # POSSIBLE SOURCE OF ERROR
    
    NuD = (f/8*(ReD-1000)*Pr)/(1+12.7*(f/8)**.5*(Pr**(2/3)-1))
    
    ho_guess = 12000  #W/m^2 K
    Ts_guess = 420  #K
    guesses = [ho_guess,Ts_guess]
    ho, Ts = fsolve(find_ho,guesses, args=(Tsat, rho_w, mu_w, cond_w, h_fg, g, L_pipe, Ao, Tavg, Pr, q))
    
    hi = find_hi(Tavg, NuD, L_tot, N_tubes, cond_w)
    UA = (1/hi/Ai + Rfoul/Ai + np.log(Do/Di)/2/np.pi/k_pipe/L_tot/N_tubes + 1/ho/Ao)**-1
    q = UA * F * dTlm_F(Tsat,To,Tin)
    
    eq1 = n_dot*Cp_w*(To - Tin) - q

    return eq1


To_guess = Tout
To = fsolve(find_To,To_guess)[0]
print(To-273.15)




