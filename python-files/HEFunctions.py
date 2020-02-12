#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

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
    mu = np.exp(A + B/Tavg + C*np.log(Tavg) + D*Tavg**E)
    return mu

def conductivity(T):
    A = -.432
    B = .0057255
    C = -8.078e-6
    D = 1.861e-9
    E = 0
    cond = A + B*T + C*T**2 + D*T**3 + E*T**4
    return cond
