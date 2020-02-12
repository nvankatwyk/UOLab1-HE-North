#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import cantera as cant

## Constants ##
Do = .25 #in
N_tubes = 56
kr = .015
L_tot = 14 #in
L = 4 #in
Ao = np.pi*Do*N_tubes*L_tot
w_tube = .022
Di = .25 - 2*w_tube
Ai = np.pi*N_tubes*L_tot
Kp = 16.3
Cp_w_in_tube = 75.364
rho_w = 54240
Tsat25 = 403.57
Tsat35 = 411.266
E_rough = .045
h_fg25 = 39159
h_fg35 = 38746
