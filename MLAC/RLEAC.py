# -*- coding: utf-8 -*-
"""
Reinforcement Learning Enhanced Analytic Continuation
RLEAC

@author: Junyin Zhang (Discussed with Bo Zhang)
"""

import numpy as np
from matplotlib import pyplot as plt
TempSlice = 16
Beta = 2.0
D_Tau = Beta / TempSlice
N_epoch = int(1e4)
N_warm = int(1e5)
step = 1e-2
E_min = -20
E_max = 20
dE = 0.2
E_range = np.arange(E_min,E_max+dE,dE)
G_k = np.zeros([TempSlice,1]);
Gauge = np.zeros([TempSlice,len(E_range)])
Theta = 1e6;
diff = 1.0
#######################################################
A = np.zeros([len(E_range),1])
for index in range(len(E_range)):
    omega = E_range[index]
    A[index] = np.exp(-0.2*omega**2)
    
A = A/(np.sum(A)*dE)
D_Tau = Beta/TempSlice
A_ori = A*1.0
for time_index in range(TempSlice):
    for E_index in range(len(E_range)):
        tau = (time_index-1)*D_Tau
        omega = E_range[E_index]
        Gauge[time_index,E_index] = dE * np.exp(-omega*(tau - Beta/2))/(2*np.cosh(Beta*omega/2))

G_k = np.dot(Gauge,A)
for E_index in range(len(E_range)):
    A[E_index] = np.random.rand()

A = A/(np.sum(A)*dE)
############################## Calculation Process ############################3
A_collect = np.zeros([len(E_range),N_epoch])
for epoch_index in range(N_warm):
    if np.mod(epoch_index,10)==0:
        print("Warm Ratio = %f, diff=%e\n"%(epoch_index/N_warm,diff))
    for A_index in range(len(E_range)):
        A_new = A*1.0
        A_new[A_index] = A_new[A_index]+ step*(np.random.rand()-0.5)
        A_new = A_new / (np.sum(A_new)*dE)
        diff = np.linalg.norm(np.dot(Gauge,A) - G_k)
        value_old = np.sum((np.dot(Gauge,A) - G_k)**2)
        value_new = np.sum((np.dot(Gauge,A_new) - G_k)**2)
        step = np.sqrt(diff)
        if np.random.rand()<np.exp(Theta*(value_old - value_new)) and A_new[A_index]>=0:
            A = A_new*1.0
            
for epoch_index in range(N_epoch):
    if np.mod(epoch_index,10)==0:
        print("MC Ratio = %f, diff=%e\n"%(epoch_index/N_warm,diff))
    for A_index in range(len(E_range)):
        A_new = A*1.0
        A_new[A_index] = A_new[A_index]+ step*(np.random.rand()-0.5)
        A_new = A_new / (np.sum(A_new)*dE)
        diff = np.linalg.norm(np.dot(Gauge,A) - G_k)
        value_old = np.sum((np.dot(Gauge,A) - G_k)**2)
        value_new = np.sum((np.dot(Gauge,A_new) - G_k)**2)
        step = np.sqrt(diff)
        if np.random.rand()<np.exp(Theta*(value_old - value_new)) and A_new[A_index]>=0:
            A = A_new*1.0
    A_collect[:,epoch_index] = A.reshape(len(E_range))

A_final = np.zeros([len(E_range),1])
for index in range(len(E_range)):
    A_final[index] = np.mean(A_collect[index,:]);
   

plt.plot(E_range,A_final)
plt.plot(E_range,A_ori)