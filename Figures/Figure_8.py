#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate

matplotlib.rcParams['text.usetex'] = True

###################  Figure 8 a #########################

A=np.genfromtxt('../Figures/data/Figure_8_all_shells.txt',usecols=[0,1,2,3,4])
v=A[:,0]
l=A[:,1]
b=A[:,2]
ratio=A[:,3]
S=A[:,4]

B=np.genfromtxt('../Figures/data/Figure_8_all_primed.txt',usecols=[0,1,2,3,4])
ratio_primed=B[:,3]
S_primed=B[:,4]

C=np.genfromtxt('../Figures/data/Figure_8_all_unprimed.txt',usecols=[0,1,2,3,4])
ratio_unprimed=C[:,3]
S_unprimed=C[:,4]


fit_both=np.polyfit(v,ratio,5)
fit_both_plot=np.poly1d(fit_both)

fit_primed=np.polyfit(v,ratio_primed,5)
fit_primed_plot=np.poly1d(fit_primed)

fit_unprimed=np.polyfit(v,ratio_unprimed,5)
fit_unprimed_plot=np.poly1d(fit_unprimed)


plt.figure()

plt.plot(v,fit_primed_plot(v),'-.',color='red',label='primed')

plt.plot(v,fit_both_plot(v),color='black',label='primed+unprimed')

plt.plot(v,fit_unprimed_plot(v),'--',color='blue',label='unprimed')

xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
ylabel(r"$v_{der}/v_{true}$",fontsize=14)
plt.xlim([90,410])
#plt.ylim([1,5])

plt.legend()



plt.savefig("../Figures/Figures/Figure_vder.eps")


###################  Figure 8 b #########################

plt.figure()



fit_both=np.polyfit(v,S,5)
fit_both_plot=np.poly1d(fit_both)

fit_primed=np.polyfit(v,S_primed,5)
fit_primed_plot=np.poly1d(fit_primed)

fit_unprimed=np.polyfit(v,S_unprimed,5)
fit_unprimed_plot=np.poly1d(fit_unprimed)




plt.plot(v,fit_primed_plot(v),'-.',color='red',label='primed')


plt.plot(v,fit_both_plot(v),color='black',label='primed+unprimed')

plt.plot(v,fit_unprimed_plot(v),'--',color='blue',label='unprimed')

xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
ylabel(r"$S/(n-2)$",fontsize=14)
plt.xlim([90,410])
#plt.ylim([1,5])

plt.legend()


plt.savefig("../Figures/Figures/Figure_Sfit.eps")