#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
matplotlib.rcParams['text.usetex'] = True

A=np.genfromtxt('../Figures/data/Figure_4.txt',usecols=[0,1,2,3])
v=A[:,0]
l=A[:,1]
b=A[:,2]
chi_a=A[:,3]

plt.figure()

one_sigma=(np.min(chi_a)+5.89)



fit=np.polyfit(v,chi_a,5)
fit_plot_chi_a=np.poly1d(fit)



plt.plot(v,fit_plot_chi_a(v),color='black')

#####  find ln B confidence intervals



A=np.genfromtxt('../Figures/data/Figure_4_lnB.txt',usecols=[0,1,2,3])
v=A[:,0]
l=A[:,1]
b=A[:,2]
ln_B=A[:,3]


f_lnB=interpolate.interp1d(ln_B[0:size(v)/2],v[0:size(v)/2],kind='linear')

lnB_1=f_lnB(1)
f_lnB=interpolate.interp1d(ln_B[size(v)/2:size(v)],v[size(v)/2:size(v)],kind='linear')

lnB_2=f_lnB(1)

print 740.6-lnB_1, lnB_2-740.6

plt.vlines(lnB_1,min(chi_a)*0.8,max(chi_a)*1.1,linestyles='dashed',color='black')
plt.text(lnB_1-35,27,'ln$B = 1$',rotation=90)


plt.vlines(lnB_2,min(chi_a)*0.8,max(chi_a)*1.1,linestyles='dashed',color='black')
plt.text(lnB_2-35,27,'ln$B = 1$',rotation=90)





xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
ylabel(r"Minimum $\chi^2_a$ across entire sky",fontsize=14)
plt.xlim([0,1500])
plt.ylim([min(chi_a)*0.8,max(chi_a)*1.1])



plt.savefig("../Figures/Figures/Chi_a_min_v.eps")

A=np.genfromtxt('../Figures/data/Figure_5.txt',usecols=[0,1,2,3])
v=A[:,0]
l=A[:,1]
b=A[:,2]
chi_a=A[:,3]



fit=np.polyfit(v,chi_a,5)
fit_plot_chi_a=np.poly1d(fit)


plt.figure()

plt.plot(v,fit_plot_chi_a(v),color='black')
#plt.plot(v,2.3*ones(size(v)),'--',color='black')

xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
ylabel(r"Minimum $\chi^2_b/\nu_b$ across entire sky",fontsize=14)
plt.xlim([0,800])
plt.ylim([0.62,0.71])



plt.savefig("../Figures/Figures/xi_wrt_dir.eps")


A=np.genfromtxt('../Figures/data/Figure_6_a.txt',usecols=[0,1,2,3])
v=A[:,0]
l=A[:,1]
b=A[:,2]
chi_a=A[:,3]



fit=np.polyfit(v,chi_a,5)
fit_plot_chi_a=np.poly1d(fit)

plt.figure()

plt.plot(v,fit_plot_chi_a(v),color='black')
A=np.genfromtxt('../Figures/data/Figure_5_lnB.txt',usecols=[0,1,2,3])
v=A[:,0]
l=A[:,1]
b=A[:,2]
ln_B=A[:,3]


f_lnB=interpolate.interp1d(ln_B,v,kind='linear')

lnB_1=f_lnB(3)

print lnB_1

plt.vlines(lnB_1,min(chi_a)*0.8,250,linestyles='dashed',color='black')
plt.text(lnB_1-35,140,'ln$B = 3$',rotation=90)
lnB_1=f_lnB(5)

print lnB_1

plt.vlines(lnB_1,min(chi_a)*0.8,250,linestyles='dashed',color='black')
plt.text(lnB_1-35,140,'ln$B = 5$',rotation=90)




xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
ylabel(r"$\chi^2_a$ in direction of minimum $\chi^2_b$",fontsize=14)
plt.xlim([0,1000])
plt.ylim([min(chi_a)*0.8,250])




plt.savefig("../Figures/Figures/chi_wrt_xi_dir.eps")


A=np.genfromtxt('../Figures/data/Figure_6_b.txt',usecols=[0,1,2,3])
v=A[:,0]
l=A[:,1]
b=A[:,2]
chi_a=A[:,3]

fit=np.polyfit(v,chi_a,5)
fit_plot_chi_a=np.poly1d(fit)

plt.figure()

plt.plot(v,fit_plot_chi_a(v),color='black')
#plt.plot(v,2.3*ones(size(v)),'--',color='black')

xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
ylabel(r"$\chi^2_b/\nu_b$ in direction of minimum $\chi^2_a$",fontsize=14)
plt.xlim([0,1000])
plt.ylim([0.62,0.76])



plt.savefig("../Figures/Figures/xi_wrt_chi_direction.eps")



