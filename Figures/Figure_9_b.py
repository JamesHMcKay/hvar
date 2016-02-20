#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate


matplotlib.rcParams['text.usetex'] = True

#def lamcutoff(LB):
#    return -0.065/(1-0.01*log(246.221/LB))

def boost_offset(r,a,b):
    y=a*(r**b)
    return y

A=np.genfromtxt('../Figures/data/Figure_9_a.txt',usecols=[0,1,2,3])
r2_p=A[0:10,0]
Del_H_p=A[0:10,1]
sig_r2_p=A[0:10,2]
sig_Del_H_p=A[0:10,3]

r2=A[11:21,0]
Del_H=A[11:21,1]
sig_r2=A[11:21,2]
sig_Del_H=A[11:21,3]



B=np.genfromtxt('../Figures/data/Figure_9_details.txt',usecols=[0,1])
b=B[:,0]
a=B[:,1]


r2_plot=np.linspace(0.00001,20000,2000)

Del_H_plot=zeros(size(r2_plot))
for i in range(0,2000):
    Del_H_plot[i]=boost_offset(r2_plot[i],a[0],b[0])

Del_H_plot_primed=zeros(size(r2_plot))

for i in range(0,2000):
    Del_H_plot_primed[i]=boost_offset(r2_plot[i],a[1],b[1])
Del_H_plot_unprimed=zeros(size(r2_plot))

for i in range(0,2000):
    Del_H_plot_unprimed[i]=boost_offset(r2_plot[i],a[2],b[2])

plt.figure()

plt.plot(r2_plot,Del_H_plot,color='black')

plt.plot(r2_plot,Del_H_plot_primed,'--',color='black')

plt.plot(r2_plot,Del_H_plot_unprimed,'--',color='black')

plt.plot(r2_plot,zeros(size(r2_plot)),color='black')
plt.errorbar(r2,Del_H,sig_Del_H,sig_r2,fmt='o',color='black')
plt.errorbar(r2_p,Del_H_p,sig_Del_H_p,sig_r2_p,fmt='o',color='black',markerfacecolor="none")

r=np.linspace(0,5e3,5)

plt.plot(r,zeros(size(r)),color='black')

xlabel(r"$<r_s^2>$ ($h^{-1}$ Mpc )$^2$",fontsize=14)
ylabel(r"$H_{LG}-H_{X}$ ($h$ km/sec/Mpc)",fontsize=14)
plt.xlim([0,5e3])
plt.ylim([-4,10])

plt.savefig("../Figures/Figures/boost_offset_122_from_LG1.eps")

B=np.genfromtxt('../Figures/data/Figure_9_b.txt',usecols=[0,1,2])
r_p=B[0:10,0]
dH_p=B[0:10,1]
sig_dH_p=B[0:10,2]

r_u=B[11:21,0]
dH_u=B[11:21,1]
sig_dH_u=B[11:21,2]

plt.figure()

r=np.linspace(0,150,5)

plt.errorbar(r_p,dH_p,xerr=0,yerr=sig_dH_p,fmt='o',color='black')
plt.errorbar(r_u,dH_u,xerr=0,yerr=sig_dH_u,fmt='o',color='black',markerfacecolor="none")
plt.plot(r,zeros(size(r)),color='black')
plt.xlim([0,140])
plt.ylim([-0.05,0.25])

xlabel(r"$<r_s>$ ($h^{-1}$ Mpc )",fontsize=14)
ylabel(r"$\delta H_{LG}$",fontsize=14)


plt.savefig("../Figures/Figures/hvar_LG1.eps")

C=np.genfromtxt('../Figures/data/Figure_9_c.txt',usecols=[0,1,2])
r_p=C[0:10,0]
dH_p=C[0:10,1]
sig_dH_p=C[0:10,2]

r_u=C[11:21,0]
dH_u=C[11:21,1]
sig_dH_u=C[11:21,2]

plt.figure()

plt.errorbar(r_p,dH_p,xerr=0,yerr=sig_dH_p,fmt='o',color='black')
plt.errorbar(r_u,dH_u,xerr=0,yerr=sig_dH_u,fmt='o',color='black',markerfacecolor="none")
plt.plot(r,zeros(size(r)),color='black')
plt.xlim([0,140])
plt.ylim([-0.05,0.25])

xlabel(r"$<r_s>$ ($h^{-1}$ Mpc )",fontsize=14)
ylabel(r"$\delta H_{X}$",fontsize=14)

plt.savefig("../Figures/Figures/hvar_Xboost1.eps")







