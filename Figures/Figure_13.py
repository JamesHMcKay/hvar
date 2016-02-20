#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
#def lamcutoff(LB):
#    return -0.065/(1-0.01*log(246.221/LB))

def boost_offset(r,a,b):
    y=a*(r**b)
    return y

A=np.genfromtxt('../Figures/data/Figure_13.txt',usecols=[0,1,2,3])
r2_p=A[0:10,0]
Del_H_p=-A[0:10,1]
sig_r2_p=A[0:10,2]
sig_Del_H_p=A[0:10,3]

r2=A[11:21,0]
Del_H=-A[11:21,1]
sig_r2=A[11:21,2]
sig_Del_H=A[11:21,3]



B=np.genfromtxt('../Figures/data/Figure_13_details.txt',usecols=[0,1])
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

xlabel(r"$<r_s^2>$ ($h^{-1}$ Mpc )$^2$",fontsize=14)
ylabel(r"$H_{CMB}-H_{LG}$ ($h$ km/sec/Mpc)",fontsize=14)
plt.xlim([0,20e3])
plt.ylim([-5,21])



plt.savefig("../Figures/Figures/Figure_13.eps")


