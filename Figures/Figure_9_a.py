#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate


A=np.genfromtxt('../Figures/data/Figure_8_all_shells.txt',usecols=[0,1,2,3,4,5,6])
v=A[:,0]
l=A[:,1]
b=A[:,2]
ratio=A[:,3]
S=A[:,4]
sig_l=A[:,5]
sig_b=A[:,6]

finv=interpolate.interp1d(ratio,v,kind='linear')

g_l=interpolate.interp1d(v,l,kind='linear')
g_b=interpolate.interp1d(v,b,kind='linear')

g_sig_l=interpolate.interp1d(v,sig_l,kind='linear')
g_sig_b=interpolate.interp1d(v,sig_b,kind='linear')




v_bestfit=finv(1)

l_bestfit=g_l(v_bestfit)
b_bestfit=g_b(v_bestfit)
sig_l_bestfit=g_sig_l(v_bestfit)
sig_b_bestfit=g_sig_b(v_bestfit)


print 'The best fit values consistent with v_der/v_true=1 and b=-1 are a boost from the LG with:'
print 'magnitude = ', v_bestfit,  ' l = ' , l_bestfit , '+/-' , sig_l_bestfit, 'b = ', b_bestfit , '+/-' , sig_b_bestfit


file = open('../Figures/data/Figure_8_best_fit.txt','w')
s=str(v_bestfit)
s2=str(l_bestfit)
s3=str(b_bestfit)
file.write(s + " " + s2 + " " + s3)
file.close()

