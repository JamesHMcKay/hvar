from mpl_toolkits.basemap import Basemap, shiftgrid
from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Hershey-Gothic-English']})
rc('font',**{'family':'serif','serif':['Times-New-Roman']})
rc('text', usetex=True)
rc('text',fontsize=10)
#import hwtCMB
#import hwtSUNLG
import sys
import math
#from sympy.solvers import nsolve
#from sympy import Symbol, sqrt
from scipy import interpolate
from matplotlib.ticker import FormatStrFormatter

Tcmb=2726.0 # mK
c=299792.458 # km/s

def weighted_avg_and_std(values, weights):
    """
    Returns the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    av = average(values, weights=weights)
    var = dot(weights, (values-av)**2)/weights.sum()  
    return (av, math.sqrt(var))

#hwtCMB.hvar()
#hwtSUNLG.hvar()

clevs=80 # Set number of H contours, keep this, adjusts the smoothness

cmb_instr='../Figures/data/Figure_23_a_chi2.txt'
cmb_hdat = loadtxt(cmb_instr)
cmb_hdat=array([float(s) for s in cmb_hdat[:,2]])

## Plot CMB temperature dipole
cmb_hdatt = transpose(reshape(cmb_hdat,(360, 180)))


cmb2_instr='../Figures/data/Figure_23_a.txt'
cmb2_hdat = loadtxt(cmb2_instr)
cmb2_hdat=array([float(s) for s in cmb2_hdat[:,2]])

## Plot CMB temperature dipole
cmb2_hdatt = transpose(reshape(cmb2_hdat,(360, 180)))





lons=linspace(360.,0.,num=360)
lats=linspace(90.,-90.,num=180)
clons, clats = meshgrid(lons, lats)
# create new figure
cmbfig=plt.figure(1)

#other option for projection is 'hammer'
m = Basemap(resolution='c',projection='moll',lon_0=180)
ax = cmbfig.add_axes([0.1,0.1,0.8,0.8])
cmbx, cmby = m(180+360-96.44,29.3) 
p3=m.plot(cmbx,cmby,'kx',color='white',markersize=20,mew=2)









x, y = m(clons, clats) 
#CSl = m.contour(x,y,hdat,clevs,linewidths=0.5,colors='k')


# The following sets up the colorbar down the bottom of the figure
# setup colorbar axes instance.
cmap=plt.cm.jet
#CSf = m.contourf(x,y,cmb_hdatt,clevs,edgecolors='none',cmap=cmap)
clevs=5
levels=[1,3,5]
#CSl = m.contour(x,y,cmb_hdatt,levels,linewidths=1,colors='k')
clevs=15
#Good colormaps: Greys (bw), summer, PuBu, YlGn
CSf = m.contourf(x,y,cmb_hdatt,clevs,edgecolors='none',drawedges=False,cmap=plt.cm.PuBu)


CSl2 = m.contour(x,y,cmb2_hdatt,levels,edgecolors='none',drawedges=False,linewidths=0.5,colors='k')






#CS4=plt.contour(x,y,cmb_hdatt,levels,colors=('k'),linewidths=(0.5,))

pos = ax.get_position()
l, b, w, h = pos.bounds
#cax = plt.axes([l+w+0.05, b, 0.025, h]) # setup colorbar axes
cax = plt.axes([l, b, w,0.025])
l_f = FormatStrFormatter('%f2.1')
cbar = plt.colorbar(CSf, drawedges=False, cax=cax, orientation='horizontal') # draw colorbar (the thing down the bottom)
cbar.ax.yaxis.set_major_formatter(l_f)
cbar.set_label(r'$\chi^2_a$',fontsize=20)
cbar.ax.tick_params(labelsize=20)


		
   
#plt.clabel(CSl,inline=1,fontsize=9,fmt='%1.2f')
plt.axes(ax)  # make the original axes current again

# Draw the graticule
#m.drawmapboundary(linewidth=0.1)

#draw parallels and meridians.
parallels = arange(-60.,90,30.)
m.drawparallels(parallels,linewidth=.5,dashes=[1,1])
meridians = arange(-360,360.,30.)
m.drawmeridians(meridians,linewidth=.5,dashes=[1,1])
mm, mp = m([0,30,60,90,120,150,180,210,240,270,300,330],[0,0,0,0,0,0,0,0,0,0,0,0])
#mlabels=[360,330,300,270,240,210,180,150,120,90,60,30]
#for name,xpt,ypt in zip(mlabels,mm,mp):
#	plt.text(xpt+250000,ypt+150000,name)
#pm, pp = m([180,180,180,180],[-60,-30,30,60])	
#plabels=['-60','-30','+30','+60']
#for name,xpt,ypt in zip(plabels,pm,pp):
#	plt.text(xpt+250000,ypt+100000,name)


#plt.title('Chi squared with respect to direction for boosts of $740$km.s$^{-1}$ from LG frame')
plt.clabel(CSl2,inline=1,fontsize=9,fmt='%.2f',manual=False)



cmbfig.savefig('../Figures/Figures/Figure_23_a.eps',dpi=600)

