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

cmb_instr='../Figures/data/Figure_3_a.txt'
cmb_hdat = loadtxt(cmb_instr)
cmb_hdat=array([float(s) for s in cmb_hdat[:,2]])

##bdata='boost_data_LG_740.dat'
##bdata_hdat=loadtxt(bdata)
##bdata_hdat=array([float(s) for s in bdata_hdat[:,2]])

##bdata_hdat= transpose(reshape(bdata_hdat,(360, 180)))





## Plot CMB temperature dipole
cmb_hdatt = transpose(reshape(cmb_hdat,(360, 180)))
lons=linspace(360.,0.,num=360)
lats=linspace(90.,-90.,num=180)
clons, clats = meshgrid(lons, lats)
# create new figure
cmbfig=plt.figure(1)

#other option for projection is 'hammer'
m = Basemap(resolution='c',projection='moll',lon_0=180)
ax = cmbfig.add_axes([0.1,0.1,0.8,0.8])

# Virgo, Fornax, Shapley, GA, Perseus in galactic coords
# SC: (l,b)=(310,30) (Watkins), b=9?
#vx, vy = m([-284.,-238.,-307.,-320.,-150.4],[74.,-54.,30.,10.,-13.4]) 
#p2=m.plot(vx,vy,'o',mfc='none',mec='r')
# Direction of sun wrt LG

# The following plots points and labels on the map

cmbx, cmby = m(180+360-96.44,29.3) 
p3=m.plot(cmbx,cmby,'kx',color='white',markersize=20,mew=2)
#plt.text(cmbx+100000,cmby+400000,'CMB')   #Plots text
# Centre of local void (Iwata et al. 2011 after Tully et al 2008)
#cmbx, cmby = m(-66.,-10.) 
#p5=m.plot(cmbx,cmby,'kx')
#plt.text(cmbx+500000,cmby+100000,'LV')
# 50 /h Mpc Bulk flow direction (Watkins et al)
bf50x, bf50y = m(-287.,8.)
#p4=m.plot(bf50x,bf50y,'ro')
#tx, ty = m(-90.,-45.)
#p4=m.plot(tx,ty,'ro')
tx1, ty1 = m(360-75.48 ,-4.69)
#p4=m.plot(tx1,ty1,'+',mfc='none')
#plt.text(tx1+400000,ty1+400000,'Minimum $\chi^2$')








# Regulus (star in constellation Leo)
#bf50x, bf50y = m(-225.98,49.1)
#p4=m.plot(bf50x,bf50y,'*r')
#plt.text(bf50x+50000,bf50y+50000,'Reg')
# Cluster labels: 
#labels=['VC','Fo','SC','GA','Pers']
#for name,xpt,ypt in zip(labels,vx,vy):
#	plt.text(xpt+50000,ypt+50000,name)
# shift data so lons go from -360 to 0 instead of 0 to 360.
#hdat,lons = shiftgrid(0.,hdat,lons,start=False)
# make a filled contour plot.
# Transform the meshgrid to the new projection
x, y = m(clons, clats) 
#CSl = m.contour(x,y,hdat,clevs,linewidths=0.5,colors='k')


# The following sets up the colorbar down the bottom of the figure
# setup colorbar axes instance.
cmap=plt.cm.jet
#CSf = m.contourf(x,y,cmb_hdatt,clevs,edgecolors='none',cmap=cmap)
clevs=5
levels=[2.26,3.97]
CSl = m.contour(x,y,cmb_hdatt,levels,linewidths=1,colors='k')
clevs=15
#Good colormaps: Greys (bw), summer, PuBu, YlGn
CSf = m.contourf(x,y,cmb_hdatt,clevs,edgecolors='none',cmap=plt.cm.PuBu)





#CS4=plt.contour(x,y,cmb_hdatt,levels,colors=('k'),linewidths=(0.5,))

pos = ax.get_position()
l, b, w, h = pos.bounds
#cax = plt.axes([l+w+0.05, b, 0.025, h]) # setup colorbar axes
cax = plt.axes([l, b, w,0.025])
l_f = FormatStrFormatter('%f2.1')
cbar = plt.colorbar(CSf, drawedges=False, cax=cax, orientation='horizontal') # draw colorbar (the thing down the bottom)
cbar.ax.yaxis.set_major_formatter(l_f)
cbar.set_label(r'$\chi_a^2$',fontsize=20)
cbar.ax.tick_params(labelsize=20)






		
   
#plt.clabel(CSl,inline=1,fontsize=9,fmt='%1.2f')
plt.axes(ax)  # make the original axes current again

# Draw the graticule
m.drawmapboundary(linewidth=.5)
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
plt.clabel(CSl,inline=1,fontsize=9,fmt='%.2f',manual=False)

cmbfig.savefig('../Figures/Figures/Figure_3_a.eps',dpi=600)

