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

cmb_instr='../Figures/data/Figure_15_b.txt'
#sunlg_instr='sunLGDipole.dat'
cmb_hdat = loadtxt(cmb_instr)
cmb_hdat=array([float(s) for s in cmb_hdat[:,2]])
#cmb_hmax=cmb_hdat.max()
#cmb_hmin=cmb_hdat.min()	
#cmb_hdiff=cmb_hmax-cmb_hmin
#print "CMB: ",cmb_hmax,cmb_hmin,cmb_hdiff
#sunlg_hdat = loadtxt(sunlg_instr)
#sunlg_hdat=array([float(s) for s in sunlg_hdat[:,2]])
#sunlg_hmax=sunlg_hdat.max()
#sunlg_hmin=sunlg_hdat.min()	
#sunlg_hdiff=sunlg_hmax-sunlg_hmin
#print "Sun-LG: ",sunlg_hmax,sunlg_hmin,sunlg_hdiff

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

cmbx, cmby = m(360-89,-28) 
p3=m.plot(cmbx,cmby,'+',color='white',markersize=20,mew=2)
#plt.text(cmbx+800000,cmby+900000,'LG')   #Plots text
# Centre of local void (Iwata et al. 2011 after Tully et al 2008)
#cmbx, cmby = m(-66.,-10.) 
#p5=m.plot(cmbx,cmby,'kx')
#plt.text(cmbx+500000,cmby+100000,'LV')
# 50 /h Mpc Bulk flow direction (Watkins et al)
bf50x, bf50y = m(-287.,8.)
#p4=m.plot(bf50x,bf50y,'ro')
#tx, ty = m(-90.,-45.)
#p4=m.plot(tx,ty,'ro')
tx1, ty1 = m(360-81.526 ,-19.04) 
#coordinates for minimum chi^2 in 634 from CMB
#p4=m.plot(tx1,ty1,'+',mfc='none',markersize=10)
#plt.text(tx1+400000,ty1+400000,'Minimum $\chi^2$')

tx1, ty1 = m(360-60 ,12.6574)
tx1, ty1 = m(360-81.526 ,-19.04) #coordinates for minimum chi^2 in 634 from CMB
p7=m.plot(tx1,ty1,'kx',mfc='black',markersize=20,mew=2)




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
levels=[0.2,1,1.6,1.8]
CSl = m.contour(x,y,cmb_hdatt,levels,linewidths=0.5,colors='k')
clevs=150

CSf = m.contourf(x,y,cmb_hdatt,clevs,edgecolors='none',cmap=plt.cm.summer)





#CS4=plt.contour(x,y,cmb_hdatt,levels,colors=('k'),linewidths=(0.5,))

pos = ax.get_position()
l, b, w, h = pos.bounds
#cax = plt.axes([l+w+0.05, b, 0.025, h]) # setup colorbar axes
cax = plt.axes([l, b, w,0.025])
l_f = FormatStrFormatter('%f2.1')
cbar = plt.colorbar(CSf,drawedges=False,ticks=[0.1, 1,1.9], cax=cax, orientation='horizontal') # draw colorbar (the thing down the bottom)
cbar.ax.yaxis.set_major_formatter(l_f)
plt.clim(0,2)

cbar.ax.set_xticklabels(['0.1', '1', '1.9'])

#cbar.set_label(r'$f_p$ $=$ $|p+1|$ when $A\geq 0$ and $2-|p+1|$ when $A<0$',fontsize=20)
#cbar.ax.tick_params(labelsize=20)
#plt.clim(0,2)
#cbar.set_label(r'$f_p$ $=$ $|p+1|$ when $A\geq 0$ and $2-|p+1|$ when $A<0$',fontsize=20)
#cbar.ax.tick_params([0, 1, 2],labelsize=20)




plt.axes(ax)  # make the original axes current againplt.clim(0,2)
cbar.set_label(r'$f_p$ $=$ $|p+1|$ when $A\geq 0$ and $2-|p+1|$ when $A<0$',fontsize=20)
cbar.ax.tick_params(labelsize=20)



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


#plt.title('Best fit parameters to a systematic boost offset of $200$km.s$^{-1}$ from LG frame')
#plt.clabel(CSl,inline=1,fontsize=9,fmt='%.2f',manual=True)

cmbfig.savefig('../Figures/Figures/Figure_15_b.eps',dpi=600)

