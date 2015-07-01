from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4

#! COBALT DATA
nc = netCDF4.Dataset("../Data/NC/ocean_cobalt_biomass_100.200601-210012.nlgz_100.nc")
Lon = np.asarray(nc.variables['geolon_t'])
Lat = np.asarray(nc.variables['geolat_t'])
Zoo = np.asarray(nc.variables['nlgz_100'])[1,:,:]
ID = np.where(Zoo!=Zoo.min())

#! Fish DATA
nc_pisc = netCDF4.Dataset("../Data/NC/Data_forecast_pristine_pisc.nc")
nc_plan = netCDF4.Dataset("../Data/NC/Data_forecast_pristine_plan.nc")
nc_detr = netCDF4.Dataset("../Data/NC/Data_forecast_pristine_detr.nc")
nc_bent = netCDF4.Dataset("../Data/NC/Data_forecast_pristine_bent.nc")

PISC = np.asarray(nc_pisc.variables['biomass'])
PLAN = np.asarray(nc_plan.variables['biomass'])
DETR = np.asarray(nc_detr.variables['biomass'])
BENT = np.asarray(nc_bent.variables['biomass'])

#! Slice timeseries for one loc
ID = 41093 # Iberian
#! ID = 28073 # Humboldt

bio_pi = PISC[:,:,ID]
bio_pl = PLAN[:,:,ID]
bio_de = DETR[:,:,ID]
bio_be = BENT[:,:,ID]

##### PLOT
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
import matplotlib
font = {'family' : 'serif',
            'weight' : 'normal',
            'size'   : 15,
            'usetex' : True}


fig = plt.figure(figsize=(7,10.5))
gs = gridspec.GridSpec(3, 2,width_ratios=[1,1])
FS = 17 # fontsize
X = np.arange(2006,2100,1./12.)

ax0 = plt.subplot(gs[0,:])
mu = np.mean(bio_pi,1)
for i in np.arange(0,bio_pi.shape[1]):
    plt.plot(X,bio_pi[:,i],'b',alpha=0.2)
plt.plot(X,mu,'r',lw=2)
#plt.xlabel(r'Time',fontsize=FS)
plt.setp( plt.gca().get_xticklabels(), visible=False)
plt.ylabel(r"Pisc ($kg$ $m^{-2}$)",fontsize=FS)
plt.gca().yaxis.labelpad = 0

ax1 = plt.subplot(gs[1,:])
mu = np.mean(bio_pl,1)
for i in np.arange(0,bio_pl.shape[1]):
    plt.plot(X,bio_pl[:,i],'g',alpha=0.2)
plt.plot(X,mu,color=[1.,0.6,0.3],lw=2)
#plt.xlabel(r'Time',fontsize=FS)
plt.setp( plt.gca().get_xticklabels(), visible=False)
plt.ylabel(r'Plan ($kg$ $m^{-2}$)',fontsize=FS)
plt.gca().yaxis.labelpad = 0

ax2 = plt.subplot(gs[2,:])
mu = np.mean(bio_de,1)
for i in np.arange(0,bio_de.shape[1]):
    plt.plot(X,bio_de[:,i],'k',alpha=0.2)
plt.plot(X,mu,color='k',lw=2)
plt.ylabel(r'Detr ($kg$ $m^{-2}$)',fontsize=FS)
plt.gca().yaxis.labelpad = 0
plt.xlabel(r'Time',fontsize=FS)

#plt.tight_layout(h_pad=1.0)
gs.tight_layout(fig)
plt.savefig('./PDF/Fig_timeseries_oneloc.pdf',dpi=200,bbox='tight')








