from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from netCDF4 import Dataset
from scipy.interpolate import griddata
from scipy import ndimage
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.patches as patches


#! Define plotting function
def plot_abs(X,Y,z,ID):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    #! Put Z onto a grid
    Z = np.zeros((200,360))
    Z[ID] = z

    #! Interpolate to regular grid
    x = np.arange(-280,81,0.5)
    y = np.arange(-81,91,0.5)
    X2,Y2 = np.meshgrid(x,y);
    Z2 = griddata((X.flatten(),Y.flatten()), Z.flatten(), (X2,Y2), method='nearest')

    #! Plot
    im = m.pcolormesh(X2,Y2,Z2,cmap=plt.cm.jet,latlon=True,
                      vmin=Z.min(),vmax=np.percentile(Z,95))
    m.fillcontinents(color=[.6,.6,.6],lake_color=[0,0,0])

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    plt.colorbar(im, cax=cax)


#! Base Map
m = Basemap(llcrnrlon=0,llcrnrlat=-80,urcrnrlon=360,urcrnrlat=80,projection='mill')

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

#! PLOT
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
import matplotlib as mpl
font = {'family' : 'serif','weight' : 'normal', 'size'   : 15}

##! Data to plot -- PISC
#! Average over ten years
bio_pi_early = PISC[-1,:,:].T

#! Average over large size classes
bio_pi_early = np.median(bio_pi_early[:,5:-1],1)

#! Plot
fig = plt.figure(figsize=(8.,6.5))
gs = gridspec.GridSpec(1, 2,width_ratios=[1,1])
ax0 = plt.subplot(gs[0,:])
plot_abs(Lon,Lat,bio_pi_early,ID)
fig.savefig('./PNG/Fig_forecast_pristine_abs_pisc.png',dpi=600,bbox_inches='tight')


##! Data to plot -- PLAN
#! Average over ten years
bio_pl_early = PLAN[-1,:,:].T

#! Average over large size classes
bio_pl_early = np.median(bio_pl_early[:,5:-1],1)

#! Plot
fig = plt.figure(figsize=(8.,6.5))
gs = gridspec.GridSpec(1, 2,width_ratios=[1,1])
ax0 = plt.subplot(gs[0,:])
plot_abs(Lon,Lat,bio_pl_early,ID)
fig.savefig('./PNG/Fig_forecast_pristine_abs_plan.png',dpi=600,bbox_inches='tight')


##! Data to plot -- DETR
#! Average over ten years
bio_de_early = DETR[-1,:,:].T

#! Average over large size classes
bio_de_early = np.median(bio_de_early[:,5:-1],1)

#! Plot
fig = plt.figure(figsize=(8.,6.5))
gs = gridspec.GridSpec(1, 2,width_ratios=[1,1])
ax0 = plt.subplot(gs[0,:])
plot_abs(Lon,Lat,bio_de_early,ID)

#plt.tight_layout(h_pad=0.5)
fig.savefig('./PNG/Fig_forecast_pristine_abs_detr.png',dpi=600,bbox_inches='tight')






