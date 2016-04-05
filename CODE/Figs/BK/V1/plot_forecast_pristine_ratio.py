from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from netCDF4 import Dataset
from scipy import ndimage
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.patches as patches


def plot_ratio(X,Y,z,ID,clim):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    #! Put Z onto a grid
    Z = np.zeros((200,360))
    Z[ID] = z

    #! INameNamonterpolate to regular grid
    x = np.arange(-280,81,0.5)
    y = np.arange(-81,91,0.5)
    X2,Y2 = np.meshgrid(x,y);
    Z2 = griddata((X.flatten(),Y.flatten()), Z.flatten(), (X2,Y2), method='nearest')

    #! Plot
    im = m.pcolormesh(X2,Y2,Z2,cmap=plt.cm.seismic,latlon=True,
                      vmin=clim[0],vmax=clim[1])
    m.fillcontinents(color=[.6,.6,.6],lake_color=[.4,.4,.4])

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
bio_early = np.median(PISC[0:12*10,:,:],0).T
bio_late  = np.median(PISC[-(12*10):-1,:,:],0).T
bio_early = np.mean(bio_early[:,5:-1],1)
bio_late = np.mean(bio_late[:,5:-1],1)
bio_pi_ratio = (bio_late - bio_early) / bio_early

bio_early = np.median(PLAN[0:12*10,:,:],0).T
bio_late  = np.median(PLAN[-(12*10):-1,:,:],0).T
bio_early = np.median(bio_early[:,5:-1],1)
bio_late = np.median(bio_late[:,5:-1],1)
bio_pl_ratio = (bio_late - bio_early) / bio_early

bio_early = np.median(DETR[0:12*10,:,:],0).T
bio_late  = np.median(DETR[-(12*10):-1,:,:],0).T
bio_early = np.median(bio_early[:,5:-1],1)
bio_late = np.median(bio_late[:,5:-1],1)
bio_de_ratio = (bio_late - bio_early) / bio_early


fig = plt.figure(figsize=(6.,10))
gs = gridspec.GridSpec(3, 2,width_ratios=[1,1])
ax0 = plt.subplot(gs[0,:])
plot_ratio(Lon,Lat,bio_pi_ratio,ID,[-0.2,0.2])
ax0 = plt.subplot(gs[1,:])
plot_ratio(Lon,Lat,bio_pl_ratio,ID,[-1,1])
ax0 = plt.subplot(gs[2,:])
plot_ratio(Lon,Lat,bio_de_ratio,ID)
plot_abs(Lon,Lat,bio_de_ratio,ID)

plt.tight_layout(h_pad=0.5)
fig.savefig('./PNG/Fig_forecast_pristine_ratio_pisc.png',dpi=600,bbox_inches='tight')



#
#bio_pl = PLAN[-1,:,:].T
#bio_de = PISC[-1,:,:].T
#bio_be = PISC[-1,:,:].T
#
##! FISH MAPS
#PISC_early = np.zeros((200,360))
#PISC_late  = np.zeros((200,360))
#PISC_ratio = np.zeros((200,360))
#
#PLAN_M = np.zeros((200,360))
#DETR_M = np.zeros((200,360))
#BENT_M = np.zeros((200,360))
#
#PISC_early[ID] = np.sum(bio_pi,1)
#PISC_early[ID] = np.sum(bio_pi,1)
#PISC_ratio[ID] = np.mean(bio_pi_ratio,1)
#
#PLAN_M[ID] = np.sum(bio_pl,1)
#DETR_M[ID] = np.sum(bio_de,1)
#BENT_M[ID] = bio_be[:,0]
#
##! Plot
#plot_globe(Lon,Lat,PISC_M,"Fig_fish_forecast_pisc.png")
#plot_globe(Lon,Lat,PISC_M,"Fig_fish_forecast_pisc.png")
#plot_ratio(Lon,Lat,PISC_ratio,"Fig_fish_forecast_pisc.png")
#
#
#plot_globe(Lon,Lat,PLAN_M,"Fig_fish_forecast_plan.png")
#plot_globe(Lon,Lat,DETR_M,"Fig_fish_forecast_detr.png")
#plot_globe(Lon,Lat,BENT_M,"Fig_fish_forecast_bent.png")
#



plt.show()
