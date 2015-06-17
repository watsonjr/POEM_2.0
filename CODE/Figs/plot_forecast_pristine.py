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
def plot_abs(X,Y,Z,Name):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    #! Interpolate to regular grid
    x = np.arange(-280,81,0.5)
    y = np.arange(-81,91,0.5)
    X2,Y2 = np.meshgrid(x,y);
    Z2 = griddata((X.flatten(),Y.flatten()), Z.flatten(), (X2,Y2), method='nearest')

    #! Plot
    fig = plt.figure(figsize=(15,6))
    im = m.pcolormesh(X2,Y2,Z2,cmap=plt.cm.jet,latlon=True,
                      vmin=Z.min(),vmax=np.percentile(Z,100))
    m.fillcontinents(color=[.6,.6,.6],lake_color=[0,0,0])

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    plt.colorbar(im, cax=cax)
    plt.title(Name)

    plt.savefig("./PNG/"+Name,dpi=300)

def plot_ratio(X,Y,Z,Name):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    #! Interpolate to regular grid
    x = np.arange(-280,81,0.5)
    y = np.arange(-81,91,0.5)
    X2,Y2 = np.meshgrid(x,y);
    Z2 = griddata((X.flatten(),Y.flatten()), Z.flatten(), (X2,Y2), method='nearest')

    #! Plot
    fig = plt.figure(figsize=(15,6))
    im = m.pcolormesh(X2,Y2,np.log10(Z2),cmap=plt.cm.seismic,latlon=True,
                      vmin=-.1,vmax=0.1)
    m.fillcontinents(color=[.6,.6,.6],lake_color=[.4,.4,.4])

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    plt.colorbar(im, cax=cax)
    plt.title(Name)

    plt.savefig("./PNG/"+Name,dpi=300)



#! Base Map
#! use this alot: http://earthpy.org/interpolation_between_grids_with_basemap.html
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

#! Data to plot
bio_pi_early = np.mean(PISC[0:12,:,:],0).T
bio_pi_late  = np.mean(PISC[-13:-1,:,:],0).T
bio_pi_ratio = bio_pi_late / bio_pi_early

bio_pl = PLAN[-1,:,:].T
bio_de = PISC[-1,:,:].T
bio_be = PISC[-1,:,:].T

#! FISH MAPS
PISC_early = np.zeros((200,360))
PISC_late  = np.zeros((200,360))
PISC_ratio = np.zeros((200,360))

PLAN_M = np.zeros((200,360))
DETR_M = np.zeros((200,360))
BENT_M = np.zeros((200,360))

PISC_early[ID] = np.sum(bio_pi,1)
PISC_early[ID] = np.sum(bio_pi,1)
PISC_ratio[ID] = np.mean(bio_pi_ratio,1)

PLAN_M[ID] = np.sum(bio_pl,1)
DETR_M[ID] = np.sum(bio_de,1)
BENT_M[ID] = bio_be[:,0]

#! Plot
plot_globe(Lon,Lat,PISC_M,"Fig_fish_forecast_pisc.png")
plot_globe(Lon,Lat,PISC_M,"Fig_fish_forecast_pisc.png")
plot_ratio(Lon,Lat,PISC_ratio,"Fig_fish_forecast_pisc.png")


plot_globe(Lon,Lat,PLAN_M,"Fig_fish_forecast_plan.png")
plot_globe(Lon,Lat,DETR_M,"Fig_fish_forecast_detr.png")
plot_globe(Lon,Lat,BENT_M,"Fig_fish_forecast_bent.png")




plt.show()
