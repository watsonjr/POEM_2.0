from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from scipy.interpolate import griddata
from scipy.io import netcdf
from scipy import ndimage
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.patches as patches


#! Define plotting function
def plot_globe(X,Y,Z,Name):
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
PISC = np.load("../Data/NPZ/Spinup_pristine_PISC.npy")
PLAN = np.load("../Data/NPZ/Spinup_pristine_PLAN.npy")
DETR = np.load("../Data/NPZ/Spinup_pristine_DETR.npy")
BENT = np.load("../Data/NPZ/Spinup_pristine_BENT.npy")

#! FISH MAPS
PISC_M = np.zeros((200,360))
PLAN_M = np.zeros((200,360))
DETR_M = np.zeros((200,360))
BENT_M = np.zeros((200,360))

PISC_M[ID] = np.sum(PISC,1)
PLAN_M[ID] = np.sum(PLAN,1)
DETR_M[ID] = np.sum(DETR,1)
BENT_M[ID] = BENT[:,0]

#! Plot
plot_globe(Lon,Lat,PISC_M,"Fig_fish_pisc.png")
plot_globe(Lon,Lat,PLAN_M,"Fig_fish_plan.png")
plot_globe(Lon,Lat,DETR_M,"Fig_fish_detr.png")
plot_globe(Lon,Lat,BENT_M,"Fig_fish_bent.png")
plt.show()
