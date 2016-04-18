from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from scipy.interpolate import griddata
from scipy.io import netcdf
from scipy import ndimage



#! Define plotting function
def plot_arctic_temp(X,Y,Z,Name):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    #! high res grid
    fc = netcdf.netcdf_file('../Data/curv_grid.nc')
    lon_curv = fc.variables['xc'][0,:,:]
    lat_curv = fc.variables['yc'][0,:,:]

    #! Interpolate NaNs out
    I = np.where(Z == Z.min())
    J = np.where(Z != Z.min())
    Zn = griddata((X[J],Y[J]),Z[J],(X[I],Y[I]), method='nearest')
    Z[I] = Zn;

    #! Interpolate to regular grid
    x = np.arange(-280,81,1)
    y = np.arange(-81,91,1)
    X2,Y2 = np.meshgrid(x,y);
    Z2 = griddata((X.flatten(),Y.flatten()), Z.flatten(), (X2,Y2), method='nearest')

    #! Basemap wants Lons to be in the -180:180 format
    X2[np.where(X2<-180)] = X2[np.where(X2<-180)]+360.

    #! The order of both lat and lon should be increasing
    x1 = X2[:,0:100]
    x2 = X2[:,100:361]
    X2 = np.hstack((x2,x1))

    z1 = Z2[:,0:100]
    z2 = Z2[:,100:361]
    Z2 = np.hstack((z2,z1))

    #! Isolate to the Arctic
    lat_min = 50;
    I = np.where(Y2 >= lat_min)
    X2 = X2[I]
    Y2 = Y2[I]
    Z2 = Z2[I]

    #! interpolate to higher resolution for plotting
    Z_curv = griddata((X2,Y2),Z2, (lon_curv, lat_curv), method='nearest')
    Z_smth = ndimage.filters.gaussian_filter(Z_curv, 3, mode='nearest')

    #! Plot
    fig = plt.figure(figsize=(10,10))
    im = m.pcolormesh(lon_curv,lat_curv,Z_smth,cmap=plt.cm.jet,latlon=True,vmin=-2,vmax=8)
    m.fillcontinents(color=[.6,.6,.6],lake_color=[0,0,0])

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    plt.colorbar(im, cax=cax)

    plt.savefig("./PNG/"+Name)


def plot_arctic_temp_diff(X,Y,Z,Name):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    #! high res grid
    fc = netcdf.netcdf_file('../Data/curv_grid.nc')
    lon_curv = fc.variables['xc'][0,:,:]
    lat_curv = fc.variables['yc'][0,:,:]

    #! Interpolate NaNs out
    I = np.where(Z == Z.min())
    J = np.where(Z != Z.min())
    Zn = griddata((X[J],Y[J]),Z[J],(X[I],Y[I]), method='nearest')
    Z[I] = Zn;

    #! Interpolate to regular grid
    x = np.arange(-280,81,1)
    y = np.arange(-81,91,1)
    X2,Y2 = np.meshgrid(x,y);
    Z2 = griddata((X.flatten(),Y.flatten()), Z.flatten(), (X2,Y2), method='nearest')

    #! Basemap wants Lons to be in the -180:180 format
    X2[np.where(X2<-180)] = X2[np.where(X2<-180)]+360.

    #! The order of both lat and lon should be increasing
    x1 = X2[:,0:100]
    x2 = X2[:,100:361]
    X2 = np.hstack((x2,x1))

    z1 = Z2[:,0:100]
    z2 = Z2[:,100:361]
    Z2 = np.hstack((z2,z1))

    #! Isolate to the Arctic
    lat_min = 50;
    I = np.where(Y2 >= lat_min)
    X2 = X2[I]
    Y2 = Y2[I]
    Z2 = Z2[I]

    #! interpolate to higher resolution for plotting
    Z_curv = griddata((X2,Y2),Z2, (lon_curv, lat_curv), method='nearest')
    Z_smth = ndimage.filters.gaussian_filter(Z_curv, 3, mode='nearest')

    #! Plot
    fig = plt.figure(figsize=(10,10))
    im = m.pcolormesh(lon_curv,lat_curv,Z_smth,cmap=plt.cm.seismic,latlon=True,vmin=-8.,vmax=8.)
    m.fillcontinents(color=[.6,.6,.6],lake_color=[0,0,0])

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    plt.colorbar(im, cax=cax)

    plt.savefig("./PNG/"+Name)



#! Base Map
#! use this alot: http://earthpy.org/interpolation_between_grids_with_basemap.html
m = Basemap(projection='npstere',boundinglat=55,lon_0=0,resolution='l')

#! high res grid
fc = netcdf.netcdf_file('../Data/curv_grid.nc')
lon_curv = fc.variables['xc'][0,:,:]
lat_curv = fc.variables['yc'][0,:,:]

#! COBALT DATA
nc = netCDF4.Dataset("../Data/GCM/ocean.200601-210012.temp_100_avg.nc")
Lon = np.asarray(nc.variables['GEOLON_T'])
Lat = np.asarray(nc.variables['GEOLAT_T'])
Temp = np.asarray(nc.variables['TEMP_100']) - 273

#! Isolate early zoo mean
Temp_fut = Temp[1140-120:1140,:,:];
Temp_fut = np.mean(Temp_fut,0)
Temp_now = Temp[0:120,:,:];
Temp_now = np.mean(Temp_now,0)
Temp_diff = Temp_fut - Temp_now

#! Plot
plot_arctic_temp(Lon,Lat,Temp_fut,"Fig_temp_future.png")
plot_arctic_temp(Lon,Lat,Temp_now,"Fig_temp_now.png")
plot_arctic_temp_diff(Lon,Lat,Temp_diff,"Fig_temp_diff.png")

##! Interpolate NaNs out
#I = np.where(Zoo == Zoo.min())
#J = np.where(Zoo != Zoo.min())
#Zn = griddata((Lon[J],Lat[J]),Zoo[J],(Lon[I],Lat[I]), method='nearest')
#Zoo[I] = Zn;
#
##! Interpolate to regular grid
#x = np.arange(-280,81,1)
#y = np.arange(-81,91,1) 
#X,Y = meshgrid(x,y);
#Z = griddata((Lon.flatten(),Lat.flatten()), Zoo.flatten(), (X,Y), method='nearest')
#
##! Basemap wants Lons to be in the -180:180 format
#X[np.where(X<-180)] = X[np.where(X<-180)]+360.
#
##! The order of both lat and lon should be increasing
#x1 = X[:,0:100]
#x2 = X[:,100:361]
#X = np.hstack((x2,x1))
#
#z1 = Z[:,0:100]
#z2 = Z[:,100:361]
#Z = np.hstack((z2,z1))
#
##! Isolate to the Arctic
#lat_min = 50;
#I = np.where(Y >= lat_min)
#X = X[I]
#Y = Y[I]
#Z = Z[I]
#
##! interpolate to higher resolution for plotting
#Z_curv = griddata((X,Y),Z, (lon_curv, lat_curv), method='nearest')
#Z_smth = ndimage.filters.gaussian_filter(Z_curv, 3, mode='nearest')
#
##! Plot
#fig = plt.figure(figsize=(10,10))
#m.pcolormesh(lon_curv,lat_curv,Z_smth,cmap=plt.cm.jet,latlon=True,vmin=0,vmax=0.013)
#plt.clim([0,0.013])
#m.fillcontinents(color=[.6,.6,.6],lake_color=[0,0,0])

