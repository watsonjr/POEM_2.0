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
import shapefile



###! EEZ DATA
def plot_eez():
    sf = shapefile.Reader("../Data/EEZ/World_EEZ_v8_2014.prj")
    recs    = sf.records()
    cns     = []
    for nshp in xrange(len(recs)):
        cns.append(recs[nshp][1]);
    cns = np.array(cns) # EEZ names

    #! Arctic EEZs
    ID = []
    ID = np.hstack((ID,np.where(cns == 'Swedish Exclusive Economic Zone')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'Icelandic Exclusive Economic Zone')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'United States Exclusive Exclusive' 
                                        + ' Economic Zone (Alaska)')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'Russian Exclusive Economic Zone')[0][0]))
    #ID = np.hstack((ID,np.where(cns == 'Russia-Japan conflict zone')[0][0]))
    #ID = np.hstack((ID,np.where(cns == 'Japanese Exclusive Economic Zone')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'Canadian Exclusive Economic Zone')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'United Kingdom Exclusive Economic Zone')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'Irish Exclusive Economic Zone')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'Faeroe Islands Exclusive Economic Zone')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'Norwegian Exclusive Economic Zone')[0][0]))
    ID = np.hstack((ID,np.where(cns == 'Greenlandic Exclusive Economic Zone')[0][0]))

    #! pull out shapes for these eezs
    shapes  = sf.shapes() # EEZ polygons
    Nshp    = len(ID)
    ax = plt.gca()

    for nshp in xrange(Nshp):
        ptchs   = []
        pts     = np.array(shapes[int(ID[nshp])].points)
        px      = pts[:,0]; py      = pts[:,1];
        x,y     = m(px,py);
        pts     = np.vstack([x,y]).T
        prt     = shapes[int(ID[nshp])].parts
        par     = list(prt) + [pts.shape[0]]
        for pij in xrange(len(prt)):
            ptchs.append(Polygon(pts[par[pij]:par[pij+1]]))
            ax.add_collection(PatchCollection(ptchs,facecolor='none',\
                    edgecolor='k', linewidths=1.))


#! Define plotting function
def plot_arctic_zoo(X,Y,Z,Name):
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
    im = m.pcolormesh(lon_curv,lat_curv,Z_smth,cmap=plt.cm.jet,latlon=True,vmin=0,vmax=0.013)
    plot_eez() # note, has to come before fillcontinents
    m.fillcontinents(color=[.6,.6,.6],lake_color=[0,0,0])


    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    plt.colorbar(im, cax=cax)

    plt.savefig("./PDF/"+Name,dpi=300)


def plot_arctic_diff(X,Y,Z,Name):
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
    im = m.pcolormesh(lon_curv,lat_curv,Z_smth,cmap=plt.cm.seismic,latlon=True,vmin=0.,vmax=2.)
    plot_eez() # note, has to come before fillcontinents
    m.fillcontinents(color=[.6,.6,.6],lake_color=[0,0,0])

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    plt.colorbar(im, cax=cax)

    plt.savefig("./PDF/"+Name,dpi=300)



#! Base Map
#! use this alot: http://earthpy.org/interpolation_between_grids_with_basemap.html
m = Basemap(projection='npstere',boundinglat=55,lon_0=0,resolution='l')

#! high res grid
fc = netcdf.netcdf_file('../Data/curv_grid.nc')
lon_curv = fc.variables['xc'][0,:,:]
lat_curv = fc.variables['yc'][0,:,:]

#! COBALT DATA
nc = netCDF4.Dataset("../Data/ocean_cobalt_biomass_100.200601-210012.nlgz_100.nc")
Lon = np.asarray(nc.variables['geolon_t'])
Lat = np.asarray(nc.variables['geolat_t'])
Zoo = np.asarray(nc.variables['nlgz_100'])

#! Isolate early zoo mean
Zoo_fut = Zoo[1140-120:1140,:,:];
Zoo_fut = np.mean(Zoo_fut,0)
Zoo_now = Zoo[0:120,:,:];
Zoo_now = np.mean(Zoo_now,0)
Zoo_diff = Zoo_fut / Zoo_now

#! Plot
plot_arctic_zoo(Lon,Lat,Zoo_fut,"Fig_zoo_future.png")
plot_arctic_zoo(Lon,Lat,Zoo_now,"Fig_zoo_now.png")
plot_arctic_diff(Lon,Lat,Zoo_diff,"Fig_zoo_diff.png")
