from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from scipy.interpolate import griddata
from scipy.io import netcdf
from scipy import ndimage
import shapefile

#! COBALT DATA
nc = netCDF4.Dataset("../Data/NC/ocean_cobalt_biomass_100.200601-210012.nlgz_100.nc")
Lon = np.asarray(nc.variables['geolon_t'])
Lat = np.asarray(nc.variables['geolat_t'])
Zoo = np.asarray(nc.variables['nlgz_100'])[1,:,:]
ID = np.where(Zoo!=Zoo.min())

#! Fish DATA
nc_pisc = netCDF4.Dataset("../Data/NC/Data_forecast_pristine_pisc.nc")

###! EEZ DATA
sf = shapefile.Reader("../Data/EEZ/World_EEZ_v8_2014.prj")
recs    = sf.records()
cns     = []
for nshp in xrange(Nshp):
    cns.append(recs[nshp][1]);
cns = np.array(cns) # EEZ names

#! Arctic EEZs 
ID = []
ID = np.hstack((ID,np.where(cns == 'United States Exclusive Exclusive Economic Zone (Alaska)')[0][0]))
ID = np.hstack((ID,np.where(cns == 'Norwegian Exclusive Economic Zone')[0][0]))
shapes  = sf.shapes() # EEZ polygons

#! Get IDs of COBALT gridcells in EEZs



#! pull out shapes for these eezs

    for nshp in xrange(Nshp):
        ptchs   = []
        pts     = np.array(shapes[ID[nshp]].points)
        px      = pts[:,0]; py      = pts[:,1];
        x,y     = m(px,py);
        pts     = np.vstack([x,y]).T
        prt     = shapes[ID[nshp]].parts
        par     = list(prt) + [pts.shape[0]]
        for pij in xrange(len(prt)):
         ptchs.append(Polygon(pts[par[pij]:par[pij+1]]))
        ax.add_collection(PatchCollection(ptchs,facecolor=cccol_ff[nshp,:],\
                edgecolor='k', linewidths=.5))


