
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from scipy.interpolate import griddata
from scipy.io import netcdf
from scipy import ndimage
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable


#! COBALT DATA
nc = netCDF4.Dataset("../Data/NC/ocean_cobalt_biomass_100.200601-210012.nlgz_100.nc")
Lon = np.asarray(nc.variables['geolon_t'])
Lat = np.asarray(nc.variables['geolat_t'])
Zoo = np.asarray(nc.variables['nlgz_100'])[1,:,:]
ID = np.where(Zoo!=Zoo.min())
JD = np.where(Zoo == Zoo.min())

#! Fish DATA
nc_sml_f = netCDF4.Dataset("../Data/NC/Data_spinup_pristine_sml_f.nc")
nc_sml_p = netCDF4.Dataset("../Data/NC/Data_spinup_pristine_sml_p.nc")
nc_sml_d = netCDF4.Dataset("../Data/NC/Data_spinup_pristine_sml_d.nc")
nc_med_f = netCDF4.Dataset("../Data/NC/Data_spinup_pristine_med_f.nc")
nc_med_p = netCDF4.Dataset("../Data/NC/Data_spinup_pristine_med_p.nc")
nc_med_d = netCDF4.Dataset("../Data/NC/Data_spinup_pristine_med_d.nc")
nc_lrg_p = netCDF4.Dataset("../Data/NC/Data_spinup_pristine_lrg_p.nc")

Sml_f = np.asarray(nc_sml_f.variables['biomass'])
Sml_p = np.asarray(nc_sml_p.variables['biomass'])
Sml_d = np.asarray(nc_sml_d.variables['biomass'])
Med_f = np.asarray(nc_med_f.variables['biomass'])
Med_p = np.asarray(nc_med_p.variables['biomass'])
Med_d = np.asarray(nc_med_d.variables['biomass'])
Lrg_p = np.asarray(nc_lrg_p.variables['biomass'])

#! FISH MAPS
Sml_f_M = np.zeros((200,360))
Sml_p_M = np.zeros((200,360))
Sml_d_M = np.zeros((200,360))
Med_f_M = np.zeros((200,360))
Med_p_M = np.zeros((200,360))
Med_d_M = np.zeros((200,360))
Lrg_p_M = np.zeros((200,360))

Sml_f_M[ID] = Sml_f[0,:]
Sml_p_M[ID] = Sml_p[0,:]
Sml_d_M[ID] = Sml_d[0,:]
Med_f_M[ID] = Med_f[0,:]
Med_p_M[ID] = Med_p[0,:]
Med_d_M[ID] = Med_d[0,:]
Lrg_p_M[ID] = Lrg_p[0,:]


#! PLOT
plt.figure()
Z = Lrg_p_M
Z[JD] = np.nan
im = plt.imshow(np.flipud(Z),cmap="jet")
plt.title('Large Piscivore')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.1)
plt.colorbar(im, cax=cax)
plt.savefig("./PNG/Fig_lrg_p.png",dpi=300)

#! PLOT
plt.figure()
Z = Med_p_M
Z[JD] = np.nan
im = plt.imshow(np.flipud(Z),cmap="jet")
plt.title('Medium Piscivore')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.1)
plt.colorbar(im, cax=cax)
plt.savefig("./PNG/Fig_med_p.png",dpi=300)

#! PLOT
plt.figure()
Z = Sml_p_M
Z[JD] = np.nan
im = plt.imshow(np.flipud(Z),cmap="jet")
plt.title('Small Piscivore')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.1)
plt.colorbar(im, cax=cax)
plt.savefig("./PNG/Fig_sml_p.png",dpi=300)

#! PLOT
plt.figure()
Z = Med_f_M
Z[JD] = np.nan
im = plt.imshow(np.flipud(Z),cmap="jet")
plt.title('Medium Planktivore')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.1)
plt.colorbar(im, cax=cax)
plt.savefig("./PNG/Fig_med_f.png",dpi=300)

#! PLOT
plt.figure()
Z = Sml_f_M
Z[JD] = np.nan
im = plt.imshow(np.flipud(Z),cmap="jet")
plt.title('Small Planktivore')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.1)
plt.colorbar(im, cax=cax)
plt.savefig("./PNG/Fig_sml_f.png",dpi=300)

#! PLOT
plt.figure()
Z = Med_d_M
Z[JD] = np.nan
im = plt.imshow(np.flipud(Z),cmap="jet")
plt.title('Medium Planktivore')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.1)
plt.colorbar(im, cax=cax)
plt.savefig("./PNG/Fig_med_d.png",dpi=300)

#! PLOT
plt.figure()
Z = Sml_d_M
Z[JD] = np.nan
im = plt.imshow(np.flipud(Z),cmap="jet")
plt.title('Small Detrivore')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.1)
plt.colorbar(im, cax=cax)
plt.savefig("./PNG/Fig_sml_d.png",dpi=300)


plt.close('all')





