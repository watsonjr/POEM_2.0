

###### OFFLINE DATA REQUIRED BY POEM 2.0
#! Zm: medium zooplankton biomass (g m-2)
#! Zl: large zooplankton biomass (g m-2)
#! dZm: medium zooplankton mortality rate (g m-2 day-1)
#! dZl: large zooplankton mortality rate (g m-2 day-1)
#! dDet: detrital flux to the benthos (g m-2 day-1)
#! Tp: pelagic temperature averaged over the top 200m (deg C) 
#! Tb: bottom temperature (deg C)
#! Time and Lon and Lat

#! 1) data to be interpolated to daily resolution,
#! and saved in monthly chunks
#! 2) daily average values will be saved too as a separate data file


###### LOAD COBALT DATA
#! ncinfo to look into netCDF files, ncread to load
#! e.g. ncinfo("./GCM/ocean.200601-210012.temp_100_avg.nc");
using NetCDF, HDF5, JLD

##! Time and Grid data
#! Time = days since 2006-01-01
#! GRD  = Lon,Lat of ESM grid cell centroids
TIME = ncread("./GCM/ocean_cobalt_biomass_100.200601-210012.nmdz_100.nc",
	"average_T1"); # timeLAT = ncread("./GCM/grid_spec.nc","geolat_t"); # lat
LON = ncread("./GCM/grid_spec.nc","geolon_t"); # lon
LAT = ncread("./GCM/grid_spec.nc","geolat_t"); # lon
Z   = ncread("./GCM/grid_spec.nc","ht"); # depth

## Pressure from depth (1 atm per 10m depth)
Pr = (Z / 10) * 1013.25

# index of water cells
ID  = find(Z.>0);
LON = LON[ID];
LAT = LAT[ID];
Z   = Z[ID];
Pr  = Pr[ID];

#xy  = zeros(360,200); # coordinates of index
#xy[GRD_ID] = [1:length(GRD_ID)];

#! save
save("./JLD/Data_grid.jld", "TIME",TIME,"LAT",LAT,"LON",LON,"Z",Z,
	"Pr",Pr,"ID",ID,"N",length(ID));















