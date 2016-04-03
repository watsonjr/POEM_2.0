

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
TIME = ncread("./GCM/Forecast/ocean_cobalt_biomass_100.200601-210012.nlgz_100.nc",
	"average_T1"); # timeLAT = ncread("./GCM/grid_spec.nc","geolat_t"); # lat
LON  = ncread("./GCM/Forecast/grid_spec.nc", "geolon_t"); # lon
LAT  = ncread("./GCM/Forecast/grid_spec.nc", "geolat_t"); # lon
Z    = ncread("./GCM/Forecast/grid_spec.nc", "ht"); # depth
zt   = ncread("./GCM/Forecast/grid_spec.nc", "zt"); # depth levels
dx   = ncread("./GCM/Forecast/grid_spec.nc", "dxt");
dy   = ncread("./GCM/Forecast/grid_spec.nc", "dyt");
kmt  = ncread("./GCM/Forecast/grid_spec.nc", "kmt");
dxtn = ncread("./GCM/Forecast/grid_spec.nc", "dxtn");
dyte = ncread("./GCM/Forecast/grid_spec.nc", "dyte");
dat  = dx.*dy # area in m
datr = 1.0./dat+eps(Float64)

ni, nj = size(LON);
nk     = length(zt);
lmask  = zeros(ni,nj,nk);
# Water mask
for k=1:nk
	 for j=1:nj
			for i=1:ni
				 if (kmt[i,j] >= k)
						lmask[i,j,k] = 1.0;
				 else
						lmask[i,j,k] = 0.0;
				 end
			end
	 end
end

## Pressure from depth (1 atm per 10m depth)
Pr = (Z / 10) * 1013.25

#xy  = zeros(360,200); # coordinates of index
#xy[GRD_ID] = [1:length(GRD_ID)];

#! save
save("./Data_grid_cp_2D.jld", "TIME",TIME,"LAT",LAT,"LON",LON,"Z",Z,"AREA",dat,
	"Pr",Pr,"Nlon",ni,"Nlat",nj,"Ndep",nk,"dxtn",dxtn,"dyte",dyte,"datr",datr,"lmask",lmask);
