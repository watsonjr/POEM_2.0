

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

########! Time and Grid data
#! Time = days since 2006-01-01
#! GRD  = Lon,Lat of ESM grid cell centroids
#TIME = ncread("./NC/ocean_cobalt_biomass_100.200601-210012.nlgz_100.nc","average_T1");
#LON = ncread("./grid_spec.nc","geolon_t"); # lon
#LAT = ncread("./grid_spec.nc","geolat_t"); # lon
#Z   = abs(ncread("./grid_spec.nc","ht")); # depth in meters
#dx = ncread("./grid_spec.nc","dxt") # width in meters
#dy = ncread(".//grid_spec.nc","dyt") # height in meters
#AREA = dx.*dy # area in m

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

## Pressure from depth (1 atm per 10m depth)
Pr = (Z / 10) * 1013.25

#######! Flip dims into map like matrix
#! Daily data is NOT flipped
#! Do NOT flipped so grid locations are consistent
#LON = flipdim(LON',1)
#LAT = flipdim(LAT',1)
#Z = flipdim(Z',1)
#dx = flipdim(dx',1)
#dy = flipdim(dy',1)
#dat = flipdim(dat',1)
#kmt  = flipdim(kmt',1)
#dxtn = flipdim(dxtn',1)
#dyte = flipdim(dyte',1)
#datr = flipdim(datr',1)

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


##########! ID of cardinal cells (land = Nan)
ID = collect(1:length(Z))
ID = reshape(ID,size(Z))
id = find(Z.>0)
IND = zeros(size(Z))
IND[id] = collect(1:length(id))
SUB = Array(Any,(length(id)))


#! retain only water cells
ID  = find(Z.>0);
LON = LON[ID] ;
LAT = LAT[ID];
Z   = Z[ID] ;
Pr  = Pr[ID];
DX = dx[ID] ;
DY = dy[ID];
AREA  = dat[ID];
dxtn  = dxtn[ID];
dyte  = dyte[ID];
datr  = datr[ID];
lmask = lmask[ID];

#! save
save("./Data_grid_forecast_NOTflipped.jld", "TIME",TIME,"LAT",LAT,"LON",LON,"Z",Z,"AREA",AREA,
		"Pr",Pr,"ID",ID,"N",length(ID),"dxtn",dxtn,"dyte",dyte,"datr",datr,"lmask",lmask);
