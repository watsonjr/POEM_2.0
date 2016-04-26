

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

TIME = ncread("./GCM/Hindcast/ocean_cobalt_biomass_100.186101-200512.nlgz_100.nc",
	"average_T1"); # timeLAT = ncread("./GCM/grid_spec.nc","geolat_t"); # lat
LON  = ncread("./GCM/Hindcast/grid_spec.nc", "geolon_t"); # lon
LAT  = ncread("./GCM/Hindcast/grid_spec.nc", "geolat_t"); # lon
Z    = ncread("./GCM/Hindcast/grid_spec.nc", "ht"); # depth
zt   = ncread("./GCM/Hindcast/grid_spec.nc", "zt"); # depth levels
dx   = ncread("./GCM/Hindcast/grid_spec.nc", "dxt");
dy   = ncread("./GCM/Hindcast/grid_spec.nc", "dyt");
kmt  = ncread("./GCM/Hindcast/grid_spec.nc", "kmt");
dxtn = ncread("./GCM/Hindcast/grid_spec.nc", "dxtn");
dyte = ncread("./GCM/Hindcast/grid_spec.nc", "dyte");
dat  = dx.*dy # area in m
datr = 1.0./dat+eps(Float64)

## Pressure from depth (1 atm per 10m depth)
Pr = (Z / 10) * 1013.25

#######! Flip dims into map like matrix
LON = flipdim(LON',1)
LAT = flipdim(LAT',1)
Z = flipdim(Z',1)
dx = flipdim(dx',1)
dy = flipdim(dy',1)
dat = flipdim(dat',1)
kmt  = flipdim(kmt',1)
dxtn = flipdim(dxtn',1)
dyte = flipdim(dyte',1)
datr = flipdim(datr',1)

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


for i in collect(1:length(id))
	#! indexes
	sub = ind2sub(ID, id[i])
	up = collect(sub) ; down = collect(sub) ; left = collect(sub) ; right = collect(sub)
	up[2] = up[2] + 1
	down[2] = down[2] - 1
	right[1] = right[1] + 1
	left[1] = left[1] - 1

	#! boundaries
	if left[1]==0
		left[1] = 360
	end
	if right[1]==361
		right[1] = 1
	end
	if down[2]==201
		down[2] = 200
	end
	if up[2]==0
		up[2] = 1
	end

	#! get indexes of cardinal cells accounting for land
	id_m = IND[sub[1],sub[2]]
	if Z[up[1],up[2]] == 0.0
		id_up = IND[sub[1],sub[2]]
	else
		id_up = IND[up[1],up[2]]
	end
	if Z[down[1],down[2]] == 0.0
		id_dw = IND[sub[1],sub[2]]
	else
		id_dw = IND[down[1],down[2]]
	end
	if Z[left[1],left[2]] == 0.0
		id_lf = IND[sub[1],sub[2]]
	else
		id_lf = IND[left[1],left[2]]
	end
	if Z[right[1],right[2]] == 0.0
		id_rt = IND[sub[1],sub[2]]
	else
		id_rt = IND[right[1],right[2]]
	end

	#! store
	SUB[i] = round(Int64,[id_m,id_up,id_dw,id_lf,id_rt])
end


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
save("./Data_grid_hindcast_NOTflipped.jld", "TIME",TIME,"LAT",LAT,"LON",LON,"Z",Z,"AREA",AREA,
		"Pr",Pr,"ID",ID,"N",length(ID),"dxtn",dxtn,"dyte",dyte,"datr",datr,"lmask",lmask,
		"Neigh",SUB);
