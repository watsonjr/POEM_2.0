

###### OFFLINE DATA REQUIRED BY POEM 2.0
#! Time and Lon and Lat

#! 1) data to be interpolated to daily resolution,
#! and saved in monthly chunks
#! 2) daily average values will be saved too as a separate data file


###### LOAD COBALT DATA
#! ncinfo to look into netCDF files, ncread to load
#! e.g. ncinfo("./GCM/ocean.200601-210012.temp_100_avg.nc");
using NetCDF, HDF5, JLD, Grid, MAT

#! Time
TIME = ncread("./GCM/Hindcast/ocean_cobalt_biomass_100.186101-200512.nmdz_100.nc",
    "average_T1"); # time
TIME = TIME[1:240]

#! Ocean currents
U = ncread("./GCM/Hindcast/feb152013_run25_ocean.198801-200712_u200_v200.nc","Ut_200");
V = ncread("./GCM/Hindcast/feb152013_run25_ocean.198801-200712_u200_v200.nc","Vt_200");

###### INTERPOLATE DATA TO SIZE-BASED MODEL TIME SCALES
#! Save in annual chunks (365 days)
#! load grid data (for pressure to calc bottom temp)

#! index of water cells
#Use Zm to be consistent with ID
Zm=ncread("./GCM/Hindcast/ocean_cobalt_biomass_100.186101-200512.nmdz_100.nc",
	"nmdz_100");
WID = find(Zm[:,:,1] .!= -1.0e10); # spatial index of water cells
NID = length(WID); # number of water cells

#! months in a year
lstd = Int(TIME[length(TIME)])+31
id1 = collect(0:365:(lstd-1))
id2 = collect(365:365:(lstd))
ID  = [id1 id2];
nyr = Int(lstd/365)
#! pull out annual information
#! *60 *60 *24 --> per day (if flux)
for i = 1:nyr
	#id = float64(ID[i,:]) #float64(x::AbstractArray) is deprecated, use map(Float64,x) instead
  id = map(Float64,ID[i,:])
	I = find(id[1].<=TIME.<=id[2])

	#! pull raw data out
	time = TIME[I];
	u = U[:,:,I];
  v = V[:,:,I];

	#! setup POEM data files
	D_u = zeros(NID,365);
  D_v = zeros(NID,365);

  #! NaN velocities = -999
  u[find(u.==minimum(u))] = 0.0
  v[find(v.==minimum(v))] = 0.0

	#! interpolate to daily resolution
	for j = 1:NID
		#! indexes
		m,n = ind2sub((360,200),WID[j]); # spatial index of water cell

    #! v currents
		Y = zeros(size(time))
		Y[:] = v[m,n,:];
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi= collect(time[1]:1:time[end]);
		Yi = yi[Xi[1:end-1]];
		D_v[j,1:length(Yi)] = (Yi *60 *60 *24); #from (m/s) to (m d-1)

		#! u currents
		Y = zeros(size(time))
		Y[:] = u[m,n,:];
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi= collect(time[1]:1:time[end]);
		Yi = yi[Xi[1:end-1]];
		D_u[j,1:length(Yi)] = Yi *60 *60 *24; #from (m/s) to (m d-1)

	end

	#! save
	println(i)
	ti = string(1000000+i); di = "./JLD/Data_hindcast_vel200_";
	save(string(di,ti[2:end],".jld"),"U",D_u,"V",D_v);


end
