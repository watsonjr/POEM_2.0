### NETCDF EXAMPLE
using NetCDF
include("toa.jl")
# Define longitudes and latitudes, day and timesteps
lat=[-89:89]
lon=[0:359]
day=1
tim=[0:23]
# Create radiation array
rad = float64([g_pot(x2,x1,day,x3) for x1=lon, x2=lat, x3=tim])
# Define some attributes of the variable (optionlal)
varatts = @Compat.Dict("longname" => "Radiation at the top of the atmosphere",
           "units"    => "W/m^2")
lonatts = @Compat.Dict("longname" => "Longitude",
           "units"    => "degrees east")
latatts = @Compat.Dict("longname" => "Latitude",
           "units"    => "degrees north")
timatts = @Compat.Dict("longname" => "Time",
           "units"    => "hours since 01-01-2000 00:00:00")
#Here we start by defining the dimensions. This is done by creating NcDim objects:
latdim = NcDim("lat",lat,latatts)
londim = NcDim("lon",lon,lonatts)
timdim = NcDim("time",tim,timatts)
# Then we create an NcVar object, the data type is defined by the corresponding julia type:
radvar = NcVar("rad",[londim,latdim,timdim],varatts,Float32)
# Now we can finally create the netcdf-file and get a file handler in return:
isfile("radiation2.nc") ? rm("radiation2.nc") : nothing
nc = NetCDF.create("radiation2.nc",radvar)
# Writing data to the file is done using putvar
NetCDF.putvar(nc, "rad", rad )
# And we close the file
NetCDF.close( nc )
###
