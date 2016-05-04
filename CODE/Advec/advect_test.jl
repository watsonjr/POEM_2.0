# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT
include("Parameters.jl")
include("sub_init.jl")
include("sub_routines.jl")
include("sub_functions.jl")
include("Experiments.jl")

make_parameters(0)

GRD = load("./Data/Data_grid_cp_2D.jld")
COBALT = load("./Data/JLD/Data_hindcast2D_000002.jld"); # 1980

### ADVECTION ###
function sub_advection(GRD,Bio_in,U,V)
	Tfield = Bio_in;
	# dtime
	dtime = 1.0/12.0

	isd = 1
	jsd = 1
	ni, nj = size(GRD["LON"]);
	ied = ni
	jed = nj

	# time loop
	for time=1:nt
	  horz_advect_tracer(time, Adv_vel, Tfield, Ttendency)
	  for j=jsd:jed
	    for i=isd:ied
	      Tfield(i,j,time+1) =  (Tfield(i,j,time) + dtime*Ttendency(i,j,time))
	    end
	  end
	end

	# return
	return Tfield
end

function horz_advect_tracer(time, U, V, field, tendency)
	isd = 1
	jsd = 1
	ied = ni
	jed = nj
	t = time

	for j=jsd:jed
		for i=isd:ied
			wrk1(i,j) = 0.0
		end
	end

  wrk1(isd:ied,jsd:jed) = -horz_advect_tracer_upwind(U(:,:,t), V(:,:,t), field(:,:,t))

	for j=jsd:jed
		for i=isd:ied
				tendency(i,j,t) = tendency(i,j,t) + wrk1(i,j)
		end
	end

	# return
	return tendency

end #function horz_advect_tracer

function horz_advect_tracer_upwind(uvel, vvel, Tracer_field)
	isd = 1
	jsd = 1
	ied = ni
	jed = nj

	# i-flux
	for j=jsd:jed
		for i=isd-1:ied
			velocity = 0.5*uvel(i,j)
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			fe(i,j)  = GRD["dyte"](i,j)*(upos*Tracer_field(i,j) + uneg*Tracer_field(i+1,j)) &
			*GRD["lmask"](i,j)*GRD["lmask"](i+1,j)
		end
	end
	# j-flux
	for j=jsd-1:jed
		for i=isd:ied
			velocity = 0.5*vvel(i,j)
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			fn(i,j)  = GRD["dxtn"](i,j)*(upos*Tracer_field(i,j) + uneg*Tracer_field(i,j+1)) &
			*GRD["lmask"](i,j)*GRD["lmask"](i,j+1)
		end
	end
	# combined
	for j=jsd:jed
		for i=isd:ied
			horz_advect_tracer_upwind(i,j) = &
			GRD["lmask"](i,j)*(fe(i,j)-fe(i-1,j)+fn(i,j)-fn(i,j-1))*GRD["datr"](i,j)
		end
	end
