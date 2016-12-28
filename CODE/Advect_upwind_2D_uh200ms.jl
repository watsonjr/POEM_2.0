###! DIFFUSION ###-------------------------------
#function sub_diffuse(Bio_in,Ns,A,DX,DY)
#	Bio_out = similar(Bio_in)
#
#	for j=jsd:jed
#		for i=isd:ied
#			# Prefactors
#			Alpha = A[i,j] .* DT ./ (DX[i,j]^2)
#			Beta  = A[i,j] .* DT ./ (DY[i,j]^2)
#
#			# calculation (2D Forward in Time Centered Scheme)
#			DIFF =  (Alpha*(bio(i,j-1)+bio(i,j+1))) + (Beta*(bio(i-1,j)+bio(i+1,j))) + (bm*(1-(2*Alpha)-(2*Beta)))
#
#			# Update biomass
#			if DIFF < 0
#				bio_out[i,j] = 0
#			else
#				bio_out[i,j] = DIFF
#			end
#		end
#	end
#
#	Bio_out[I] = bio_out
#
#	return Bio_out
#end


### ADVECTION ###-------------------------------
function sub_advection(GRD,Bio_in,U,V,ni,nj,dep)
	# ntime = time steps in a day
	# dtime = # seconds in ntime
	dtime = 60.0*60.0 #60.0*60.0 = 24 hours
	ntime = (60.0*60.0*24.0) / dtime
	nt = Int(ntime)
	# biol concentration
	Tfield = zeros(Float64,360,200,nt);
	#Tfield = zeros(Float64,360,200,2);
	#Tfield[:,:,1] = Bio_in/ntime;
	Tfield[:,:,1] = Bio_in;
	Ttendency = zeros(Float64,360,200,nt);
  # grid size
	isd = 1
	jsd = 2 #ignore j=1 b/c land (Antarctica)
	ied = ni
	jed = nj

	# time loop
	for time=1:nt-1
		t = time
		#println(t)
		wrk1 = zeros(Float64,ni,nj);
	  wrk1 = -horz_advect_tracer_upwind(U,V,Tfield[:,:,t],ni,nj,dep)
		for j=jsd:jed
			for i=1:ied
				Ttendency[i,j,t] = Ttendency[i,j,t] + wrk1[i,j]
			end
		end
		for j=jsd:jed
	    for i=1:ied
				Tfield[i,j,time+1]=(Tfield[i,j,time] + dtime.*Ttendency[i,j,time])
	    end
	  end
	end

	# return
	#return Tfield[:,:,nt] * ntime
	return Tfield[:,:,nt]
	#return Tfield[:,:,2]
end


function horz_advect_tracer_upwind(uvel,vvel,Tracer_field,ni,nj,dep)
	isd = 1
	jsd = 2 #ignore j=1 b/c land (Antarctica)
	ied = ni
	jed = nj

	fe = zeros(Float64,ni,nj);
	fn = zeros(Float64,ni,nj);
	upwind = zeros(Float64,ni,nj);

	# i-flux
	for j=jsd:jed
		for i=isd:ied #i=isd-1:ied

			velocity = 0.5*uvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			if (i == ied)
				if (dep[i,j] > 0.0 && dep[isd,j] > 0.0)
					fe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer_field[i,j]./dep[i,j] + uneg.*Tracer_field[isd,j]./dep[isd,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][isd,j,1]
				else
					fe[i,j]  = 0.0
				end
			else
				if (dep[i,j] > 0.0 && dep[i+1,j] > 0.0)
					fe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer_field[i,j]./dep[i,j] + uneg.*Tracer_field[i+1,j]./dep[i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i+1,j,1]
				else
					fe[i,j]  = 0.0
				end
			end
		end
	end

	# j-flux
	for j=jsd:jed #j=jsd-1:jed
		for i=isd:ied

			velocity = 0.5*vvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			if (j < jed)
				if (dep[i,j] > 0.0 && dep[i,j+1] > 0.0)
					fn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer_field[i,j]./dep[i,j] + uneg.*Tracer_field[i,j+1]./dep[i,j+1]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i,j+1,1]
				else
					fn[i,j]  = 0.0
				end
			else
				if (dep[i,j] > 0.0 && dep[ni-i+1,j] > 0.0)
					fn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer_field[i,j]./dep[i,j] + uneg.*Tracer_field[ni-i+1,j]./dep[ni-i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][ni-i+1,j,1]
				else
					fn[i,j]  = 0.0
				end
			end
		end
	end

	# combined
	for j=jsd:jed
		for i=isd:ied
			if (j > 1)
				if (i > 1)
					upwind[i,j] = GRD["lmask"][i,j,1].*(fe[i,j]-fe[i-1,j]+fn[i,j]-fn[i,j-1]).*GRD["datr"][i,j]
				else
					upwind[i,j] = GRD["lmask"][i,j,1].*(fe[i,j]-fe[ied,j]+fn[i,j]-fn[i,j-1]).*GRD["datr"][i,j]
				end
			# no need for else at South Pole
			end
		end
	end

	# return
	return upwind
end
