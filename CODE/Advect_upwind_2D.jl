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
function sub_advection(GRD,Bio_in,U,V,ni,nj)
	# time steps within one day
	ntime = 24.0
	nt = Int(ntime)
	dtime = 1.0/ntime
	# biol concentration
	Tfield = zeros(Float64,360,200,nt);
	Tfield[:,:,1] = Bio_in;
	Ttendency = zeros(Float64,360,200,nt);
  # grid size
	isd = 1
	jsd = 1
	ied = ni
	jed = nj

	# time loop
	for time=1:nt-1
		t = time
		wrk1 = zeros(Float64,ni,nj);
	  wrk1 = -horz_advect_tracer_upwind(U,V,Tfield[:,:,t],ni,nj)
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
	return Tfield[:,:,nt]
end


function horz_advect_tracer_upwind(uvel,vvel,Tracer_field,ni,nj)
	isd = 1
	jsd = 1
	ied = ni
	jed = nj

	fe = zeros(Float64,ni,nj);
	fn = zeros(Float64,ni,nj);
	upwind = zeros(Float64,ni,nj);

	# i-flux
	for j=jsd:jed
		for i=isd:ied #i=isd-1:ied
			# Reflective BC
			# #west
			# if (i > 1)
			# 	if (uvel[i,j].<0.0 + GRD["lmask"][i-1,j,1].==0.0)
			# 		uvel[i,j] *= -1.0
			# 	end
			# else
			# 	if (uvel[i,j].<0.0 + GRD["lmask"][ied,j,1].==0.0)
			# 		uvel[i,j] *= -1.0
			# 	end
			# end
			# #east
			# if (i == ied)
			# 	if (uvel[i,j].>0.0 + GRD["lmask"][isd,j,1].==0.0)
			# 		uvel[i,j] *= -1.0
			# 	end
			# else
			# 	if (uvel[i,j].>0.0 + GRD["lmask"][i+1,j,1].==0.0)
			# 		uvel[i,j] *= -1.0
			# 	end
			# end

			velocity = 0.5*uvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			if (i == ied)
				fe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer_field[i,j] + uneg.*Tracer_field[isd,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][isd,j,1]
			else
				fe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer_field[i,j] + uneg.*Tracer_field[i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i+1,j,1]
			end
		end
	end

	# j-flux
	for j=jsd:jed #j=jsd-1:jed
		for i=isd:ied
			# # Reflective BC
			# #south
			# if (j > 1)
			# 	if (vvel[i,j].<0.0 + GRD["lmask"][i,j-1,1].==0.0)
			# 		vvel[i,j] *= -1.0
			# 	end
			# #don't need else at j=1 (South Pole) because all land
			# end
			# #north
			# if (j < jed)
			# 	if (vvel[i,j].>0.0 + GRD["lmask"][i,j+1,1].==0.0)
			# 		vvel[i,j] *= -1.0
			# 	end
			# else #N Pole keeps same lat, changes lon
			# 	if (vvel[i,j].>0.0 + GRD["lmask"][ni-i+1,j,1].==0.0)
			# 		vvel[i,j] *= -1.0
			# 	end
			# end

			velocity = 0.5*vvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			if (j < jed)
				fn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer_field[i,j] + uneg.*Tracer_field[i,j+1]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i,j+1,1]
			else
				fn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer_field[i,j] + uneg.*Tracer_field[ni-i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][ni-i+1,j,1]
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
