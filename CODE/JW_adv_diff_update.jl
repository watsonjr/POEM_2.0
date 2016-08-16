###! Diffusion
function sub_diffuse(Bio_in,Ns,A,DX,DY)
	Bio_out = similar(Bio_in)
	for I = 1:NX
		# Biomasses in Neighboring cells
		bio = Bio_in[GRD["Neigh"][I]]
		bio_out = zeros(Ns)
		#!
		bm = bio[1]
		bu = bio[2]
		bd = bio[3]
		bl = bio[4]
		br = bio[5]

		# Prefactors
		Alpha = A * DT / (DX[I]^2)
		Beta  = A * DT / (DY[I]^2)

		# calculation (2D Forward in Time Centered Scheme)
		DIFF =  (Alpha*(bu+bd)) + (Beta*(bl+br)) + (bm*(1-(2*Alpha)-(2*Beta)))

		# Update biomass
		if DIFF < 0
			bio_out = 0
		else
			bio_out = DIFF
		end
		Bio_out[I] = bio_out
	end
	return Bio_out
end


### ADVECTION ###
function sub_advection(Bio_in,Nu,U,V,DX,DY,Tp,Tb,tdif,wgt)
	#Bio_in = .bio; Nu = .nu; DX = GRD["dxtn"]; DY = GRD["dyte"]; U = ENVR.U; V = ENVR.V;

	#! Calc swimming speed (m/d)
	T = (Tp.*tdif) + (Tb.*(1.0-tdif))
	Q = ((3.9*wgt.^0.13 * exp(0.149*T)) /100*60*60*24)

	#! Fix velocities - ALREADY FIXED
	#U[find(U.==minimum(U))] = 0.0
	#V[find(V.==minimum(V))] = 0.0

	#! Find direction to swim in
	#Nu_n = Nu[GRD["Neigh"]]
	#data = zeros(NX,Ns,5)
	#for I = 1:NX
	#	for J = 1:5
	#		data[I,:,J] = Nu_n[I][J]
	#	end
	#end
	#M = maximum(data, 3)
	#IND = ismax = data .== M
	#KK = find(KK,3)
	KK = zeros(Int64,NX)
	LID = zeros(NX,5)
	for I = 1:NX
		NID = GRD["Neigh"][I]
		Nu_n = Nu[NID]
		LID[I,:] = NID - NID[1]
		KK[I] = findmax([Nu_n[1],Nu_n[2],Nu_n[3],Nu_n[4],Nu_n[5]])[2]
	end
	#LID[find(LID)] = 1
	#LID[:,1] = 1


	#! Calculate effective speed
	u = zeros(NX)
	v = zeros(NX)
	I1 = find(KK .== 1)
	I2 = find(KK .== 2)
	I3 = find(KK .== 3)
	I4 = find(KK .== 4)
	I5 = find(KK .== 5)

	u[I1] = U[I1]
	v[I1] = V[I1]

	u[I2] = U[I2]
	v[I2] = V[I2] + Q[I2]

	u[I3] = U[I3]
	v[I3] = V[I3] - Q[I3]

	u[I4] = U[I4] - Q[I4]
	v[I4] = V[I4]

	u[I5] = U[I5] + Q[I5]
	v[I5] = V[I5]

	#! Land
	for J = 1:Ns
		I1 = find(v.>0. + LID[:,2].==0.)
		I2 = find(v.<0. + LID[:,3].==0.)
		I3 = find(u.<0. + LID[:,4].==0.)
		I4 = find(u.>0. + LID[:,5].==0.)

		v[I1] *= -1
		v[I2] *= -1
		u[I3] *= -1
		u[I4] *= -1
	end

	#! U advection
	Bio_minus = similar(Bio_in)
	for I = 1:NX
		Bio_minus[I] = 0.0
	end
	for I = 1:NX
		bio = Bio_in[GRD["Neigh"][I]]
		UU = u[I] * DT / DX[I]
		BB = bio[1] - ((UU/2)*(bio[3]-bio[2])) -
							(((UU^2)/2)*(bio[3]-(2*bio[1])+bio[2]))
		if BB < 0
			Bio_minus[I] = 0
		else
			Bio_minus[I] = BB
		end
	end

	#! V advection
	Bio_plus = similar(Bio_in)
	for I = 1:NX
		Bio_plus[I] = 0.0
	end
	for I = 1:NX
		bio = Bio_minus[GRD2["Neigh"][I]]
		VV = v[I] * DT / DY[I]
		BB = bio[1] - ((VV/2)*(bio[5]-bio[4])) -
							(((VV^2)/2)*(bio[5]-(2*bio[1])+bio[4]))
		if BB < 0
			Bio_plus[I] = 0
		else
			Bio_plus[I] = BB
		end
	end

	# return
	return Bio_plus
end
