###! Diffusion
function sub_diffuse(Bio_in,Ns,A,DX,DY)
	Bio_out = similar(Bio_in)
	for I = 1:NX
		# Biomasses in Neighboring cells
		bio = Bio_in[GRD["Neigh"][I]]
		bio_out = zeros(Ns)

		for J = 1:Ns

			#!
			bm = bio[1][J]
			bu = bio[2][J]
			bd = bio[3][J]
			bl = bio[4][J]
			br = bio[5][J]

			# Prefactors
			Alpha = A[J] * DT / (DX[I]^2)
			Beta  = A[J] * DT / (DY[I]^2)

			# calculation (2D Forward in Time Centered Scheme)
			DIFF =  (Alpha*(bu+bd)) + (Beta*(bl+br)) + (bm*(1-(2*Alpha)-(2*Beta)))

			# Update biomass
			if DIFF < 0
				bio_out[J] = 0
			else
				bio_out[J] = DIFF
			end
		end
		Bio_out[I] = bio_out
	end
	return Bio_out
end


### ADVECTION ###
function sub_advection(Bio_in,Nu,Ns,Q,U,V,DX,DY)
	#Bio_in = PISC.bio; Nu = PISC.nu; DX = GRD["DX"]; DY = GRD["DY"]; U = ENVR.U; V = ENVR.V; Ns = PI_N; Q = PI_U;

	#! Fix velocities
	U[find(U.==minimum(U))] = 0.0
	V[find(V.==minimum(V))] = 0.0

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
	KK = zeros(Int64,NX,Ns)
	Nu_n = Nu[GRD["Neigh"]]
	LID = zeros(NX,5)
	for I = 1:NX
		NID = GRD["Neigh"][I]
		LID[I,:] = NID - NID[1]
		nu_n = Nu_n[I]
		for J = 1:Ns
			KK[I,J] = findmax([nu_n[1][J],nu_n[2][J],nu_n[3][J],nu_n[4][J],nu_n[5][J]])[2]
		end
	end
	#LID[find(LID)] = 1
	#LID[:,1] = 1


	#! Calculate effective speed
	u = zeros(NX,Ns)
	v = zeros(NX,Ns)
	for J = 1:Ns
		I1 = find(KK[:,J] .== 1)
		I2 = find(KK[:,J] .== 2)
		I3 = find(KK[:,J] .== 3)
		I4 = find(KK[:,J] .== 4)
		I5 = find(KK[:,J] .== 5)

		u[I1,J] = U[I1]
		v[I1,J] = V[I1]

		u[I2,J] = U[I2]
		v[I2,J] = V[I2] + Q[J]

		u[I3,J] = U[I3]
		v[I3,J] = V[I3] - Q[J]

		u[I4,J] = U[I4] - Q[J]
		v[I4,J] = V[I4]

		u[I5,J] = U[I5] + Q[J]
		v[I5,J] = V[I5]
	end

	#! Land
	for J = 1:Ns
		I1 = find(v[:,J].>0. + LID[:,2].==0.)
		I2 = find(v[:,J].<0. + LID[:,3].==0.)
		I3 = find(u[:,J].<0. + LID[:,4].==0.)
		I4 = find(u[:,J].>0. + LID[:,5].==0.)

		v[I1,J] *= -1
		v[I2,J] *= -1
		u[I3,J] *= -1
		u[I4,J] *= -1
	end

	#! U advection
	Bio_minus = similar(Bio_in)
	for I = 1:NX
		Bio_minus[I] = zeros(Ns)
	end
	for I = 1:NX
		bio = Bio_in[GRD["Neigh"][I]]
		for J = 1:Ns
			UU = u[I,J] * DT / DX[I]
			BB = bio[1][J] - ((UU/2)*(bio[3][J]-bio[2][J])) -
								(((UU^2)/2)*(bio[3][J]-(2*bio[1][J])+bio[2][J]))
			if BB < 0
				Bio_minus[I][J] = 0
			else
				Bio_minus[I][J] = BB
			end
		end
	end

	#! V advection
	Bio_plus = similar(Bio_in)
	for I = 1:NX
		Bio_plus[I] = zeros(Ns)
	end
	for I = 1:NX
		bio = Bio_minus[GRD["Neigh"][I]]
		for J = 1:Ns
			VV = v[I,J] * DT / DY[I]
			BB = bio[1][J] - ((VV/2)*(bio[5][J]-bio[4][J])) -
								(((VV^2)/2)*(bio[5][J]-(2*bio[1][J])+bio[4][J]))
			if BB < 0
				Bio_plus[I][J] = 0
			else
				Bio_plus[I][J] = BB
			end

		end
	end

	# return
	return Bio_plus
end
