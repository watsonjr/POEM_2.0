! PIECES OF THE MOM5 OCEAN_TRACER_ADVECT ROUTINE
! Tracer_field = tracer concentration
! Grd%dyte = cell length on east side
! Grd%dxtn = cell length on north side
! Grd%dat = cell area
! Grd%datr = 1/cell area
! Grd%tmask = land mask
! Adv_vel%uhrho_et = eastward velocity * h * rho
! Adv_vel%vhrho_nt = northward velocity * h * rho
! Adv_vel%wrho_bt = vertical velocity * rho

! Horizontal

Tracer%wrk1(isc:iec,jsc:jec,:) =  &
-horz_advect_tracer_upwind(Adv_vel, Tracer%field(:,:,:,taum1))

function horz_advect_tracer_upwind(Adv_vel, Tracer_field)
  do k=1,nk

     ! i-flux
     do j=jsc,jec
        do i=isc-1,iec
           velocity = 0.5*Adv_vel%uhrho_et(i,j,k)
           upos     = velocity + abs(velocity)
           uneg     = velocity - abs(velocity)
           fe(i,j)  = Grd%dyte(i,j)*(upos*Tracer_field(i,j,k) + uneg*Tracer_field(i+1,j,k)) &
                      *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
           flux_x(i,j,k) = fe(i,j)
        enddo
     enddo

     ! j-flux
     do j=jsc-1,jec
        do i=isc,iec
           velocity = 0.5*Adv_vel%vhrho_nt(i,j,k)
           upos     = velocity + abs(velocity)
           uneg     = velocity - abs(velocity)
           fn(i,j)  = Grd%dxtn(i,j)*(upos*Tracer_field(i,j,k) + uneg*Tracer_field(i,j+1,k)) &
                      *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
           flux_y(i,j,k) = fn(i,j)
        enddo
     enddo

    ! together
     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_upwind(i,j,k) = &
            Grd%tmask(i,j,k)*(fe(i,j)-fe(i-1,j)+fn(i,j)-fn(i,j-1))*Grd%datr(i,j)
        enddo
     enddo
  enddo
end function horz_advect_tracer_upwind

do k=1,nk
   do j=jsc,jec
      do i=isc,iec
         Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
      enddo
   enddo
enddo




! Vertical
Tracer%wrk1(isc:iec,jsc:jec,:) = &
-vert_advect_tracer_upwind(Adv_vel, Tracer%field(:,:,:,taum1))

function vert_advect_tracer_upwind(Adv_vel, Tracer_field)
  ft1  = 0.0
  do k=1,nk
    kp1 = min(k+1,nk)
    do j=jsc,jec
      do i=isc,iec
        velocity = 0.5*Adv_vel%wrho_bt(i,j,k)
        wpos     = velocity + abs(velocity)
        wneg     = velocity - abs(velocity)
        ft2(i,j) = (wneg*Tracer_field(i,j,k) + wpos*Tracer_field(i,j,kp1)) &
                   *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
        flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
        vert_advect_tracer_upwind(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        ft1(i,j) = ft2(i,j)
      enddo
    enddo
  enddo
end function vert_advect_tracer_upwind


do k=1,nk
   do j=jsc,jec
      do i=isc,iec
         Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
      enddo
   enddo
enddo
