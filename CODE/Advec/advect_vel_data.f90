program advect_vel_data

use netcdf
implicit none

!------------------------- DECLARATIONS -------------------------
integer :: ncidu, ncidv, LonDimID, LatDimID
integer :: uVarID, vVarID
integer :: isd, ied, jsd, jed, tsd, ted
integer :: i, j, k
integer, parameter :: ni = 360, nj = 200, nt = 1140, nk = 50
real, dimension(ni,nj,nt) :: uValues  ! lon-width of T-cells at grid point (m)
real, dimension(ni,nj,nt) :: vValues  ! lat-width of T-cells at grid point (m)
character(len=256) :: u_file = '/Volumes/GFDL/GCM_DATA/Forecast/ocean.200601-210012.u_surf.nc'
character(len=256) :: v_file = '/Volumes/GFDL/GCM_DATA/Forecast/ocean.200601-210012.v_surf.nc'

type :: ocean_adv_vel_type
  real, dimension(ni,nj,nk,nt)   :: uhrho_et  !rho_dzu * advect vel (kg/(m*s)) on i-face of T-cell
  real, dimension(ni,nj,nk,nt)   :: vhrho_nt  !rho*dzu * advect vel (kg/(m*s)) on j-face of T-cell
  real, dimension(ni,nj,0:nk,nt) :: wrho_bt   !rho * vertical advect vel (kg/(m^2*s)) on T-bottom
end type ocean_adv_vel_type

type(ocean_adv_vel_type) :: Adv_vel


!------------------------- NETCDF -------------------------
  ! U VELOCITY
    !open netCDF dataset
    call check(nf90_open(u_file, nf90_nowrite, ncid=ncidu))
    ! get dimension IDs
    call check(nf90_inq_dimid(ncidu, "xt_ocean", LonDimID))
    call check(nf90_inq_dimid(ncidu, "yt_ocean", LatDimID))
    ! get variable IDs
    call check(nf90_inq_varid(ncidu, "u_surf", uVarID))
    ! get values of variables
    call check(nf90_get_var(ncidu, uVarId, uValues))
    ! close netCDF dataset
    call check(nf90_close(ncidu))

    ! V VELOCITY
      !open netCDF dataset
      call check(nf90_open(v_file, nf90_nowrite, ncid=ncidv))
      ! get variable IDs
      call check(nf90_inq_varid(ncidv, "v_surf", vVarID))
      ! get values of variables
      call check(nf90_get_var(ncidv, vVarId, vValues))
      ! close netCDF dataset
      call check(nf90_close(ncidv))

!------------------------- VEL COMPONENTS -------------------------
    isd = 1; ied=ni; jsd=1; jed=nj; tsd = 1; ted = nt;

    ! U
    Adv_vel%uhrho_et(isd:ied,jsd:jed,1,tsd:ted) = uValues
    ! V
    Adv_vel%vhrho_nt(isd:ied,jsd:jed,1,tsd:ted) = vValues


!------------------------- WRITE DATA -------------------------
    write(*,*) 'u(150,80,1,1:2) : ', Adv_vel%uhrho_et(150,80,1,1:2)
    write(*,*) 'v(150,80,1,1:2) : ', Adv_vel%vhrho_nt(150,80,1,1:2)
    write(*,*) 'u(150,80,2,1) : ', Adv_vel%uhrho_et(150,80,2,1)
    write(*,*) 'v(150,80,2,1) : ', Adv_vel%vhrho_nt(150,80,2,1)

  contains
    subroutine check(status)
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
    end subroutine check

end program advect_vel_data
