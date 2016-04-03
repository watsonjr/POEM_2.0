program advect_grid_data

use netcdf
implicit none

!------------------------- DECLARATIONS -------------------------
integer :: ncid, LonDimID, LatDimID
integer :: dxtVarID, dytVarID, kmtVarID, dxtnVarID, dyteVarID
integer :: isd, ied, jsd, jed
integer :: i, j, k
integer, parameter :: ni = 360, nj = 200, nk = 50
!real :: dxtValues, dytValues, kmtValues, dxtnValues, dyteValues
real, dimension(ni,nj) :: dxtValues  ! lon-width of T-cells at grid point (m)
real, dimension(ni,nj) :: dytValues  ! lat-width of T-cells at grid point (m)
real, dimension(ni,nj) :: datValues  ! area of T-cells (m^2)
real, dimension(ni,nj) :: dxtnValues ! long-width of north face of T-cells (m)
real, dimension(ni,nj) :: dyteValues ! lat-width of east face of T-cells (m)
real, dimension(ni,nj) :: datrValues ! 1/[area of T-cells (m^2)]
real, dimension(ni,nj) :: kmtValues
character(len=256) :: grid_file = '/Volumes/GFDL/GCM_DATA/Forecast/grid_spec.nc'
character(len=256) :: grid_version = 'VERSION_0'
character(len=256) :: x_boundary_type
character(len=256) :: y_boundary_type

type :: ocean_grid_type
  character(len=32) :: name
  ! geometry and topology and rotation
  logical                   :: cyclic_x ! true if cyclic in the i dir
  logical                   :: cyclic_y ! true if cyclic in the j dir
  logical                   :: tripolar ! folded at "top" Arctic row
  integer                   :: nk       ! # of vertical grid points
  integer                   :: ni, nj   ! # of global horizontal pts
  integer, dimension(ni,nj) :: kmt      ! # of t-levels
  real, dimension(ni,nj)    :: dxt      ! lon-width of T-cells at gridpt(m)
  real, dimension(ni,nj)    :: dyt      ! lat-width of T-cells at gridpt(m)
  real, dimension(ni,nj)    :: dat      ! area of T-cells (m^2)
  real, dimension(ni,nj)    :: dxtn     ! lon-width of N face of T-cells(m)
  real, dimension(ni,nj)    :: dyte     ! lat-width of E face of T-cells(m)
  real, dimension(ni,nj)    :: datr     ! 1/[area of T-cells (m^2)]
  ! land/sea masks
  real, dimension(ni,nj,nk) :: tmask !land/sea mask for T cells
end type ocean_grid_type

!type(ocean_grid_type), intent(inout) :: Grid
!type(ocean_grid_type), intent(out) :: Grid
type(ocean_grid_type) :: Grid


!------------------------- NETCDF -------------------------
    !open netCDF dataset
    call check(nf90_open(grid_file, nf90_nowrite, ncid=ncid))
    ! get dimension IDs
    call check(nf90_inq_dimid(ncid, "gridlon_t", LonDimID))
    call check(nf90_inq_dimid(ncid, "gridlat_t", LatDimID))
    ! get variable IDs
    call check(nf90_inq_varid(ncid, "dxt", dxtVarID))
    call check(nf90_inq_varid(ncid, "dyt", dytVarID))
    call check(nf90_inq_varid(ncid, "kmt", kmtVarID))
    call check(nf90_inq_varid(ncid, "dxtn", dxtnVarID))
    call check(nf90_inq_varid(ncid, "dyte", dyteVarID))
    ! get attribute values
    call check(nf90_get_att(ncid, nf90_global, "y_boundary_type", y_boundary_type))
    call check(nf90_get_att(ncid, nf90_global, "x_boundary_type", x_boundary_type))
    ! get values of variables
    call check(nf90_get_var(ncid, dxtVarId, dxtValues))
    call check(nf90_get_var(ncid, dytVarId, dytValues))
    call check(nf90_get_var(ncid, dxtnVarId, dxtnValues))
    call check(nf90_get_var(ncid, dyteVarId, dyteValues))
    call check(nf90_get_var(ncid, kmtVarId, kmtValues))
    ! close netCDF dataset
    call check(nf90_close(ncid))

!------------------------- GRID COMPONENTS -------------------------
    isd = 1; ied=ni; jsd=1; jed=nj;

! Horizontal
!allocate (Grid%dxt(ni,nj))
!allocate (Grid%dyt(ni,nj))
!allocate (Grid%dat(ni,nj))
!allocate (Grid%dxtn(ni,nj))
!allocate (Grid%dyte(ni,nj))
!allocate (Grid%datr(ni,nj))

    Grid%ni = ni
    Grid%nj = nj
    Grid%nk = nk

    Grid%cyclic_x=.false.; Grid%cyclic_y=.false.; Grid%tripolar=.false.

      if(x_boundary_type == 'cyclic') Grid%cyclic_x = .true.

      if(y_boundary_type == 'cyclic') then
          Grid%cyclic_y = .true.
      else if(y_boundary_type == 'fold_north_edge') then
          Grid%tripolar = .true.
      endif

    if(Grid%cyclic_x) then
      write(*,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): x_boundary_type is cyclic'
    else
      write(*,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): x_boundary_type is solid_walls'
    endif

    if(Grid%tripolar) then
      write(*,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is fold_north_edge'
    else if(Grid%cyclic_y) then
      write(*,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is cyclic'
    else
      write(*,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is solid_walls'
    endif

    if(Grid%tripolar) then
      write (*,'(1x,/a)')  ' ==> Note: Energy conversion errors are nontrivial when using tripolar=.true.'
      write (*,'(7x,a)')   'The cause is related to the need to update redundantly computed information'
      write (*,'(7x,a)')   'across the Arctic bipolar fold in a bit-wise exact manner for terms contributing'
      write (*,'(7x,a)')   'to the energy conversion analysis.  The extra code and mpp calls have not been'
      write (*,'(7x,a)')   'implemented.'
    endif

    Grid%dxt(isd:ied,jsd:jed) = dxtValues
    Grid%dyt(isd:ied,jsd:jed) = dxtValues
    Grid%dxtn(isd:ied,jsd:jed) = dxtnValues
    Grid%dyte(isd:ied,jsd:jed) = dyteValues
    Grid%dat(:,:) = Grid%dxt(:,:)*Grid%dyt(:,:)
    Grid%datr(:,:) = 1.0/(Grid%dat(:,:)+epsilon(1.0))


! Vertical
!allocate (Grid%kmt(ni,nj))
!allocate (Grid%tmask(ni,nj,nk))

  Grid%kmt(isd:ied,jsd:jed) = kmtValues

    ! construct T cell and U cell land/sea masks
    do k=1,nk
       do j=jsd,jed
          do i=isd,ied
             if (Grid%kmt(i,j) .ge. k) then
                Grid%tmask(i,j,k) = 1.0
             else
                Grid%tmask(i,j,k) = 0.0
             endif
          enddo
       enddo
    enddo

!------------------------- WRITE DATA -------------------------
! HOW DO I SAVE GRID SO THAT IT CAN BE PASSED TO ADVECT
! INSTEAD OF CALCULATING EVERY TIME?
    write(*,*) 'dat(10,10) : ', Grid%dat(10,10)
    write(*,*) 'datr(10,10) : ', Grid%datr(10,10)
    write(*,*) 'dxtn(111,111) : ', Grid%dxtn(111,111)
    write(*,*) 'dyte(15,15) : ', Grid%dyte(15,15)
    write(*,*) 'tmask(150,80,5) : ', Grid%tmask(150,80,5)


  contains
    subroutine check(status)
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
    end subroutine check

end program advect_grid_data
