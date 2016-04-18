program surf_advect_upwind_io
  !
  !<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Matt Harrison
  !</CONTACT>
  !
  !<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> S. M. Griffies
  !</CONTACT>
  !
  !<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> John Dunne
  !</CONTACT>
  !
  !<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Alistair Adcroft
  !</CONTACT>
  !
  !<OVERVIEW>
  ! This module computes thickness weighted tracer advection tendencies.
  !</OVERVIEW>
  !
  !<DESCRIPTION>
  ! This module computes tracer advection using an upwind advection schemes.
  ! Grid and velocity data are read from netcdfs
  ! Tracer concentration data are random
  ! Only surface advection for 12 months is tested
  !</DESCRIPTION>
  !
  !<INFO>
  !
  !<REFERENCE>
  ! S.-J. Lin
  ! A "Vertically Lagrangian" Finite-Volume Dynamical Core for Global Models
  ! Month. Weather Rev. (2004) 132, 2293-2307
  ! (Appendix B)
  !</REFERENCE>
  !
  !<REFERENCE>
  ! H.T. Huynh
  ! Schemes and Constraints for advection
  ! 15th Intern. Conf. on Numeric. Meth. in Fluid Mech., Springer (1997)
  !</REFERENCE>
  !
  !<REFERENCE>
  ! Colella P. and P.R. Woodward
  ! The piedewise parabloic method (PPM) for gasdynamical simulations
  ! J. Comput. Phys. (1984) 54, 174-201
  !</REFERENCE>
  !
  !<REFERENCE>
  ! A. Suresh and H.T. Huynh
  ! Accurate Monotonicity-Preserving Schemes with Runge-Kutta Time Splitting
  ! J. Comput. Phys. (1997) 136, 83-99
  !</REFERENCE>
  !
  !<REFERENCE>
  ! V. Daru and  C. Tenaud
  ! High order one-step monotonicity-preserving schemes for unsteady
  ! compressible flow calculations.
  ! J. Comp. Phys. (2004) 193, 563-594
  !</REFERENCE>
  !
  !<REFERENCE>
  ! R.C. Easter
  ! Two modified versions of Botts positive-definite numerical advection scheme.
  ! Month. Weath. Rev. (1993) 121, 297-304
  !</REFERENCE>
  !
  !<REFERENCE>
  ! Prather, M. J.,"Numerical Advection by Conservation of Second-Order Moments"
  ! JGR, Vol 91, NO. D6, p 6671-6681, May 20, 1986
  !</REFERENCE>
  !
  !<REFERENCE>
  ! Merryfield and Holloway (2003), "Application of an accurate advection
  ! algorithm to sea-ice modelling". Ocean Modelling, Vol 5, p 1-15.
  !</REFERENCE>
  !
  !<REFERENCE>
  ! Hundsdorder and Trompert (1994), "Method of lines and
  ! direct disdretization: a comparison for linear
  ! advection", Applied Numerical Mathematics,
  ! pages 469--490.
  !</REFERENCE>
  !
  !<REFERENCE>
  ! Sweby (1984): "High-resolution schemes using flux
  ! limiters for hyperbolic conservation laws",
  ! SIAM Journal of Numerical Analysis, vol. 21
  ! pages 995-1011.
  !</REFERENCE>
  !
  !</INFO>
  !
  !<NOTE>
  ! Contains a version of quicker with MOM3 masking.
  !</NOTE>
  !
  !<NAMELIST NAME="ocean_tracer_advect_nml">
  !  <DATA NAME="limit_with_upwind" TYPE="logical">
  !  If true, will compute tracer fluxes entering a cell using upwind
  !  if the tracer value is outside a specified range. Implemented
  !  only for quick at this time. This is an ad hoc and incomplete attempt
  !  to maintain monotonicity with the quicker scheme.
  !  </DATA>
  !  <DATA NAME="advect_sweby_all" TYPE="logical">
  !  For running all tracers with sweby, thereby utilizing a bitwise same
  !  routine that reorganizes loops and can be faster for certain configurations.
  !  Default advect_sweby_all=.false.
  !  </DATA>
  !  <DATA NAME="zero_tracer_advect_horz" TYPE="logical">
  !  For debugging.  Set to .true. to turn off horizontal advection.
  !  </DATA>
  !  <DATA NAME="zero_tracer_advect_vert" TYPE="logical">
  !  For debugging.  Set to .true. to turn off vertical advection.
  !  </DATA>
  !  <DATA NAME="psom_limit_prather" TYPE="logical">
  !  For running with the original Prather limiter for the PSOM scheme.
  !  The limiter is positive definite, but not monotonic.  This limiter
  !  is NOT recommended for most applications.  The default is
  !  psom_limit_prather=.false., since we prefer to use the limiter
  !  from Merryfield and Holloway (2003).
  !  </DATA>
  !
  !  <DATA NAME="debug_this_module" TYPE="logical">
  !  For debugging
  !  </DATA>
  !
  !  <DATA NAME="write_a_restart" TYPE="logical">
  !  Set true to write a restart.  False setting only for rare
  !  cases where wish to benchmark model without measuring the cost
  !  of writing restarts and associated chksums.
  !  Default is write_a_restart=.true.
  !  </DATA>
  !
  !  <DATA NAME="read_basin_mask" TYPE="logical">
  !  For reading in a mask that selects regions of the domain
  !  for performing gyre and overturning diagnostics.
  !  The basin-mask convention used at GFDL has
  !  Southern=1.0,Atlantic=2.0,Pacific=3.0,Arctic=4.0,Indian=5.0
  !  Default read_basin_mask=.false., whereby basin_mask
  !  is set to tmask(k=1).
  !  </DATA>
  !
  !</NAMELIST>

use netcdf
implicit none

! FROM OCEAN_PARAMETERS_MOD
! parameters for tracer advection
integer, parameter :: ADVECT_UPWIND            = 1
integer, parameter :: ADVECT_QUICKER           = 5

! FROM OCEAN_TYPES_MOD
integer, parameter :: ni = 360, nj = 200, nk = 50, nt = 1, dt = 1440
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
  real, dimension(ni,nj)    :: dat      ! area of T-cells (m^1)
  real, dimension(ni,nj)    :: dxtn     ! lon-width of N face of T-cells(m)
  real, dimension(ni,nj)    :: dyte     ! lat-width of E face of T-cells(m)
  real, dimension(ni,nj)    :: datr     ! 1/[area of T-cells (m^2)]
  ! land/sea masks
  real, dimension(ni,nj)    :: tmask    !land/sea mask for T cells
end type ocean_grid_type

type :: ocean_adv_vel_type
  real, dimension(ni,nj,nt)   :: uhrho_et  !rho_dzu * advect vel (kg/(m*s)) on i-face of T-cell
  real, dimension(ni,nj,nt)   :: vhrho_nt  !rho*dzu * advect vel (kg/(m*s)) on j-face of T-cell
  real, dimension(ni,nj,nt)   :: wrho_bt   !rho * vertical advect vel (kg/(m^2*s)) on T-bottom
end type ocean_adv_vel_type

type :: ocean_prog_tracer_type
  character(len=32) :: name
  integer :: horz_advect_scheme=1  ! id for horizontal advection scheme
  integer :: vert_advect_scheme=1  ! id for vertical advection scheme
  real, dimension(ni,nj,nt*dt)   :: field        !tracer concentration at 1 time
  real, dimension(ni,nj)         :: wrk1         !work array
  real, dimension(ni,nj,nt*dt)   :: tendency     !advection tendency
end type ocean_prog_tracer_type

! FROM OCEAN_TRACER_ADVECT_MOD
character(len=256) :: version='CVS $Id: ocean_tracer_advect.F90,v 1.1.2.4.18.1 2013/03/06 18:19:34 Niki.Zadeh Exp $'
character(len=256) :: tagname='Tag $Name: siena_201303 $'
type(ocean_grid_type)  , pointer :: Grd =>NULL()
logical :: limit_with_upwind       = .false.
logical :: advect_sweby_all        = .false.
logical :: zero_tracer_advect_horz = .false.
logical :: zero_tracer_advect_vert = .false.
logical :: write_a_restart         = .false.

type(ocean_grid_type),    target, save :: Grid
type(ocean_adv_vel_type), target, save :: Adv_vel
real,                     target, save :: tmin, tmax
type(ocean_prog_tracer_type)  :: T_prog !If we pass all fish types at once, find how to make array
type(ocean_prog_tracer_type)  :: Tracer
real                          :: dtime
logical                       :: store_flux
integer                       :: i, j, time, d, count
integer                       :: isd, ied, jsd, jed, tsd, ted
real,dimension(ni,nj,nt)      :: tempk
real,dimension(ni,nj)         :: temp

character(len=256) :: grid_filename = '/Volumes/GFDL/GCM_DATA/Forecast/grid_spec.nc'
character(len=256) :: u_fn = '/Volumes/GFDL/GCM_DATA/Forecast/ocean.200601-210012.u_surf.nc'
character(len=256) :: v_fn = '/Volumes/GFDL/GCM_DATA/Forecast/ocean.200601-210012.v_surf.nc'

! Grid
call grid_data(Grid, grid_filename)

! Adv_vel
call adv_vel_data(Adv_vel, u_fn, v_fn)

! T_prog
T_prog%name = "PL"
T_prog%horz_advect_scheme = 1
T_prog%vert_advect_scheme = 1

! Tracer
call random_number(tempk)
call random_number(temp)
Tracer%field(:,:,1) = 0.0
Tracer%field(:,94:99,1) = 1.0e9
Tracer%field(:,:,1) = Tracer%field(:,:,1) * Grid%tmask
Tracer%wrk1 = floor(temp)
Tracer%tendency(:,:,1) = floor(temp)

! dtime
dtime = 1.0/1440.0
count = 0

! store_flux
store_flux = .false.

isd = 1
jsd = 1
ied = Grid%ni
jed = Grid%nj

!time=1

call ocean_tracer_advect_init (Grid, T_prog)

! time loop
do time=1,nt
  do d=1,dt
    count = count + 1
    !write(*,*) 'count : ', count
    call horz_advect_tracer(time, count, Adv_vel, Tracer, store_flux)
    do j=jsd,jed
      do i=isd,ied
        Tracer%field(i,j,count+1) =  (Tracer%field(i,j,count) &
        + dtime*Tracer%tendency(i,j,count))
      enddo
    enddo
  enddo
enddo

call get_tracer_stats(Tracer%field, Grid%tmask, tmin, tmax)
!write(*,*) 'Tracer min : ', tmin
!write(*,*) 'Tracer max : ', tmax
! EQUAL -INF AND INF AFTER 2 TIME STEPS, MEANS THE MODEL BLOWS UP
! TRACER CONCENTRATION SHOULD NOT BE ABLE TO GO NEG
! HOW FIX?

contains

!#######################################################################
! <SUBROUTINE NAME="check">
!
!   <OVERVIEW>
!     Checks status of netcdf i/o calls
!   </OVERVIEW>
!
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="grid_data">
!
!   <OVERVIEW>
!     Gets gridspec variables needed for advection routine
!   </OVERVIEW>
!   <DESCRIPTION>
!     Gets gridspec variables needed for advection routine
!   </DESCRIPTION>
!   <TEMPLATE>
!    call grid_data ( filename )
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character(len=*)">
!    string that contains the gridspec.nc filename
!   </IN>
!
    subroutine grid_data (Grid, filename)
      integer :: ncid, LonDimID, LatDimID
      integer :: dxtVarID, dytVarID, kmtVarID, dxtnVarID, dyteVarID
      integer :: isd, ied, jsd, jed
      integer :: i, j, k
      integer, parameter :: ni = 360, nj = 200, nk = 50
      !real :: dxtValues, dytValues, kmtValues, dxtnValues, dyteValues
      real, dimension(ni,nj) :: dxtValues  ! lon-width of T-cells at grid point (m)
      real, dimension(ni,nj) :: dytValues  ! lat-width of T-cells at grid point (m)
      real, dimension(ni,nj) :: dxtnValues ! long-width of north face of T-cells (m)
      real, dimension(ni,nj) :: dyteValues ! lat-width of east face of T-cells (m)
      real, dimension(ni,nj) :: kmtValues
      character(len=256) :: grid_file
      character(len=256) :: filename
      character(len=256) :: x_boundary_type
      character(len=256) :: y_boundary_type
      type(ocean_grid_type), intent(inout) :: Grid
      !type(ocean_grid_type) :: Grid

      grid_file = filename
      write(*,'(1x,a)') &
      trim(grid_file) //' used for grid specifications.'

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
            write(*,*) '==>Note from ocean_grid: x_boundary_type is cyclic'
          else
            write(*,*) '==>Note from ocean_grid: x_boundary_type is solid_walls'
          endif

          if(Grid%tripolar) then
            write(*,*) '==>Note from ocean_grid: y_boundary_type is fold_north_edge'
          else if(Grid%cyclic_y) then
            write(*,*) '==>Note from ocean_grid: y_boundary_type is cyclic'
          else
            write(*,*) '==>Note from ocean_grid: y_boundary_type is solid_walls'
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
        Grid%kmt(isd:ied,jsd:jed) = nint(kmtValues)
        ! construct T cell and U cell land/sea masks
        k=1
       do j=jsd,jed
          do i=isd,ied
             if (Grid%kmt(i,j) .ge. k) then
                Grid%tmask(i,j) = 1.0
             else
                Grid%tmask(i,j) = 0.0
             endif
          enddo
       enddo

    end subroutine grid_data
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="adv_vel_data">

!   <OVERVIEW>
!     Gets gridspec variables needed for advection routine
!   </OVERVIEW>
!   <DESCRIPTION>
!     Gets gridspec variables needed for advection routine
!   </DESCRIPTION>
!   <TEMPLATE>
!    call grid_data ( filename )
!   </TEMPLATE>
!   <IN NAME="filename" TYPE="character(len=*)">
!    string that contains the gridspec.nc filename
!   </IN>
!
  subroutine adv_vel_data (Adv_vel, u_filename, v_filename)
    integer :: ncidu, ncidv, uVarID, vVarID
    integer :: isd, ied, jsd, jed, tsd, ted
    integer, parameter :: ni = 360, nj = 200, nt=1
    real, dimension(ni,nj,nt) :: uValues  ! lon-width of T-cells at grid point (m)
    real, dimension(ni,nj,nt) :: vValues  ! lat-width of T-cells at grid point (m)
    character(len=256) :: u_file
    character(len=256) :: v_file
    character(len=256) :: u_filename
    character(len=256) :: v_filename

    type(ocean_adv_vel_type), intent(inout) :: Adv_vel

    u_file = u_filename
    v_file = v_filename
    !------------------------- NETCDF -------------------------
    ! U VELOCITY
    !open netCDF dataset
    call check(nf90_open(u_file, nf90_nowrite, ncid=ncidu))
    ! get variable IDs
    call check(nf90_inq_varid(ncidu, "u_surf", uVarID))
    ! get values of variables
    call check(nf90_get_var(ncidu, uVarId, uValues, start = (/ 1, 1, 1 /), &
    count = (/ ni, nj, nt /)))
    ! close netCDF dataset
    call check(nf90_close(ncidu))

    ! V VELOCITY
    !open netCDF dataset
    call check(nf90_open(v_file, nf90_nowrite, ncid=ncidv))
    ! get variable IDs
    call check(nf90_inq_varid(ncidv, "v_surf", vVarID))
    ! get values of variables
    call check(nf90_get_var(ncidv, vVarId, vValues, start = (/ 1, 1, 1 /), &
    count = (/ ni, nj, nt /)))
    ! close netCDF dataset
    call check(nf90_close(ncidv))

    !------------------------- VEL COMPONENTS -------------------------
    isd = 1; ied=ni; jsd=1; jed=nj; tsd = 1; ted = nt;

    ! Horizontal
    Adv_vel%uhrho_et(isd:ied,jsd:jed,tsd:ted) = uValues
    Adv_vel%vhrho_nt(isd:ied,jsd:jed,tsd:ted) = vValues
    ! Veritcal
    Adv_vel%wrho_bt(isd:ied,jsd:jed,tsd:ted) = 0.0

end subroutine adv_vel_data
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="write_version_number">

!   <OVERVIEW>
!     Prints to the screen (or a specified unit) the (cvs) version id string and
!     (cvs) tag name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Prints to the screen or a specified unit the (cvs) version id string
!      and (cvs) tag name.
!   </DESCRIPTION>
!   <TEMPLATE>
!    call write_version_number ( version [, tag, unit] )
!   </TEMPLATE>
!   <IN NAME="version" TYPE="character(len=*)">
!    string that contains routine name and version number.
!   </IN>
!   <IN NAME="tag" TYPE="character(len=*)">
!    The tag/name string, this is usually the Name string
!    returned by CVS when checking out the code.
!   </IN>
!   <IN NAME="unit" TYPE="integer">
!    The Fortran unit number of an open formatted file.
!   </IN>
! prints module version number to the log file of specified unit number

  subroutine write_version_number (version, tag)

    !   in:  version = string that contains routine name and version number
    !
    !   optional in:
    !        tag = cvs tag name that code was checked out with

    character(len=*), intent(in) :: version
    character(len=*), intent(in), optional :: tag

    if (present(tag)) then
      write (*,'(/,80("="),/(a))') trim(version), trim(tag)
    else
      write (*,'(/,80("="),/(a))') trim(version)
    endif

  end subroutine write_version_number
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_advect_init">
!
! <DESCRIPTION>
! Initialize the tracer advection module.
! </DESCRIPTION>
!
  subroutine ocean_tracer_advect_init (Grid, T_prog)

    type(ocean_grid_type),        intent(in), target   :: Grid
    type(ocean_prog_tracer_type), intent(inout)        :: T_prog

    call write_version_number( version, tagname )


    Grd => Grid

    ! Use this if we give all fish types at once
    ! If just pass one-by-one, then not necessary
    !num_prog_tracers = size(T_prog(:))

    if(.not. write_a_restart) then
      write(*,'(a)') '==>Note: running ocean_tracer_advect with write_a_restart=.false.'
      write(*,'(a)') '   Will NOT write restart file, so cannot restart if using PSOM advection.'
    endif

    if(zero_tracer_advect_horz) then
      write(*,*)'==>WARNING: have turned off horizontal tracer advection. Unrealistic simulation.'
    endif
    if(zero_tracer_advect_vert) then
      write(*,*)'==>WARNING: have turned off vertical tracer advection. Unrealistic simulation.'
    endif

    if(.not. advect_sweby_all) then
      write(*,'(a)') ' '
      write(*,'(a)') ' From ocean_tracer_advect_init: SUMMARY OF TRACER ADVECTION SCHEMES'

      ! Use this if we give all fish types at once
      !do n=1,num_prog_tracers
      !  if(T_prog(n)%horz_advect_scheme==ADVECT_UPWIND) then
      !    write(*,'(1x,a)') &
      !    trim(T_prog(n)%name) //'is using first order upwind for horz advection.'
      !  endif
      !  if(T_prog(n)%vert_advect_scheme==ADVECT_UPWIND) then
      !    write(*,'(1x,a)') &
      !    trim(T_prog(n)%name) //'is using first order upwind for vert advection.'
      !  endif
      !  if(T_prog(n)%horz_advect_scheme==ADVECT_QUICKER) then
      !    write(*,'(1x,a)') &
      !    trim(T_prog(n)%name) //' is using Quicker for horz advection.'
      !    if(limit_with_upwind) then
      !      write(*,'(1x,a)')'limit_with_upwind reverts quicker to upwind if tracer outside limits'
      !    endif
      !  endif
      !  if(T_prog(n)%vert_advect_scheme==ADVECT_QUICKER) then
      !    write(*,'(1x,a)') &
      !    trim(T_prog(n)%name) //' is using Quicker for vert advection.'
      !    if(limit_with_upwind) then
      !      write(*,'(1x,a)')'limit_with_upwind reverts quicker to upwind if tracer outside limits'
      !    endif
      !  endif
      !enddo

      if(T_prog%horz_advect_scheme==ADVECT_UPWIND) then
        write(*,'(1x,a)') &
        trim(T_prog%name) //' is using first order upwind for horz advection.'
      endif
      if(T_prog%vert_advect_scheme==ADVECT_UPWIND) then
        write(*,'(1x,a)') &
        trim(T_prog%name) //' is using first order upwind for vert advection.'
      endif
      if(T_prog%horz_advect_scheme==ADVECT_QUICKER) then
        write(*,'(1x,a)') &
        trim(T_prog%name) //' is using Quicker for horz advection.'
        if(limit_with_upwind) then
          write(*,'(1x,a)')'limit_with_upwind reverts quicker to upwind if tracer outside limits'
        endif
      endif
      if(T_prog%vert_advect_scheme==ADVECT_QUICKER) then
        write(*,'(1x,a)') &
        trim(T_prog%name) //' is using Quicker for vert advection.'
        if(limit_with_upwind) then
          write(*,'(1x,a)')'limit_with_upwind reverts quicker to upwind if tracer outside limits'
        endif
      endif
      write(*,'(a)') ' '

    endif !(.not. advect_sweby_all)

  end subroutine ocean_tracer_advect_init
! </SUBROUTINE> NAME="ocean_tracer_advect_init"

!#######################################################################
! <SUBROUTINE NAME="horz_advect_tracer">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers
! </DESCRIPTION>
!
  subroutine horz_advect_tracer(time, count, Adv_vel, Tracer, store_flux)

    integer, parameter :: ni = 360, nj = 200, nt = 1
    type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
    type(ocean_prog_tracer_type), intent(inout) :: Tracer
    integer,                      intent(in)    :: time  !day/month of model run
    integer,                      intent(in)    :: count
    logical,            optional, intent(in)    :: store_flux
    integer                                     :: i, j, t, c
    logical                                     :: store
    integer                                     :: isd, ied, jsd, jed

    isd = 1
    jsd = 1
    ied = Grd%ni
    jed = Grd%nj
    t = time
    c = count

    if(zero_tracer_advect_horz) return

    store = .FALSE.
    if(present(store_flux)) store = store_flux

    if(.not. advect_sweby_all) then

      do j=jsd,jed
        do i=isd,ied
          Tracer%wrk1(i,j) = 0.0
        enddo
      enddo

      select case (Tracer%horz_advect_scheme)
      case (ADVECT_UPWIND)
        Tracer%wrk1(isd:ied,jsd:jed) =  &
        -horz_advect_tracer_upwind(Adv_vel%uhrho_et(:,:,t), Adv_vel%vhrho_nt(:,:,t), &
        Tracer%field(:,:,c))
      case default
        write(*,*)'==>Error from ocean_tracer_advect_mod (horz_advect_tracer): chose invalid horz advection scheme'
      end select

      do j=jsd,jed
        do i=isd,ied
            Tracer%tendency(i,j,c) = Tracer%tendency(i,j,c) + Tracer%wrk1(i,j)
        enddo
      enddo

    endif   ! endif for (.not. advect_sweby_all)

  end subroutine horz_advect_tracer
! </SUBROUTINE> NAME="horz_advect_tracer"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_upwind">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from first order upwind.
! This scheme is positive definite but very diffusive.
! </DESCRIPTION>
!
  function horz_advect_tracer_upwind(uvel, vvel, Tracer_field)

    integer, parameter :: ni = 360, nj = 200
    real, dimension(ni,nj), intent(in)    :: uvel
    real, dimension(ni,nj), intent(in)    :: vvel
    real, dimension(ni,nj), intent(in)    :: Tracer_field
    real, dimension(ni,nj)                :: horz_advect_tracer_upwind
    real, dimension(ni,nj)                :: fe, fn
    real                                  :: velocity, upos, uneg
    integer                               :: i, j
    integer                               :: isd, ied, jsd, jed

    isd = 1
    jsd = 1
    ied = Grd%ni
    jed = Grd%nj

    ! i-flux
    do j=jsd,jed
      do i=isd-1,ied
        velocity = 0.5*uvel(i,j)
        upos     = velocity + abs(velocity)
        uneg     = velocity - abs(velocity)
        fe(i,j)  = Grd%dyte(i,j)*(upos*Tracer_field(i,j) + uneg*Tracer_field(i+1,j)) &
        *Grd%tmask(i,j)*Grd%tmask(i+1,j)
      enddo
    enddo
    ! j-flux
    do j=jsd-1,jed
      do i=isd,ied
        velocity = 0.5*vvel(i,j)
        upos     = velocity + abs(velocity)
        uneg     = velocity - abs(velocity)
        fn(i,j)  = Grd%dxtn(i,j)*(upos*Tracer_field(i,j) + uneg*Tracer_field(i,j+1)) &
        *Grd%tmask(i,j)*Grd%tmask(i,j+1)
      enddo
    enddo
    ! combined
    do j=jsd,jed
      do i=isd,ied
        horz_advect_tracer_upwind(i,j) = &
        Grd%tmask(i,j)*(fe(i,j)-fe(i-1,j)+fn(i,j)-fn(i,j-1))*Grd%datr(i,j)
      enddo
    enddo

  end function horz_advect_tracer_upwind
! </FUNCTION> NAME="horz_advect_tracer_upwind"


!#######################################################################
! <SUBROUTINE NAME="get_tracer_stats">
!
! <DESCRIPTION>
! Compute the upper/lower values of a 3D field, returning values via arguments
! </DESCRIPTION>
  subroutine get_tracer_stats(A,Mask,tmin,tmax)
    integer, parameter :: ni = 360, nj = 200, nt = 1, dt=1440
    real, dimension(ni,nj,nt*dt), intent(in) :: A
    real, dimension(ni,nj),       intent(in) :: Mask
    real, intent(out)                        :: tmin,tmax
    integer                                  :: i,j,k
    integer                                  :: isd, ied, jsd, jed

    tmin=1.e30; tmax=-1.e30
    isd = 1
    jsd = 1
    ied = Grd%ni
    jed = Grd%nj

    do k=1,nt*dt
      do j=jsd,jed
        do i=isd,ied
          if (Mask(i,j)>0.) then
            tmin=min(tmin,A(i,j,k))
            tmax=max(tmax,A(i,j,k))
          endif
        enddo
      enddo
      write(*,*) 'Count : ', k
      write(*,*) 'Tracer min : ', tmin
      write(*,*) 'Tracer max : ', tmax
      if (tmax > 1.0e30) stop
    enddo
  end subroutine get_tracer_stats
! </SUBROUTINE> NAME="get_tracer_stats"



end program surf_advect_upwind_io
