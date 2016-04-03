module ocean_tracer_advect_mod
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
  ! This module computes thickness weighted tracer advection tendencies
  ! using a variety of advection schemes.
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
  ! The piecewise parabloic method (PPM) for gasdynamical simulations
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
  ! direct discretization: a comparison for linear
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




  ! ========================================================================================================================
  ! FROM FMS_MOD
  !use FMS_mod,    only: write_version_number, error_mesg
  !use FMS_mod,    only: close_file, read_data
  !-----------------------------------------------------------------------
  !
  !         A collection of commonly used routines.
  !
  !  The routines are primarily I/O related, however, there also
  !  exists several simple miscellaneous utility routines.
  !
  !-----------------------------------------------------------------------
  !
  !  read_data           Reads distributed data from a single threaded file.
  !
  !  write_version_number  Prints to the log file (or a specified unit)
  !                        the (cvs) version id string and (cvs) tag name.
  !
  !  error_mesg          Print notes, warnings and error messages,
  !                      terminates program for error messages.
  !                      (use error levels NOTE,WARNING,FATAL)
  !
  !  close_file          Closes a file that was opened using
  !                      open_namelist_file, open_restart_file, or
  !                      open_ieee32_file.
  !
  !-----------------------------------------------------------------------
  ! routines for opening/closing specific types of file
  private :: open_namelist_file, close_file

  ! miscellaneous i/o routines
  private :: write_version_number, error_mesg

  logical :: module_is_initialized = .FALSE.

  integer, parameter :: NOTE=0, WARNING=1, FATAL=2


contains

  ! <SUBROUTINE NAME="read_data">
  !<DESCRIPTION>
  ! This routine performs reading "fieldname" stored in "filename". The data values of fieldname
  ! will be stored in "data" at the end of this routine. For fieldname with multiple timelevel
  ! just repeat the routine with explicit timelevel in each call.
  !</DESCRIPTION>
  !   <TEMPLATE>
  ! call read_data(filename,fieldname,data,domain,timelevel)
  !   </TEMPLATE>
  !   <IN NAME="filename" TYPE="character" DIM="(*)">
  !    File name
  !   </IN>
  !   <IN NAME="fieldname" TYPE="character" DIM="(*)">
  !    Field  name
  !   </IN>
  !   <IN NAME="domain"  TYPE="domain, optional">
  !   domain of fieldname
  !   </IN>
  !   <IN NAME="timelevel" TYPE="integer, optional">
  !     time level of fieldname
  !   </IN>
  !   <OUT NAME="data"  TYPE="real">
  !   array containing data of fieldname
  !   </OUT>
  !=====================================================================================
  subroutine read_data(filename,fieldname,data)
    character(len=*),                  intent(in) :: filename, fieldname
    real, dimension(:,:,:),         intent(inout) :: data ! 3 dimensional data

    character(len=256)            :: fname
    integer                       :: unit, siz_in(4)
    integer                       :: file_index  ! index of the opened file in array files
    integer                       :: tlev=1
    integer                       :: index_field ! position of the fieldname in the list of variables
    integer                       :: cxsize, cysize
    integer                       :: dxsize, dysize
    integer                       :: gxsize, gysize
    integer                       :: ishift, jshift
    logical                       :: found_file

    found_file = get_file_name(filename, fname)
    if(.not.found_file) call error_mesg(FATAL, 'fms_io_mod(read_data_3d_new): file ' //trim(filename)// &
    '(with the consideration of tile number) and corresponding distributed file are not found')
    call get_file_unit(fname, unit, file_index)

    siz_in(3) = size(data,3)
    gxsize = size(data,1)
    gysize = size(data,2)

    tlev = 1

    d_ptr =>NULL()

    return
  end subroutine read_data
  ! </SUBROUTINE>


  !#######################################################################
  ! <SUBROUTINE NAME="error_mesg">

  !   <OVERVIEW>
  !     Print notes, warnings and error messages; terminates program for warning
  !     and error messages. (use error levels NOTE,WARNING,FATAL, see example below)
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Print notes, warnings and error messages; and terminates the program for
  !     error messages. This routine is a wrapper around error_mesg, and is provided
  !     for backward compatibility. This module also publishes error_mesg,
  !      <B>users should try to use the error_mesg interface</B>.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call error_mesg ( routine, message, level )
  !   </TEMPLATE>

  !   <IN NAME="routine"  TYPE="character" >
  !     Routine name where the warning or error has occurred.
  !   </IN>
  !   <IN NAME="message"  TYPE="character" >
  !     Warning or error message to be printed.
  !   </IN>
  !   <IN NAME="level"  TYPE="integer" >
  !     Level of severity; set to NOTE, WARNING, or FATAL Termination always occurs
  !     for FATAL, never for NOTE, and is settable for WARNING (see namelist).
  !   </IN>
  !   <NOTE>
  !
  !     Examples:
  !     <PRE>
  !        use fms_mod, only: error_mesg, FATAL, NOTE

  !        call error_mesg ('fms_mod', 'initialization not called', FATAL)
  !        call error_mesg ('fms_mod', 'fms_mod message', NOTE)
  !     </PRE>
  !   </NOTE>
  ! wrapper for the mpp error handler
  ! users should try to use the error_mesg interface

  subroutine error_mesg (level, message, routine)
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: routine
    integer,          intent(in) :: level
    character(len=512)           :: text
    logical                      :: opened
    integer                      :: istat, errunit, outunit

    !  input:
    !      routine   name of the calling routine (character string)
    !      message   message written to output   (character string)
    !      level     set to NOTE, MESSAGE, or FATAL (integer)

    select case( level )
    case(NOTE)
      text = 'NOTE'         !just FYI
    case(WARNING)
      text = 'WARNING'      !probable error
    case(FATAL)
      text = 'FATAL'        !fatal error
    case default
      text = 'WARNING: non-existent errortype (must be NOTE|WARNING|FATAL)'
    end select

    if( PRESENT(message) )text = trim(text)//': '//trim(message)

    errunit = stderr()
    outunit = stdout()

    select case( level )
    case(NOTE)
      write( outunit,'(a)' )trim(text)
    case default
      write( errunit,'(/a/)' )trim(text)
      write( outunit,'(/a/)' )trim(text)
      if( level.EQ.FATAL .OR. warnings_are_fatal )then
        call FLUSH(outunit)
        call ABORT() !automatically calls traceback on Cray systems
      end if
    end select

    return

  end subroutine error_mesg
  ! </SUBROUTINE>


  !#######################################################################
  ! <SUBROUTINE NAME="write_version_number">

  !   <OVERVIEW>
  !     Prints to the log file (or a specified unit) the (cvs) version id string and
  !     (cvs) tag name.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     Prints to the log file (stdlog) or a specified unit the (cvs) version id string
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
  !    The Fortran unit number of an open formatted file. If this unit number
  !    is not supplied the log file unit number is used (stdlog).
  !   </IN>
  ! prints module version number to the log file of specified unit number

  subroutine write_version_number (version, tag, unit)

    !   in:  version = string that contains routine name and version number
    !
    !   optional in:
    !        tag = cvs tag name that code was checked out with
    !        unit    = alternate unit number to direct output
    !                  (default: unit=stdlog)

    character(len=*), intent(in) :: version
    character(len=*), intent(in), optional :: tag
    integer,          intent(in), optional :: unit

    integer :: logunit

    logunit = stdlog()
    if (present(unit)) then
      logunit = unit
    else
      return
    endif

    if (present(tag)) then
      write (logunit,'(/,80("="),/(a))') trim(version), trim(tag)
    else
      write (logunit,'(/,80("="),/(a))') trim(version)
    endif

  end subroutine write_version_number
  ! </SUBROUTINE>

  !#####################################################################
  ! <FUNCTION NAME="stdout">
  !  <OVERVIEW>
  !    Standard fortran unit numbers.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !    This function returns the current  standard fortran unit numbers for output.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   stdout()
  !  </TEMPLATE>
  ! </FUNCTION>
  function stdout()
    integer :: stdout
    stdout = out_unit
    return
  end function stdout
  ! <FUNCTION>

  ! ========================================================================================================================




  ! ========================================================================================================================
  !GRID INFO FROM OCEAN_GRIDS.F90
  ! grid file name
  character(len=256) :: ocean_hgrid        ! will be set in set_ocean_grid_size
  character(len=256) :: ocean_vgrid = 'INPUT/ocean_vgrid.nc'

  !#######################################################################
  ! <SUBROUTINE NAME="set_ocean_grid_size">
  !
  ! <DESCRIPTION>
  ! Set the ocean grid size.  Model expects the grid specification file
  ! to be called grid_spec.nc.
  ! </DESCRIPTION>
  !
  subroutine set_ocean_grid_size(Grid, grid_file, grid_name)

    type(ocean_grid_type), intent(inout)        :: Grid
    character(len=*),      intent(in), optional :: grid_file
    character(len=*),      intent(in), optional :: grid_name

    integer                                   :: siz(4)
    integer                                   :: m
    integer                                   :: ntiles, ncontacts
    integer, dimension(2)                     :: tile1, tile2
    integer, dimension(2)                     :: istart1, iend1, jstart1, jend1
    integer, dimension(2)                     :: istart2, iend2, jstart2, jend2
    character(len=256)                        :: grd_file, ocean_mosaic, attvalue
    real                                      :: f_plane_latitude

    integer :: stdoutunit
    stdoutunit=stdout()

    Grid%ni=0 ; Grid%nj=0 ; Grid%nk=0
    Grid%mosaic = .false.
    grd_file = "INPUT/grid_spec.nc"
    if(present(grid_file)) grd_file = grid_file
    if (.not. file_exist(grd_file) ) then
      call error_mesg(FATAL, '==> Error from ocean_grids_mod(set_ocean_grid_size): '// &
      'grid specification file '//grd_file//' does not exist')
    endif
    !
    !  Determine if the grid is mosaic file
    !
    if(field_exist(grd_file, 'ocn_mosaic_file') .or. field_exist(grd_file, 'gridfiles') ) then ! read from mosaic file
      write(stdoutunit,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from mosaic version grid'
      grid_version = VERSION_2
      Grid%mosaic = .true.
      if( field_exist(grd_file, 'ocn_mosaic_file') ) then ! coupler mosaic
        call read_data(grd_file, "ocn_mosaic_file", ocean_mosaic)
        ocean_mosaic = "INPUT/"//trim(ocean_mosaic)
      else
        ocean_mosaic = trim(grd_file)
      end if
      ntiles = get_mosaic_ntiles(ocean_mosaic)
      if(ntiles .NE. 1) call error_mesg(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
      'ntiles should be 1 for ocean mosaic, contact developer')
      call read_data(ocean_mosaic, "gridfiles", ocean_hgrid)
      ocean_hgrid = 'INPUT/'//trim(ocean_hgrid)
      call field_size(ocean_hgrid, 'x', siz)

      if(mod(siz(1),2) .NE. 1) call error_mesg(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
      'x-size of x in file '//trim(ocean_hgrid)//' should be 2*ni+1')
      if(mod(siz(2),2) .NE. 1) call error_mesg(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
      'y-size of x in file '//trim(ocean_hgrid)//' should be 2*nj+1')
      Grid%ni = siz(1)/2
      Grid%nj = siz(2)/2
      call field_size(ocean_vgrid , "zeta", siz)
      if(mod(siz(1),2) .NE. 1) call error_mesg(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
      'size of dimension zeta in file '//trim(ocean_vgrid)//' should be 2*nk+1')
      Grid%nk  = siz(1)/2
    else  if(field_exist(grd_file, 'x_T')) then
      ocean_hgrid = grd_file
      write(stdoutunit,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from new version grid'
      grid_version = VERSION_1
      call field_size( ocean_hgrid, 'x_T', siz)
      Grid%ni = siz(1)
      Grid%nj = siz(2)
      call field_size( ocean_hgrid, "zt", siz)
      Grid%nk  = siz(1)
    else if(field_exist(grd_file, 'geolon_t')) then
      ocean_hgrid = grd_file
      write(stdoutunit,*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from old version grid'
      grid_version = VERSION_0
      call field_size( ocean_hgrid, 'geolon_t', siz)
      Grid%ni = siz(1)
      Grid%nj = siz(2)
      call field_size( ocean_hgrid, "zt", siz)
      Grid%nk  = siz(1)
    else
      call error_mesg(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
      'x_T, geolon_t, ocn_mosaic_file, gridfiles does not exist in file ' //trim(grd_file))
    endif

    if (Grid%ni == 0 .or. Grid%nj == 0 .or. Grid%nk == 0) then
      write(stdoutunit,*) '==>Error reading grid information from ',trim(grd_file),'. Make sure file exists'
      call error_mesg(FATAL,'==>Error reading grid information from grid file.  Are you sure file exists?')
    endif

    Grid%cyclic_x=.false.;Grid%cyclic_y=.false.;Grid%tripolar=.false.
    Grid%f_plane=.false.; Grid%beta_plane=.false.
    Grid%f_plane_latitude = -999.0

    if(grid_version == VERSION_2) then
      !z1l: f_plane, beta_plane area not supported in mosaic grid. Need to think about to implement this.
      if(field_exist(ocean_mosaic, "contacts") ) then
        ncontacts = get_mosaic_ncontacts(ocean_mosaic)
        if(ncontacts < 1) call error_mesg(FATAL,'==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
        'number of contacts should be larger than 0 when field contacts exist in file '//trim(ocean_mosaic) )
        if(ncontacts > 2) call error_mesg(FATAL,'==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
        'number of contacts should be no larger than 2')
        call get_mosaic_contact( ocean_mosaic, tile1(1:ncontacts), tile2(1:ncontacts),           &
        istart1(1:ncontacts), iend1(1:ncontacts), jstart1(1:ncontacts), jend1(1:ncontacts), &
        istart2(1:ncontacts), iend2(1:ncontacts), jstart2(1:ncontacts), jend2(1:ncontacts)  )
        do m = 1, ncontacts
          if(istart1(m) == iend1(m) ) then  ! x-direction contact, only cyclic condition
            if(istart2(m) .NE. iend2(m) ) call error_mesg(FATAL,  &
            "==>Error from ocean_grids_mod(set_ocean_grid_size): only cyclic condition is allowed for x-boundary")
            Grid%cyclic_x = .true.
          else if( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic or folded-north
            if(jstart2(m) .NE. jend2(m) ) call error_mesg(FATAL,  &
            "==>Error from ocean_grids_mod(set_ocean_grid_size): "//&
            "only cyclic/folded-north condition is allowed for y-boundary")
            if( jstart1(m) == jstart2(m) ) then ! folded north
              Grid%tripolar = .true.
            else
              Grid%cyclic_y = .true.
            endif
          else
            call error_mesg(FATAL,  &
            "==>Error from ocean_grids_mod(set_ocean_grid_size): invalid boundary contact")
          end if
        end do
      end if
    else
      if( get_global_att_value(ocean_hgrid, "x_boundary_type", attvalue) ) then
        if(attvalue == 'cyclic') Grid%cyclic_x = .true.
      end if
      if( get_global_att_value(ocean_hgrid, "y_boundary_type", attvalue) ) then
        if(attvalue == 'cyclic') then
          Grid%cyclic_y = .true.
        else if(attvalue == 'fold_north_edge') then
          Grid%tripolar = .true.
        end if
      end if
      if(get_global_att_value(ocean_hgrid, "f_plane", attvalue) ) then
        if(attvalue == 'y') Grid%f_plane = .true.
      end if
      if(get_global_att_value(ocean_hgrid, "beta_plane", attvalue) ) then
        if(attvalue == 'y') Grid%beta_plane = .true.
      end if
      if( get_global_att_value(ocean_hgrid, "f_plane_latitude", f_plane_latitude) ) then
        Grid%f_plane_latitude = f_plane_latitude
      end if
    end if

    if(Grid%cyclic_x) then
      call error_mesg(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): x_boundary_type is cyclic')
    else
      call error_mesg(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): x_boundary_type is solid_walls')
    end if

    if(Grid%tripolar) then
      call error_mesg(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is fold_north_edge')
    else if(Grid%cyclic_y) then
      call error_mesg(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is cyclic')
    else
      call error_mesg(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is solid_walls')
    endif

    if(Grid%f_plane .or. Grid%beta_plane) then
      if(abs(Grid%f_plane_latitude) > 90.0 ) then
        write(stdoutunit,*)  "==>Error from ocean_grids_mod(set_ocean_grid_size): "// &
        "f_plane_latitude must be between -90 and 90 degrees. Please check your grid file."
        call error_mesg(FATAL, "==>Error from ocean_grids_mod(set_ocean_grid_size): "// &
        "f_plane_latitude must be between -90 and 90 degrees. Please check your grid file.")
      end if
      write(stdoutunit,'(a,f8.2)') '==>f-plane/beta-plane latitude = ',Grid%f_plane_latitude
    end if

    if(Grid%f_plane) then
      write(stdoutunit,*) '==>NOTE: Setting geometry and Coriolis according to f-plane'
    endif
    if(Grid%beta_plane) then
      write(stdoutunit,*) '==>NOTE: Setting geometry and Coriolis according to beta-plane'
    endif

    if(Grid%tripolar) then
      write (stdoutunit,'(1x,/a)')  ' ==> Note: Energy conversion errors are nontrivial when using tripolar=.true.'
      write (stdoutunit,'(7x,a)')   'The cause is related to the need to update redundantly computed information'
      write (stdoutunit,'(7x,a)')   'across the Arctic bipolar fold in a bit-wise exact manner for terms contributing'
      write (stdoutunit,'(7x,a)')   'to the energy conversion analysis.  The extra code and mpp calls have not been'
      write (stdoutunit,'(7x,a)')   'implemented.'
    endif

    if (PRESENT(grid_name)) then
      Grid%name = grid_name
    else
      Grid%name = 'ocean'
    endif

  end subroutine set_ocean_grid_size
  ! </SUBROUTINE> NAME="set_ocean_grid_size"


  !#######################################################################
  ! <SUBROUTINE NAME="set_ocean_hgrid_arrays">
  !
  ! <DESCRIPTION>
  ! Define horizontal (and some vertical) grid arrays.
  !
  !---------------------------------------------------------------------------------------------------------------------
  ! Grid%       grid_spec      grid_spec     grid_spec                       Description
  !  var        field          field         field
  !             VERSION_0      VERSION_1     VERSION_2(mosaic)
  !---------------------------------------------------------------------------------------------------------------------
  !                                                    ocean_vgrid.nc
  !                                                    k=1,nk
  ! zt          zt             zt            zeta(2k-1)
  ! zw          zw             zb            zeta(2k)
  !
  !                                                    ocean_hgrid.nc
  !                                                    i=1,ni
  !                                                    j=1,nj
  ! grid_x_t    gridlon_t      grid_x_T      x(2i  ,2)
  ! grid_x_u    gridlon_vert_t grid_x_C      x(2i+1,1)
  ! grid_y_t    gridlat_t      grid_y_T      y(ni/4,2j)
  ! grid_y_u    gridlat_vert_t grid_y_C      y(ni/4,2j+1)
  !
  !T
  ! xt(i,j)     geolon_t(i,j)  x_T(i,j)      x(2i,2j)
  ! yt          geolat_t       y_T           y(2i,2j)
  ! dtw         dtw            ds_01_11_T    dx(2i-1,2j)                 distance to western face of t cell
  ! dte         dte            ds_11_21_T    dx(2i,2j)                   distance to eastern face of t cell
  ! dts         dts            ds_10_11_T    dy(2i,2j-1)                 distance to southern face of t cell
  ! dtn         dtn            ds_11_12_T    dy(2i,2j)                   distance to northern face of t cell
  !
  ! dxt         dxt            ds_01_21_T    dx(2i,2j)    +dx(2i-1,2j)   width of t cell
  ! dxtn        dxtn           ds_02_22_T    dx(2i-1,2j+1)+dx(2i,2j+1)   width of northern face of t cell
  ! dxte        dxte           ds_00_20_C    dx(2i,2j)    +dx(2i+1,2j)   distance to adjacent t cell to the east!
  ! dyt         dyt            ds_10_12_T    dy(2i,2j)    +dy(2i,2j-1)   height of t cell
  ! dytn        dytn           ds_00_02_C    dy(2i,2j)    +dy(2i,2j+1)   distance to adjacent t cell to the north!
  ! dyte        dyte           ds_20_22_T    dy(2i+1,2j-1)+dy(2i+1,2j)   height of eastern face of t cell
  !
  !C
  ! NOTE: The "first" (I,J) C-cell is the one shifted NE of the "first" (I,J) T-cell
  !
  !
  ! xu          geolon_c       x_C           x(2i+1,2j+1)
  ! yu          geolat_c       y_c           y(2i+1,2j+1)
  ! dxu         dxu            ds_01_21_C    dx(2i+1,2j+1)+dx(2i,2j+1)   width of u cell
  ! dxun        dxun           ds_02_22_C    dx(2i,2j+2)+dx(2i+1,2j+2)   width of northern face of u cell
  ! dyu         dyu            ds_10_12_C    dy(2i+1,2j+1)+dy(2i+1,2j)   height of u cell
  ! dyue        dyue           ds_20_22_C    dy(2i+2,2j)+dy(2i+2,2j+1)   height of eastern face of u cell
  !
  ! dyun        dyun           ds_11_12_C    dy(2i+1,2j+1)+dy(2i+1,2j+2) distance to adjacent u cell to the north
  !                            +ds_10_11_C(i,j+1)                         satisfies sum rule dyte(i,j)=dyun(i,j-1)
  ! dxue        dxue           ds_11_21_C    dx(2i+1,2j+1)+dx(2i+2,2j+1) distance to adjacent u cell to the east!
  !                            +ds_01_11_C(i+1,j)
  !
  ! duw         duw            ds_01_11_C    dx(2i,2j+1)                 distance to western face of u cell
  ! due         due            ds_11_21_C    dx(2i+1,2j+1)               distance to eastern face of u cell
  ! dus         dus            ds_10_11_C    dy(2i+1,2j)                 distance to southern face of u cell
  ! dun         dun            ds_11_12_C    dy(2i+1,2j+1)               distance to northern face of u cell
  !
  ! sin_rot     sin_rot        angle_C       sin(angle_dx(2*i+1,2*j+1)   sin of rotation angle at corner cell centers
  ! cos_rot     cos_rot        angle_C       cos(angle_dx(2*i+1,2*j+1)   cos of rotation angle at corner cell centers
  !
  !Following are the available fields in mosaic files
  !--------------------------------------------------------
  !Mosaic file     fields
  !--------------------------------------------------------
  !ocean_hgrid.nc  x, y, dx, dy, angle_dx, area
  !ocean_vgrid.nc  zeta
  !topog.nc        depth
  !
  ! </DESCRIPTION>
  subroutine set_ocean_hgrid_arrays(Grid)

    type(ocean_grid_type),   intent(inout) :: Grid

    real, dimension(:),   allocatable          :: data
    real, dimension(:,:), allocatable          :: tmp
    real, dimension(:,:), allocatable          :: tmp_local
    real, dimension(:,:), allocatable          :: tmp1_local
    real, dimension(:,:), allocatable          :: tmp2_local
    real, dimension(:,:), allocatable          :: tmpx
    real, dimension(:,:), allocatable          :: tmpy
    integer, dimension(:),allocatable          :: xext, yext

    integer :: isd2, ied2, jsd2, jed2
    integer :: isc2, iec2, jsc2, jec2
    integer :: isg2, ieg2, jsg2, jeg2
    integer :: isc1, iec1, jsc1, jec1
    integer :: i, j, k, lon_Tripol
    integer :: ioff, joff
    integer :: start(4), nread(4)

    integer :: stdoutunit
    stdoutunit=stdout()

    ! set the grid points coordinates (degrees) and grid spacing (degrees)

    #ifndef MOM_STATIC_ARRAYS

    ni = Grid%ni;nj = Grid%nj; nk = Grid%nk

    ! WHAT ARE THE VECTOR NOTATIONS? ARE THEY LEFT OVER FROM DOMAINS???
    ! integer :: isc, iec, jsc, jec        ! computational domain indices
    ! integer :: isd, ied, jsd, jed        ! local indices, consistent with domain2d
    ! integer :: isg, ieg, jsg, jeg        ! global indices
    ! integer :: isa, iea, jsa, jea        ! active indices (for minimizing comm2d calls)

    allocate (Grid%xt(isd:ied,jsd:jed))
    allocate (Grid%yt(isd:ied,jsd:jed))
    allocate (Grid%xu(isd:ied,jsd:jed))
    allocate (Grid%yu(isd:ied,jsd:jed))
    allocate (Grid%grid_x_t(ni))
    allocate (Grid%grid_y_t(nj))
    allocate (Grid%grid_x_u(ni))
    allocate (Grid%grid_y_u(nj))

    allocate (Grid%zt(nk))
    allocate (Grid%zw(nk))

    allocate (Grid%dat(isd:ied,jsd:jed))
    allocate (Grid%dxtn(isd:ied,jsd:jed))
    allocate (Grid%dyte(isd:ied,jsd:jed))
    allocate (Grid%datr(isd:ied,jsd:jed))

    #endif

    !--- initialize grid data

    Grid%dxtn=1.0;  Grid%dxte=1.0;
    Grid%dxt=epsln; Grid%dxu=epsln ;
    Grid%zw=0.0;

    select case( grid_version )
    case( VERSION_0 )
       call read_data(ocean_hgrid, "zw", Grid%zw)

    case( VERSION_2 )
       allocate(data(2*nk+1) )
       call read_data(ocean_vgrid, "zeta", data)
       do k=1,nk
          Grid%zt(k) = data(2*k)
       enddo
       deallocate(data)

       !--- The following adjust is for the consideration of static memory.
       ioff = isc1 - isc; joff = jsc1 - jsc
       isc2 = isc2 - 2*ioff; iec2 = iec2 - 2*ioff
       jsc2 = jsc2 - 2*joff; jec2 = jec2 - 2*joff
       isd2 = isd2 - 2*ioff; ied2 = ied2 - 2*ioff
       jsd2 = jsd2 - 2*joff; jed2 = jed2 - 2*joff


       start = 1; nread = 1
       start(2) = 2; nread(1) = 2*ni+1; nread(2) = 2
       allocate(tmpx(2*ni+1,2), tmpy(2,2*nj+1) )
       call read_data(ocean_hgrid, "x", tmpx, start, nread)

       do i = 1, ni
          Grid%grid_x_t(i) = tmpx(2*i,  1)
          Grid%grid_x_u(i) = tmpx(2*i+1,2)
       enddo
       lon_Tripol = ni/4 ! 90 for 1 degree grid
       start = 1; nread = 1
       start(1) = 2*lon_Tripol+1; nread(1) = 2; nread(2) = 2*nj+1
       call read_data(ocean_hgrid, "y", tmpy, start, nread)

       deallocate(tmpx, tmpy)
       allocate(tmpx(isc2:iec2+1, jsc2:jec2+1))
       allocate(tmpy(isc2:iec2+1, jsc2:jec2+1))
       call read_data(ocean_hgrid, "x", tmpx)
       call read_data(ocean_hgrid, "y", tmpy)
    end select


    ! --- Grid%dxt
    select case(grid_version)
    case(VERSION_0)
       call read_data(ocean_hgrid, 'dxt', Grid%dxt)
    case(VERSION_1)
       call read_data(ocean_hgrid, 'ds_01_21_T', Grid%dxt)
    case(VERSION_2)
       Grid%dxt(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2,2*jsc:2*jec:2) + tmpx(2*isc:2*iec:2,2*jsc:2*jec:2)
    end select

    !--- Grid%dyt
    select case(grid_version)
    case(VERSION_0)
       call read_data(ocean_hgrid, 'dyt', Grid%dyt)
    case(VERSION_1)
       call read_data(ocean_hgrid, 'ds_10_12_T', Grid%dyt)
    case(VERSION_2)
       Grid%dyt(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2,2*jsc-1:2*jec-1:2) + tmpy(2*isc:2*iec:2,2*jsc:2*jec:2)
    end select


    ! --- Grid%dat
    Grid%dat(:,:)  = Grid%dxt(:,:)*Grid%dyt(:,:)

    ! --- Grid%dxtn
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dxtn', Grid%dxtn)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_02_22_T', Grid%dxtn)
    case(VERSION_2)
      Grid%dxtn(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2,2*jsc+1:2*jec+1:2) + tmpx(2*isc:2*iec:2,2*jsc+1:2*jec+1:2)
    end select

    ! --- Grid%dyte
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dyte', Grid%dyte)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_20_22_T', Grid%dyte)
    case(VERSION_2)
      Grid%dyte(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2,2*jsc-1:2*jec-1:2) + tmpy(2*isc+1:2*iec+1:2,2*jsc:2*jec:2)
    end select

    ! set reciprocals and related quantities
    do j=jsd,jed
      do i=isd,ied
        Grid%datr(i,j)       = 1.0/(Grid%dat(i,j)+epsln)
      enddo
    enddo

    ! set quantities which account for rotation of basis vectors
    allocate(tmp_local(isd:ied,jsd:jed))
    tmp_local=1.0

    deallocate(tmp_local, tmp1_local, tmp2_local)

  end subroutine set_ocean_hgrid_arrays
  ! </SUBROUTINE> NAME="set_ocean_hgrid_arrays"


  !#######################################################################
  ! <SUBROUTINE NAME="set_ocean_vgrid_arrays">
  !
  ! <DESCRIPTION>
  ! Compute vertical (and some horizontal) grids for ocean model.
  ! Also compute axes information for diagnostic manager.
  ! </DESCRIPTION>
  !
  subroutine set_ocean_vgrid_arrays (Grid, grid_file, vert_coordinate_type)

    type(ocean_grid_type),   intent(inout)        :: Grid
    character(len=*),        intent(in), optional :: grid_file
    integer,                 intent(in)           :: vert_coordinate_type

    character(len=128)  :: grd_file
    character(len=128)  :: ocean_topog = "INPUT/topog.nc"

    grd_file = "INPUT/grid_spec.nc"
    if(present(grid_file)) grd_file = grid_file

    real    :: ht_fax, ht_fay, ht_jp1_fax, ht_ip1_fay
    real    :: tcella_sphere, ucella_sphere
    real    :: wet_t_points, wet_u_points

    integer :: i, j, k, kzt, kzu, kmin

    logical :: have_obc=.false.
    logical :: obc=.false.

    real    :: sigma_sign, zwnk, zwnk_r

    integer :: stdoutunit
    stdoutunit=stdout()


    if (PRESENT(obc)) have_obc = obc

    #ifndef MOM_STATIC_ARRAYS
      allocate (Grid%tmask(isd:ied,jsd:jed,nk))
      allocate (Grid%dst(nk))
    #endif

    ! set depth arrays
    Grid%dzt(1)         = Grid%zw(1)
    do k=2,nk-1
       Grid%dzt(k)      = Grid%zw(k)-Grid%zw(k-1)
    enddo
    Grid%dzt(nk)        = Grid%zw(nk)-Grid%zw(nk-1)

    Grid%kmt=0;Grid%kmu=0;Grid%ht=0;Grid%hu=0
    if(field_exist(grd_file, 'depth_t') ) then ! new grid file
       call read_data(grd_file, "num_levels", Grid%kmt(isc:iec,jsc:jec))
       if(field_exist(grd_file, 'depth_c')) then
          if(.not. field_exist(grd_file, 'num_levels_c')) call error_mesg(FATAL, &
             'ocean_topog_mod: depth_c exist but num_levels_c does not exist in file '//trim(grd_file) )
       elseif(field_exist(grd_file, 'num_levels_c')) then
          call error_mesg(FATAL, &
             'ocean_topog_mod: num_levels_c exist but depth_c does not exist in file '//trim(grd_file) )
       endif
    elseif(field_exist(grd_file, 'ht') ) then ! old grid file
       call read_data(grd_file, "kmt", Grid%kmt(isc:iec,jsc:jec))
    elseif(file_exist(ocean_topog)) then
       call read_data(ocean_topog, 'depth', Grid%ht(isc:iec,jsc:jec))
       !--- calculate kmt based on ht and zw.
       do j=jsc,jec
          do i=isc,iec
             if(Grid%ht(i,j) <= 0.0) then
                Grid%kmt(i,j) = 0
             else
                Grid%kmt(i,j) = nearest_index(Grid%ht(i,j), Grid%zw)
                if( Grid%ht(i,j) >= Grid%zw(Grid%kmt(i,j)) + min_thickness) then
                   if(Grid%kmt(i,j) < nk) then
                      Grid%kmt(i,j) = Grid%kmt(i,j) + 1
                   endif
                endif
             endif
          enddo
       enddo
    else
       call error_mesg(FATAL, 'ocean_topog_mod: depth_t and ht do not exist in file '//trim(grd_file)// &
                      ', also file '//trim(ocean_topog)//' does not exist')
    endif

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

    allocate(rho0_profile(nk))
    rho0_profile(:) = rho0
    if(read_rho0_profile) then
      write(stdoutunit,'(a)') '==>Warning: ocean_grids_mod: read rho0 profile to set dst grid spacing.'
      write(stdoutunit,'(a)') '            This option is experimental, and NOT generally supported.'
      if(vert_coordinate_class==DEPTH_BASED) then
        write(stdoutunit,'(a)') '   Since using DEPTH_BASED vertical coordinates, rho0_profile = rho0.'
      endif
      if(vert_coordinate_class==PRESSURE_BASED) then
        call read_data('INPUT/rho0_profile.nc','rho0_profile', rho0_profile)
      endif
    endif

    ! compute static s-grid information
    if(vert_coordinate == GEOPOTENTIAL .or. vert_coordinate == ZSTAR) then
      do k=1,nk
        Grid%dst(k)   = Grid%dzt(k)
      enddo
    elseif(vert_coordinate == PRESSURE .or. vert_coordinate == PSTAR) then
      do k=1,nk
        Grid%dst(k)   = -rho0_profile(k)*grav*Grid%dzt(k)
      enddo
    elseif(vert_coordinate == ZSIGMA .or. vert_coordinate == PSIGMA) then
       ! assume initial bottom pressure = rho0*g*ht, which means that
       ! PSIGMA and ZSIGMA initialization are same, to within sign.
       ! -1 <= ZSIGMA <= 0, so dst and dsw are positive
       !  0 <= PSIGMA <= 1, so dst and dsw are negative
       ! Recall s-coordinate is dimensionless when use ZSIGMA or PSIGMA.
       if(vert_coordinate==ZSIGMA) then
           sigma_sign = -1.0
       elseif(vert_coordinate==PSIGMA) then
           sigma_sign =  1.0
       endif

       zwnk_r = 1.0/(epsln + zwnk)

       do k=1,nk
          Grid%sw(k) = sigma_sign*Grid%zw(k)*zwnk_r
       enddo

       ! set s-grid increments from st and sw
       Grid%dst(1) = - Grid%sw(1)
       do k=2,nk-1
          Grid%dst(k) = -(Grid%sw(k)  - Grid%sw(k-1))
       enddo
       Grid%dst(nk) = -(Grid%sw(nk)-Grid%sw(nk-1))
       Grid%dsw(nk) = -(Grid%sw(nk)-Grid%st(nk))

    endif

  end subroutine set_ocean_vgrid_arrays
  ! </SUBROUTINE> NAME="set_ocean_vgrid_arrays"




  !#######################################################################
  ! <SUBROUTINE NAME="ocean_thickness_init">
  !
  ! <DESCRIPTION>
  ! Initialize the thickness type.
  !
  ! For pressure-based vertical coordinates, this initialization here
  ! is preliminary.
  ! </DESCRIPTION>
  !
  subroutine ocean_thickness_init  (Time, Time_steps, Grid, Ext_mode, Thickness, &
                                   ver_coordinate, ver_coordinate_class, ver_coordinate_type, dtimein)

    type(ocean_time_type),          intent(in)           :: Time
    type(ocean_time_steps_type),    intent(in)           :: Time_steps
    type(ocean_grid_type),          intent(inout)        :: Grid
    type(ocean_external_mode_type), intent(in)           :: Ext_mode
    type(ocean_thickness_type),     intent(inout)        :: Thickness
    integer,                        intent(in)           :: ver_coordinate
    integer,                        intent(in)           :: ver_coordinate_class
    integer,                        intent(in)           :: ver_coordinate_type
    real,                           intent(in)           :: dtimein

    integer :: ioun, io_status, ierr
    integer :: i,j,k
    logical :: error_flag=.false.

    integer :: stdoutunit,stdlogunit
    stdoutunit=stdout();stdlogunit=stdlog()

    grav_r    = 1.0/grav
    grav_rho0 = grav*rho0
    dtime     = dtimein

    vert_coordinate       = ver_coordinate
    vert_coordinate_class = ver_coordinate_class
    vert_coordinate_type  = ver_coordinate_type

    if(vert_coordinate_class==DEPTH_BASED) then
      convert_factor=1.0
    elseif(vert_coordinate_class==PRESSURE_BASED) then
      convert_factor=c2dbars
    endif

    if(enforce_positive_dzt) then
      write(stdoutunit,'(/a)') '==>Warning: running with enforce_positive_dzt=.true. '
      write(stdoutunit,'(a)')  '            This option artifically truncates cell size to enforce dzt > 0.'
      write(stdoutunit,'(a,f10.4/)') '      Minimum dzt is set to ', thickness_dzt_min
    endif

    if(vert_coordinate == ZSIGMA .or. vert_coordinate == PSIGMA) then
       write(stdoutunit,'(a,f12.4,a)') &
       '==>To remove division by zero with sigma-coord, add fictitious water of depth(m) ' &
       ,depth_min_for_sigma, ' over land. This water does not contribute to budgets.'
      if(Grid%tripolar) then
         call error_mesg(WARNING, &
         '==>Warning ocean_thickness_mod: zsigma and psigma may have problems with tripolar. Testing incomplete.')
      endif
    endif

    if(vert_coordinate == GEOPOTENTIAL .and. linear_free_surface) then
       write(stdoutunit,'(a)') &
       '==>Running GEOPOTENTIAL model with linear_free_surface, so cell thicknesses are constant in time.'
       call error_mesg(WARNING, &
         '==>Running GEOPOTENTIAL model with linear_free_surface, so cell thicknesses are constant in time.')
    endif

    if(thickness_method=='energetic') then
        Thickness%method = ENERGETIC
        write(stdoutunit,'(a)') &
        '==>Note: running ocean_thickness with thickness_method=energetic.'
    elseif(thickness_method=='finitevolume') then
        Thickness%method = FINITEVOLUME
        write(stdoutunit,'(a)') &
        '==>Warning: ocean_thickness w/ thickness_method=finitevolume is experimental. Testing incomplete.'
         call error_mesg(WARNING, &
         '==>Warning ocean_thickness w/ thickness_method=finitevolume is experimental. Testing incomplete.')
    endif

  #ifndef MOM_STATIC_ARRAYS
    nk = Grid%nk

    allocate (Thickness%rho_dzt(isd:ied,jsd:jed,nk,3) )

  #endif

    ! set rescale_rho0_mask=rescale_rho0_value in those regions where we modify
    ! the rho0 value by a fraction. Otherwise rho0 is unscaled.
    allocate (rescale_rho0_mask(isd:ied,jsd:jed) )
    rescale_rho0_mask(:,:) = 1.0
    if(read_rescale_rho0_mask .and. vert_coordinate_class==PRESSURE_BASED) then

        allocate (data(isd:ied,jsd:jed) )
        data(:,:) = 0.0
        call read_data('INPUT/basin_mask','basin_mask',data)

        ! the following is specific to the mask used at GFDL
        if(rescale_rho0_mask_gfdl) then
            write(stdoutunit,'(a,f12.6)') &
            '==>Note: running ocean_thickness with rho0 in selected basin modified to rho0(kg/m3) = ',&
            rho0*rescale_rho0_value
            do j=jsc,jec
               do i=isc,iec
                  if(data(i,j)==rescale_rho0_basin_label) then
                      rescale_rho0_mask(i,j)=rescale_rho0_value
                  else
                      rescale_rho0_mask(i,j)=1.0
                  endif
               enddo
            enddo
        endif
    endif

    ! rho0_profile for setting dst with pressure based vertical coordinate models
    allocate(rho0_profile(nk))
    rho0_profile(:) = rho0
    if(read_rho0_profile) then

        write(stdoutunit,'(a)') &
         '==>Warning: ocean_thickness_mod: read rho0 profile to set dst w/ pressure models.'
        write(stdoutunit,'(a)') &
        '             This option is experimental, and NOT generally supported.'
        if(vert_coordinate_class==DEPTH_BASED) then
            write(stdoutunit,'(a)') &
                 '   Since using DEPTH_BASED vertical coordinates, rho0_profile = rho0.'
        endif

        if(vert_coordinate_class==PRESSURE_BASED) then
            call read_data('INPUT/rho0_profile.nc','rho0_profile', rho0_profile)
            write(stdoutunit,'(a)') 'rho0_profile is used to define pressure grid and pressure gradients.'
            do k=1,nk
               write(stdoutunit,'(a,i4,a,e22.12)') 'rho0_profile(',k,') = ',rho0_profile(k)
               if(rho0_profile(k) <= 0.0) then
                   call error_mesg(FATAL, &
                   '==>ocean_thickness_mod: rho0_profile must be > 0.0, even in rock, since we divide by rho0_profile.')
               endif
            enddo
        endif

    endif

    if(pbot0_simple) then
       write(stdoutunit,'(/a)') '==>Warning from ocean_thickness_mod: pbot0_simple=.true.'
    endif

    ! initialize vertical grid arrays
    call thickness_initialize(Grid, Thickness)

  end subroutine ocean_thickness_init
  ! </SUBROUTINE> NAME="ocean_thickness_init"



  !#######################################################################
  ! <SUBROUTINE NAME="thickness_initialize">
  !
  ! <DESCRIPTION>
  ! Initialize vertical thicknesses of grid cells.
  ! For Boussinesq models, this code assumes the
  ! surface heights eta_t and eta_u are zero.
  !
  ! The values here are relevant for the time=0 initialization
  ! of the model. Some time independent arrays are also set here,
  ! but they are over-written if there is a restart file.
  !
  ! For pressure based vertical coordinates, the results here
  ! assume density = rho0_profile(k). The values are readjusted in
  ! ocean_thickness_init_adjust after we have determined the initial
  ! in situ density.  This readjustment may involve enforcing an
  ! initial eta field of zero value, or the default which allows
  ! for eta to initially be nonzero.
  !
  ! </DESCRIPTION>
  !
  subroutine thickness_initialize(Grid, Thickness)

   type(ocean_grid_type),      intent(inout) :: Grid
   type(ocean_thickness_type), intent(inout) :: Thickness

   real                              :: sigma_sign, density_tmp
   integer                           :: i,j,k,n,kb,kmin

   ! initialization to zero
   Thickness%sea_lev(:,:) = 0.0
   do k=1,nk
      do j=jsd,jed
         do i=isd,ied
            Thickness%rho_dzt_tendency(i,j,k) = 0.0
         enddo
      enddo
   enddo

   ! vertical increments for non-terrain following coordinates
   if(vert_coordinate_type /= TERRAIN_FOLLOWING) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                Thickness%dzt(i,j,k)       = Grid%dzt(k)
                Thickness%depth_zwt(i,j,k) = Grid%zw(k)
             enddo
          enddo
       enddo

       ! modifications for partial step representation of bottom topography
       if(full_step_topography) then
           do j=jsd,jed
              do i=isd,ied
                 kb = Grid%kmt(i,j)
                 if (kb > 1) then
                     Grid%ht(i,j) = Thickness%depth_zwt(i,j,kb)
                 endif
              enddo
           enddo
       else
           do j=jsd,jed
              do i=isd,ied
                 kb = Grid%kmt(i,j)
                 if (kb > 1) then
                     Thickness%depth_zwt(i,j,kb)= Grid%ht(i,j)
                     Thickness%dzt(i,j,kb)      = &
                          Thickness%depth_zwt(i,j,kb)   - Thickness%depth_zwt(i,j,kb-1)
                 endif
              enddo
           enddo
       endif

   ! vertical increments for terrain following coordinates
   else

       ! -1 <= ZSIGMA <= 0, so dst and dsw are positive
       !  0 <= PSIGMA <= 1, so dst and dsw are negative
       ! Recall s-coordinate is dimensionless when use ZSIGMA or PSIGMA.
       if(vert_coordinate==ZSIGMA) then
           sigma_sign =-1.0
       elseif(vert_coordinate==PSIGMA) then
           sigma_sign = 1.0
       endif

       ! set depth_st and depth_swt according to Grid arrays
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                Thickness%dst(i,j,k)       = Grid%dst(k)
             enddo
          enddo
       enddo

       ! now get the z-increments
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                Thickness%dzt_dst(i,j,k) = -sigma_sign*max(Grid%ht(i,j),depth_min_for_sigma)
                Thickness%dzt(i,j,k)     = Thickness%dzt_dst(i,j,k)*Thickness%dst(i,j,k)
             enddo
          enddo
       enddo


    endif  ! endif for ZSIGMA and PSIGMA

   ! compute rho*dz values with rho0_profile*rescale_rho0_mask for initialization
   Thickness%mass_u(:,:,:) = 0.0
   do n=1,3
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Thickness%rho_dzt(i,j,k,n) = rho0_profile(k)*rescale_rho0_mask(i,j)*Thickness%dzt(i,j,k)
            enddo
         enddo
      enddo
   enddo

   end subroutine thickness_initialize
   ! </SUBROUTINE> NAME="thickness_initialize"



  !#############################################################################
  ! This routine will get the actual file name, as well as if read_dist is true or false.
  ! return true if such file exist and return false if not.
  function get_file_name(orig_file, actual_file)
    character(len=*),                 intent(in) :: orig_file
    character(len=*),                intent(out) :: actual_file
    logical                                      :: get_file_name

    logical                       :: fexist
    character(len=256)            :: fname

    fexist          = .false.
    get_file_name   = .false.

    !--- The file maybe not netcdf file, we just check the original file.
    if(index(orig_file, '.nc', back=.true.) == 0) then
      inquire (file=trim(orig_file), exist=fexist)
      if(fexist) then
        actual_file = orig_file
        get_file_name = .true.
        return
      endif
    endif

    if(fexist) then
      get_file_name = .true.
      return
    endif

    inquire (file=trim(actual_file), exist=fexist)
    if(fexist) then
      get_file_name = .true.
      return
    endif

  end function get_file_name

  !#############################################################################
  subroutine get_file_unit(filename, unit, index_file)
    character(len=*),         intent(in) :: filename
    integer,                 intent(out) :: unit, index_file

    logical  :: file_opened
    integer  :: i

    ! Need to check if filename has been opened or not
    file_opened=.false.
    do i=1,num_files_r
      if (files_read(i)%name == trim(filename))  then
        index_file = i
        unit = files_read(index_file)%unit
        return
      endif
    enddo

    ! need to open the file now
    ! Increase num_files_r and set file_type
    if(num_files_r == max_files_r) &  ! need to have bigger max_files_r
    call error_mesg(FATAL,'fms_io(get_file_unit): max_files_r exceeded, increase it via fms_io_nml')
    num_files_r=num_files_r + 1

    files_read(num_files_r)%name = trim(filename)
    allocate(files_read(num_files_r)%var (max_fields) )
    files_read(num_files_r)%nvar = 0
    index_file = num_files_r
    files_read(index_file)%unit = unit

  end subroutine get_file_unit

  ! OR
  !#######################################################################
  ! This function is only for global meta and vgrid. So get_file_unit should
  ! always linked to file opened_file(0)
  ! <FUNCTION NAME="get_file_unit">
  !   <OVERVIEW>
  !    returns the io unit corresponding to filename.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !    If the file filename is already open, return the io unit of this
  !    opened file. Otherwise will open the file and return the io unit.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     get_file_unit(filename)
  !   </TEMPLATE>
  !   <IN NAME="filename" TYPE= "character(len=*)">
  !    The name of the grid file to be generated.
  !   </IN>

  function get_file_unit(filename)
    character(len=*), intent(in) :: filename
    integer                      :: get_file_unit

    if(trim(opened_files(0)) == trim(filename) ) then
          get_file_unit = files_unit(0)
          return
    endif

    !--- if file is not opened, open the file
    ! REPLACE WITH REGULAR METHOD TO OPEN NETCDF FILES
    call mpp_open(get_file_unit, trim(filename), MPP_OVERWR,MPP_NETCDF, &
          threading=MPP_SINGLE,  fileset=MPP_SINGLE )

   files_unit(num_files) = get_file_unit
   opened_files(num_files) = trim(filename)

  end function get_file_unit
  ! </FUNCTION>

  !#######################################################################
  ! <FUNCTION NAME="field_exist">

  !   <OVERVIEW>
  !     check if a given field name exists in a given file name.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     check if a given field name exists in a given file name.
  !     If the field_name string has zero length or the
  !     first character is blank return a false result.
  !     if the file file_name don't exist, return a false result.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     field_exist ( file_name, field_name )
  !   </TEMPLATE>

  !   <IN NAME="file_name"  TYPE="character" >
  !     A file name (or path name) that is checked for existence.
  !   </IN>
  !   <IN NAME="field_name"  TYPE="character" >
  !     A field name that is checked for existence.
  !   </IN>
  !   <OUT NAME=""  TYPE="logical" >
  !     This function returns a logical result.  If field exists in the
  !     file file_name, the result is true, otherwise false is returned.
  !     If the length of character string "field_name" is zero or the first
  !     character is blank, then the returned value will be false.
  !     if the file file_name don't exist, return a false result.
  !   </OUT>

  function field_exist (file_name, field_name)
    character(len=*),                 intent(in) :: file_name
    character(len=*),                 intent(in) :: field_name
    logical                      :: field_exist
    integer                      :: unit, ndim, nvar, natt, ntime, i, nfile
    character(len=64)            :: name
    type(fieldtype), allocatable :: fields(:)
    logical                      :: file_exist
    character(len=256)           :: fname

    field_exist = .false.
    if (len_trim(field_name) == 0) return
    if (field_name(1:1) == ' ')    return

    file_exist=get_file_name(file_name, fname)
    if(file_exist) then
      call get_file_unit(fname, unit, nfile)
      allocate(fields(nvar))

      do i=1, nvar
        if(lowercase(trim(name)) == lowercase(trim(field_name))) field_exist = .true.
      enddo
      deallocate(fields)
    endif
    if(field_exist) return
    file_exist =  get_file_name(file_name, fname, read_dist)
    if(file_exist) then, io_domain_exist)
      allocate(fields(nvar))
      do i=1, nvar
        if(lowercase(trim(name)) == lowercase(trim(field_name))) field_exist = .true.
      enddo
      deallocate(fields)
    endif

    return

  end function field_exist
  ! </FUNCTION>

  ! ========================================================================================================================






  ! ========================================================================================================================
  ! FROM OCEAN_PARAMETERS_MOD
  !use ocean_parameters_mod,  only: ADVECT_UPWIND, ADVECT_QUICKER
  !use ocean_parameters_mod,  only: TWO_LEVEL, missing_value

  ! some numerical constants
  real, parameter, private  :: missing_value=-1.e20
  ! parameters for choosing the time tendency calculation
  integer, parameter, private :: TWO_LEVEL   = 2
  ! parameters for tracer advection
  integer, parameter, private :: ADVECT_UPWIND            = 1
  integer, parameter, private :: ADVECT_QUICKER           = 5


  ! ========================================================================================================================




  ! ========================================================================================================================
  ! FROM OCEAN_TYPES_MOD
  !use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type, ocean_density_type, ocean_time_type
  !use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_thickness_type, ocean_adv_vel_type

  type time_type
     private
     integer:: seconds
     integer:: days
     integer:: ticks
     integer:: dummy ! added as a workaround bug on IRIX64 (AP)
  end type time_type

  #ifdef MOM_STATIC_ARRAYS
  !########################
  type, private :: ocean_thickness_type
    integer :: method       ! energetic or finite volume
    real, dimension(isd:ied,jsd:jed,nk,3) :: rho_dzt   ! rho(kg/m^3)*thickness (m) of T cell at 3 times
    real, dimension(isd:ied,jsd:jed,nk)   :: rho_dztr  ! 1.0/(rho*dzt) at time taup1
    real, dimension(isd:ied,jsd:jed,nk)   :: rho_dzt_tendency ! rho_dzt tendency (kg/m^3)*(m/s)
    real, dimension(isd:ied,jsd:jed)      :: sea_lev     ! eta_t + patm/(rho0*grav) - eta_geoid - eta_tide (m) at time taup1 for coupler
    real, dimension(isd:ied,jsd:jed,nk)   :: dzt         ! thickness (m) of T cell at time tau/taup1
    real, dimension(isd:ied,jsd:jed,nk)   :: dzu         ! thickness (m) of U cell at time tau/taup1
    real, dimension(isd:ied,jsd:jed,0:nk) :: dzwu        ! vertical distance (m) between U points at tau/taup1
    real, dimension(isd:ied,jsd:jed,nk)   :: depth_zwt   ! vertical distance (m) from column top to T-bottom
    real, dimension(isd:ied,jsd:jed,nk)   :: depth_st    ! s-distance to T-point
    real, dimension(isd:ied,jsd:jed,nk)   :: depth_swt   ! s-distance to T-bottom
    real, dimension(isd:ied,jsd:jed,nk)   :: dst         ! vertical increment of s-coordinate on T-cell
    real, dimension(isd:ied,jsd:jed,nk)   :: dzt_dst     ! specific thickness (metre/s-coordinate) on T-cell
    ! The following arrays are never allocated with MOM_STATIC_ARRAYS,
    ! but they need to be declared allocatable in order to compile.
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztT     _NULL ! rho(kg/m^3)*thickness (m) of T cell at 3 times (total)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzuL     _NULL ! L system contribution to rho_dzuT (3 time levels)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzuT     _NULL ! rho (kg/m^3) * thickness (m) of U cell at 3 times (total)
    real, dimension(:,:,:),   _ALLOCATABLE :: dztL         _NULL ! L system contribution to dztT
    real, dimension(:,:,:,:), _ALLOCATABLE :: dztT         _NULL ! thickness (m) of T cell at time tau/taup1
    real, dimension(:,:,:),   _ALLOCATABLE :: dzuL         _NULL ! L system contribution to dzuT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzuT         _NULL ! thickness (m) of U cell at time tau/taup1
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwtL        _NULL ! L system contribution to dzwtT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwtT        _NULL ! vertical distance (m) between T points at tau/taup1
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwuL        _NULL ! L system contribution to dzwuT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwuT        _NULL ! vertical distance (m) between U points at tau/taup1
    real, dimension(:,:,:),   _ALLOCATABLE :: geodepth_zwt _NULL ! vert distance (m) from z=0 to bottom of T-cell
  end type ocean_thickness_type

  type, private :: ocean_grid_type
    character(len=32) :: name
    ! geometry and topology and rotation
    logical                           :: cyclic_x         ! true if domain is cyclic in the i direction
    logical                           :: cyclic_y         ! true if domain is cyclic in the j direction
    logical                           :: tripolar         ! folded connectivity at "top" row in bipolar Arctic
    logical                           :: mosaic           ! true when using a mosaic grid
    logical                           :: beta_plane       ! beta plane Cartesian
    logical                           :: f_plane          ! f-plane Cartesian
    real                              :: f_plane_latitude ! latitude where f_plane is centered
    real, dimension(isd:ied,jsd:jed)  :: f                ! coriolis parameter at u-cell points (sec^-1)
    real, dimension(isd:ied,jsd:jed)  :: fstar            ! horizontal coriolis parameter at u-cell points (sec^-1)
    real, dimension(isd:ied,jsd:jed)  :: beta             ! df/dy at u-cell points (1/(sec*m))
    real, dimension(isd:ied,jsd:jed)  :: beta_eff         ! beta plus topographic beta at u-cell (1/(sec*m))
    ! vertical grid information (time independent)
    integer                             :: nk     ! number of vertical grid points
    integer, dimension(isd:ied,jsd:jed) :: kmt    ! number of t-levels
    integer, dimension(isd:ied,jsd:jed) :: kmu    ! number of u-levels
    real, dimension(nk)                 :: zt     ! full cell depth (m) from surface to level k T-cell
    real, dimension(nk)                 :: zw     ! full cell depth (m) from surface to bottom of k T-cell
    real, dimension(nk)                 :: dzt    ! initial vertical resolution of T or U grid cells (m)
    real, dimension(0:nk)               :: dzw    ! initial vertical resolution of W grid cells (m)
    real, dimension(0:nk)               :: dzwr   ! reciprocal of dzw (W cell vertical resolution)
    real, dimension(nk)                 :: st     ! full cell s-depth from surface to level k T-cell
    real, dimension(nk)                 :: sw     ! full cell s-depth from surface to bottom of k T-cell
    real, dimension(nk)                 :: dst    ! initial vertical s-resolution of T or U grid cells (m)
    real, dimension(0:nk)               :: dsw    ! initial vertical s-resolution of W grid cells (m)
    real, dimension(nk,0:1)             :: fracdz ! fractional distance between grid point and cell top/bot
    real, dimension(isd:ied,jsd:jed)    :: ht     ! depth to bottom of ocean (m) on t-cells from z=0
    real, dimension(isd:ied,jsd:jed)    :: htr    ! inverse depth to bottom of ocean (m^-1) on t-cells
    real, dimension(isd:ied,jsd:jed)    :: hu     ! depth to bottom of ocean (m) on u-cells from z=0
    real, dimension(isd:ied,jsd:jed)    :: dht_dx ! d(ht)/dx on u-cells (m/m)
    real, dimension(isd:ied,jsd:jed)    :: dht_dy ! d(ht)/dy on u-cells (m/m)
    real, dimension(isd:ied,jsd:jed)    :: gradH  ! sqrt(dht_dx**2+dht_dyx**2) on u-cells (m/m)
    ! horizontal grid information (time independent)
    integer                          :: ni, nj     ! number of global points in the two horizontal directions
    real, dimension(isd:ied,jsd:jed) :: xt         ! longitude of the T grid points in degrees.
    real, dimension(isd:ied,jsd:jed) :: xu         ! longitude of the U grid points in degrees.
    real, dimension(isd:ied,jsd:jed) :: yt         ! latitude of the T grid points in degrees.
    real, dimension(isd:ied,jsd:jed) :: yu         ! latitude of the U grid points in degrees.
    real, dimension(isd:ied,jsd:jed) :: dxt        ! longitudinal width of T-cells at grid point (m)
    real, dimension(isd:ied,jsd:jed) :: dxu        ! longitudinal width of U-cells at grid point (m)
    real, dimension(isd:ied,jsd:jed) :: dyt        ! latitudinal width of T-cells at grid point (m)
    real, dimension(isd:ied,jsd:jed) :: dyu        ! latitudinal width of U-cells at grid point (m)
    real, dimension(isd:ied,jsd:jed) :: dau        ! area of U-cells (m^2)
    real, dimension(isd:ied,jsd:jed) :: dat        ! area of T-cells (m^2)
    real, dimension(isd:ied,jsd:jed) :: dxtn       ! long-width of north face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dyte       ! lat-width of east face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dts        ! width from grid point to south face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dxtr       ! 1/dxt
    real, dimension(isd:ied,jsd:jed) :: dxur       ! 1/dxu
    real, dimension(isd:ied,jsd:jed) :: dytr       ! 1/dyt
    real, dimension(isd:ied,jsd:jed) :: dyur       ! 1/dyu
    real, dimension(isd:ied,jsd:jed) :: datr       ! 1/[area of T-cells (m^2)]
    real, dimension(isd:ied,jsd:jed) :: dyter      ! 1/dyte
    real, dimension(isd:ied,jsd:jed) :: dxtnr      ! 1/dxtn
    real, dimension(isd:ied,jsd:jed) :: dxt_dxtnr  ! dxt/dxtn
    real, dimension(isd:ied,jsd:jed) :: dxtn_dxtr  ! dxtn/dxt
    real, dimension(isd:ied,jsd:jed) :: dyt_dyter  ! dyt/dyte
    real, dimension(isd:ied,jsd:jed) :: dyte_dytr  ! dyte/dyt
    ! land/sea masks
    real, dimension(isd:ied,jsd:jed,nk)   :: tmask       ! land/sea mask for T cells based on s-coordinate
  end type ocean_grid_type

  type, private :: ocean_time_type
    type(time_type) :: model_time     ! fms variable
    type(time_type) :: Time_step      ! time step for tracers (and for ocean model)
    type(time_type) :: Time_init      ! time of initial conditions
    integer :: calendar               ! calendar type defined by time_manager_mod
    logical :: init                   ! true at beginning of run (initial condition time)
    integer :: itt                    ! timestep counter measured relative to time at restart
    integer :: itt0                   ! timestep counter measured relative to initial condition time
    integer :: taum1, tau, taup1      ! time level indices
    integer :: tau_m2, tau_m1, tau_m0 ! time level indices for Adams-Bashforth velocity advection and coriolis
  end type ocean_time_type

  type, private :: ocean_adv_vel_type
    real, dimension(isd:ied,jsd:jed,nk)   :: uhrho_et   ! rho_dzu * advect vel (kg/(m*s)) on i-face of T-cell
    real, dimension(isd:ied,jsd:jed,nk)   :: vhrho_nt   ! rho*dzu * advect vel (kg/(m*s)) on j-face of T-cell
    real, dimension(isd:ied,jsd:jed,0:nk) :: wrho_bt    ! rho * vertical advect vel (kg/(m^2*s)) on T-bottom
  end type ocean_adv_vel_type

  type, private :: ocean_prog_tracer_type
    character(len=32)  :: name
    character(len=32)  :: units
    character(len=32)  :: type
    character(len=128) :: longname
    logical :: use_only_advection     ! for testing purposes, evolve using ONLY advection
    logical :: neutral_physics_limit  ! neutral physics reduce to horz diffusion if tracer out of bounds
    logical :: complete               ! to determine if ready to do mpp updates
    integer :: sfc_flux_id=-1         ! index for time_interp_external
    integer :: horz_advect_scheme=1  ! id for horizontal advection scheme
    integer :: vert_advect_scheme=1  ! id for vertical advection scheme
    integer :: ppm_hlimiter=1          ! Limiter for use with PPM in horizontal
    integer :: ppm_vlimiter=1          ! Limiter for use with PPM in vertical
    integer :: mdt_scheme=4            ! Version of Multi-Dim. Modified Daru & Tenaud (MDMDT)
    type(obc_flux), _ALLOCATABLE, dimension(:) :: otf   _NULL ! flux through open boundaries, allocate nobc
    real, dimension(isd:ied,jsd:jed,nk,3) :: field            ! tracer concentration at 3 time levels
    real, dimension(isd:ied,jsd:jed,nk)   :: th_tendency      ! thickness weighted tracer tendency
    real, dimension(isd:ied,jsd:jed,nk)   :: tendency         ! for diagnostics: tendency concentration [concentration/sec]
    real, dimension(isd:ied,jsd:jed,nk)   :: wrk1             ! work array
    real, dimension(isd:ied,jsd:jed,nk)   :: K33_implicit     ! m^2/sec vert-diffusivity from neutral diffusion
    real                                  :: conversion         ! conversion of dimensions
    real                                  :: offset             ! offset in dimensions (e.g., Celsius to Kelvin)
    real                                  :: min_tracer         ! min acceptable value--model stopped if less
    real                                  :: max_tracer         ! max acceptable value--model stopped if greater
    real                                  :: min_range          ! min value used for calls to diagnostic manager
    real                                  :: max_range          ! max value used for calls to diagnostic manager
    real                                  :: min_tracer_limit   ! min value used to limit quicker and neutral fluxes
    real                                  :: max_tracer_limit   ! max value used to limit quicker and neutral fluxes
    real                                  :: min_flux_range     ! min and max values used for flux diagnostics
    real                                  :: max_flux_range     ! min and max values used for flux diagnostics
    real                                  :: const_init_value   ! value used for constant tracer init
    logical                               :: const_init_tracer  ! false (default) if the tracer must exist in the restart file
    !   otherwise will initialize with const_init_value
    character(len=32)                     :: flux_units         ! units for the tracer flux
  end type ocean_prog_tracer_type

  #else
  !#########################
  ! not MOM_STATIC_ARRAYS
  type, private :: ocean_thickness_type
    integer :: method       ! energetic or finite volume
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzt   _NULL ! E system contribution to rho_dztT (3 time levels)
    real, dimension(:,:,:),   _ALLOCATABLE :: rho_dztr  _NULL ! 1.0/(rho*dzt) at time taup1 (E system)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztT  _NULL ! rho(kg/m^3)*thickness (m) of T cell at 3 times (total)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzu   _NULL ! E system contribution to rho_dzuT (3 time levels)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzuL  _NULL ! L system contribution to rho_dzuT (3 time levels)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzuT  _NULL ! rho (kg/m^3) * thickness (m) of U cell at 3 times (total)
    real, dimension(:,:,:),   _ALLOCATABLE :: rho_dzt_tendency _NULL ! rho_dzt tendency (kg/m^3)*(m/s)
    real, dimension(:,:),     _ALLOCATABLE :: sea_lev _NULL ! eta_t + patm/(rho0*grav) - eta_geoid - eta_tide (m) at time taup1 for coupler
    real, dimension(:,:,:),   _ALLOCATABLE :: dzt    _NULL ! E system contribution to dztT
    real, dimension(:,:,:,:), _ALLOCATABLE :: dzten  _NULL ! E system contribution to dzt at east/north face of T-cell
    real, dimension(:,:,:),   _ALLOCATABLE :: dztL   _NULL ! L system contribution to dztT
    real, dimension(:,:,:,:), _ALLOCATABLE :: dztT   _NULL ! thickness (m) of T cell at time tau/taup1
    real, dimension(:,:,:),   _ALLOCATABLE :: dzu    _NULL ! E system contribution to dzuT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzuL   _NULL ! L system contribution to dzuT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzuT   _NULL ! thickness (m) of U cell at time tau/taup1
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwtL  _NULL ! L system contribution to dzwtT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwtT  _NULL ! vertical distance (m) between T points at tau/taup1
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwu   _NULL ! E system contribution to dzwuT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwuL  _NULL ! L system contribution to dzwuT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwuT  _NULL ! vertical distance (m) between U points at tau/taup1
    real, dimension(:,:,:),   _ALLOCATABLE :: geodepth_zt  _NULL ! vert distance (m) from z=0 to T-point
    real, dimension(:,:,:),   _ALLOCATABLE :: geodepth_zwt _NULL ! vert distance (m) from z=0 to bottom of T-cell
    real, dimension(:,:,:),   _ALLOCATABLE :: depth_zt     _NULL ! vert distance (m) from column top to T-point
    real, dimension(:,:,:),   _ALLOCATABLE :: depth_zwt    _NULL ! vert distance (m) from column top to T-bottom
    real, dimension(:,:,:),   _ALLOCATABLE :: depth_zu     _NULL ! vert distance (m) from column top to U-point
    real, dimension(:,:,:),   _ALLOCATABLE :: depth_zwu    _NULL ! vert distance (m) from column top to U-bottom
    real, dimension(:,:,:),   _ALLOCATABLE :: depth_st    _NULL ! s-distance to T cell grid point
    real, dimension(:,:,:),   _ALLOCATABLE :: depth_swt   _NULL ! s-distance to T cell grid bottom
    real, dimension(:,:,:),   _ALLOCATABLE :: dst         _NULL ! s-increment of T-cell at time tau
    real, dimension(:,:,:),   _ALLOCATABLE :: dswt        _NULL ! s-increment of T-point at time tau
    real, dimension(:,:,:),   _ALLOCATABLE :: dzt_dst     _NULL ! T-cell specific thickness (m/s-coordinate)
  end type ocean_thickness_type


  type, private :: ocean_grid_type
    character(len=32) :: name
    ! geometry and topology and rotation
    logical                            :: cyclic_x          ! true if domain is cyclic in the i direction
    logical                            :: cyclic_y          ! true if domain is cyclic in the j direction
    logical                            :: tripolar          ! folded connectivity at "top" row w/i bipolar Arctic
    logical                            :: mosaic            ! true when using a mosaic grid
    logical                            :: beta_plane        ! beta plane Cartesian
    logical                            :: f_plane           ! f-plane Cartesian
    real                               :: f_plane_latitude  ! latitude where f_plane is centered
    real, dimension(:,:), _ALLOCATABLE :: f        _NULL ! coriolis parameter at u-cell points (sec^-1)
    real, dimension(:,:), _ALLOCATABLE :: fstar    _NULL ! horizontal coriolis parameter at u-cell points (sec^-1)
    real, dimension(:,:), _ALLOCATABLE :: beta     _NULL ! df/dy at u-cell points (1/(sec*m))
    real, dimension(:,:), _ALLOCATABLE :: beta_eff _NULL ! df/dy plus topographic beta at u-cell (1/(sec*m))
    ! vertical grid information (time independent)
    integer                               :: nk           ! number of vertical grid points
    integer, dimension(:,:), _ALLOCATABLE :: kmt    _NULL ! number of t-levels
    integer, dimension(:,:), _ALLOCATABLE :: kmu    _NULL ! number of u-levels
    real,    dimension(:),   _ALLOCATABLE :: zt     _NULL ! distance from surface to grid point in level k (m)
    real,    dimension(:),   _ALLOCATABLE :: zw     _NULL ! distance from surface down to bottom of level k (m)
    real,    dimension(:),   _ALLOCATABLE :: dzt    _NULL ! initial vertical resolution of T or U grid cells (m)
    real,    dimension(:),   _ALLOCATABLE :: dzw    _NULL ! initial vertical resolution of W grid cells (m)
    real,    dimension(:),   _ALLOCATABLE :: dzwr   _NULL ! reciprocal of dzw (W cell vertical resolution)
    real,    dimension(:),   _ALLOCATABLE :: st     _NULL ! s-distance from surface to grid point in level k
    real,    dimension(:),   _ALLOCATABLE :: sw     _NULL ! s-distance from surface down to bottom of level k
    real,    dimension(:),   _ALLOCATABLE :: dst    _NULL ! initial s-vertical resolution of T or U grid cells
    real,    dimension(:),   _ALLOCATABLE :: dsw    _NULL ! initial s-vertical resolution of W grid cells
    real,    dimension(:,:), _ALLOCATABLE :: ht     _NULL ! depth to bottom of ocean (m) on t-cells
    ! horizontal grid information (time independent)
    integer                            :: ni, nj           ! global points in the two horizontal directions
    real, dimension(:,:), _ALLOCATABLE :: xt         _NULL ! longitude of the T grid points in degrees
    real, dimension(:,:), _ALLOCATABLE :: xu         _NULL ! longitude of the U grid points in degrees
    real, dimension(:,:), _ALLOCATABLE :: yt         _NULL ! latitude of the T grid points in degrees
    real, dimension(:,:), _ALLOCATABLE :: yu         _NULL ! latitude of the U grid points in degrees
    real, dimension(:,:), _ALLOCATABLE :: dxt        _NULL ! longitudinal width of T-cells at grid point (m)
    real, dimension(:,:), _ALLOCATABLE :: dyt        _NULL ! latitudinal width of T-cells at grid point (m)
    real, dimension(:,:), _ALLOCATABLE :: dat        _NULL ! area of T-cells (m^2)
    real, dimension(:,:), _ALLOCATABLE :: dxtn       _NULL ! i-width of north face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dyte       _NULL ! j-width of east  face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dtn        _NULL ! width from grid point to north face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dts        _NULL ! width from grid point to south face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dtw        _NULL ! width from grid point to west  face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dte        _NULL ! width from grid point to east  face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dxtr       _NULL ! 1/dxt
    real, dimension(:,:), _ALLOCATABLE :: dytr       _NULL ! 1/dyt
    real, dimension(:,:), _ALLOCATABLE :: datr       _NULL ! 1/[area of T-cells (m^2)]
    real, dimension(:,:), _ALLOCATABLE :: dyter      _NULL ! 1/dyte
    real, dimension(:,:), _ALLOCATABLE :: dxtnr      _NULL ! 1/dxtn
    real, dimension(:,:), _ALLOCATABLE :: dxt_dxtnr  _NULL ! dxt/dxtn
    real, dimension(:,:), _ALLOCATABLE :: dxtn_dxtr  _NULL ! dxtn/dxt
    real, dimension(:,:), _ALLOCATABLE :: dyt_dyter  _NULL ! dyt/dyte
    real, dimension(:,:), _ALLOCATABLE :: dyte_dytr  _NULL ! dyte/dyt
    ! land/sea masks
    real, dimension(:,:,:),   _ALLOCATABLE :: tmask       _NULL ! land/sea mask for T cells based on s-coordinate
  end type ocean_grid_type

  type, private :: ocean_time_type
    type(time_type) :: model_time     ! fms variable
    type(time_type) :: time_step      ! ocean tracer timestep
    type(time_type) :: Time_init      ! time of initial conditions
    integer :: calendar               ! calendar type defined by time_manager_mod
    logical :: init                   ! true at beginning of run (initial condition time)
    integer :: itt                    ! timestep counter measured relative to time at restart
    integer :: itt0                   ! timestep counter measured relative to initial condition time
    integer :: taum1, tau, taup1      ! time level indices
    integer :: tau_m2, tau_m1, tau_m0 ! time level indices for Adams-Bashforth velocity advection and coriolis
  end type ocean_time_type

  type, private :: ocean_adv_vel_type
    real, _ALLOCATABLE, dimension(:,:,:)  :: uhrho_et   _NULL ! rho_dzu weight advect vel (kg/(m*s)) on T-cell i-face
    real, _ALLOCATABLE, dimension(:,:,:)  :: vhrho_nt   _NULL ! rho_dzu weight advect vel (kg/(m*s)) on T-cell j-face
    real, _ALLOCATABLE, dimension(:,:,:)  :: wrho_bt    _NULL ! rho weight (kg/(m^2*s)) vert advect vel on T-bottom
  end type ocean_adv_vel_type

  type, private :: ocean_prog_tracer_type
    character(len=32)  :: name
    character(len=32)  :: units
    character(len=32)  :: type
    character(len=128) :: longname
    logical :: use_only_advection      ! for testing purposes, evolve using ONLY advection
    logical :: neutral_physics_limit   ! revert neutral physics to horz diffusion where tracer out of bounds
    logical :: complete                ! to determine if ready to do mpp updates
    integer :: sfc_flux_id=-1          ! index for time_interp_external
    integer :: horz_advect_scheme=1   ! id for horizontal advection scheme
    integer :: vert_advect_scheme=1   ! id for vertical advection scheme
    integer :: id_obc                  ! id to identify tracer in OBC-subroutines
    integer :: ppm_hlimiter=1          ! Limiter for use with PPM in horizontal
    integer :: ppm_vlimiter=1          ! Limiter for use with PPM in vertical
    integer :: mdt_scheme=4            ! Version of Multi-Dim. Modified Daru & Tenaud (MDMDT)
    real, _ALLOCATABLE, dimension(:,:,:,:) :: field          _NULL ! tracer concentration at 3 time levels (E system)
    real, _ALLOCATABLE, dimension(:,:,:)   :: th_tendency    _NULL ! thickness weighted tracer tendency
    real, _ALLOCATABLE, dimension(:,:,:)   :: tendency       _NULL ! for diagnostics: tendency concentration [concentration/sec]
    real, _ALLOCATABLE, dimension(:,:,:)   :: wrk1                _NULL ! work array
    real, _ALLOCATABLE, dimension(:,:,:)   :: K33_implicit        _NULL ! m^2/sec vert-diffusivity from neutral diffusion
    real, _ALLOCATABLE, dimension(:,:)     :: flux_int            _NULL ! integrated sfc tracer flux for diagnostics
    real                              :: conversion            ! conversion between  dimensions
    real                              :: offset                ! offset in dimensions (e.g., Celsius to Kelvin)
    real                              :: min_tracer            ! min acceptable value used for error checking
    real                              :: max_tracer            ! max acceptable value used for error checking
    real                              :: min_range             ! min value used for calls to diagnostic manager
    real                              :: max_range             ! max value used for calls to diagnostic manager
    real                              :: min_tracer_limit      ! min value used to limit quicker & neutral fluxes
    real                              :: max_tracer_limit      ! max value used to limit quicker & neutral fluxes
    real                              :: min_flux_range        ! min and max values used for flux diagnostics
    real                              :: max_flux_range        ! min and max values used for flux diagnostics
    real                              :: const_init_value      ! value used to initialize constant tracer
    logical                           :: const_init_tracer     ! false (default) if the tracer must exist in the restart file
    !   otherwise will initialize with const_init_value
    character(len=32)                 :: flux_units            ! units for the tracer flux
  end type ocean_prog_tracer_type

  #endif
  !############################################################################################
  ! end of STATIC_MEMORY

  ! ========================================================================================================================






  ! =========================================================================================================================
  ! FROM OCEAN_TRACER_ADVECT_MOD
  implicit none

  private
  integer :: now
  integer :: id_clock_up_horz
  integer :: id_clock_up_vert
  integer :: id_clock_quick_horz
  integer :: id_clock_quick_vert
  integer :: id_clock_adv_diss

  #include <ocean_memory.h>

  #ifdef MOM_STATIC_ARRAYS

  real, dimension(isd:ied,jsd:jed,nk) :: flux_x
  real, dimension(isd:ied,jsd:jed,nk) :: flux_y
  real, dimension(isd:ied,jsd:jed,nk) :: flux_z

  real, dimension(isd:ied,jsd:jed,nk) :: advect_tendency

  #else

  real, dimension(:,:,:), allocatable :: flux_x
  real, dimension(:,:,:), allocatable :: flux_y
  real, dimension(:,:,:), allocatable :: flux_z

  real, dimension(:,:,:), allocatable :: advect_tendency

  #endif

  !work array on neutral density space
  integer :: neutralrho_nk
  real, dimension(:,:,:),   allocatable :: nrho_work

  integer :: num_prog_tracers = 0

  character(len=256) :: version='CVS $Id: ocean_tracer_advect.F90,v 1.1.2.4.18.1 2013/03/06 18:19:34 Niki.Zadeh Exp $'
  character(len=256) :: tagname='Tag $Name: siena_201303 $'


  type(ocean_grid_type)  , pointer :: Grd =>NULL()

  real, parameter :: a4=7.0/12.0, b4=-1.0/12.0

  logical :: used

  integer, dimension(:), allocatable :: id_tracer_advection
  integer, dimension(:), allocatable :: id_tracer2_advection
  integer, dimension(:), allocatable :: id_tracer_adv_diss

  integer, dimension(:), allocatable :: id_vert_advect
  integer, dimension(:), allocatable :: id_horz_advect
  integer, dimension(:), allocatable :: id_advection_x
  integer, dimension(:), allocatable :: id_advection_y
  integer, dimension(:), allocatable :: id_advection_z
  integer, dimension(:), allocatable :: id_xflux_adv
  integer, dimension(:), allocatable :: id_yflux_adv
  integer, dimension(:), allocatable :: id_zflux_adv
  integer, dimension(:), allocatable :: id_xflux_adv_int_z
  integer, dimension(:), allocatable :: id_yflux_adv_int_z

  private horz_advect_tracer
  private horz_advect_tracer_upwind

  private vert_advect_tracer
  private vert_advect_tracer_upwind

  private advection_diag_init

  private compute_adv_diss

  private ocean_tracer_advect_init
  private ocean_tracer_advect_end
  private ocean_tracer_advect_restart

  logical :: module_is_initialized   = .false.
  logical :: have_obc                = .false.

  logical :: limit_with_upwind       = .false.
  logical :: debug_this_module       = .false.
  logical :: advect_sweby_all        = .false.
  logical :: zero_tracer_advect_horz = .false.
  logical :: zero_tracer_advect_vert = .false.
  logical :: write_a_restart         = .false.
  logical :: async_domain_update     = .false.

  namelist /ocean_tracer_advect_nml/ debug_this_module, limit_with_upwind, advect_sweby_all,   &
  zero_tracer_advect_horz, zero_tracer_advect_vert,         &
  write_a_restart, read_basin_mask,     &
  async_domain_update

contains


  !#######################################################################
  ! <SUBROUTINE NAME="ocean_tracer_advect_init">
  !
  ! <DESCRIPTION>
  ! Initialize the tracer advection module.
  ! </DESCRIPTION>
  !
  subroutine ocean_tracer_advect_init (Grid, Time, T_prog, dtime)

    type(ocean_grid_type),        intent(in), target   :: Grid
    type(ocean_time_type),        intent(in)           :: Time
    type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
    real,                         intent(in)           :: dtime

    integer :: n
    integer :: ioun, io_status, ierr

    integer :: stdoutunit,stdlogunit
    stdoutunit=stdout();stdlogunit=stdlog()

    write( stdlogunit,'(/a/)') trim(version)

    call write_version_number( version, tagname )


    #ifndef MOM_STATIC_ARRAYS
    nk  = Grid%nk
    allocate(flux_x(isd:ied,jsd:jed,nk))
    allocate(flux_y(isd:ied,jsd:jed,nk))
    allocate(flux_z(isd:ied,jsd:jed,nk))
    allocate(advect_tendency(isd:ied,jsd:jed,nk))
    #endif

    Grd => Grid

    flux_x = 0.0
    flux_y = 0.0
    flux_z = 0.0

    advect_tendency(:,:,:)     = 0.0

    num_prog_tracers = size(T_prog(:))


    if(.not. write_a_restart) then
      write(stdoutunit,'(a)') '==>Note: running ocean_tracer_advect with write_a_restart=.false.'
      write(stdoutunit,'(a)') '   Will NOT write restart file, so cannot restart if using PSOM advection.'
    endif

    if(zero_tracer_advect_horz) then
      write(stdoutunit,*)'==>WARNING: have turned off horizontal tracer advection. Unrealistic simulation.'
    endif
    if(zero_tracer_advect_vert) then
      write(stdoutunit,*)'==>WARNING: have turned off vertical tracer advection. Unrealistic simulation.'
    endif

    if(.not. advect_sweby_all) then
      write(stdoutunit,'(a)') ' '
      write(stdoutunit,'(a)') ' From ocean_tracer_advect_init: SUMMARY OF TRACER ADVECTION SCHEMES'

      do n=1,num_prog_tracers

        if(T_prog(n)%horz_advect_scheme==ADVECT_UPWIND) then
          write(stdoutunit,'(1x,a)') &
          trim(T_prog(n)%name) //'is using first order upwind for horz advection.'
        endif
        if(T_prog(n)%vert_advect_scheme==ADVECT_UPWIND) then
          write(stdoutunit,'(1x,a)') &
          trim(T_prog(n)%name) //'is using first order upwind for vert advection.'
        endif

        if(T_prog(n)%horz_advect_scheme==ADVECT_QUICKER) then
          write(stdoutunit,'(1x,a)') &
          trim(T_prog(n)%name) //' is using Quicker for horz advection.'
          if(limit_with_upwind) then
            write(stdoutunit,'(1x,a)')'limit_with_upwind reverts quicker to upwind if tracer outside limits'
          endif
        endif
        if(T_prog(n)%vert_advect_scheme==ADVECT_QUICKER) then
          write(stdoutunit,'(1x,a)') &
          trim(T_prog(n)%name) //' is using Quicker for vert advection.'
          if(limit_with_upwind) then
            write(stdoutunit,'(1x,a)')'limit_with_upwind reverts quicker to upwind if tracer outside limits'
          endif
        endif

      enddo
      write(stdoutunit,'(a)') ' '

    endif !(.not. advect_sweby_all)

  end subroutine ocean_tracer_advect_init
  ! </SUBROUTINE> NAME="ocean_tracer_advect_init"


!#######################################################################
! VELOCITY READ IN FROM COBALT MAY NOT BE THE CORRECT VELOCITY

! <SUBROUTINE NAME="read_advect_velocity">
!
! <DESCRIPTION>
! For reading in the advection velocity components.  Assume that
! the advection velocity components read in from a file are in units
! of meter/sec and placed on the T-cell faces, as in a C-grid ocean model.
!
! This routine assumes that the read-in velocity components already
! have the proper masking.
!
! The main application of this routine is for developing idealized
! test cases for tracer advection.
! </DESCRIPTION>
!
subroutine read_advect_velocity(Time, Thickness, Adv_vel)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_adv_vel_type),   intent(inout) :: Adv_vel

  integer            :: k
  character(len=128) :: filename
if(horz_grid == MOM_BGRID) then

   do k=1,nk
      Adv_vel%uhrho_et(:,:,k) = rho0*ue(:,:,k)*BAY(Thickness%dzu(:,:,k))
      Adv_vel%vhrho_nt(:,:,k) = rho0*vn(:,:,k)*BAX(Thickness%dzu(:,:,k))
      Adv_vel%wrho_bt(:,:,k)  = rho0*wb(:,:,k)
      Adv_vel%uhrho_eu(:,:,k) = REMAP_ET_TO_EU(Adv_vel%uhrho_et(:,:,k))
      Adv_vel%vhrho_nu(:,:,k) = REMAP_NT_TO_NU(Adv_vel%vhrho_nt(:,:,k))
      Adv_vel%wrho_bu(:,:,k)  = REMAP_BT_TO_BU(Adv_vel%wrho_bt(:,:,k))
   enddo

else

   ! U-cell advective components are not used for Cgrid
   do k=1,nk
      Adv_vel%uhrho_et(:,:,k) = rho0*ue(:,:,k)*Thickness%dzten(:,:,k,1)
      Adv_vel%vhrho_nt(:,:,k) = rho0*vn(:,:,k)*Thickness%dzten(:,:,k,2)
      Adv_vel%wrho_bt(:,:,k)  = rho0*wb(:,:,k)
   enddo

endif
end subroutine read_advect_velocity
! </SUBROUTINE> NAME="read_advect_velocity"
!#######################################################################

  !#######################################################################
  ! <SUBROUTINE NAME="horz_advect_tracer">
  !
  ! <DESCRIPTION>
  ! Compute horizontal advection of tracers
  ! </DESCRIPTION>
  !
  subroutine horz_advect_tracer(Time, Adv_vel, Thickness, T_prog, Tracer, ntracer, dtime, store_flux)

    type(ocean_time_type),        intent(in)    :: Time
    type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
    type(ocean_thickness_type),   intent(in)    :: Thickness
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
    type(ocean_prog_tracer_type), intent(inout) :: Tracer
    integer,                      intent(in)    :: ntracer
    real,                         intent(in)    :: dtime
    logical,            optional, intent(in)    :: store_flux

    real,dimension(isd:ied,jsd:jed)             :: tmp_flux
    integer                                     :: i, j, k
    integer                                     :: taum1, tau
    logical                                     :: store

    if(zero_tracer_advect_horz) return

    store = .TRUE.
    if(present(store_flux)) store = store_flux

    taum1 = Time%taum1
    tau   = Time%tau

    if(.not. advect_sweby_all) then

      do k=1,nk
        do j=jsd,jed
          do i=isd,ied
            Tracer%wrk1(i,j,k) = 0.0
          enddo
        enddo
      enddo

      select case (Tracer%horz_advect_scheme)

      case (ADVECT_UPWIND)
        Tracer%wrk1(isc:iec,jsc:jec,:) =  &
        -horz_advect_tracer_upwind(Adv_vel, Tracer%field(:,:,:,taum1))


      case default
        write(stdoutunit,*)'==>Error from ocean_tracer_advect_mod (horz_advect_tracer): chose invalid horz advection scheme')
      end select


      do k=1,nk
        do j=jsc,jec
          do i=isc,iec
            Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
          enddo
        enddo
      enddo


      ! fill some diagnostic fields
      advect_tendency(:,:,:) = 0.0
      do k=1,nk
        do j=jsc,jec
          do i=isc,iec
            advect_tendency(i,j,k) = Tracer%wrk1(i,j,k)
          enddo
        enddo
      enddo

      if (id_xflux_adv_int_z(ntracer) > 0) then
        tmp_flux(:,:) = 0.0
        do k=1,nk
          tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) +  flux_x(isc:iec,jsc:jec,k)
        enddo
      endif
      if (id_yflux_adv_int_z(ntracer) > 0) then
        tmp_flux(:,:) = 0.0
        do k=1,nk
          tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) +  flux_y(isc:iec,jsc:jec,k)
        enddo
      endif

    endif   ! endif for (.not. advect_sweby_all)


  end subroutine horz_advect_tracer
  ! </SUBROUTINE> NAME="horz_advect_tracer"


  !#######################################################################
  ! <SUBROUTINE NAME="vert_advect_tracer">
  !
  ! <DESCRIPTION>
  ! Compute vertical advection of tracers for case
  ! with advect_sweby_all=.false.
  ! </DESCRIPTION>
  !
  subroutine vert_advect_tracer(Time, Adv_vel, Thickness, T_prog, Tracer, ntracer, dtime)

    type(ocean_time_type),        intent(in)    :: Time
    type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
    type(ocean_thickness_type),   intent(in)    :: Thickness
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
    type(ocean_prog_tracer_type), intent(inout) :: Tracer
    integer,                      intent(in)    :: ntracer
    real,                         intent(in)    :: dtime

    integer :: tau, taum1
    integer :: i,j,k
    real    :: temporary

    if(zero_tracer_advect_vert) return

    tau   = Time%tau
    taum1 = Time%taum1

    if(.not. advect_sweby_all) then

      do k=1,nk
        do j=jsd,jed
          do i=isd,ied
            Tracer%wrk1(i,j,k) = 0.0
          enddo
        enddo
      enddo

      select case (Tracer%vert_advect_scheme)

      case (ADVECT_UPWIND)
        Tracer%wrk1(isc:iec,jsc:jec,:) = &
        -vert_advect_tracer_upwind(Adv_vel, Tracer%field(:,:,:,taum1))

      case default
        write(stdoutunit,*)'==>Error from ocean_tracer_advect_mod (vert_advect_tracer): invalid advection scheme chosen')
      end select

      do k=1,nk
        do j=jsc,jec
          do i=isc,iec
            Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
          enddo
        enddo
      enddo

      ! for diagnostics, here and in compute_adv_diss
      do k=1,nk
        do j=jsc,jec
          do i=isc,iec
            advect_tendency(i,j,k) = advect_tendency(i,j,k) + Tracer%wrk1(i,j,k)
          enddo
        enddo
      enddo

    endif  ! endif for (.not. advect_sweby_all)


    if(id_tracer_adv_diss(ntracer) > 0) then
      call compute_adv_diss(Time, Adv_vel, Thickness, T_prog(:), Tracer, ntracer, dtime)
    endif


  end subroutine vert_advect_tracer
  ! </SUBROUTINE> NAME="vert_advect_tracer"


  !#######################################################################
  ! <FUNCTION NAME="horz_advect_tracer_upwind">
  !
  ! <DESCRIPTION>
  ! Compute horizontal advection of tracers from first order upwind.
  ! This scheme is positive definite but very diffusive.
  ! </DESCRIPTION>
  !
  function horz_advect_tracer_upwind(Adv_vel, Tracer_field)

    type(ocean_adv_vel_type),     intent(in) :: Adv_vel
    real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field

    real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_upwind

    real, dimension(isd:ied,jsd:jed) :: fe, fn
    real                             :: velocity, upos, uneg
    integer                          :: i, j, k

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

      do j=jsc,jec
        do i=isc,iec
          horz_advect_tracer_upwind(i,j,k) = &
          Grd%tmask(i,j,k)*(fe(i,j)-fe(i-1,j)+fn(i,j)-fn(i,j-1))*Grd%datr(i,j)
        enddo
      enddo

    enddo

  end function horz_advect_tracer_upwind
  ! </FUNCTION> NAME="horz_advect_tracer_upwind"


  !#######################################################################
  ! <FUNCTION NAME="vert_advect_tracer_upwind">
  !
  ! <DESCRIPTION>
  ! Compute vertical advection of tracers from first order upwind.
  ! This scheme is positive definite, but very diffusive.
  ! </DESCRIPTION>
  !
  function vert_advect_tracer_upwind(Adv_vel, Tracer_field)

    type(ocean_adv_vel_type),     intent(in) :: Adv_vel
    real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field

    real,dimension(isc:iec,jsc:jec,nk) :: vert_advect_tracer_upwind

    real,dimension(isc:iec,jsc:jec)    :: ft1, ft2
    real                               :: velocity, wpos, wneg
    integer                            :: k, i, j, kp1

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
  ! </FUNCTION> NAME="vert_advect_tracer_upwind"



  !#######################################################################
  ! <SUBROUTINE NAME="compute_adv_diss">
  !
  ! <DESCRIPTION>
  !
  ! Compute the dissipation due to advection truncation errors.
  ! This diagnostic requires computation of advection operator acting
  ! on the squared tracer concentration.
  !
  ! NOTE: This scheme isolates the dissipation from trucation errors
  ! in advection ONLY for the following vertical coordinates:
  ! 1/ geopotential: for all k-levels, except for k=1
  ! (due to undulating surface height)
  ! 2/ pressure: for all k-levels, except for k=kmt
  ! (due to undulating bottom pressure)
  !
  ! NOTE: For the Quicker advection scheme, we assume the preferred
  ! two_level time scheme is used here, so that taum1=tau.
  !
  ! NOTE: If PSOM is used for temp or salt, then we MUST also enable
  ! a new passive tracer in the field table, with name
  ! passive_temp_sq and passive_salt_sq.  This extra tracer is
  ! used for diagnostics alone, and is required due to the extra
  ! moment fields used for computing the PSOM tendency.
  ! If PSOM is used for another tracer besides temp or salt, then
  ! some extra code needs to be written inside ocean_passive.F90,
  ! emulating the work done for temp and salt.
  !
  ! </DESCRIPTION>
  !
  subroutine compute_adv_diss(Time, Adv_vel, Thickness, T_prog, Tracer, ntracer, dtime)

    type(ocean_time_type),        intent(in)    :: Time
    type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
    type(ocean_thickness_type),   intent(in)    :: Thickness
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
    type(ocean_prog_tracer_type), intent(inout) :: Tracer
    integer,                      intent(in)    :: ntracer
    real,                         intent(in)    :: dtime

    integer :: tau, taup1
    integer :: i,j,k
    logical :: use_psom=.false.
    real    :: dtimer, term1, term2

    tau    = Time%tau
    taup1  = Time%taup1
    dtimer = 1.0/dtime

    wrk1 = 0.0
    wrk2 = 0.0
    wrk3 = 0.0
    wrk4 = 0.0

    do k=1,nk
      do j=jsd,jed
        do i=isd,ied
          wrk1(i,j,k) = (Tracer%field(i,j,k,tau))**2
        enddo
      enddo
    enddo

    ! advection operator acting on squared tracer concentration
    select case (Tracer%horz_advect_scheme)

    case (ADVECT_UPWIND)
      wrk2(isc:iec,jsc:jec,:) =  &
      -horz_advect_tracer_upwind(Adv_vel, wrk1)

    end select

    select case (Tracer%vert_advect_scheme)

    case (ADVECT_UPWIND)
      wrk3(isc:iec,jsc:jec,:) = &
      -vert_advect_tracer_upwind(Adv_vel, wrk1)

    end select

    ! wrk1 is now the total tendency for the squared tracer concentration
    wrk1(:,:,:) = 0.0
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          wrk1(i,j,k) = wrk2(i,j,k) + wrk3(i,j,k)
        enddo
      enddo
    enddo

    ! compute the dissipation from advection truncation errors
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          term1 = advect_tendency(i,j,k) &
          *(2.0*Thickness%rho_dzt(i,j,k,tau)*Tracer%field(i,j,k,tau) + dtime*advect_tendency(i,j,k))
          term2 = -Thickness%rho_dzt(i,j,k,taup1)*wrk1(i,j,k)
          wrk4(i,j,k) =  -(Tracer%conversion**2)*dtimer*(term1+term2)
        enddo
      enddo
    enddo

    ! HOW DO YOU USE THIS INFO? CURRENTLY, WRK4 GETS SENT SOMEWHERE
    ! SHOULD IT GET STORED AS TENDANCY OR SOMETHING?

    used = send_data(id_tracer_adv_diss(ntracer), wrk4(:,:,:), &
    Time%model_time, rmask=Grd%tmask(:,:,:),            &
    is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if(id_tracer2_advection(ntracer) > 0) then
      used = send_data(id_tracer2_advection(ntracer), wrk1(:,:,:)*Tracer%conversion**2, &
      Time%model_time, rmask=Grd%tmask(:,:,:),                                     &
      is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif

  end subroutine compute_adv_diss
  ! </SUBROUTINE> NAME="compute_adv_diss"


  !#######################################################################
  ! <SUBROUTINE NAME="get_tracer_stats">
  !
  ! <DESCRIPTION>
  ! Compute the upper/lower values of a 3D field, returning values via arguments
  ! </DESCRIPTION>
  subroutine get_tracer_stats(A,Mask,tmin,tmax)
    real, dimension(isc:iec,jsc:jec,nk), intent(in) :: A
    real, dimension(isc:iec,jsc:jec,nk), intent(in) :: Mask
    real, intent(out) :: tmin,tmax
    integer :: i,j,k
    tmin=1.e30; tmax=-1.e30
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          if (Mask(i,j,k)>0.) then
            tmin=min(tmin,A(i,j,k))
            tmax=max(tmax,A(i,j,k))
          endif
        enddo
      enddo
    enddo
  end subroutine get_tracer_stats
  ! </SUBROUTINE> NAME="get_tracer_stats"


  !#######################################################################
  ! <SUBROUTINE NAME="tracer_stats">
  !
  ! <DESCRIPTION>
  ! Check the upper/lower values of a 3D field fall between speficied bounds
  ! reporting points that fall outside. Bring model down uncleanly if bounds are
  ! exceeded.
  !
  ! NOTE: This is a debugging tool and not for normal use.
  !
  ! </DESCRIPTION>
  subroutine tracer_stats(A,Mask,tmin0,tmax0,label)
    real, dimension(isc:iec,jsc:jec,nk), intent(in) :: A
    real, dimension(isc:iec,jsc:jec,nk), intent(in) :: Mask
    real, intent(in) :: tmin0,tmax0
    character(len=*), intent(in) :: label
    integer :: i,j,k
    real :: tmin,tmax
    tmin=1.e30; tmax=-1.e30

    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          if (tmask_mdppm(i,j,k)>0.) then
            tmin=min(tmin,A(i,j,k))
            tmax=max(tmax,A(i,j,k))
            if (j==1.and.A(i,j,k)<tmin0*0.99999999999999) write(0,*) 'under i,k=',i,k,A(i,j,k)
            if (j==1.and.A(i,j,k)>tmax0*1.00000000000001) write(0,*) 'over i,k=',i,k,A(i,j,k)
          endif
        enddo
      enddo
    enddo
    write(0,*) 't_stats: ',label,' min=',tmin,' max=',tmax
    if ((tmin<tmin0*0.99999999999999).or.(tmax>tmax0*1.00000000000001)) then
      write(0,*) 'tmin0=',tmin0,'tmax0=',tmax0
      write(0,*) 'tmin =',tmin, 'tmax =',tmax
      !   stop 'Overshoots!'
    endif
    if ((tmin<tmin0*0.9999999999999).or.(tmax>tmax0*1.0000000000001)) then
      stop 'Overshoots!'
    endif
  end subroutine
  ! </SUBROUTINE> NAME="tracer_stats"




  !#######################################################################
  ! <SUBROUTINE NAME="update_advection_only">
  !
  ! <DESCRIPTION>
  !
  ! Redo tracer updates for those that use only advection--nothing else.
  ! This method is useful for testing advection schemes.
  !
  ! T_prog(n)%use_only_advection==.true. ignores all boundary forcing
  ! and sources, so if T_prog(n)%stf or pme, rivers, sources
  ! are nonzero, tracer diagnostics will spuriously indicate
  ! non-conservation.
  !
  ! Assume for these tests that
  !  (1) vertical advection is done fully explictly in time
  !  (2) pme, rivers, stf, btf, and other sources are zero
  !  (3) do not use advect_sweby_all
  !
  ! </DESCRIPTION>
  !
  subroutine update_advection_only(Time, Adv_vel, Dens, Thickness, T_prog, n)

    type(ocean_time_type),        intent(in)    :: Time
    type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
    type(ocean_density_type),     intent(in)    :: Dens
    type(ocean_thickness_type),   intent(in)    :: Thickness
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
    integer,                      intent(in)    :: n
    integer                                     :: i, j, k
    integer                                     :: taup1, taum1

    taup1 = Time%taup1
    taum1 = Time%taum1

    T_prog(n)%th_tendency = 0.0

    call horz_advect_tracer(Time, Adv_vel, Thickness, Dens, &
         T_prog(1:num_prog_tracers), T_prog(n), n, dtime, store_flux=.FALSE.)
    call vert_advect_tracer(Time, Adv_vel, Dens, Thickness, &
         T_prog(1:num_prog_tracers), T_prog(n), n, dtime)

    do k=1,nk
       do j=jsc,jec
          do i=isc,iec
             T_prog(n)%field(i,j,k,taup1) =                                    &
                  (Thickness%rho_dzt(i,j,k,taum1)*T_prog(n)%field(i,j,k,taum1) &
                 + dtime*T_prog(n)%th_tendency(i,j,k)) * Thickness%rho_dztr(i,j,k)
          enddo
       enddo
    enddo

  end subroutine update_advection_only
  ! </SUBROUTINE> NAME="update_advection_only">




end module ocean_tracer_advect_mod
