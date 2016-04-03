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

  !#######################################################################

  ! <SUBROUTINE NAME="fms_init">

  !   <OVERVIEW>
  !     Initializes the FMS module and also calls the initialization routines for all
  !     modules in the MPP package. Will be called automatically if the user does
  !     not call it.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !      Initialization routine for the fms module. It also calls initialization routines
  !      for the mpp, mpp_domains, and mpp_io modules. Although this routine
  !      will be called automatically by other fms_mod routines, users should
  !      explicitly call fms_init. If this routine is called more than once it will
  !      return silently. There are no arguments.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call fms_init ( )
  !   </TEMPLATE>


  !   <ERROR MSG="invalid entry for namelist variable warning_level" STATUS="FATAL">
  !     The namelist variable warning_level must be either 'fatal' or 'warning'
  !     (case-insensitive).
  !   </ERROR>
  !   <ERROR MSG="invalid entry for namelist variable clock_grain" STATUS="FATAL">
  !     The namelist variable clock_grain must be one of the following values:
  !     'NONE', 'COMPONENT', 'SUBCOMPONENT', 'MODULE_DRIVER', 'MODULE', 'ROUTINE',
  !     'LOOP', or 'INFRA' (case-insensitive).
  !   </ERROR>

  ! initializes the fms module/package
  ! also calls mpp initialization routines and reads fms namelist

  subroutine fms_init (localcomm )
    integer, intent(in), optional :: localcomm
    integer :: unit, ierr, io

    if (module_is_initialized) return    ! return silently if already called
    module_is_initialized = .true.

    !---- set severity level for warnings ----
    select case( trim(lowercase(warning_level)) )
    case default
      call error_mesg ( 'fms_init',  &
      'invalid entry for namelist variable warning_level', FATAL )
    end select

    !--- write version info and namelist to logfile ---
    call write_version_number (version, tagname)

  end subroutine fms_init
  ! </SUBROUTINE>


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
    logical,                optional,  intent(in) :: scalar_or_1d
    logical,                optional,  intent(in) :: is_compressed

    character(len=256)            :: fname
    integer                       :: unit, siz_in(4)
    integer                       :: file_index  ! index of the opened file in array files
    integer                       :: tlev=1
    integer                       :: index_field ! position of the fieldname in the list of variables
    integer                       :: cxsize, cysize
    integer                       :: dxsize, dysize
    integer                       :: gxsize, gysize
    integer                       :: ishift, jshift
    logical                       :: is_scalar_or_1d = .false.
    logical                       :: is_compressed = .false.
    logical                       :: found_file

    ! read disttributed files is used when reading restart files that are NOT mppnccombined. In this
    ! case PE 0 will read file_res.nc.0000, PE 1 will read file_res.nc.0001 and so forth.
    !
    ! namelist to be used with read_dist_files: threading_read=multi,
    ! threading_write=multi, fileset_write=multi.

    ! Initialize files to default values
    if(.not.module_is_initialized) call error_mesg(FATAL,'fms_io(read_data_3d_new):  module not initialized')

    scalar_or_1d = .false.
    if(present(scalar_or_1d)) is_scalar_or_1d = scalar_or_1d

    compressed = .false.
    if(present(is_compressed)) compressed = is_compressed

    found_file = get_file_name(filename, fname, read_dist)
    if(.not.found_file) call error_mesg(FATAL, 'fms_io_mod(read_data_3d_new): file ' //trim(filename)// &
    '(with the consideration of tile number) and corresponding distributed file are not found')
    call get_file_unit(fname, unit, file_index, read_dist)

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
    character(len=*), intent(in) :: routine, message
    integer,          intent(in) :: level
    character(len=512)           :: text
    logical                      :: opened
    integer                      :: istat, errunit, outunit

    !  input:
    !      routine   name of the calling routine (character string)
    !      message   message written to output   (character string)
    !      level     set to NOTE, MESSAGE, or FATAL (integer)

    if (.not.module_is_initialized) call fms_init ( )

    if( .NOT.module_is_initialized )call ABORT()

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

    if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)

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

    error_state = level
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

    if (.not.module_is_initialized) call fms_init ( )

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
    call get_local_indices(isd, ied, jsd, jed, isc, iec, jsc, jec)
    call get_global_indices(isg, ieg, jsg, jeg)

    ni = Grid%ni;nj = Grid%nj; nk = Grid%nk

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

    allocate(Grid%phit(isd:ied,jsd:jed))
    allocate(Grid%phiu(isd:ied,jsd:jed))
    allocate(Grid%h1t(isd:ied,jsd:jed))
    allocate(Grid%h2t(isd:ied,jsd:jed))
    allocate(Grid%h1u(isd:ied,jsd:jed))
    allocate(Grid%h2u(isd:ied,jsd:jed))

    allocate (Grid%dxt(isd:ied,jsd:jed))
    allocate (Grid%dxu(isd:ied,jsd:jed))
    allocate (Grid%dyt(isd:ied,jsd:jed))
    allocate (Grid%dyu(isd:ied,jsd:jed))
    allocate (Grid%dat(isd:ied,jsd:jed))
    allocate (Grid%dat_frac(isd:ied,jsd:jed))
    allocate (Grid%dau(isd:ied,jsd:jed))

    allocate (Grid%dxtn(isd:ied,jsd:jed))
    allocate (Grid%dytn(isd:ied,jsd:jed))
    allocate (Grid%dxte(isd:ied,jsd:jed))
    allocate (Grid%dyte(isd:ied,jsd:jed))
    allocate (Grid%dxun(isd:ied,jsd:jed))
    allocate (Grid%dyun(isd:ied,jsd:jed))
    allocate (Grid%dxue(isd:ied,jsd:jed))
    allocate (Grid%dyue(isd:ied,jsd:jed))

    allocate (Grid%dxtr(isd:ied,jsd:jed))
    allocate (Grid%dxur(isd:ied,jsd:jed))
    allocate (Grid%dytr(isd:ied,jsd:jed))
    allocate (Grid%dyur(isd:ied,jsd:jed))
    allocate (Grid%dxuer(isd:ied,jsd:jed))
    allocate (Grid%dyuer(isd:ied,jsd:jed))
    allocate (Grid%dxunr(isd:ied,jsd:jed))
    allocate (Grid%dyunr(isd:ied,jsd:jed))
    allocate (Grid%dxter(isd:ied,jsd:jed))
    allocate (Grid%dyter(isd:ied,jsd:jed))
    allocate (Grid%dxtnr(isd:ied,jsd:jed))
    allocate (Grid%dytnr(isd:ied,jsd:jed))
    allocate (Grid%dyue_dxuer(isd:ied,jsd:jed))
    allocate (Grid%dxun_dyunr(isd:ied,jsd:jed))
    allocate (Grid%datr(isd:ied,jsd:jed))
    allocate (Grid%daur(isd:ied,jsd:jed))
    allocate (Grid%dater(isd:ied,jsd:jed))
    allocate (Grid%datnr(isd:ied,jsd:jed))

    allocate (Grid%dxt_dxter(isd:ied,jsd:jed))
    allocate (Grid%dxte_dxtr(isd:ied,jsd:jed))
    allocate (Grid%dxt_dxtnr(isd:ied,jsd:jed))
    allocate (Grid%dxtn_dxtr(isd:ied,jsd:jed))
    allocate (Grid%dyt_dyter(isd:ied,jsd:jed))
    allocate (Grid%dyte_dytr(isd:ied,jsd:jed))
    allocate (Grid%dyt_dytnr(isd:ied,jsd:jed))
    allocate (Grid%dytn_dytr(isd:ied,jsd:jed))

    allocate (Grid%dh1dy(isd:ied,jsd:jed))
    allocate (Grid%dh2dx(isd:ied,jsd:jed))

    allocate ( Grid%duw(isd:ied,jsd:jed) )
    allocate ( Grid%due(isd:ied,jsd:jed) )
    allocate ( Grid%dus(isd:ied,jsd:jed) )
    allocate ( Grid%dun(isd:ied,jsd:jed) )
    allocate ( Grid%dtw(isd:ied,jsd:jed) )
    allocate ( Grid%dte(isd:ied,jsd:jed) )
    allocate ( Grid%dts(isd:ied,jsd:jed) )
    allocate ( Grid%dtn(isd:ied,jsd:jed) )

    allocate(Grid%sin_rot(isd:ied,jsd:jed))
    allocate(Grid%cos_rot(isd:ied,jsd:jed))

    allocate (Grid%obc_tmask(isd:ied,jsd:jed))
    allocate (Grid%obc_umask(isd:ied,jsd:jed))

    #endif

    !--- initialize grid data

    Grid%xt=0.0;    Grid%yt=0.0;     Grid%grid_x_t=0.0; Grid%grid_y_t=0.0
    Grid%xu=0.0;    Grid%yu=0.0;     Grid%grid_x_u=0.0; Grid%grid_y_u=0.0
    Grid%dxtn=1.0;  Grid%dytn=1.0;   Grid%dxte=1.0;     Grid%dyte=1.0
    Grid%dxun=1.0;  Grid%dyun=1.0;   Grid%dxue=1.0;     Grid%dyue=1.0
    Grid%duw=0.0 ;  Grid%due=0.0 ;   Grid%dus=0.0 ;     Grid%dun=0.0
    Grid%dtw=0.0 ;  Grid%dte=0.0 ;   Grid%dts=0.0;      Grid%dtn=0.0
    Grid%zt=0.0;    Grid%zw=0.0;     Grid%sin_rot=0.0;  Grid%cos_rot=0.0
    Grid%dxt=epsln; Grid%dxu=epsln ; Grid%dyt=epsln ;   Grid%dyu=epsln
    Grid%dh1dy=0.0; Grid%dh2dx=0.0

    select case( grid_version )
    case( VERSION_0 )
      call read_data(ocean_hgrid, "zt",             Grid%zt)
      call read_data(ocean_hgrid, "zw",             Grid%zw)
      call read_data(ocean_hgrid, "gridlon_t",      Grid%grid_x_t)
      call read_data(ocean_hgrid, "gridlat_t",      Grid%grid_y_t)
      allocate(data(ni+1))
      call read_data(ocean_hgrid, "gridlon_vert_t", data)
      Grid%grid_x_u = data(2:ni+1)
      deallocate(data)
      allocate(data(nj+1))
      call read_data(ocean_hgrid, "gridlat_vert_t", data)
      Grid%grid_y_u = data(2:nj+1)
      deallocate(data)
    case( VERSION_1 )
      call read_data(ocean_hgrid, "zt",             Grid%zt)
      call read_data(ocean_hgrid, "zb",             Grid%zw)
      call read_data(ocean_hgrid, "grid_x_T",      Grid%grid_x_t)
      call read_data(ocean_hgrid, "grid_y_T",      Grid%grid_y_t)
      call read_data(ocean_hgrid, "grid_x_C",      Grid%grid_x_u)
      call read_data(ocean_hgrid, "grid_y_C",      Grid%grid_y_u)
    case( VERSION_2 )
      allocate(data(2*nk+1) )
      call read_data(ocean_vgrid, "zeta", data)
      do k=1,nk
        Grid%zt(k) = data(2*k)
        Grid%zw(k) = data(2*k+1)
      enddo
      deallocate(data)

      !--- The following adjust is for the consideration of static memory.
      ioff = isc1 - isc; joff = jsc1 - jsc
      isc2 = isc2 - 2*ioff; iec2 = iec2 - 2*ioff
      jsc2 = jsc2 - 2*joff; jec2 = jec2 - 2*joff
      isd2 = isd2 - 2*ioff; ied2 = ied2 - 2*ioff
      jsd2 = jsd2 - 2*joff; jed2 = jed2 - 2*joff

      if(isc2 .NE. 2*isc-1 .OR. iec2 .NE. 2*iec .OR. jsc2 .NE. 2*jsc-1 .OR. jec2 .NE. 2*jec ) then
        call error_mesg(FATAL, 'ocean_grids_mod (set_ocean_hgrid_arrays): supergrid compute domain is not set properly')
      endif
      if(isd2 .NE. 2*isd-1 .OR. ied2 .NE. 2*ied .OR. jsd2 .NE. 2*jsd-1 .OR. jed2 .NE. 2*jed ) then
        print*, isd, ied, jsd, jed, isd2, ied2, jsd2, jed2
        call error_mesg(FATAL, 'ocean_grids_mod (set_ocean_hgrid_arrays): supergrid data domain is not set properly')
      endif

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

      do j = 1, nj
        Grid%grid_y_t(j) = tmpy(1 , 2*j)
        Grid%grid_y_u(j) = tmpy(1 , 2*j+1)
      enddo
      deallocate(tmpx, tmpy)
      allocate(tmpx(isc2:iec2+1, jsc2:jec2+1))
      allocate(tmpy(isc2:iec2+1, jsc2:jec2+1))
    end select

    !--- Grid%xt
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'geolon_t', Grid%xt)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'x_T', Grid%xt)
    case(VERSION_2)
      Grid%xt(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2,2*jsc:2*jec:2)
    end select

    !--- Grid%yt
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'geolat_t', Grid%yt)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'y_T', Grid%yt)
    case(VERSION_2)
      Grid%yt(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2,2*jsc:2*jec:2)
    end select

    !--- Grid%xu
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'geolon_c', Grid%xu)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'x_C', Grid%xu)
    case(VERSION_2)
      Grid%xu(isc:iec,jsc:jec) = tmpx(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2)
    end select

    !--- Grid%yu
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'geolat_c', Grid%yu)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'y_C', Grid%yu)
    case(VERSION_2)
      Grid%yu(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2)
    end select

    if (Grid%beta_plane .or. Grid%f_plane) then
      Grid%phiu(:,:) = Grid%f_plane_latitude/radian
      Grid%phit(:,:) = Grid%phiu(:,:)
    else
      Grid%phit(:,:) = Grid%yt(:,:)/radian
      Grid%phiu(:,:) = Grid%yu(:,:)/radian
    endif

    Grid%h1t(:,:) = radius*cos(Grid%phit(:,:))
    where (cos(Grid%phit) == 0.0) Grid%h1t = radius*abs(epsln)
    Grid%h1u(:,:) = radius*cos(Grid%phiu(:,:))
    where (cos(Grid%phiu) == 0.0) Grid%h1u = radius*abs(epsln)
    Grid%h2t(:,:) = radius
    Grid%h2u(:,:) = radius

    ! set cell widths (meters) and area (meters^2)
    if(grid_version == VERSION_2) then
      deallocate(tmpx, tmpy)
      allocate(tmpx(isd2:ied2,jsd2:jed2+1), tmpy(isd2:ied2+1,jsd2:jed2))
      tmpx = 0.0
      tmpy = 0.0
      call read_data(ocean_hgrid, "dx", tmpx, position=NORTH)
      call read_data(ocean_hgrid, "dy", tmpy, position=EAST)
    end if

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

    ! --- Grid%dxu
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dxu', Grid%dxu)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_01_21_C', Grid%dxu)
    case(VERSION_2)
      Grid%dxu(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2,2*jsc+1:2*jec+1:2) + tmpx(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2)
    end select

    !--- Grid%dyu
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dyu', Grid%dyu)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_10_12_C', Grid%dyu)
    case(VERSION_2)
      Grid%dyu(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2,2*jsc:2*jec:2) + tmpy(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2)
    end select

    Grid%dat(:,:)  = Grid%dxt(:,:)*Grid%dyt(:,:)
    Grid%dau(:,:)  = Grid%dxu(:,:)*Grid%dyu(:,:)

    ! set lengths at edges of grid cells (meters)
    allocate(tmp1_local(isd:ied,jsd:jed), tmp2_local(isd:ied,jsd:jed))

    ! --- Grid%dxtn
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dxtn', Grid%dxtn)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_02_22_T', Grid%dxtn)
    case(VERSION_2)
      Grid%dxtn(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2,2*jsc+1:2*jec+1:2) + tmpx(2*isc:2*iec:2,2*jsc+1:2*jec+1:2)
    end select

    ! --- Grid%dytn
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dytn', Grid%dytn)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_00_02_C', Grid%dytn)
    case(VERSION_2)
      Grid%dytn(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2,2*jsc:2*jec:2) + tmpy(2*isc:2*iec:2,2*jsc+1:2*jec+1:2)
    end select

    ! --- Grid%dxte
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dxte', Grid%dxte)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_00_20_C', Grid%dxte)
    case(VERSION_2)
      Grid%dxte(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2,2*jsc:2*jec:2) + tmpx(2*isc+1:2*iec+1:2,2*jsc:2*jec:2)
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

    ! --- Grid%dxun
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dxun', Grid%dxun)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_02_22_C', Grid%dxun)
    case(VERSION_2)
      Grid%dxun(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2,2*jsc+2:2*jec+2:2) + tmpx(2*isc+1:2*iec+1:2,2*jsc+2:2*jec+2:2)
    end select

    ! --- Grid%dyun
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dyun', Grid%dyun)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_10_11_C', tmp2_local)
      call read_data(ocean_hgrid, 'ds_11_12_C', tmp1_local)
      Grid%dyun(isc:iec,jsc:jec) = tmp1_local(isc:iec,jsc:jec)+tmp2_local(isc:iec,jsc+1:jec+1)
    case(VERSION_2)
      Grid%dyun(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2) + tmpy(2*isc+1:2*iec+1:2,2*jsc+2:2*jec+2:2)
    end select

    ! --- Grid%dxue
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dxue', Grid%dxue)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_01_11_C', tmp2_local)
      call read_data(ocean_hgrid, 'ds_11_21_C', tmp1_local)
      Grid%dxue(isc:iec,jsc:jec) = tmp1_local(isc:iec,jsc:jec)+tmp2_local(isc+1:iec+1,jsc:jec)
    case(VERSION_2)
      Grid%dxue(isc:iec,jsc:jec) = tmpx(2*isc+1:2*iec+1:2,2*jsc+1:2*jec+1:2) + tmpx(2*isc+2:2*iec+2:2,2*jsc+1:2*jec+1:2)
    end select

    ! --- Grid%dyue
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dyue', Grid%dyue)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_20_22_C', Grid%dyue)
    case(VERSION_2)
      Grid%dyue(isc:iec,jsc:jec) = tmpy(2*isc+2:2*iec+2:2,2*jsc:2*jec:2) + tmpy(2*isc+2:2*iec+2:2,2*jsc+1:2*jec+1:2)
    end select

    ! set reciprocals and related quantities

    do j=jsd,jed
      do i=isd,ied
        Grid%dxtr(i,j)       = 1.0/(Grid%dxt(i,j)+epsln)
        Grid%dxur(i,j)       = 1.0/(Grid%dxu(i,j)+epsln)
        Grid%dytr(i,j)       = 1.0/(Grid%dyt(i,j)+epsln)
        Grid%dyur(i,j)       = 1.0/(Grid%dyu(i,j)+epsln)
        Grid%dxuer(i,j)      = 1.0/(Grid%dxue(i,j)+epsln)
        Grid%dyuer(i,j)      = 1.0/(Grid%dyue(i,j)+epsln)
        Grid%dxunr(i,j)      = 1.0/(Grid%dxun(i,j)+epsln)
        Grid%dyunr(i,j)      = 1.0/(Grid%dyun(i,j)+epsln)
        Grid%dxter(i,j)      = 1.0/(Grid%dxte(i,j)+epsln)
        Grid%dyter(i,j)      = 1.0/(Grid%dyte(i,j)+epsln)
        Grid%dxtnr(i,j)      = 1.0/(Grid%dxtn(i,j)+epsln)
        Grid%dytnr(i,j)      = 1.0/(Grid%dytn(i,j)+epsln)
        Grid%dyue_dxuer(i,j) = Grid%dyue(i,j)/(Grid%dxue(i,j)+epsln)
        Grid%dxun_dyunr(i,j) = Grid%dxun(i,j)/(Grid%dyun(i,j)+epsln)
        Grid%datr(i,j)       = 1.0/(Grid%dat(i,j)+epsln)
        Grid%daur(i,j)       = 1.0/(Grid%dau(i,j)+epsln)
        Grid%dater(i,j)      = 1.0/(Grid%dxte(i,j)*Grid%dyte(i,j)+epsln)
        Grid%datnr(i,j)      = 1.0/(Grid%dxtn(i,j)*Grid%dytn(i,j)+epsln)

        Grid%dxt_dxter(i,j)  = Grid%dxt(i,j)/(Grid%dxte(i,j)+epsln)
        Grid%dxte_dxtr(i,j)  = Grid%dxte(i,j)/(Grid%dxt(i,j)+epsln)
        Grid%dxt_dxtnr(i,j)  = Grid%dxt(i,j)/(Grid%dxtn(i,j)+epsln)
        Grid%dxtn_dxtr(i,j)  = Grid%dxtn(i,j)/(Grid%dxt(i,j)+epsln)
        Grid%dyt_dyter(i,j)  = Grid%dyt(i,j)/(Grid%dyte(i,j)+epsln)
        Grid%dyte_dytr(i,j)  = Grid%dyte(i,j)/(Grid%dyt(i,j)+epsln)
        Grid%dyt_dytnr(i,j)  = Grid%dyt(i,j)/(Grid%dytn(i,j)+epsln)
        Grid%dytn_dytr(i,j)  = Grid%dytn(i,j)/(Grid%dyt(i,j)+epsln)

      enddo
    enddo

    ! set quantities which account for rotation of basis vectors
    allocate(tmp_local(isd:ied,jsd:jed))
    tmp_local=1.0

    ! expanded version of dh2dx(:,:) = BDX_EU(tmp(:,:))
    do j=jsc,jec
      do i=isc,iec
        Grid%dh2dx(i,j) = (Grid%dyue(i,j)*tmp_local(i,j) - Grid%dyue(i-1,j)*tmp_local(i-1,j))*Grid%daur(i,j)
      enddo
    enddo
    if (.not. Grid%tripolar .and. jec+joff == nj) then
      do i=isc,iec
        Grid%dh2dx(i,jec) = 0.0
      enddo
    endif

    ! expanded version of dh1dy(:,:) = BDY_NU(tmp(:,:))
    do j=jsc,jec
      do i=isc,iec
        Grid%dh1dy(i,j) = (Grid%dxun(i,j)*tmp_local(i,j) - Grid%dxun(i,j-1)*tmp_local(i,j-1))*Grid%daur(i,j)
      enddo
    enddo
    if (.not. Grid%tripolar .and. jec+joff == nj) then
      do i=isc,iec
        Grid%dh1dy(i,jec) = 0.0
      enddo
    endif

    ! build distances from grid point to cell faces {meters}
    ! --- Grid%duw
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'duw', Grid%duw)
      call read_data(ocean_hgrid, 'due', tmp_local)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_01_11_C', Grid%duw)
      call read_data(ocean_hgrid, 'ds_11_21_C', tmp_local)
    case(VERSION_2)
      Grid%duw(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2, 2*jsc+1:2*jec+1:2)
      tmp_local(isc:iec,jsc:jec) = tmpx(2*isc+1:2*iec+1:2, 2*jsc+1:2*jec+1:2)
    end select

    ! --- Grid%due
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'due', Grid%due)
      call read_data(ocean_hgrid, 'duw', tmp_local)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_11_21_C', Grid%due)
      call read_data(ocean_hgrid, 'ds_01_11_C', tmp_local)
    case(VERSION_2)
      Grid%due(isc:iec,jsc:jec) = tmpx(2*isc+1:2*iec+1:2, 2*jsc+1:2*jec+1:2)
      tmp_local(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2, 2*jsc+1:2*jec+1:2)
    end select

    ! --- Grid%dus
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dus', Grid%dus(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'dun', tmp_local)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_10_11_C', Grid%dus(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'ds_11_12_C', tmp_local)
    case(VERSION_2)
      Grid%dus(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2, 2*jsc:2*jec:2)
      tmp_local(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2, 2*jsc+1:2*jec+1:2)
    end select

    ! --- Grid%dun
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dun', Grid%dun(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'dus', tmp_local)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_11_12_C', Grid%dun(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'ds_10_11_C', tmp_local)
    case(VERSION_2)
      Grid%dun(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2, 2*jsc+1:2*jec+1:2)
      tmp_local(isc:iec,jsc:jec) = tmpy(2*isc+1:2*iec+1:2, 2*jsc:2*jec:2)
    end select

    ! --- Grid%dtw
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dtw', Grid%dtw(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'dte', tmp_local)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_01_11_T', Grid%dtw(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'ds_11_21_T', tmp_local)
    case(VERSION_2)
      Grid%dtw(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2, 2*jsc:2*jec:2)
      tmp_local(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2, 2*jsc:2*jec:2)
    end select

    ! --- Grid%dte
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dte', Grid%dte(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'dtw', tmp_local)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_11_21_T', Grid%dte(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'ds_01_11_T', tmp_local)
    case(VERSION_2)
      Grid%dte(isc:iec,jsc:jec) = tmpx(2*isc:2*iec:2, 2*jsc:2*jec:2)
      tmp_local(isc:iec,jsc:jec) = tmpx(2*isc-1:2*iec-1:2, 2*jsc:2*jec:2)
    end select

    ! --- Grid%dts
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dts', Grid%dts(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'dtn', tmp_local)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_10_11_T', Grid%dts(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'ds_11_12_T', tmp_local)
    case(VERSION_2)
      Grid%dts(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2, 2*jsc-1:2*jec-1:2)
      tmp_local(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2, 2*jsc:2*jec:2)
    end select

    ! --- Grid%dtn
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'dtn', Grid%dtn(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'dts', tmp_local)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'ds_11_12_T', Grid%dtn(isc:iec,jsc:jec))
      call read_data(ocean_hgrid, 'ds_10_11_T', tmp_local)
    case(VERSION_2)
      Grid%dtn(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2, 2*jsc:2*jec:2)
      tmp_local(isc:iec,jsc:jec) = tmpy(2*isc:2*iec:2, 2*jsc-1:2*jec-1:2)
      deallocate(tmpx, tmpy)
    end select

    ! calculate rotation angles on velocity points
    select case(grid_version)
    case(VERSION_0)
      call read_data(ocean_hgrid, 'sin_rot', Grid%sin_rot)
      call read_data(ocean_hgrid, 'cos_rot', Grid%cos_rot)
    case(VERSION_1)
      call read_data(ocean_hgrid, 'angle_C', tmp_local)
      tmp_local     = tmp_local*deg_to_rad
      Grid%sin_rot(isc:iec,jsc:jec) = sin(tmp_local(isc:iec,jsc:jec))
      Grid%cos_rot(isc:iec,jsc:jec) = cos(tmp_local(isc:iec,jsc:jec))
    case(VERSION_2)
      allocate(tmp(isc2:iec2+1,jsc2:jec2+1) )
      call read_data(ocean_hgrid, "angle_dx", tmp, position=CORNER)
      tmp = tmp*deg_to_rad
      do j = jsc, jec
        do i = isc, iec
          Grid%sin_rot(i,j) = sin(tmp(2*i+1,2*j+1) )
          Grid%cos_rot(i,j) = cos(tmp(2*i+1,2*j+1) )
        end do
      end do
      deallocate(tmp)
    end select


    ! ensure that all grid factors are set properly at the southern boundary.
    ! many grid factors with index j=0 have no affect within the computational domain
    ! Even so, for completeness they are explicity set below.
    tmp_local=1.0
    if (jsc+joff == 1) then

      call error_mesg(NOTE, '==>Note from ocean_grids_mod (set_ocean_hgrid_arrays): altering U-grid arrays at j=0')

      ! phiu(:,0), h1u(:,0), and dxu(:,0) are within the computational domain and are
      ! not arbitrary. their values follow directly from grid_spec definitions.
      ! dyu(:,0) is arbitrary and has been set equal to dyu(:,1) by update_boundaries.
      ! However, a symmetric ocean domain requires that dyu be symmetric so dyu(:,0) = dyu(:,nj)

      if (Grid%beta_plane .or. Grid%f_plane) then
        Grid%phiu(:,0) = Grid%f_plane_latitude/radian
      else
        Grid%phiu(:,0) = Grid%yu(:,1)/radian - Grid%dyte(:,1)/radius
      endif
      Grid%phiu(:,0) = max(min(Grid%phiu(:,0),90.0/radian),-90.0/radian)
      where (cos(Grid%phiu(:,0)) == 0.0)
        Grid%h1u(:,0) = radius*abs(epsln)
      elsewhere
        Grid%h1u(:,0) = radius*cos(Grid%phiu(:,0))
      end where

      Grid%dxu(:,0)  = (Grid%dxu(:,1)/Grid%h1u(:,1))*Grid%h1u(:,0)
      Grid%dau(:,0)  = Grid%dxu(:,0)*Grid%dyu(:,0)

      Grid%dxur(:,0) = 1.0/(Grid%dxu(:,0)+epsln)
      Grid%dyur(:,0) = 1.0/(Grid%dyu(:,0)+epsln)
      Grid%daur(:,0) = 1.0/(Grid%dau(:,0)+epsln)

      Grid%dun(:,0)  = Grid%dyte(:,1) - Grid%dus(:,1)
      Grid%dus(:,0)  = Grid%dau(:,0)/(Grid%dxu(:,0)+epsln) - Grid%dun(:,0)
      Grid%due(:,0)  = (Grid%due(:,1)/Grid%h1u(:,1))*Grid%h1u(:,0)
      Grid%duw(:,0)  = Grid%dxu(:,0) - Grid%due(:,0)

      Grid%dxue(:,0)  = (Grid%dxue(:,1)/Grid%h1u(:,1))*Grid%h1u(:,0)
      Grid%dxun(:,0)  = (Grid%dxun(:,1)/Grid%h1t(:,2))*Grid%h1t(:,1)
      Grid%dxuer(:,0) = 1.0/(Grid%dxue(:,0)+epsln)
      Grid%dxunr(:,0) = 1.0/(Grid%dxun(:,0)+epsln)

      Grid%dyun(:,0)  = Grid%dyte(:,1)
      Grid%dyunr(:,0) = 1.0/(Grid%dyun(:,0)+epsln)
      Grid%dxun_dyunr(:,0) = Grid%dxun(:,0)/(Grid%dyun(:,0)+epsln)

      ! not correct, but any use should be zeroed out by the land mask at j=0
      Grid%dyue(:,0)  = Grid%dyu(:,0)
      Grid%dyuer(:,0) = 1.0/(Grid%dyue(:,0)+epsln)
      Grid%dyue_dxuer(:,0) = Grid%dyue(:,0)/(Grid%dxue(:,0)+epsln)

      Grid%dh2dx(:,0) = 0.0
      Grid%dh1dy(:,1) = (Grid%dxun(:,1)*tmp_local(:,1) - Grid%dxun(:,0)*tmp_local(:,0))*Grid%daur(:,1)
      Grid%dh1dy(:,0) = 0.0

    endif

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
  subroutine set_ocean_vgrid_arrays (Grid)

    type(ocean_grid_type),   intent(inout)        :: Grid

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
    allocate (Grid%tcella(nk))
    allocate (Grid%ucella(nk))
    allocate (Grid%mask(isd:ied,jsd:jed,nk))
    allocate (Grid%tmask(isd:ied,jsd:jed,nk))
    allocate (Grid%umask(isd:ied,jsd:jed,nk))
    allocate (Grid%tmasken(isd:ied,jsd:jed,nk,2))
    allocate (Grid%tmask_depth(isd:ied,jsd:jed,nk))
    allocate (Grid%umask_depth(isd:ied,jsd:jed,nk))

    allocate (Grid%dht_dx(isd:ied,jsd:jed))
    allocate (Grid%dht_dy(isd:ied,jsd:jed))
    allocate (Grid%gradH(isd:ied,jsd:jed))

    allocate (Grid%dzt(nk))
    allocate (Grid%dztlo(nk))
    allocate (Grid%dztup(nk))
    allocate (Grid%dzw(0:nk))

    allocate (Grid%st(nk))
    allocate (Grid%sw(nk))
    allocate (Grid%dst(nk))
    allocate (Grid%dstlo(nk))
    allocate (Grid%dstup(nk))
    allocate (Grid%dsw(0:nk))

    allocate (Grid%fracdz(nk,0:1))
    allocate (Grid%dzwr(0:nk))

    #endif

    ! set depth arrays
    Grid%dzw(0)         = Grid%zt(1)
    Grid%dzw(1)         = Grid%zt(2)-Grid%zt(1)
    Grid%dzwr(0)        = 1.0/Grid%dzw(0)
    Grid%dzwr(1)        = 1.0/Grid%dzw(1)
    Grid%dzt(1)         = Grid%zw(1)
    Grid%fracdz(1,0)    = Grid%zt(1)/Grid%dzt(1)
    Grid%fracdz(1,1)    = (Grid%zw(1) - Grid%zt(1))/Grid%dzt(1)
    do k=2,nk-1
      Grid%dzt(k)      = Grid%zw(k)-Grid%zw(k-1)
      Grid%dzw(k)      = Grid%zt(k+1)-Grid%zt(k)
      Grid%dzwr(k)     = 1.0/Grid%dzw(k)
      Grid%fracdz(k,0) = (Grid%zt(k) - Grid%zw(k-1))/Grid%dzt(k)
      Grid%fracdz(k,1) = (Grid%zw(k) - Grid%zt(k))/Grid%dzt(k)
    enddo
    Grid%dzt(nk)        = Grid%zw(nk)-Grid%zw(nk-1)
    Grid%dzw(nk)        = Grid%zw(nk) - Grid%zt(nk)
    Grid%dzwr(nk)       = 1.0/Grid%dzw(nk)
    Grid%fracdz(nk,0)   = (Grid%zt(nk) - Grid%zw(nk-1))/Grid%dzt(nk)
    Grid%fracdz(nk,1)   = (Grid%zw(nk) - Grid%zt(nk))/Grid%dzt(nk)

    ! half-thicknesses
    k=1
    Grid%dztlo(k) = Grid%zw(k)-Grid%zt(k)
    Grid%dztup(k) = Grid%dzw(k-1)
    do k=2,nk
      Grid%dztlo(k) = Grid%zw(k)-Grid%zt(k)
      Grid%dztup(k) = Grid%zt(k)-Grid%zw(k-1)
    enddo


    ! construct T cell and U cell land/sea masks
    do k=1,nk
      do j=jsd,jed
        do i=isd,ied
          if (Grid%kmt(i,j) .ge. k) then
            Grid%tmask(i,j,k) = 1.0
          else
            Grid%tmask(i,j,k) = 0.0
          endif
          if (Grid%kmu(i,j) .ge. k) then
            Grid%umask(i,j,k) = 1.0
          else
            Grid%umask(i,j,k) = 0.0
          endif
        enddo
      enddo
    enddo

    if(horz_grid == MOM_CGRID) then
      Grid%mask(:,:,:) = Grid%tmask(:,:,:)
      do k=1,nk
        do j=jsd,jec
          do i=isd,iec
            Grid%tmasken(i,j,k,1) = min(Grid%tmask(i,j,k),Grid%tmask(i+1,j,k))
            Grid%tmasken(i,j,k,2) = min(Grid%tmask(i,j,k),Grid%tmask(i,j+1,k))
          enddo
        enddo
      enddo
    else
      Grid%mask(:,:,:) = Grid%umask(:,:,:)
      do k=1,nk
        do j=jsd,jed
          do i=isd,ied
            Grid%tmasken(i,j,k,1) = Grid%umask(i,j,k)
            Grid%tmasken(i,j,k,2) = Grid%umask(i,j,k)
          enddo
        enddo
      enddo
    endif


    ! Compute masks as if depth or pressure were the vertical coordinates.
    ! This new mask is used for case when have vert_coordinate=ZSIGMA or PSIGMA,
    ! in which case tmask and umask are 1 for all k=1,nk except where ht=0.0.
    ! We need tmask_depth and umask_depth for vertical remapping in case when
    ! remap fields into depth or pressure space from terrain following s-space.
    Grid%tmask_depth(:,:,:) = 0.0
    Grid%umask_depth(:,:,:) = 0.0

    k=1
    do j=jsd,jed
      do i=isd,ied
        if (Grid%ht(i,j) > 0.0) then
          Grid%tmask_depth(i,j,k) = 1.0
        endif
        if (Grid%hu(i,j) > 0.0) then
          Grid%umask_depth(i,j,k) = 1.0
        endif
      enddo
    enddo
    do k=1,nk-1
      do j=jsd,jed
        do i=isd,ied
          if (Grid%ht(i,j) >= Grid%zw(k)) then
            Grid%tmask_depth(i,j,k+1) = 1.0
          endif
          if (Grid%hu(i,j) >= Grid%zw(k)) then
            Grid%umask_depth(i,j,k+1) = 1.0
          endif
        enddo
      enddo
    enddo

    ! compute surface area and volume of ocean (T cells and U cells)
    Grid%tcellv    = 0.0 !total ocean volume on T cells (assuming eta=0)
    Grid%ucellv    = 0.0 !total ocean volume on U cells (assuming eta=0)
    Grid%tcella(:) = 0.0 !total ocean surface area on T cells in level k
    Grid%ucella(:) = 0.0 !total ocean surface area on U cells in level k

    ! exclude points at open boundaries from diagnostics
    do j=jsc,jec
      do i=isc,iec
        kzt = Grid%kmt(i,j)
        if(have_obc) kzt = kzt * int(Grid%obc_tmask(i,j))
        if (kzt .gt. 0) then
          do k=1,kzt
            Grid%tcella(k) = Grid%tcella(k) + Grid%dat(i,j)
          enddo
          Grid%tcellv = Grid%tcellv + Grid%dat(i,j)*Grid%ht(i,j)
        endif
      enddo
    enddo

    do j=jsc,jec
      do i=isc,iec
        kzu = Grid%kmu(i,j)
        if(have_obc) kzu = kzu * int(Grid%obc_umask(i,j))
        if (kzu .gt. 0) then
          do k=1,kzu
            Grid%ucella(k) = Grid%ucella(k) + Grid%dau(i,j)
          enddo
          Grid%ucellv = Grid%ucellv + Grid%dau(i,j)*Grid%hu(i,j)
        endif
      enddo
    enddo

    ! compute area of full domain, including land regions
    tcella_sphere=0.0
    ucella_sphere=0.0
    do j=jsc,jec
      do i=isc,iec
        tcella_sphere = tcella_sphere + Grid%dat(i,j)
        ucella_sphere = ucella_sphere + Grid%dau(i,j)
      enddo
    enddo

    ! information about model grid size
    Grid%total_t_points = Grid%ni*Grid%nj*Grid%nk
    if(global_sum_flag==BITWISE_EXACT_SUM) then
      Grid%wet_t_points = nint(Grid%tmask(:,:,:), BITWISE_EXACT_SUM))
      Grid%wet_u_points = nint(Grid%umask(:,:,:), BITWISE_EXACT_SUM))
    else
      wet_t_points = real(Grid%tmask(:,:,:), NON_BITWISE_EXACT_SUM), KIND=4)
      wet_u_points = real(Grid%umask(:,:,:), NON_BITWISE_EXACT_SUM), KIND=4)
      Grid%wet_t_points = nint(wet_t_points)
      Grid%wet_u_points = nint(wet_u_points)
    endif

    ! compute bitwise reproducible surface areas
    if(global_sum_flag==BITWISE_EXACT_SUM) then
      if(have_obc) then
        Grid%tcellsurf = &
        mpp_global_sum(Grid%dat(:,:)*Grid%tmask(:,:,1)*Grid%obc_tmask(:,:), BITWISE_EXACT_SUM)
        Grid%ucellsurf = &
        mpp_global_sum(Grid%dau(:,:)*Grid%umask(:,:,1)*Grid%obc_umask(:,:), BITWISE_EXACT_SUM)
      else
        Grid%tcellsurf = &
        mpp_global_sum(Grid%dat(:,:)*Grid%tmask(:,:,1), BITWISE_EXACT_SUM)
        Grid%ucellsurf = &
        mpp_global_sum(Grid%dau(:,:)*Grid%umask(:,:,1), BITWISE_EXACT_SUM)
      endif
    else
      if(have_obc) then
        Grid%tcellsurf = &
        real(mpp_global_sum(Grid%dat(:,:)*Grid%tmask(:,:,1)*Grid%obc_tmask(:,:), NON_BITWISE_EXACT_SUM), KIND=4)
        Grid%ucellsurf = &
        real(mpp_global_sum(Grid%dau(:,:)*Grid%umask(:,:,1)*Grid%obc_umask(:,:), NON_BITWISE_EXACT_SUM), KIND=4)
      else
        Grid%tcellsurf = &
        real(mpp_global_sum(Grid%dat(:,:)*Grid%tmask(:,:,1), NON_BITWISE_EXACT_SUM), KIND=4)
        Grid%ucellsurf = &
        real(mpp_global_sum(Grid%dau(:,:)*Grid%umask(:,:,1), NON_BITWISE_EXACT_SUM), KIND=4)
      endif
    endif


    ! compute fraction of total wet area occupied by a single tracer cell
    do j=jsc,jec
      do i=isc,iec
        Grid%dat_frac(i,j) = Grid%dat(i,j)/Grid%tcellsurf
      enddo
    enddo

    ! compute topographic slopes at U-cell. Note: cannot use
    ! operators as ocean_operators_init has not yet been called.
    do j=jsc,jec
      do i=isc,iec
        ht_fax           = 0.5*( Grid%ht(i,j) + Grid%ht(i+1,j))
        ht_fay           = 0.5*( Grid%ht(i,j) + Grid%ht(i,j+1))
        ht_jp1_fax       = 0.5*( Grid%ht(i,j+1) + Grid%ht(i+1,j+1))
        ht_ip1_fay       = 0.5*( Grid%ht(i+1,j) + Grid%ht(i+1,j+1))
        Grid%dht_dx(i,j) = (ht_ip1_fay-ht_fay)*Grid%dxur(i,j)*Grid%umask(i,j,1)
        Grid%dht_dy(i,j) = (ht_jp1_fax-ht_fax)*Grid%dyur(i,j)*Grid%umask(i,j,1)
        Grid%gradH(i,j)  = sqrt(Grid%dht_dx(i,j)**2 + Grid%dht_dy(i,j)**2)
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
        Grid%st(k)    = Grid%zt(k)
        Grid%sw(k)    = Grid%zw(k)
        Grid%dst(k)   = Grid%dzt(k)
        Grid%dstlo(k) = Grid%dztlo(k)
        Grid%dstup(k) = Grid%dztup(k)
      enddo
      do k=0,nk
        Grid%dsw(k) = Grid%dzw(k)
      enddo

    elseif(vert_coordinate == PRESSURE .or. vert_coordinate == PSTAR) then
      do k=1,nk
        Grid%st(k)    =  rho0_profile(k)*grav*Grid%zt(k)
        Grid%sw(k)    =  rho0_profile(k)*grav*Grid%zw(k)
        Grid%dst(k)   = -rho0_profile(k)*grav*Grid%dzt(k)
        Grid%dstlo(k) = -rho0_profile(k)*grav*Grid%dztlo(k)
        Grid%dstup(k) = -rho0_profile(k)*grav*Grid%dztup(k)
      enddo
      do k=0,nk
        kmin = max(1,k)
        Grid%dsw(k)   = -rho0_profile(kmin)*grav*Grid%dzw(k)
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

      zwnk   = Grid%zw(nk)
      zwnk_r = 1.0/(epsln + zwnk)

      do k=1,nk
        Grid%st(k) = sigma_sign*Grid%zt(k)*zwnk_r
        Grid%sw(k) = sigma_sign*Grid%zw(k)*zwnk_r
      enddo

      ! set s-grid increments from st and sw
      Grid%dsw(0) = - Grid%st(1)
      Grid%dsw(1) = -(Grid%st(2)-Grid%st(1))
      Grid%dst(1) = - Grid%sw(1)
      do k=2,nk-1
        Grid%dst(k) = -(Grid%sw(k)  - Grid%sw(k-1))
        Grid%dsw(k) = -(Grid%st(k+1)- Grid%st(k))
      enddo
      Grid%dst(nk) = -(Grid%sw(nk)-Grid%sw(nk-1))
      Grid%dsw(nk) = -(Grid%sw(nk)-Grid%st(nk))

      k=1
      Grid%dstlo(k) = Grid%sw(k)-Grid%st(k)
      Grid%dstup(k) = Grid%dsw(k-1)
      do k=2,nk
        Grid%dstlo(k) = Grid%sw(k)-Grid%st(k)
        Grid%dstup(k) = Grid%st(k)-Grid%sw(k-1)
      enddo

    endif

    ! set up the axis definition for variables
    call axes_info(Grid)


  end subroutine set_ocean_vgrid_arrays
  ! </SUBROUTINE> NAME="set_ocean_vgrid_arrays"

  !#######################################################################
  ! <SUBROUTINE NAME="axes_info">
  !
  ! <DESCRIPTION>
  ! Set up axes definitions.
  ! </DESCRIPTION>
  !
  subroutine axes_info(Grid)

    type(ocean_grid_type),   intent(inout) :: Grid

    integer :: id_xt, id_yt, id_xu, id_yu
    integer :: id_zt_bounds, id_zt, id_zw_bounds, id_zw
    integer :: id_st_bounds, id_st, id_sw_bounds, id_sw
    integer :: k

    real, allocatable, dimension(:) :: zt_bounds
    real, allocatable, dimension(:) :: zw_bounds
    real, allocatable, dimension(:) :: st_bounds
    real, allocatable, dimension(:) :: sw_bounds

    real :: meter2dbar
    meter2dbar = c2dbars*rho0*grav


    ! horizontal axes



    ! vertical axes (needs to be re-thought for partial cells and general vertical coordinates)

    allocate(zt_bounds(Grid%nk+1))
    allocate(zw_bounds(Grid%nk+1))
    zt_bounds(1) = Grid%zt(1)-Grid%dzt(1)/2.0
    zw_bounds(1) = Grid%zw(1)-Grid%dzw(1)/2.0
    do k=2,Grid%nk+1
      zt_bounds(k)=zt_bounds(k-1)+Grid%dzt(k-1)
      zw_bounds(k)=zw_bounds(k-1)+Grid%dzw(k-1)
    enddo

    allocate(st_bounds(Grid%nk+1))
    allocate(sw_bounds(Grid%nk+1))

    if(vert_coordinate_class==DEPTH_BASED) then
      ! for depth-based, st and sw are > 0 and dst and dsw are > 0
      st_bounds(1) = Grid%st(1)-Grid%dst(1)/2.0
      sw_bounds(1) = Grid%sw(1)-Grid%dsw(1)/2.0
      do k=2,Grid%nk+1
        st_bounds(k)=st_bounds(k-1)+Grid%dst(k-1)
        sw_bounds(k)=sw_bounds(k-1)+Grid%dsw(k-1)
      enddo
    else
      ! for pressure-based, st and sw are > 0 but dst and dsw are < 0
      st_bounds(1) = Grid%st(1)+Grid%dst(1)/2.0
      sw_bounds(1) = Grid%sw(1)+Grid%dsw(1)/2.0
      do k=2,Grid%nk+1
        st_bounds(k)=st_bounds(k-1)-Grid%dst(k-1)
        sw_bounds(k)=sw_bounds(k-1)-Grid%dsw(k-1)
      enddo
    endif


    ! attributes for variables

    Grid%tracer_axes        = (/ id_xt, id_yt, id_st /)
    Grid%tracer_axes_flux_x = (/ id_xu, id_yt, id_st /)
    Grid%tracer_axes_flux_y = (/ id_xt, id_yu, id_st /)
    Grid%tracer_axes_wt     = (/ id_xt, id_yt, id_sw /)

    Grid%tracer_axes_depth        = (/ id_xt, id_yt, id_zt /)
    Grid%tracer_axes_flux_x_depth = (/ id_xu, id_yt, id_zt /)
    Grid%tracer_axes_flux_y_depth = (/ id_xt, id_yu, id_zt /)
    Grid%tracer_axes_wt_depth     = (/ id_xt, id_yt, id_zw /)

    Grid%vel_axes_uv              = (/ id_xu, id_yu, id_st /)
    Grid%vel_axes_wu_depth        = (/ id_xu, id_yu, id_zw /)
    Grid%vel_axes_wt_depth        = (/ id_xt, id_yt, id_zw /)
    Grid%vel_axes_uv_depth        = (/ id_xu, id_yu, id_zt /)
    Grid%vel_axes_flux_x_depth    = (/ id_xt, id_yu, id_zt /)
    Grid%vel_axes_flux_y_depth    = (/ id_xu, id_yt, id_zt /)

    Grid%vel_axes_wu              = (/ id_xu, id_yu, id_sw /)
    Grid%vel_axes_wt              = (/ id_xt, id_yt, id_sw /)
    Grid%vel_axes_flux_x          = (/ id_xt, id_yu, id_st /)
    Grid%vel_axes_flux_y          = (/ id_xu, id_yt, id_st /)

    if(horz_grid == MOM_BGRID) then
      Grid%vel_axes_u = (/ id_xu, id_yu, id_st /)
      Grid%vel_axes_v = (/ id_xu, id_yu, id_st /)
    else
      Grid%vel_axes_u = (/ id_xu, id_yt, id_st /)
      Grid%vel_axes_v = (/ id_xt, id_yu, id_st /)
    endif

    deallocate(zt_bounds, zw_bounds)
    deallocate(st_bounds, sw_bounds)

  end subroutine axes_info
  ! </SUBROUTINE> NAME="axes_info"

  !#############################################################################
  ! This routine will get the actual file name, as well as if read_dist is true or false.
  ! return true if such file exist and return false if not.
  function get_file_name(orig_file, actual_file, read_dist)
    character(len=*),                 intent(in) :: orig_file
    character(len=*),                intent(out) :: actual_file
    logical,                         intent(out) :: read_dist
    logical                                      :: get_file_name

    logical                       :: fexist
    character(len=256)            :: fname

    fexist          = .false.
    read_dist       = .false.
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
  subroutine get_file_unit(filename, unit, index_file, read_dist)
    character(len=*),         intent(in) :: filename
    integer,                 intent(out) :: unit, index_file
    logical,                  intent(in) :: read_dist

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
    if (read_dist .and. thread_r == MPP_SINGLE) then
      call error_mesg(FATAL,'fms_io(get_file_unit): single-threaded read from distributed fileset not allowed' &
      //'change threading_read to MULTI')
    endif

    files_read(num_files_r)%name = trim(filename)
    allocate(files_read(num_files_r)%var (max_fields) )
    files_read(num_files_r)%nvar = 0
    index_file = num_files_r
    files_read(index_file)%unit = unit

  end subroutine get_file_unit

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
    logical                      :: file_exist, read_dist
    character(len=256)           :: fname

    field_exist = .false.
    if (len_trim(field_name) == 0) return
    if (field_name(1:1) == ' ')    return

    file_exist=get_file_name(file_name, fname, read_dist)
    if(file_exist) then
      call get_file_unit(fname, unit, nfile, read_dist)
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

  #ifdef MOM_STATIC_ARRAYS
  !########################

  type, private :: ocean_thickness_type

    integer :: method       ! energetic or finite volume

    real, dimension(isd:ied,jsd:jed,nk,3) :: rho_dzt   ! rho(kg/m^3)*thickness (m) of T cell at 3 times
    real, dimension(isd:ied,jsd:jed,nk,2) :: rho_dzten ! rho(kg/m^3)*thickness (m) at east/north face of T-cell
    real, dimension(isd:ied,jsd:jed,nk)   :: rho_dztr  ! 1.0/(rho*dzt) at time taup1
    real, dimension(isd:ied,jsd:jed,nk,3) :: rho_dzu   ! rho (kg/m^3) * thickness (m) of U cell at 3 times
    real, dimension(isd:ied,jsd:jed,nk)   :: rho_dzur  ! 1.0/(rho*dzu) at time taup1

    real, dimension(isd:ied,jsd:jed,nk)   :: rho_dzt_tendency ! rho_dzt tendency (kg/m^3)*(m/s)

    real, dimension(isd:ied,jsd:jed)      :: sea_lev     ! eta_t + patm/(rho0*grav) - eta_geoid - eta_tide (m) at time taup1 for coupler
    real, dimension(isd:ied,jsd:jed,nk)   :: dzt         ! thickness (m) of T cell at time tau/taup1
    real, dimension(isd:ied,jsd:jed,nk,2) :: dzten       ! thickness (m) of east/north face of T cell at time tau/taup1
    real, dimension(isd:ied,jsd:jed,nk)   :: dzu         ! thickness (m) of U cell at time tau/taup1
    real, dimension(isd:ied,jsd:jed,0:nk) :: dzwt        ! vertical distance (m) between T points at tau/taup1
    real, dimension(isd:ied,jsd:jed,0:nk) :: dzwu        ! vertical distance (m) between U points at tau/taup1

    real, dimension(isd:ied,jsd:jed,nk)   :: dztup       ! distance (m) from T-cell point to top of T-cell
    real, dimension(isd:ied,jsd:jed,nk)   :: dztlo       ! distance (m) from T-cell point to bottom of T-cell

    real, dimension(isd:ied,jsd:jed,nk)   :: geodepth_zt ! vertical distance (m) from z=0 to T-point
    real, dimension(isd:ied,jsd:jed,nk)   :: depth_zt    ! vertical distance (m) from column top to T-point
    real, dimension(isd:ied,jsd:jed,nk)   :: depth_zwt   ! vertical distance (m) from column top to T-bottom
    real, dimension(isd:ied,jsd:jed,nk)   :: depth_zu    ! vertical distance (m) from column top to U-point
    real, dimension(isd:ied,jsd:jed,nk)   :: depth_zwu   ! vertical distance (m) from column top to U-bottom

    real, dimension(isd:ied,jsd:jed,nk)   :: depth_st    ! s-distance to T-point
    real, dimension(isd:ied,jsd:jed,nk)   :: depth_swt   ! s-distance to T-bottom
    real, dimension(isd:ied,jsd:jed,nk)   :: dst         ! vertical increment of s-coordinate on T-cell
    real, dimension(isd:ied,jsd:jed,0:nk) :: dswt        ! s-coordinate T-point increment at time tau
    real, dimension(isd:ied,jsd:jed,nk)   :: dzt_dst     ! specific thickness (metre/s-coordinate) on T-cell

    real, dimension(isd:ied,jsd:jed,nk)   :: dstup       ! s-distance (s-units) from T-cell point to top of T-cell
    real, dimension(isd:ied,jsd:jed,nk)   :: dstlo       ! s-distance (s-units) from T-cell point to bottom of T-cell

    real, dimension(isd:ied,jsd:jed)      :: pbot0       ! reference bottom pressure (Pa)
    real, dimension(isd:ied,jsd:jed)      :: pbot0r      ! inverse of reference bottom pressure (1/Pa)

    real, dimension(isd:ied,jsd:jed,nk)   :: mass_source ! mass source (kg/m^3)*(m/sec)
    real, dimension(isd:ied,jsd:jed,3)    :: mass_u      ! mass per area (kg/m^2) in a velocity column
    real, dimension(isd:ied,jsd:jed,2)    :: mass_en     ! vertical sum rho_dzten (kg/m^2)
    real, dimension(isd:ied,jsd:jed,3)    :: thicku      ! thickness on U-cell (metre) from z=eta_u to z=-H.
    real, dimension(isd:ied,jsd:jed,2)    :: thicken     ! sum of dzten (metre)


    ! The following arrays are never allocated with MOM_STATIC_ARRAYS,
    ! but they need to be declared allocatable in order to compile.
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztL     _NULL ! L system contribution to rho_dztT (3 time levels)
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
    real, dimension(:,:,:),   _ALLOCATABLE :: dztupL       _NULL ! L system contribution to dztupT
    real, dimension(:,:,:),   _ALLOCATABLE :: dztupT       _NULL ! distance (m) from T-cell point to top of T-cell
    real, dimension(:,:,:),   _ALLOCATABLE :: dztloL       _NULL ! L system contribution to dztloT
    real, dimension(:,:,:),   _ALLOCATABLE :: dztloT       _NULL ! distance (m) from T-cell point to bottom of T-cell
    real, dimension(:,:,:),   _ALLOCATABLE :: mass_uT      _NULL ! mass per area (kg/m^2) in a velocity column
    real, dimension(:,:,:),   _ALLOCATABLE :: geodepth_zwt _NULL ! vert distance (m) from z=0 to bottom of T-cell
    real, dimension(:,:),     _ALLOCATABLE :: blob_source _NULL  ! mass source (sink) for the L system (E system)

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
    real, dimension(nk)                 :: dztlo  ! distance (m) from T-cell point to bottom of T-cell
    real, dimension(nk)                 :: dztup  ! distance (m) from T-cell point to top of T-cell
    real, dimension(0:nk)               :: dzw    ! initial vertical resolution of W grid cells (m)
    real, dimension(0:nk)               :: dzwr   ! reciprocal of dzw (W cell vertical resolution)

    real, dimension(nk)                 :: st     ! full cell s-depth from surface to level k T-cell
    real, dimension(nk)                 :: sw     ! full cell s-depth from surface to bottom of k T-cell
    real, dimension(nk)                 :: dst    ! initial vertical s-resolution of T or U grid cells (m)
    real, dimension(nk)                 :: dstlo  ! s-distance (s-units) from T-cell point to bottom of T-cell
    real, dimension(nk)                 :: dstup  ! s-distance (s-units) from T-cell point to top of T-cell
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
    real, dimension(isd:ied,jsd:jed) :: phiu       ! latitude of U grid point in radians
    real, dimension(isd:ied,jsd:jed) :: phit       ! latitude of T grid point in radians
    real, dimension(isd:ied,jsd:jed) :: h1t        ! metric factors in i grid direction
    real, dimension(isd:ied,jsd:jed) :: h1u        ! metric factors in j grid jirection
    real, dimension(isd:ied,jsd:jed) :: h2t        ! metric factors in i grid direction
    real, dimension(isd:ied,jsd:jed) :: h2u        ! metric factors in j grid jirection
    real, dimension(isd:ied,jsd:jed) :: dh2dx      ! (1/delta_y)*d(delta_y)/dx (1/m)
    real, dimension(isd:ied,jsd:jed) :: dh1dy      ! (1/delta_x)*d(delta_x)/dy (1/m)
    real, dimension(isd:ied,jsd:jed) :: dxt        ! longitudinal width of T-cells at grid point (m)
    real, dimension(isd:ied,jsd:jed) :: dxu        ! longitudinal width of U-cells at grid point (m)
    real, dimension(isd:ied,jsd:jed) :: dyt        ! latitudinal width of T-cells at grid point (m)
    real, dimension(isd:ied,jsd:jed) :: dyu        ! latitudinal width of U-cells at grid point (m)
    real, dimension(isd:ied,jsd:jed) :: dau        ! area of U-cells (m^2)
    real, dimension(isd:ied,jsd:jed) :: dat        ! area of T-cells (m^2)
    real, dimension(isd:ied,jsd:jed) :: dat_frac   ! fraction of total wet area occuped by a T-cell (dimensionless)
    real, dimension(isd:ied,jsd:jed) :: dxte       ! long-width between grid points at i+1 and i in T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dxtn       ! long-width of north face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dyte       ! lat-width of east face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dytn       ! lat-width between grid points at j+1 and j in T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dxue       ! long-width between grid points at i+1 and i in U-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dxun       ! long-width of north face of U-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dyue       ! lat-width of east face of U-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dyun       ! lat-width between grid points at j+1 and j in U-cells (m)
    real, dimension(isd:ied,jsd:jed) :: datnr      ! reciprocal area at north face of T-cell
    real, dimension(isd:ied,jsd:jed) :: dater      ! reciprocal area at east face of T-cell
    real, dimension(isd:ied,jsd:jed) :: dun        ! width from grid point to north face of U-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dus        ! width from grid point to south face of U-cells (m)
    real, dimension(isd:ied,jsd:jed) :: duw        ! width from grid point to west  face of U-cells (m)
    real, dimension(isd:ied,jsd:jed) :: due        ! width from grid point to east  face of U-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dtn        ! width from grid point to north face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dts        ! width from grid point to south face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dtw        ! width from grid point to west  face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dte        ! width from grid point to east  face of T-cells (m)
    real, dimension(isd:ied,jsd:jed) :: dxtr       ! 1/dxt
    real, dimension(isd:ied,jsd:jed) :: dxur       ! 1/dxu
    real, dimension(isd:ied,jsd:jed) :: dytr       ! 1/dyt
    real, dimension(isd:ied,jsd:jed) :: dyur       ! 1/dyu
    real, dimension(isd:ied,jsd:jed) :: daur       ! 1/[area of U-cells (m^2)]
    real, dimension(isd:ied,jsd:jed) :: datr       ! 1/[area of T-cells (m^2)]
    real, dimension(isd:ied,jsd:jed) :: dxter      ! 1/dxte
    real, dimension(isd:ied,jsd:jed) :: dyter      ! 1/dyte
    real, dimension(isd:ied,jsd:jed) :: dxtnr      ! 1/dxtn
    real, dimension(isd:ied,jsd:jed) :: dytnr      ! 1/dytn
    real, dimension(isd:ied,jsd:jed) :: dxuer      ! 1/dxue
    real, dimension(isd:ied,jsd:jed) :: dxunr      ! 1/dxun
    real, dimension(isd:ied,jsd:jed) :: dyuer      ! 1/dyue
    real, dimension(isd:ied,jsd:jed) :: dyunr      ! 1/dyun
    real, dimension(isd:ied,jsd:jed) :: dyue_dxuer ! dyue/dxue
    real, dimension(isd:ied,jsd:jed) :: dxun_dyunr ! dxun/dyun
    real, dimension(isd:ied,jsd:jed) :: dxt_dxter  ! dxt/dxte
    real, dimension(isd:ied,jsd:jed) :: dxte_dxtr  ! dxte/dxt
    real, dimension(isd:ied,jsd:jed) :: dxt_dxtnr  ! dxt/dxtn
    real, dimension(isd:ied,jsd:jed) :: dxtn_dxtr  ! dxtn/dxt
    real, dimension(isd:ied,jsd:jed) :: dyt_dyter  ! dyt/dyte
    real, dimension(isd:ied,jsd:jed) :: dyte_dytr  ! dyte/dyt
    real, dimension(isd:ied,jsd:jed) :: dyt_dytnr  ! dyt/dytn
    real, dimension(isd:ied,jsd:jed) :: dytn_dytr  ! dytn/dyt

    ! land/sea masks
    real, dimension(isd:ied,jsd:jed,nk)   :: mask        ! land/sea tmask or umask depending on C/B grids
    real, dimension(isd:ied,jsd:jed,nk)   :: tmask       ! land/sea mask for T cells based on s-coordinate
    real, dimension(isd:ied,jsd:jed,nk)   :: umask       ! land/sea mask for B-grid U cells based on s-coordinate
    real, dimension(isd:ied,jsd:jed,nk,2) :: tmasken     ! land/sea mask for east/north face of t-cell based on s-coordinate
    real, dimension(isd:ied,jsd:jed,nk)   :: tmask_depth ! based on depth-based vert_coordinate
    real, dimension(isd:ied,jsd:jed,nk)   :: umask_depth ! based on depth-based vert_coordinate

    ! grid areas and volumes
    real                 :: tcellv    ! initial T cell volume m^3 (entire resting ocean)
    real                 :: ucellv    ! initial U cell volume m^3 (entire resting ocean)
    real                 :: tcellsurf ! T cell surface area (k=1) in bitwise reproducible form
    real                 :: ucellsurf ! U cell surface area (k=1) in bitwise reproducible form
    real, dimension(nk)  :: tcella    ! T cell surface area m^2 (entire ocean)
    real, dimension(nk)  :: ucella    ! U cell surface area m^2 (entire ocean)

    ! model grid information
    integer :: wet_t_points   ! total number of wet tracer points
    integer :: wet_u_points   ! total number of wet B-grid velocity points
    integer :: total_t_points ! total number of wet or dry tracer points

    ! sine and cosine of rotation angles (clockwise) of velocity for tripolar
    real, dimension(isd:ied,jsd:jed) :: sin_rot
    real, dimension(isd:ied,jsd:jed) :: cos_rot

    ! 1-d grid coordinates for COARDS NetCDF files
    real, dimension(ni) :: grid_x_t
    real, dimension(nj) :: grid_y_t
    real, dimension(ni) :: grid_x_u
    real, dimension(nj) :: grid_y_u

    ! axes id for diagnostic manager
    integer, dimension(3)  :: tracer_axes
    integer, dimension(3)  :: vel_axes_uv
    integer, dimension(3)  :: vel_axes_u
    integer, dimension(3)  :: vel_axes_v
    integer, dimension(3)  :: vel_axes_wu
    integer, dimension(3)  :: vel_axes_wt
    integer, dimension(3)  :: tracer_axes_wt
    integer, dimension(3)  :: tracer_axes_flux_x
    integer, dimension(3)  :: tracer_axes_flux_y
    integer, dimension(3)  :: vel_axes_flux_x
    integer, dimension(3)  :: vel_axes_flux_y

    ! axes id for diagnostic manager, appropriate
    ! when with to remap native vertical fields
    ! to depth or pressure levels.  Appropropriate
    ! when vert_coordinate == ZSIGMA or PSIGMA, or other
    ! vertical coordinats whose iso-surfaces are
    ! not quasi-horizontal
    integer, dimension(3)  :: tracer_axes_depth
    integer, dimension(3)  :: vel_axes_uv_depth
    integer, dimension(3)  :: vel_axes_u_depth
    integer, dimension(3)  :: vel_axes_v_depth
    integer, dimension(3)  :: vel_axes_wu_depth
    integer, dimension(3)  :: vel_axes_wt_depth
    integer, dimension(3)  :: tracer_axes_wt_depth
    integer, dimension(3)  :: tracer_axes_flux_x_depth
    integer, dimension(3)  :: tracer_axes_flux_y_depth
    integer, dimension(3)  :: vel_axes_flux_x_depth
    integer, dimension(3)  :: vel_axes_flux_y_depth

  end type ocean_grid_type


  type, private :: ocean_domain_type
    type(domain2d) :: domain2d           ! fms variable, used by mpp routines
    integer :: isc, iec, jsc, jec        ! computational domain indices
    integer :: isd, ied, jsd, jed        ! local indices including halo, consistent with domain2d
    integer :: isg, ieg, jsg, jeg        ! global indices
    integer :: isa, iea, jsa, jea        ! active indices (for minimizing comm2d calls)
    integer :: xhalo, yhalo              ! halo sizes
    integer :: xflags, yflags, layout(2) ! options to mpp_define_domains
    integer :: ioff , joff               ! index offset for absolute indices if MOM_STATIC_ARRAYS (0 otherwise)
    integer :: io_layout(2)              ! options to define io_domain.
    integer :: x_cyclic_offset           ! offset applied to x-direction cyclic boundary condition
    integer :: y_cyclic_offset           ! offset applied to y-direction cyclic boundary condition
    logical, pointer :: maskmap(:,:) =>NULL() ! option to mpp_define_domains
  end type ocean_domain_type

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
    real, dimension(isd:ied,jsd:jed,nk)   :: uhrho_eu   ! remapped uhrho_et (kg/(m*s)) onto i-face of U-cell
    real, dimension(isd:ied,jsd:jed,nk)   :: vhrho_nu   ! remapped vhrho_nt (kg/(m*s)) onto j-face of U-cell
    real, dimension(isd:ied,jsd:jed,0:nk) :: wrho_bt    ! rho * vertical advect vel (kg/(m^2*s)) on T-bottom
    real, dimension(isd:ied,jsd:jed,0:nk) :: wrho_bu    ! remapped wrho_bt onto U-bottom
    real, dimension(isd:ied,jsd:jed,nk)   :: diverge_t  ! divergence on T-cell of horiz momentum per mass (kg/m3)*(m/s)
    real, dimension(isd:ied,jsd:jed,nk)   :: diverge_u  ! divergence on U-cell of horiz momentum per mass (kg/m3)*(m/s)
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
    integer :: horz_advect_scheme=-1  ! id for horizontal advection scheme
    integer :: vert_advect_scheme=-1  ! id for vertical advection scheme
    integer :: ppm_hlimiter=1          ! Limiter for use with PPM in horizontal
    integer :: ppm_vlimiter=1          ! Limiter for use with PPM in vertical
    integer :: mdt_scheme=4            ! Version of Multi-Dim. Modified Daru & Tenaud (MDMDT)

    type(obc_flux), _ALLOCATABLE, dimension(:) :: otf   _NULL ! flux through open boundaries, allocate nobc

    real, dimension(isd:ied,jsd:jed,nk,3) :: field            ! tracer concentration at 3 time levels
    real, dimension(isd:ied,jsd:jed,nk)   :: th_tendency      ! thickness weighted tracer tendency
    real, dimension(isd:ied,jsd:jed,nk)   :: tendency         ! for diagnostics: tendency concentration [concentration/sec]
    real, dimension(isd:ied,jsd:jed,nk)   :: source           ! tracer source [=tracer concentration per time]
    real, dimension(isd:ied,jsd:jed,nk)   :: wrk1             ! work array
    real, dimension(isd:ied,jsd:jed,nk)   :: tmask_limit      ! to limit advection &/or neutral physics fluxes
    real, dimension(isd:ied,jsd:jed,nk)   :: K33_implicit     ! m^2/sec vert-diffusivity from neutral diffusion
    real, dimension(isd:ied,jsd:jed,nk)   :: radiation        ! radiation absorbed within a cell [W/m^2]
    real, dimension(isd:ied,jsd:jed)      :: stf              ! surface tracer flux [rho*m/sec*tracer concen]
    real, dimension(isd:ied,jsd:jed)      :: btf              ! bottom tracer flux [rho*m/sec*tracer concen]
    real, dimension(isd:ied,jsd:jed)      :: tpme             ! tracer concentration in precip-evap
    real, dimension(isd:ied,jsd:jed)      :: triver           ! tracer concentration in river(=runoff+calving) water
    real, dimension(isd:ied,jsd:jed)      :: trunoff          ! tracer concentration in liquid runoff from land
    real, dimension(isd:ied,jsd:jed)      :: tcalving       ! tracer concentration in frozen runoff from land (e.g., calving ice)
    real, dimension(isd:ied,jsd:jed) :: runoff_tracer_flux  ! tracer flux in liquid runoff (e.g., kg*degC/(m^2 s) for temp)
    real, dimension(isd:ied,jsd:jed) :: calving_tracer_flux ! tracer flux in solid  runoff (e.g., kg*psu/(m^2 s)  for salinity)
    real, dimension(isd:ied,jsd:jed) :: riverdiffuse        ! sets where to enhance diff_cbt according to rivers
    real, dimension(isd:ied,jsd:jed) :: eta_smooth          ! tendency [tracer*(kg/m^3)*(m/s)] from eta_t smoother
    real, dimension(isd:ied,jsd:jed) :: pbot_smooth         ! tendency [tracer*(kg/m^3)*(m/s)] from pbot_t smoother

    ! The following arrays are never allocated with MOM_STATIC_ARRAYS,
    ! but they need to be declared allocatable in order to compile.
    real, _ALLOCATABLE, dimension(:,:,:)   :: fieldT     _NULL ! Total tracer concentration
    real, _ALLOCATABLE, dimension(:,:,:,:) :: sum_blob   _NULL ! tracer content [concentration*kg] from the L system
    real, _ALLOCATABLE, dimension(:,:,:)   :: tend_blob  _NULL ! blob contribution to dat*th_tendency

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
    character(len=128)                    :: restart_file       ! name for restart file
  end type ocean_prog_tracer_type

  #else
  !#########################
  ! not MOM_STATIC_ARRAYS

  type, private :: ocean_thickness_type
    integer :: method       ! energetic or finite volume

    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzt   _NULL ! E system contribution to rho_dztT (3 time levels)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzten _NULL ! rho*dz at east/north face of T-cell
    real, dimension(:,:,:),   _ALLOCATABLE :: rho_dztr  _NULL ! 1.0/(rho*dzt) at time taup1 (E system)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztL  _NULL ! L system contribution to rho_dztT (3 time levels)
    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dztT  _NULL ! rho(kg/m^3)*thickness (m) of T cell at 3 times (total)

    real, dimension(:,:,:,:), _ALLOCATABLE :: rho_dzu   _NULL ! E system contribution to rho_dzuT (3 time levels)
    real, dimension(:,:,:),   _ALLOCATABLE :: rho_dzur  _NULL ! 1.0/(rho*dzu) at time taup1 (E system)
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

    real, dimension(:,:,:),   _ALLOCATABLE :: dzwt   _NULL ! E system contribution to dzwtT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwtL  _NULL ! L system contribution to dzwtT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwtT  _NULL ! vertical distance (m) between T points at tau/taup1

    real, dimension(:,:,:),   _ALLOCATABLE :: dzwu   _NULL ! E system contribution to dzwuT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwuL  _NULL ! L system contribution to dzwuT
    real, dimension(:,:,:),   _ALLOCATABLE :: dzwuT  _NULL ! vertical distance (m) between U points at tau/taup1

    real, dimension(:,:,:),   _ALLOCATABLE :: dztup  _NULL ! E system contribution to dztupT
    real, dimension(:,:,:),   _ALLOCATABLE :: dztupL _NULL ! L system contribution to dztupT
    real, dimension(:,:,:),   _ALLOCATABLE :: dztupT _NULL ! distance (m) from T-cell point to top of T-cell

    real, dimension(:,:,:),   _ALLOCATABLE :: dztlo  _NULL ! E system contribution to dztloT
    real, dimension(:,:,:),   _ALLOCATABLE :: dztloL _NULL ! L system contribution to dztloT
    real, dimension(:,:,:),   _ALLOCATABLE :: dztloT _NULL ! distance (m) from T-cell point to bottom of T-cell

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

    real, dimension(:,:,:),   _ALLOCATABLE :: dstlo       _NULL ! s-distance (s-units) from T-cell point to top of T-cell
    real, dimension(:,:,:),   _ALLOCATABLE :: dstup       _NULL ! s-distance (s-units) from T-cell point to bottom of T-cell

    real, dimension(:,:),     _ALLOCATABLE :: pbot0       _NULL ! reference bottom pressure (Pa)
    real, dimension(:,:),     _ALLOCATABLE :: pbot0r      _NULL ! inverse of reference bottom pressure (1/Pa)

    real, dimension(:,:,:),   _ALLOCATABLE :: mass_source _NULL ! mass source (kg/m^3)*(m/sec)
    real, dimension(:,:),     _ALLOCATABLE :: blob_source _NULL ! mass source (sink) for the L system (E system)
    real, dimension(:,:,:),   _ALLOCATABLE :: mass_u      _NULL ! mass per area (kg/m^2) in a velocity column
    real, dimension(:,:,:),   _ALLOCATABLE :: mass_en     _NULL ! vertical sum rho_dzten (kg/m^2)
    real, dimension(:,:,:),   _ALLOCATABLE :: mass_uT     _NULL ! mass per area (kg/m^2) in a velocity column
    real, dimension(:,:,:),   _ALLOCATABLE :: thicku      _NULL ! thickness on U-cell (metre) from z=eta_u to z=-H.
    real, dimension(:,:,:),   _ALLOCATABLE :: thicken     _NULL ! sum of dzten (metre)

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
    real,    dimension(:),   _ALLOCATABLE :: dztlo  _NULL ! z-distance (m) from T-cell point to bottom of T-cell
    real,    dimension(:),   _ALLOCATABLE :: dztup  _NULL ! z-distance (m) from T-cell point to top of T-cell
    real,    dimension(:),   _ALLOCATABLE :: dzw    _NULL ! initial vertical resolution of W grid cells (m)
    real,    dimension(:),   _ALLOCATABLE :: dzwr   _NULL ! reciprocal of dzw (W cell vertical resolution)

    real,    dimension(:),   _ALLOCATABLE :: st     _NULL ! s-distance from surface to grid point in level k
    real,    dimension(:),   _ALLOCATABLE :: sw     _NULL ! s-distance from surface down to bottom of level k
    real,    dimension(:),   _ALLOCATABLE :: dst    _NULL ! initial s-vertical resolution of T or U grid cells
    real,    dimension(:),   _ALLOCATABLE :: dstlo  _NULL ! s-distance (s-units) from T-cell point to bottom of T-cell
    real,    dimension(:),   _ALLOCATABLE :: dstup  _NULL ! s-distance (s-units) from T-cell point to top of T-cell
    real,    dimension(:),   _ALLOCATABLE :: dsw    _NULL ! initial s-vertical resolution of W grid cells

    real,    dimension(:,:), _ALLOCATABLE :: fracdz _NULL ! fractional distance between grid point & cell top/bot
    real,    dimension(:,:), _ALLOCATABLE :: ht     _NULL ! depth to bottom of ocean (m) on t-cells
    real,    dimension(:,:), _ALLOCATABLE :: htr    _NULL ! inverse depth to bottom of ocean (m^-1) on t-cells
    real,    dimension(:,:), _ALLOCATABLE :: hu     _NULL ! depth to bottom of ocean (m) on u-cells
    real,    dimension(:,:), _ALLOCATABLE :: dht_dx _NULL ! d(ht)/dx on u-cells (m/m)
    real,    dimension(:,:), _ALLOCATABLE :: dht_dy _NULL ! d(ht)/dy on u-cells (m/m)
    real,    dimension(:,:), _ALLOCATABLE :: gradH  _NULL ! sqrt(dht_dx**2+dht_dyx**2) on u-cells (m/m)

    ! horizontal grid information (time independent)
    integer                            :: ni, nj           ! global points in the two horizontal directions
    real, dimension(:,:), _ALLOCATABLE :: xt         _NULL ! longitude of the T grid points in degrees
    real, dimension(:,:), _ALLOCATABLE :: xu         _NULL ! longitude of the U grid points in degrees
    real, dimension(:,:), _ALLOCATABLE :: yt         _NULL ! latitude of the T grid points in degrees
    real, dimension(:,:), _ALLOCATABLE :: yu         _NULL ! latitude of the U grid points in degrees
    real, dimension(:,:), _ALLOCATABLE :: phiu       _NULL ! latitude of U grid point in radians
    real, dimension(:,:), _ALLOCATABLE :: phit       _NULL ! latitude of T grid point in radians
    real, dimension(:,:), _ALLOCATABLE :: h1t        _NULL ! metric factors in i grid direction
    real, dimension(:,:), _ALLOCATABLE :: h1u        _NULL ! metric factors in j grid jirection
    real, dimension(:,:), _ALLOCATABLE :: h2t        _NULL ! metric factors in i grid direction
    real, dimension(:,:), _ALLOCATABLE :: h2u        _NULL ! metric factors in j grid jirection
    real, dimension(:,:), _ALLOCATABLE :: dh2dx      _NULL ! (1/delta_y)*d(delta_y)/dx (1/m)
    real, dimension(:,:), _ALLOCATABLE :: dh1dy      _NULL ! (1/delta_x)*d(delta_x)/dy (1/m)
    real, dimension(:,:), _ALLOCATABLE :: dxt        _NULL ! longitudinal width of T-cells at grid point (m)
    real, dimension(:,:), _ALLOCATABLE :: dxu        _NULL ! longitudinal width of U-cells at grid point (m)
    real, dimension(:,:), _ALLOCATABLE :: dyt        _NULL ! latitudinal width of T-cells at grid point (m)
    real, dimension(:,:), _ALLOCATABLE :: dyu        _NULL ! latitudinal width of U-cells at grid point (m)
    real, dimension(:,:), _ALLOCATABLE :: dau        _NULL ! area of U-cells (m^2)
    real, dimension(:,:), _ALLOCATABLE :: dat        _NULL ! area of T-cells (m^2)
    real, dimension(:,:), _ALLOCATABLE :: dat_frac   _NULL ! fraction of total wet area occuped by a T-cell (dimensionless)
    real, dimension(:,:), _ALLOCATABLE :: dxte       _NULL ! i-width between i+1 and i points in T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dxtn       _NULL ! i-width of north face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dyte       _NULL ! j-width of east  face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dytn       _NULL ! j-width between j+1 and j points in T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dxue       _NULL ! i-width between i+1 and i points in U-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dxun       _NULL ! i-width of north face of U-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dyue       _NULL ! j-width of east  face of U-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dyun       _NULL ! j-width between j+1 and j points in U-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: datnr      _NULL ! reciprocal area at north face of T-cell
    real, dimension(:,:), _ALLOCATABLE :: dater      _NULL ! reciprocal area at east face of T-cell
    real, dimension(:,:), _ALLOCATABLE :: dun        _NULL ! width from grid point to north face of U-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dus        _NULL ! width from grid point to south face of U-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: duw        _NULL ! width from grid point to west  face of U-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: due        _NULL ! width from grid point to east  face of U-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dtn        _NULL ! width from grid point to north face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dts        _NULL ! width from grid point to south face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dtw        _NULL ! width from grid point to west  face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dte        _NULL ! width from grid point to east  face of T-cells (m)
    real, dimension(:,:), _ALLOCATABLE :: dxtr       _NULL ! 1/dxt
    real, dimension(:,:), _ALLOCATABLE :: dxur       _NULL ! 1/dxu
    real, dimension(:,:), _ALLOCATABLE :: dytr       _NULL ! 1/dyt
    real, dimension(:,:), _ALLOCATABLE :: dyur       _NULL ! 1/dyu
    real, dimension(:,:), _ALLOCATABLE :: daur       _NULL ! 1/[area of U-cells (m^2)]
    real, dimension(:,:), _ALLOCATABLE :: datr       _NULL ! 1/[area of T-cells (m^2)]
    real, dimension(:,:), _ALLOCATABLE :: dxter      _NULL ! 1/dxte
    real, dimension(:,:), _ALLOCATABLE :: dyter      _NULL ! 1/dyte
    real, dimension(:,:), _ALLOCATABLE :: dxtnr      _NULL ! 1/dxtn
    real, dimension(:,:), _ALLOCATABLE :: dytnr      _NULL ! 1/dytn
    real, dimension(:,:), _ALLOCATABLE :: dxuer      _NULL ! 1/dxue
    real, dimension(:,:), _ALLOCATABLE :: dxunr      _NULL ! 1/dxun
    real, dimension(:,:), _ALLOCATABLE :: dyuer      _NULL ! 1/dyue
    real, dimension(:,:), _ALLOCATABLE :: dyunr      _NULL ! 1/dyun
    real, dimension(:,:), _ALLOCATABLE :: dyue_dxuer _NULL ! dyue/dxue
    real, dimension(:,:), _ALLOCATABLE :: dxun_dyunr _NULL ! dxun/dyun
    real, dimension(:,:), _ALLOCATABLE :: dxt_dxter  _NULL ! dxt/dxte
    real, dimension(:,:), _ALLOCATABLE :: dxte_dxtr  _NULL ! dxte/dxt
    real, dimension(:,:), _ALLOCATABLE :: dxt_dxtnr  _NULL ! dxt/dxtn
    real, dimension(:,:), _ALLOCATABLE :: dxtn_dxtr  _NULL ! dxtn/dxt
    real, dimension(:,:), _ALLOCATABLE :: dyt_dyter  _NULL ! dyt/dyte
    real, dimension(:,:), _ALLOCATABLE :: dyte_dytr  _NULL ! dyte/dyt
    real, dimension(:,:), _ALLOCATABLE :: dyt_dytnr  _NULL ! dyt/dytn
    real, dimension(:,:), _ALLOCATABLE :: dytn_dytr  _NULL ! dytn/dyt

    ! land/sea masks
    real, dimension(:,:,:),   _ALLOCATABLE :: mask        _NULL ! land/sea tmask or umask depending on C/B grids
    real, dimension(:,:,:),   _ALLOCATABLE :: tmask       _NULL ! land/sea mask for T cells based on s-coordinate
    real, dimension(:,:,:),   _ALLOCATABLE :: umask       _NULL ! land/sea mask for U cells based on s-coordinate
    real, dimension(:,:,:,:), _ALLOCATABLE :: tmasken     _NULL ! land/sea mask for east/north of t-cell based on s-coordinate
    real, dimension(:,:,:),   _ALLOCATABLE :: tmask_depth _NULL ! based on depth-based vert_coordinate
    real, dimension(:,:,:),   _ALLOCATABLE :: umask_depth _NULL ! based on depth-based vert_coordinate

    ! grid areas and volumes
    real                              :: tcellv          ! initial T cell volume m^3 (entire ocean)
    real                              :: ucellv          ! initial U cell volume m^3 (entire ocean)
    real                              :: tcellsurf       ! T cell surface area (k=1) (bitwise reproducible)
    real                              :: ucellsurf       ! U cell surface area (k=1) (bitwise reproducible)
    real, dimension(:), _ALLOCATABLE  :: tcella   _NULL  ! T cell surface area m^2 (entire ocean)
    real, dimension(:), _ALLOCATABLE  :: ucella   _NULL  ! U cell surface area m^2 (entire ocean)

    ! model grid information
    integer :: wet_t_points   ! total number of wet tracer points
    integer :: wet_u_points   ! total number of wet velocity points
    integer :: total_t_points ! total number of wet or dry tracer points

    ! sine and cosine of rotation angles (clockwise) of velocity for tripolar
    real, dimension(:,:), _ALLOCATABLE :: sin_rot  _NULL
    real, dimension(:,:), _ALLOCATABLE :: cos_rot  _NULL

    ! 1-d grid coordinates for COARDS NetCDF files
    real, dimension(:), _ALLOCATABLE  :: grid_x_t _NULL
    real, dimension(:), _ALLOCATABLE  :: grid_y_t _NULL
    real, dimension(:), _ALLOCATABLE  :: grid_x_u _NULL
    real, dimension(:), _ALLOCATABLE  :: grid_y_u _NULL

    ! axes id for diagnostic manager
    integer, dimension(3)  :: tracer_axes
    integer, dimension(3)  :: vel_axes_uv
    integer, dimension(3)  :: vel_axes_u
    integer, dimension(3)  :: vel_axes_v
    integer, dimension(3)  :: vel_axes_wu
    integer, dimension(3)  :: vel_axes_wt
    integer, dimension(3)  :: tracer_axes_wt
    integer, dimension(3)  :: tracer_axes_flux_x
    integer, dimension(3)  :: tracer_axes_flux_y
    integer, dimension(3)  :: vel_axes_flux_x
    integer, dimension(3)  :: vel_axes_flux_y

    ! axes id for diagnostic manager, appropriate
    ! when with to remap native vertical fields
    ! to depth or pressure levels.  Appropropriate
    ! when vert_coordinate == ZSIGMA or PSIGMA, or other
    ! vertical coordinats whose iso-surfaces are
    ! not quasi-horizontal
    integer, dimension(3)  :: tracer_axes_depth
    integer, dimension(3)  :: vel_axes_uv_depth
    integer, dimension(3)  :: vel_axes_u_depth
    integer, dimension(3)  :: vel_axes_v_depth
    integer, dimension(3)  :: vel_axes_wu_depth
    integer, dimension(3)  :: vel_axes_wt_depth
    integer, dimension(3)  :: tracer_axes_wt_depth
    integer, dimension(3)  :: tracer_axes_flux_x_depth
    integer, dimension(3)  :: tracer_axes_flux_y_depth
    integer, dimension(3)  :: vel_axes_flux_x_depth
    integer, dimension(3)  :: vel_axes_flux_y_depth


  end type ocean_grid_type

  type, private :: ocean_domain_type
    type(domain2d) :: domain2d           ! fms variable, used by mpp routines
    integer :: isc, iec, jsc, jec        ! computational domain indices
    integer :: isd, ied, jsd, jed        ! local indices, consistent with domain2d
    integer :: isg, ieg, jsg, jeg        ! global indices
    integer :: isa, iea, jsa, jea        ! active indices (for minimizing comm2d calls)
    integer :: xhalo, yhalo              ! halo sizes
    integer :: xflags, yflags, layout(2) ! options to mpp_define_domains
    integer :: ioff, joff                ! offset to get absolute indices when MOM_STATIC_ARRAYS (0 otherwise)
    integer :: io_layout(2)              ! options to define io_domain.
    integer :: x_cyclic_offset           ! offset applied to x-direction cyclic boundary condition
    integer :: y_cyclic_offset           ! offset applied to y-direction cyclic boundary condition
    logical, pointer :: maskmap(:,:) =>NULL() ! option to mpp_define_domains
  end type ocean_domain_type

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
    real, _ALLOCATABLE, dimension(:,:,:)  :: uhrho_eu   _NULL ! remapped uhrho_et (kg/(m*s)) onto i-face of U-cell
    real, _ALLOCATABLE, dimension(:,:,:)  :: vhrho_nu   _NULL ! remapped vhrho_nt (kg/(m*s)) onto j-face of U-cell
    real, _ALLOCATABLE, dimension(:,:,:)  :: wrho_bt    _NULL ! rho weight (kg/(m^2*s)) vert advect vel on T-bottom
    real, _ALLOCATABLE, dimension(:,:,:)  :: wrho_bu    _NULL ! remapped wrho_bt onto U-bottom
    real, _ALLOCATABLE, dimension(:,:,:)  :: diverge_t  _NULL ! divergence on T-cell of horiz momentum per mass (kg/m3)*(m/s)
    real, _ALLOCATABLE, dimension(:,:,:)  :: diverge_u  _NULL ! divergence on U-cell of horiz momentum per mass (kg/m3)*(m/s)
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
    integer :: horz_advect_scheme=-1   ! id for horizontal advection scheme
    integer :: vert_advect_scheme=-1   ! id for vertical advection scheme
    integer :: id_obc                  ! id to identify tracer in OBC-subroutines
    integer :: ppm_hlimiter=1          ! Limiter for use with PPM in horizontal
    integer :: ppm_vlimiter=1          ! Limiter for use with PPM in vertical
    integer :: mdt_scheme=4            ! Version of Multi-Dim. Modified Daru & Tenaud (MDMDT)

    real, _ALLOCATABLE, dimension(:,:,:,:) :: field          _NULL ! tracer concentration at 3 time levels (E system)
    real, _ALLOCATABLE, dimension(:,:,:)   :: fieldT         _NULL ! Total tracer concentration
    real, _ALLOCATABLE, dimension(:,:,:)   :: th_tendency    _NULL ! thickness weighted tracer tendency
    real, _ALLOCATABLE, dimension(:,:,:)   :: tendency       _NULL ! for diagnostics: tendency concentration [concentration/sec]
    real, _ALLOCATABLE, dimension(:,:,:)   :: source         _NULL ! tracer source [=tracer concentration per time]
    real, _ALLOCATABLE, dimension(:,:)     :: eta_smooth     _NULL ! tendency from eta_t smoother
    real, _ALLOCATABLE, dimension(:,:)     :: pbot_smooth    _NULL ! tendency from pbot_t smoother

    real, _ALLOCATABLE, dimension(:,:,:)   :: wrk1                _NULL ! work array
    real, _ALLOCATABLE, dimension(:,:,:)   :: tmask_limit         _NULL ! to limit advective and/or neutral physics flux
    real, _ALLOCATABLE, dimension(:,:,:)   :: K33_implicit        _NULL ! m^2/sec vert-diffusivity from neutral diffusion
    real, _ALLOCATABLE, dimension(:,:,:)   :: radiation           _NULL ! radiation absorbed within a cell [W/m^2]
    real, _ALLOCATABLE, dimension(:,:)     :: stf                 _NULL ! surface tracer flux [rho*m/sec*tracer concen]
    real, _ALLOCATABLE, dimension(:,:)     :: btf                 _NULL ! bottom tracer flux [rho*m/sec*tracer concen]
    real, _ALLOCATABLE, dimension(:,:)     :: tpme                _NULL ! tracer concentration in precip-evap
    real, _ALLOCATABLE, dimension(:,:)     :: triver              _NULL ! tracer concentration in river(=runoff+calving)
    real, _ALLOCATABLE, dimension(:,:)     :: trunoff             _NULL ! tracer concentration in river runoff
    real, _ALLOCATABLE, dimension(:,:)     :: tcalving            _NULL ! tracer concentration in calving lang ice
    real, _ALLOCATABLE, dimension(:,:)     :: runoff_tracer_flux  _NULL ! flux in liquid runoff (e.g., kg*degC/(m^2 s) for temp)
    real, _ALLOCATABLE, dimension(:,:)     :: calving_tracer_flux _NULL ! flux in solid runoff (e.g., kg*psu/(m^2 s) for salt)
    real, _ALLOCATABLE, dimension(:,:)     :: riverdiffuse        _NULL ! where to enhance diff_cbt according to rivers
    real, _ALLOCATABLE, dimension(:,:)     :: flux_int            _NULL ! integrated sfc tracer flux for diagnostics
    real, _ALLOCATABLE, dimension(:,:,:,:) :: sum_blob            _NULL ! tracer content [concentration*kg] from the L system
    real, _ALLOCATABLE, dimension(:,:,:)   :: tend_blob           _NULL ! blob contribution to dat*th_tendency


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
    character(len=128)                :: restart_file          ! name for restart file
  end type ocean_prog_tracer_type

  #endif
  !############################################################################################
  ! end of STATIC_MEMORY

  ! ========================================================================================================================




  ! ========================================================================================================================
  ! FROM OCEAN_WORKSPACE_MOD
  !use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4
  #include <ocean_memory.h>

  #ifdef MOM_STATIC_ARRAYS
  real, private, dimension(isd:ied,jsd:jed,nk)   :: wrk1
  real, private, dimension(isd:ied,jsd:jed,nk)   :: wrk2
  real, private, dimension(isd:ied,jsd:jed,nk)   :: wrk3
  real, private, dimension(isd:ied,jsd:jed,nk)   :: wrk4

  #else

  real, private, allocatable, dimension(:,:,:)   :: wrk1
  real, private, allocatable, dimension(:,:,:)   :: wrk2
  real, private, allocatable, dimension(:,:,:)   :: wrk3
  real, private, allocatable, dimension(:,:,:)   :: wrk4

  #endif
  !#######################################################################
  ! <SUBROUTINE NAME="ocean_workspace_init">
  !
  ! <DESCRIPTION>
  ! Initialize MOM workspace module.
  ! </DESCRIPTION>
  !
  subroutine ocean_workspace_init(Grid)

    type(ocean_grid_type), intent(in)   :: Grid

    module_is_initialized = .TRUE.

    #ifndef MOM_STATIC_ARRAYS
    allocate(wrk1(Grid%nk))
    allocate(wrk2(Grid%nk))
    allocate(wrk3(Grid%nk))
    allocate(wrk4(Grid%nk))
    #endif

  end subroutine ocean_workspace_init
  ! </SUBROUTINE> NAME="ocean_workspace_init">

  !#######################################################################
  ! <SUBROUTINE NAME="ocean_workspace_end">
  !
  ! <DESCRIPTION>
  ! End MOM workspace.
  ! </DESCRIPTION>
  !
  subroutine ocean_workspace_end()

    module_is_initialized = .FALSE.

    #ifndef MOM_STATIC_ARRAYS
    deallocate(wrk1,wrk2,wrk3,wrk4,wrk5,wrk6)
    #endif

  end subroutine ocean_workspace_end
  ! </SUBROUTINE> NAME="ocean_workspace_end">
  ! ========================================================================================================================





  ! =========================================================================================================================
  ! FROM OCEAN_TRACER_MOD
  !#######################################################################
  ! <FUNCTION NAME="ocean_prog_tracer_init">
  !
  ! <DESCRIPTION>
  ! Initialization code for prognostic tracers, returning a pointer to
  ! the T_prog array.
  ! </DESCRIPTION>
  !
  function ocean_prog_tracer_init (Grid, Thickness, Ocean_options, Time, Time_steps, &
                                   num_prog, vert_coordinate_type, cmip_units)   &
                                   result (T_prog)

    type(ocean_grid_type),       intent(in), target   :: Grid
    type(ocean_thickness_type),  intent(in)           :: Thickness
    type(ocean_options_type),    intent(inout)        :: Ocean_options
    type(ocean_time_type),       intent(in)           :: Time
    type(ocean_time_steps_type), intent(in)           :: Time_steps
    integer,                     intent(out)          :: num_prog
    integer,                     intent(in)           :: vert_coordinate_type
    logical,                     intent(in)           :: cmip_units

    ! return value
    type(ocean_prog_tracer_type), dimension(:), pointer :: T_prog

    integer               :: i, j, k, n, kb, l
    integer               :: ierr, num_diag
    integer               :: tau, taum1, taup1
    integer               :: ioun, io_status
    integer, dimension(4) :: siz
    integer               :: frazil_heating_order=0
    real                  :: fact
    character(len=32)     :: name
    character(len=128)    :: filename
    logical               :: initialize_as_a_passive_tracer=.false.
    character(len=33)     :: prog_name
    character(len=138)    :: prog_longname

    character(len=48),  parameter :: sub_name = 'ocean_prog_tracer_init'
    character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                    '(' // trim(sub_name) // '): '
    character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                   '(' // trim(sub_name) // '): '
    character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                   '(' // trim(sub_name) // '): '

    ! variables for tracer package
    integer                               :: ind
    real, dimension(2)                    :: range_array
    character(len=64)                     :: caller_str
    character(len=fm_string_len)          :: string_fm
    character(len=fm_type_name_len)       :: typ
    character(len=fm_string_len), pointer, dimension(:) :: good_list
    integer :: stdoutunit, stdlogunit
    stdoutunit=stdout();stdlogunit=stdlog()

    if (prog_module_initialized) then
      call error_mesg(FATAL, trim(error_header) // ' Prognostic tracers already initialized')
    endif

    nullify(T_prog)

    write( stdlogunit,'(/a/)') trim(version)

    dtts          = Time_steps%dtts
    dtuv          = Time_steps%dtuv
    cp_oceanr     = 1.0/cp_ocean

    if(cmip_units) then
        cmip_offset = kelvin
        temp_units='degrees K'
    else
        cmip_offset = 0.0
        temp_units='degrees C'
    endif

    if(zero_tendency) then
        call error_mesg(NOTE, trim(note_header) // ' zero_tendency=true so will not time step tracer fields.')
        Ocean_options%tracer_tendency = 'Did NOT time step prognostic tracer fields.'
    else
        Ocean_options%tracer_tendency = 'Time stepped the prognostic tracer fields.'
    endif

    if(zero_tracer_source) then
        call error_mesg(NOTE, &
        trim(note_header) // ' zero_tracer_source=true so remove T_prog%source from evolution.')
    endif

    if(.not. write_a_restart) then
      write(stdoutunit,'(a)') '==>Note: running ocean_tracer with write_a_restart=.false.'
      write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
    endif

    if(frazil_heating_before_vphysics) then
        frazil_heating_order=frazil_heating_order+1
        write(stdoutunit,'(a)') '==>Note: frazil heating called before vertical physics and before boundary fluxes.'
        write(stdoutunit,'(a)') '         This method is retained for legacy purposes: it is NOT recommended for new runs. '
        write(stdoutunit,'(a)') ' '
    endif
    if(frazil_heating_after_vphysics) then
        frazil_heating_order=frazil_heating_order+1
        write(stdoutunit,'(a)') '==>Note: frazil heating called after vertical physics and after boundary fluxes.'
        write(stdoutunit,'(a)') '         This is the recommended method. '
        write(stdoutunit,'(a)') ' '
    endif
    if(frazil_heating_order>1) then
        write(stdoutunit,'(a)') '==>Error from ocean_tracer_mod: choose just one temporal order for frazil heating.'
        write(stdoutunit,'(a)') ' '
        call error_mesg(FATAL, &
        trim(note_header) // ' can only choose one temporal order for frazil heating.')
    endif
    if(frazil_heating_order==0) then
        write(stdoutunit,'(a)') '==>Error from ocean_tracer_mod: MUST specify order for frazil heating in ocean_tracer.'
        write(stdoutunit,'(a)') ' '
        call error_mesg(FATAL, &
        trim(note_header) // '  Need to specify an order for calling frazil heating from ocean_tracer_mod.')
    endif

    ! some time step information
    if (dtts /= dtuv) then
       write (stdoutunit,'(/a)') trim(warn_header) // ' Asynchronous timesteps (dtts > dtuv) imply inaccurate transients.'
       write (stdoutunit,'(a)')  '          and total tracer (i.e. heat content) is not conserved.'
       write (stdoutunit,'(a,f5.2,a/)') '            dtts =',dtts/dtuv,' times larger than dtuv.'
    else
       call error_mesg(NOTE, trim(note_header) // ' Synchronous timesteps have been specified (dtts = dtuv).')
    endif

    tendency   = Time_steps%tendency
    dtime      = Time_steps%dtime_t
    dtimer     = 1.0/(dtime+epsln)

    tau   = Time%tau
    taum1 = Time%taum1
    taup1 = Time%taup1

    ! call routine to initialize the required tracer package
    if (otpm_set_tracer_package('required', caller=trim(mod_name)//'('//trim(sub_name)//')') .le. 0) then
      call error_mesg(FATAL, trim(error_header) // ' Could not set required packages list')
    endif

    ! call the initialization routine for the ocean tracer packages
    !call ocean_tpm_init(Grid, Time, Time_steps, Ocean_options)
    call ocean_passive_init (Domain, Grid, Ocean_options)

    ! get the number of tracers
    ! for now, we will not try to dynamically allocate the tracer arrays here
    num_prog_tracers = fm_get_length('/ocean_mod/prog_tracers')
    num_prog=num_prog_tracers

    ! allocate arrays based on the number of prognostic tracers
    allocate( T_prog            (num_prog_tracers) )
    allocate( id_eta_smooth     (num_prog_tracers) )
    allocate( id_pbot_smooth    (num_prog_tracers) )
    allocate( id_prog           (num_prog_tracers) )
    allocate( id_prog_explicit  (num_prog_tracers) )
    allocate( id_prog_rhodzt    (num_prog_tracers) )
    allocate( id_prog_int_rhodz (num_prog_tracers) )
    allocate( id_prog_on_depth  (num_prog_tracers) )
    allocate( id_tendency_conc  (num_prog_tracers) )
    allocate( id_tendency       (num_prog_tracers) )
    allocate( id_tendency_expl  (num_prog_tracers) )
    allocate( id_surf_tracer    (num_prog_tracers) )
    allocate( id_surf_tracer_sq (num_prog_tracers) )
    allocate( id_bott_tracer    (num_prog_tracers) )

    id_eta_smooth(:)     = -1
    id_pbot_smooth(:)    = -1
    id_prog(:)           = -1
    id_prog_explicit(:)  = -1
    id_prog_rhodzt(:)    = -1
    id_prog_int_rhodz(:) = -1
    id_prog_on_depth(:)  = -1
    id_tendency_conc(:)  = -1
    id_tendency(:)       = -1
    id_tendency_expl(:)  = -1
    id_surf_tracer(:)    = -1
    id_surf_tracer_sq(:) = -1
    id_bott_tracer(:)    = -1

    ! set logical to determine when to mpp_update tracers
    do n=1,num_prog_tracers-1
      T_prog(n)%complete=.false.
    enddo
    T_prog(num_prog_tracers)%complete=.true.

    !  dump the lists for the tracer packages
    write (stdoutunit,*)
    write (stdoutunit,*) 'Dumping tracer_packages tracer tree'
    if (.not. fm_dump_list('/ocean_mod/tracer_packages', recursive = .true.)) then
      call error_mesg(FATAL, trim(error_header) // ' Problem dumping tracer_packages tracer tree')
    endif

    write (stdoutunit,*)
    write (stdoutunit,*) 'Dumping prog_tracers tracer tree'
    if (.not. fm_dump_list('/ocean_mod/prog_tracers', recursive = .true.)) then
      call error_mesg(FATAL, trim(error_header) // ' Problem dumping prog_tracers tracer tree')
    endif

    write (stdoutunit,*)
    write (stdoutunit,*) 'Dumping namelists tracer tree'
    if (.not. fm_dump_list('/ocean_mod/namelists', recursive = .true.)) then
      call error_mesg(FATAL, trim(error_header) // ' Problem dumping namelists tracer tree')
    endif

    !### finished with initializing t_prog arrays
    !### finished with call to tracer startup routine


    ! set local array indices
    Grd => Grid

  #ifndef MOM_STATIC_ARRAYS
    call get_local_indices(isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk=Grd%nk
  #endif

    do n=1,num_prog_tracers
  #ifndef MOM_STATIC_ARRAYS
      allocate( T_prog(n)%field(isd:ied,jsd:jed,nk,3))
      allocate( T_prog(n)%th_tendency(isd:ied,jsd:jed,nk))
      allocate( T_prog(n)%tendency(isd:ied,jsd:jed,nk))
      allocate( T_prog(n)%source(isd:ied,jsd:jed,nk))
      allocate( T_prog(n)%eta_smooth(isd:ied,jsd:jed))
      allocate( T_prog(n)%pbot_smooth(isd:ied,jsd:jed))
      allocate( T_prog(n)%wrk1(isd:ied,jsd:jed,nk))
      allocate( T_prog(n)%tmask_limit(isd:ied,jsd:jed,nk))
      allocate( T_prog(n)%K33_implicit(isd:ied,jsd:jed,nk))
  #endif
      T_prog(n)%field(:,:,:,:)        = 0.0
      T_prog(n)%th_tendency(:,:,:)    = 0.0
      T_prog(n)%tendency(:,:,:)       = 0.0
      T_prog(n)%source(:,:,:)         = 0.0
      T_prog(n)%eta_smooth(:,:)       = 0.0
      T_prog(n)%pbot_smooth(:,:)      = 0.0
      T_prog(n)%wrk1(:,:,:)           = 0.0
      T_prog(n)%tmask_limit(:,:,:)    = 0.0
      T_prog(n)%K33_implicit(:,:,:)   = 0.0
      T_prog(n)%neutral_physics_limit = .false.
    enddo

    ! for global surface area normalization
    cellarea_r = 1.0/(epsln + Grd%tcellsurf)

    ! fill the field table entries for the prognostic tracers
    n = 0
    do while (fm_loop_over_list('/ocean_mod/prog_tracers', name, typ, ind))  !{

       if (typ .ne. 'list') then  !{
           call error_mesg(FATAL, trim(error_header) // ' ' // trim(name) // ' is not a list')
       else  !}{
           n = n + 1  ! increment the array index

           if (n .ne. ind) then  !{
               write (stdoutunit,*) trim(warn_header), ' Tracer index, ', ind,   &
                    ' does not match array index, ', n, ' for ', trim(name)
           endif  !}

           ! save the name
           T_prog(n)%name = name
           if (.not. fm_change_list('/ocean_mod/prog_tracers/' // trim(name))) then  !{
               call error_mesg(FATAL, trim(error_header) // ' Problem changing to ' // trim(name))
           endif  !}
           caller_str = 'ocean_tracer_mod(ocean_prog_tracer_init)'

           ! save the units
           T_prog(n)%units = fm_util_get_string('units', caller = caller_str, scalar = .true.)

           ! save the type
           T_prog(n)%type = fm_util_get_string('type', caller = caller_str, scalar = .true.)

           ! save the longname
           T_prog(n)%longname = fm_util_get_string('longname', caller = caller_str, scalar = .true.)

           ! save the conversion
           T_prog(n)%conversion = fm_util_get_real('conversion', caller = caller_str, scalar = .true.)

           ! save the offset
           T_prog(n)%offset = fm_util_get_real('offset', caller = caller_str, scalar = .true.)

           ! get the min and max of the tracer
           T_prog(n)%min_tracer = fm_util_get_real('min_tracer', caller = caller_str, scalar = .true.)
           T_prog(n)%max_tracer = fm_util_get_real('max_tracer', caller = caller_str, scalar = .true.)

           ! get the min and max of the range for analysis
           T_prog(n)%min_range = fm_util_get_real('min_range', caller = caller_str, scalar = .true.)
           T_prog(n)%max_range = fm_util_get_real('max_range', caller = caller_str, scalar = .true.)

           ! get the flux unit
           T_prog(n)%flux_units = fm_util_get_string('flux_units', caller = caller_str, scalar = .true.)

           ! get the min and max of the flux range for analysis
           T_prog(n)%min_flux_range = fm_util_get_real('min_flux_range', caller = caller_str, scalar = .true.)
           T_prog(n)%max_flux_range = fm_util_get_real('max_flux_range', caller = caller_str, scalar = .true.)

           ! save the restart file
           T_prog(n)%restart_file = fm_util_get_string('restart_file', caller = caller_str, scalar = .true.)

           ! save flag for whether the tracer must have a value in the restart file
           T_prog(n)%const_init_tracer = fm_util_get_logical('const_init_tracer', caller = caller_str, scalar = .true.)

           ! save value to globally initialize this tracer (optional)
           T_prog(n)%const_init_value = fm_util_get_real('const_init_value', caller = caller_str, scalar = .true.)

           ! get the horizontal-advection-scheme
           string_fm = fm_util_get_string('horizontal-advection-scheme', caller = caller_str, scalar = .true.)

           select case (trim(string_fm))
           case ('upwind')
               T_prog(n)%horz_advect_scheme = ADVECT_UPWIND
           case ('quicker')
               T_prog(n)%horz_advect_scheme = ADVECT_QUICKER
           case default
            call error_mesg(FATAL, trim(error_header) // ' Invalid horz-advect-scheme '  // trim(string_fm) // ' for ' // trim(name))
           end select

           ! get the vertical-advection-scheme
           string_fm = fm_util_get_string('vertical-advection-scheme', caller = caller_str, scalar = .true.)

           select case (trim(string_fm))
           case ('upwind')
               T_prog(n)%vert_advect_scheme = ADVECT_UPWIND
           case ('quicker')
               T_prog(n)%vert_advect_scheme = ADVECT_QUICKER
           case default
             call error_mesg(FATAL, trim(error_header) // ' Invalid vert-advect scheme ' // trim(string_fm) // ' for ' // trim(name))
           end select

           ! save the max_tracer_limit
           T_prog(n)%max_tracer_limit = fm_util_get_real('max_tracer_limit', caller = caller_str, scalar = .true.)

           ! save the min_tracer_limit
           T_prog(n)%min_tracer_limit = fm_util_get_real('min_tracer_limit', caller = caller_str, scalar = .true.)

           ! save flag for whether to transport tracer using only advection
           T_prog(n)%use_only_advection = fm_util_get_logical('use_only_advection', caller = caller_str, scalar = .true.)
           if(T_prog(n)%use_only_advection) then
               call error_mesg(NOTE, &
               trim(note_header) // ' Will evolve '// trim(T_prog(n)%name) // ' via advection alone.')
           endif

       endif  !}
    enddo  !}


    allocate(id_tmask_limit(num_prog_tracers))

  #ifdef USE_OCEAN_BGC
    !Get the %filed for "generic" tracers as it might have already been set.
    !nnz: find a way to use their already allocated field pointer directly.
    do n=1,num_prog_tracers
      if(T_prog(n)%type .eq. 'generic') then
         call ocean_generic_get_field(T_prog(n)%name,T_prog(n)%field)
      endif
    enddo
  #endif

    ! read prognostic tracer initial conditions or restarts
    write (stdoutunit,*) ' '
    write (stdoutunit,*) trim(note_header), &
     ' Reading prognostic tracer initial conditions or restarts'

    do n=1,num_prog_tracers  !{

      write (stdoutunit,*)
      write (stdoutunit,*) &
      'Initializing tracer number', n,' at time level tau. This tracer is called ',trim(T_prog(n)%name)

      if(.not. Time%init) then
          if(tendency==TWO_LEVEL) then
              write (stdoutunit,'(/a)')'Expecting only one time record from the tracer restart.'
          elseif(tendency==THREE_LEVEL) then
              write (stdoutunit,'(/a)')'Expecting two time records from the tracer restart.'
          endif
      endif

      T_prog(n)%field(:,:,:,taum1) = 0.0

      ! initialization logic
      !
      ! 1. Always call passive_tracer_init once.
      !    If Time%init=.false. then return w/o doing anything
      !    If not using passive tracers, then return w/o doing anything
      !
      ! 2. If filename exists with field in it, then fill the tracer with field in filename.
      !
      ! 3. If filename exists but there is no field in it, then fill the tracer with const_init_value
      !    if const_init_tracer is true, otherwise end with a fatal error.
      !
      ! 4. If filename does not exist, then fill the tracer with const_init_value
      !    if const_init_tracer is true, otherwise end with a fatal error.
      !

      call passive_tracer_init(Time%init, T_prog(n), initialize_as_a_passive_tracer)

      call field_size(filename,T_prog(n)%name, siz)

      ! initialize tracer at taup1 to tracer at tau
        T_prog(n)%field(:,:,:,taup1) = T_prog(n)%field(:,:,:,tau)

    enddo  !} n-loop

    write (stdoutunit,*) ' '
    write (stdoutunit,*) trim(note_header), ' finished reading prognostic tracer restarts.'

    prog_module_initialized = .true.

  end function ocean_prog_tracer_init
  ! </FUNCTION> NAME="ocean_prog_tracer_init">
 !=========================================================================================================================


 !#######################################################################
 ! <SUBROUTINE NAME="ocean_passive_init">
 !
 ! <DESCRIPTION>
 ! Initialize the indices for passive tracer fields.
 ! This routine is called by ocean_model.F90.
 ! </DESCRIPTION>
 !
 subroutine ocean_passive_init(Domain, Grid, Ocean_options, debug)

   type(ocean_domain_type),  intent(in), target   :: Domain
   type(ocean_grid_type),    intent(in), target   :: Grid
   type(ocean_options_type), intent(inout)        :: Ocean_options
   logical,                  intent(in), optional :: debug
   integer                                        :: ioun, io_status, ierr

 ! -------------------------------------------------------------------------
 ! start of tracer package material

   ! local parameters
   character(len=64), parameter :: sub_name = 'ocean_passive_init'
   character(len=256), parameter   :: error_header =                               &
      '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
   character(len=256), parameter   :: note_header =                                &
      '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   ! local variables
   integer                                             :: n
   character(len=fm_field_name_len)                    :: name
   character(len=fm_path_name_len)                     :: path_to_names
   character(len=fm_field_name_len+1)                  :: suffix
   character(len=fm_field_name_len+3)                  :: long_suffix
   character(len=256)                                  :: caller_str
   character(len=fm_string_len), pointer, dimension(:) :: good_list

   integer :: stdoutunit,stdlogunit
   stdoutunit=stdout();stdlogunit=stdlog()

   if ( module_is_initialized ) then
       call mpp_error(FATAL, '==>Error in ocean_passive_mod (ocean_passive_init): module initialized.')
   endif
   module_is_initialized=.true.

   ! initialize the ocean passive tracer package and set some defaults.
   ! each of these defaults can be altered via entries to the field_table.
   package_index = otpm_set_tracer_package(package_name,                             &
      units = 'dimensionless', min_tracer=t_min, max_tracer=t_max,                   &
      min_tracer_limit=t_min_limit, max_tracer_limit=t_max_limit,                    &
      restart_file=default_restart_file,                                             &
      caller = trim(mod_name) // '(' // trim(sub_name) // ')',                       &
      conversion=1.0, offset=0.0, min_range=-10.0, max_range=100.0,                  &
      flux_units = 'dimensionless', min_flux_range=-1.0e+16, max_flux_range=1.0e+16, &
      vert_adv_scheme='mdppm', horiz_adv_scheme='mdppm', psom_limit=.true.)

   ! check for number of entries in the field_table for passive tracers.
   path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
   instances = fm_get_length(path_to_names)
   if (instances < 0) then
     call mpp_error( &
     FATAL, trim(error_header) // ' Could not get number of passive tracer instances.')
   endif

   ! determine whether to use this module.
   write (stdoutunit,*) ' '
   if(instances==0) then
     write (stdoutunit,*) trim(note_header), ' No instances of passive tracers in field_table.'
     use_this_module=.false.
     write(stdoutunit,'(a)') ' '
     write(stdoutunit,'(a)') &
     '==>Note: NOT running with idealized passive tracers.'
     call mpp_error(NOTE, '==>Note: ocean_passive_mod: NOT using idealized passive tracer module.')
     Ocean_options%passive_tracers = 'Did NOT use the idealized passive tracer module.'
     return
   else
     write (stdoutunit,*) trim(note_header), ' ', instances, ' instances of passive tracers in field_table.'
     use_this_module=.true.
     Ocean_options%passive_tracers = 'Ran with idealized passive tracer module.'
   endif

   ! allocate the passive array
   allocate (passive(instances))

   ! loop over the names of the passive tracers, saving them into the passive array
   do n=1,instances
      if (fm_get_value(path_to_names, name, index = n)) then
          passive(n)%name = name
      else
          write (name,*) n
          call mpp_error(FATAL, trim(error_header) //        &
               ' Bad field name for index ' // trim(name))
      endif
   enddo

   ! determine the tracer name for this instance,
   ! and set the prog tracer through otpm_set_prog_tracer
   do n=1,instances
      name = passive(n)%name
      if (name(1:1) .eq. '_') then
          suffix = ' '
          long_suffix = ' '
      else
          suffix = '_' // name
          long_suffix = ' (' // trim(name) // ')'
      endif
      passive(n)%index = otpm_set_prog_tracer('passive' // trim(suffix), package_name, &
           longname = 'passive' // trim(long_suffix),                                  &
           caller = trim(mod_name)//'('//trim(sub_name)//')')
   enddo

   ! add the package name to the list of good namelists; used later for consistency check
   if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then
     call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
   endif

   ! set defaults for the available namelist parameters
   caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
   do n=1,instances
      call fm_util_start_namelist(package_name, passive(n)%name, caller=caller_str, &
                                  no_overwrite=.true., check=.true.)
      call fm_util_set_value('init_condition', 'patch')
      call fm_util_set_value('init_surface', 1030.0)
      call fm_util_set_value('init_value', 1.0)
      call fm_util_set_value('restore', .false.)
      call fm_util_end_namelist(package_name, passive(n)%name, check = .true., caller = caller_str)
   enddo

   ! check for errors in number of fields in the namelists for this package
   good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
        caller = trim(mod_name) // '(' // trim(sub_name) // ')')
   if (associated(good_list)) then
       call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
            caller = trim(mod_name) // '(' // trim(sub_name) // ')')
       deallocate(good_list)
   else
       call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
   endif

 ! end of tracer package material
 ! -------------------------------------------------------------------------

 #ifndef MOM_STATIC_ARRAYS
   call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
   nk = Grid%nk
 #endif

   Dom => Domain
   Grd => Grid

   call write_version_number()

   ! provide for namelist override of defaults
 #ifdef INTERNAL_FILE_NML
   read (input_nml_file, nml=ocean_passive_nml, iostat=io_status)
   ierr = check_nml_error(io_status,'ocean_passive_nml')
 #else
   ioun = open_namelist_file()
   read (ioun,ocean_passive_nml,IOSTAT=io_status)
   ierr = check_nml_error(io_status, 'ocean_passive_nml')
   call close_file(ioun)
 #endif
   write (stdoutunit,'(/)')
   write (stdoutunit,ocean_passive_nml)
   write (stdlogunit,ocean_passive_nml)

   if (PRESENT(debug) .and. .not. debug_this_module) then
       debug_this_module = debug
   endif
   if(debug_this_module) then
       write(stdoutunit,'(a)') '==>Note: running ocean_passive with debug_this_module=.true.'
   endif

   ! set default initial condition according to setting in ocean_passive_nml
   do n=1,instances
      passive(n)%init_condition = trim(common_init_condition)
   enddo

   ! namelist information for each passive tracer from the field_table.
   ! in particular, can override the common_init_condition so that each
   ! tracer can, for example, have a distinct initial condition.
   caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
   do n=1,instances
      call fm_util_start_namelist(package_name, passive(n)%name, caller = caller_str)
      passive(n)%init_condition =  fm_util_get_string ('init_condition', scalar = .true.)
      passive(n)%init_surface   =  fm_util_get_real('init_surface', scalar = .true.)
      passive(n)%init_value     =  fm_util_get_real('init_value', scalar = .true.)
      passive(n)%restore        =  fm_util_get_logical('restore', scalar = .true.)
      call fm_util_end_namelist(package_name, passive(n)%name, caller = caller_str)
   enddo
   do n=1,instances
      if (passive(n)%restore ) then
         name = passive(n)%name
         if (name(1:1) .eq. '_') then
            suffix = ' '
            long_suffix = ' '
         else
            suffix = '_' // name
            long_suffix = ' (' // trim(name) // ')'
         endif
         passive(n)%diag_index = otpm_set_diag_tracer('passive' // trim(suffix) // 'mask',         &
            longname = 'passive' // trim(long_suffix), restart_file='ocean_passive_masks.res.nc',  &
            min_range=-1.0, max_range=2.0, const_init_tracer=.true.,const_init_value=0.0,          &
            caller = caller_str )
         use_tracer_restore=.true.
      endif
   enddo


 end subroutine ocean_passive_init
 ! </SUBROUTINE> NAME="ocean_passive_init"






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

  !for restart file
  type(restart_file_type), save :: Adv_restart

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

  private  horz_advect_tracer
  private horz_advect_tracer_upwind

  private  vert_advect_tracer
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

    module_is_initialized = .true.

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

    if ( .not. module_is_initialized ) then
      write(stdoutunit,*)'==>Error from ocean_tracer_advect (horz_advect_tracer_upwind): ocean_tracer_advect_mod not yet initialized')
    endif

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
