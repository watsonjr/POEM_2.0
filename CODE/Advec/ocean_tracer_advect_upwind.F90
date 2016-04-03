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
!use FMS_mod,    only: check_nml_error, write_version_number, error_mesg
!use FMS_mod,    only: open_namelist_file, close_file
!-----------------------------------------------------------------------
!
!         A collection of commonly used routines.
!
!  The routines are primarily I/O related, however, there also
!  exists several simple miscellaneous utility routines.
!
!-----------------------------------------------------------------------
!
!  check_nml_error    Checks the iostat argument that is returned after
!                     reading a namelist and determines if the error
!                     code is valid.
!
!  write_version_number  Prints to the log file (or a specified unit)
!                        the (cvs) version id string and (cvs) tag name.
!
!  error_mesg          Print notes, warnings and error messages,
!                      terminates program for error messages.
!                      (use error levels NOTE,WARNING,FATAL)
!
!  open_namelist_file  Opens namelist file for reading only.
!
!  close_file          Closes a file that was opened using
!                      open_namelist_file, open_restart_file, or
!                      open_ieee32_file.
!
!-----------------------------------------------------------------------
! routines for opening/closing specific types of file
public :: open_namelist_file, close_file

! miscellaneous i/o routines
public :: check_nml_error, write_version_number, error_mesg

! Namelist read error values
  TYPE nml_errors_type
     INTEGER :: multipleNMLSinFile
     INTEGER :: badType1
     INTEGER :: badType2
     INTEGER :: missingVar
  END TYPE nml_errors_type
  TYPE(nml_errors_type), SAVE :: nml_errors


!------ namelist interface -------
!------ adjustable severity level for warnings ------

  logical           :: read_all_pe   = .true.
  character(len=16) :: clock_grain = 'NONE', clock_flags='NONE'
  character(len=8)  :: warning_level = 'warning'
  character(len=64) :: iospec_ieee32 = '-N ieee_32'
  integer           :: stack_size = 0
  integer           :: domains_stack_size = 0
  logical, public   :: print_memory_usage = .FALSE.

!------ namelist interface -------

! <NAMELIST NAME="fms_nml">
!   <DATA NAME="clock_grain"  TYPE="character"  DEFAULT="'NONE'">
!     The level of clock granularity used for performance timing sections
!     of code. Possible values in order of increasing detail are:
!     'NONE', 'COMPONENT', 'SUBCOMPONENT', 'MODULE_DRIVER', 'MODULE', 'ROUTINE',
!     'LOOP', and 'INFRA'.  Code sections are defined using routines in MPP
!     module: mpp_clock_id, mpp_clock_begin, and mpp_clock_end.
!     The fms module makes these routines public.
!     A list of timed code sections will be printed to STDOUT.
!     See the <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp.html">MPP</LINK>
!     module for more details.
!   </DATA>
!   <DATA NAME="clock_flags"  TYPE="character"  DEFAULT="'NONE'">
!     Possible values are 'NONE', 'SYNC', or 'DETAILED'.
!     SYNC will give accurate information on load balance of the clocked
!     portion of code.
!     DETAILED also turns on detailed message-passing performance diagnosis.
!     Both SYNC and DETAILED will  work correctly on innermost clock nest
!     and distort outer clocks, and possibly the overall code time.
!     See the <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp.html">MPP</LINK>
!     module for more details.
!   </DATA>
!   <DATA NAME="read_all_pe"  TYPE="logical"  DEFAULT="true">
!     Read global data on all processors extracting local part needed (TRUE) or
!     read global data on PE0 and broadcast to all PEs (FALSE).
!   </DATA>
!   <DATA NAME="warning_level"  TYPE="character"  DEFAULT="'warning'">
!     Sets the termination condition for the WARNING flag to interfaces
!     error_mesg/mpp_error. set warning_level = 'fatal' (program crashes for
!     warning messages) or 'warning' (prints warning message and continues).
!   </DATA>
!   <DATA NAME="iospec_ieee32"  TYPE="character"  DEFAULT="'-N ieee_32'">
!     iospec flag used with the open_ieee32_file interface.
!   </DATA>
!   <DATA NAME="stack_size"  TYPE="integer"  DEFAULT="0">
!     The size in words of the MPP user stack. If stack_size > 0, the following
!     MPP routine is called: call mpp_set_stack_size (stack_size). If stack_size
!     = 0 (default) then the default size set by mpp_mod is used.
!   </DATA>
!   <DATA NAME="domains_stack_size" TYPE="integer"  DEFAULT="0">
!     The size in words of the MPP_DOMAINS user stack. If
!     domains_stack_size > 0, the following MPP_DOMAINS routine is called:
!     call mpp_domains_set_stack_size (domains_stack_size). If
!     domains_stack_size = 0 (default) then the default size set by
!     mpp_domains_mod is used.
!   </DATA>
!   <DATA NAME="print_memory_usage"  TYPE="logical"  DEFAULT=".FALSE.">
!     If set to .TRUE., memory usage statistics will be printed at various
!     points in the code. It is used to study memory usage, e.g to detect
!     memory leaks.
!   </DATA>
! </NAMELIST>

  namelist /fms_nml/  read_all_pe, clock_grain, clock_flags,    &
                      warning_level, iospec_ieee32, &
                      stack_size, domains_stack_size, &
                      print_memory_usage

!   ---- private data for check_nml_error ----

   integer, private :: num_nml_error_codes, nml_error_codes(20)
   logical, private :: do_nml_error_init = .true.
   private  nml_error_init


!  ---- version number -----

  character(len=128) :: version = '$Id: fms.F90,v 19.0.6.1 2013/02/25 18:32:54 Zhi.Liang Exp $'
  character(len=128) :: tagname = '$Name: siena_201303 $'

  logical :: module_is_initialized = .FALSE.


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
!---- initialize mpp routines ----
    if(present(localcomm)) then
       call mpp_init(localcomm=localcomm)
    else
       call mpp_init()
    endif
    call mpp_domains_init
    call fms_io_init

!---- read namelist input ----

    call nml_error_init  ! first initialize namelist iostat error codes

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, fms_nml, iostat=io)
      ierr = check_nml_error(io,'fms_nml')
#else
    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
          read  (unit, nml=fms_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'fms_nml')  ! also initializes nml error codes
       enddo
 10    call mpp_close (unit)
    endif
#endif

!---- define mpp stack sizes if non-zero -----

    if (        stack_size > 0) call         mpp_set_stack_size (        stack_size)
    if (domains_stack_size > 0) call mpp_domains_set_stack_size (domains_stack_size)

!---- set severity level for warnings ----

    select case( trim(lowercase(warning_level)) )
    case( 'fatal' )
        call mpp_set_warn_level ( FATAL )
    case( 'warning' )
        call mpp_set_warn_level ( WARNING )
    case default
        call error_mesg ( 'fms_init',  &
             'invalid entry for namelist variable warning_level', FATAL )
    end select

!--- set granularity for timing code sections ---

    select case( trim(uppercase(clock_grain)) )
    case( 'NONE' )
        call mpp_clock_set_grain (0)
    case( 'COMPONENT' )
        call mpp_clock_set_grain (CLOCK_COMPONENT)
    case( 'SUBCOMPONENT' )
        call mpp_clock_set_grain (CLOCK_SUBCOMPONENT)
    case( 'MODULE_DRIVER' )
        call mpp_clock_set_grain (CLOCK_MODULE_DRIVER)
    case( 'MODULE' )
        call mpp_clock_set_grain (CLOCK_MODULE)
    case( 'ROUTINE' )
        call mpp_clock_set_grain (CLOCK_ROUTINE)
    case( 'LOOP' )
        call mpp_clock_set_grain (CLOCK_LOOP)
    case( 'INFRA' )
        call mpp_clock_set_grain (CLOCK_INFRA)
    case default
        call error_mesg ( 'fms_init',  &
             'invalid entry for namelist variable clock_grain', FATAL )
    end select
!Balaji
    select case( trim(uppercase(clock_flags)) )
    case( 'NONE' )
       clock_flag_default = 0
    case( 'SYNC' )
       clock_flag_default = MPP_CLOCK_SYNC
    case( 'DETAILED' )
       clock_flag_default = MPP_CLOCK_DETAILED
    case default
       call error_mesg ( 'fms_init',  &
            'invalid entry for namelist variable clock_flags', FATAL )
   end select

!--- write version info and namelist to logfile ---

    call write_version_number (version, tagname)
    if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog()
      write (unit, nml=fms_nml)
      write (unit,*) 'nml_error_codes=', nml_error_codes(1:num_nml_error_codes)
    endif

    call memutils_init( print_memory_usage )
    call print_memuse_stats('fms_init')

    call write_version_number (constants_version,constants_tagname)

end subroutine fms_init
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="error_mesg">

!   <OVERVIEW>
!     Print notes, warnings and error messages; terminates program for warning
!     and error messages. (use error levels NOTE,WARNING,FATAL, see example below)
!   </OVERVIEW>
!   <DESCRIPTION>
!     Print notes, warnings and error messages; and terminates the program for
!     error messages. This routine is a wrapper around mpp_error, and is provided
!     for backward compatibility. This module also publishes mpp_error,
!      <B>users should try to use the mpp_error interface</B>.
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
! users should try to use the mpp_error interface

 subroutine error_mesg (routine, message, level)
  character(len=*), intent(in) :: routine, message
  integer,          intent(in) :: level

!  input:
!      routine   name of the calling routine (character string)
!      message   message written to output   (character string)
!      level     set to NOTE, MESSAGE, or FATAL (integer)

    if (.not.module_is_initialized) call fms_init ( )
    call mpp_error ( routine, message, level )

 end subroutine error_mesg
! </SUBROUTINE>


!#######################################################################
! <FUNCTION NAME="check_nml_error">

!   <OVERVIEW>
!     Checks the iostat argument that is returned after reading a namelist
!     and determines if the error code is valid.
!   </OVERVIEW>
!   <DESCRIPTION>
!     The FMS allows multiple namelist records to reside in the same file.
!     Use this interface to check the iostat argument that is returned after
!     reading a record from the namelist file. If an invalid iostat value
!     is detected this routine will produce a fatal error. See the NOTE below.
!   </DESCRIPTION>
!   <TEMPLATE>
!     check_nml_error ( iostat, nml_name )
!   </TEMPLATE>

!   <IN NAME="iostat"  TYPE="integer" >
!     The iostat value returned when reading a namelist record.
!   </IN>
!   <IN NAME="nml_name"  TYPE="character" >
!     The name of the namelist. This name will be printed if an error is
!     encountered, otherwise the name is not used.
!   </IN>
!   <OUT NAME=""  TYPE="integer" >
!     This function returns the input iostat value (integer) if it is an
!     allowable error code. If the iostat error code is not
!     allowable, an error message is printed and the program terminated.
!   </OUT>
!   <NOTE>
!     Some compilers will return non-zero iostat values when reading through
!     files with multiple namelist. This routine
!     will try skip these errors and only terminate for true namelist errors.
!
!     Examples
!
!       The following example checks if a file exists, reads a namelist input
!       from that file, and checks for errors in that
!       namelist. When the correct namelist is read and it has no errors the
!       routine check_nml_error will return zero and the while loop will exit.
!       This code segment should be used to read namelist files.
!       <PRE>
!          integer :: unit, ierr, io
!
!          if ( file_exist('input.nml') ) then
!              unit = open_namelist_file ( )
!              ierr=1
!              do while (ierr > 0)
!                read  (unit, nml=moist_processes_nml, iostat=io)
!                ierr = check_nml_error(io,'moist_processes_nml')
!              enddo
!              call close_file (unit)
!          endif
!       </PRE>
!   </NOTE>

!   <ERROR MSG="Unknown error while reading namelist ...., (IOSTAT = ####)" STATUS="FATAL">
!     There was an error reading the namelist specified. Carefully examine all namelist and variables
!     for anything incorrect (e.g. malformed, hidden characters).
!   </ERROR>
!   <ERROR MSG="Unknown namelist, or mistyped namelist variable in namelist ...., (IOSTAT = ####)" STATUS="FATAL">
!     The name list given doesn't exist in the namelist file, or a variable in the namelist is mistyped or isn't a
!     namelist variable.
!   </ERROR>

! used to check the iostat argument that is
! returned after reading a namelist
! see the online documentation for how this routine might be used
  INTEGER FUNCTION check_nml_error(IOSTAT, NML_NAME)
    INTEGER, INTENT(in) :: IOSTAT
    CHARACTER(len=*), INTENT(in) :: NML_NAME

    CHARACTER(len=256) :: err_str

    IF ( .NOT.module_is_initialized) CALL fms_init()

    check_nml_error = IOSTAT

    ! Return on valid IOSTAT values
    IF ( IOSTAT <= 0 .OR. IOSTAT == nml_errors%multipleNMLSinFile ) RETURN

    ! Everything else is a FATAL
    IF ( mpp_pe() == mpp_root_pe() ) THEN
       IF ( (IOSTAT == nml_errors%badType1 .OR. IOSTAT == nml_errors%badType2) .OR. IOSTAT == nml_errors%missingVar ) THEN
          WRITE (err_str,*) 'Unknown namelist, or mistyped namelist variable in namelist ',TRIM(NML_NAME),', (IOSTAT = ',IOSTAT,')'
          CALL error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
          CALL mpp_sync()
       ELSE
          WRITE (err_str,*) 'Unknown error while reading namelist ',TRIM(NML_NAME),', (IOSTAT = ',IOSTAT,')'
          CALL error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
          CALL mpp_sync()
       END IF
    ELSE
       CALL mpp_sync()
    END IF
  END FUNCTION check_nml_error
! </FUNCTION>

!-----------------------------------------------------------------------
!   private routine for initializing allowable error codes

  SUBROUTINE nml_error_init
    ! Determines the IOSTAT error value for some common Namelist errors.
    ! Also checks if the compiler returns a non-zero status if there are
    ! multiple namelist records in a single file.
    INTEGER, PARAMETER :: unit_begin = 20, unit_end = 1024
    INTEGER :: fileunit, io_stat
    INTEGER, DIMENSION(4) :: nml_iostats
    LOGICAL :: opened

    ! Variables for sample namelists
    INTEGER :: i1, i2
    REAL :: r1, r2
    LOGICAL :: l1
    NAMELIST /a_nml/ i1, r1
    NAMELIST /b_nml/ i2, r2, l1
    NAMELIST /badType1_nml/ i1, r1
    NAMELIST /badType2_nml/ i1, r1
    NAMELIST /missingVar_nml/ i2, r2

    ! Initialize the sample namelist variables
    i1 = 1
    i2 = 2
    r1 = 1.0
    r2 = 2.0
    l1 = .FALSE.

    ! Create a dummy namelist file
    IF ( mpp_pe() == mpp_root_pe() ) THEN
       ! Find a free file unit for a scratch file
       file_opened: DO fileunit = unit_begin, unit_end
          INQUIRE(UNIT=fileunit, OPENED=opened)
          IF ( .NOT.opened ) EXIT file_opened
       END DO file_opened

#if defined __PGI
       OPEN (UNIT=fileunit, FILE='_read_error.nml', IOSTAT=io_stat)
#else
       OPEN (UNIT=fileunit, STATUS='SCRATCH', IOSTAT=io_stat)
#endif

       ! Write sample namelist to the SCRATCH file.
       WRITE (UNIT=fileunit, NML=a_nml, IOSTAT=io_stat)
       WRITE (UNIT=fileunit, NML=b_nml, IOSTAT=io_stat)
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&badType1_nml  i1=1, r1=''bad'' /",/)')
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&badType2_nml  i1=1, r1=.true. /",/)')
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&missingVar_nml  i2=1, r2=1.0e0, l1=.true. /",/)')

       ! Rewind for reading
       REWIND(UNIT=fileunit)

       ! Read the second namelist from the file -- check for namelist bug
       READ (UNIT=fileunit, NML=b_nml, IOSTAT=nml_iostats(1))
       REWIND(UNIT=fileunit)

       ! Read in bad type 1 --- Some compilers treat the string cast differently
       READ (UNIT=fileunit, NML=badType1_nml, IOSTAT=nml_iostats(2))
       REWIND(UNIT=fileunit)

       ! Read in bad type 2
       READ (UNIT=fileunit, NML=badType2_nml, IOSTAT=nml_iostats(3))
       REWIND(UNIT=fileunit)

       ! Read in missing variable/misstyped
       READ (UNIT=fileunit, NML=missingVar_nml, IOSTAT=nml_iostats(4))

       ! Done, close file
       CLOSE (UNIT=fileunit)

       ! Some compilers don't handle the type casting as well as we would like.
       IF ( nml_iostats(2) * nml_iostats(3) .EQ. 0 ) THEN
          IF ( nml_iostats(2) .NE. 0 .AND. nml_iostats(3) .EQ. 0 ) THEN
             nml_iostats(3) = nml_iostats(2)
          ELSE IF ( nml_iostats(2) .EQ. 0 .AND. nml_iostats(3) .NE.0 ) THEN
             nml_iostats(2) = nml_iostats(3)
          ELSE
             nml_iostats(2) = nml_iostats(4)
             nml_iostats(2) = nml_iostats(4)
          END IF
       END IF
    END IF

    ! Broadcast nml_errors
    CALL mpp_broadcast(nml_iostats,4,mpp_root_pe())
    nml_errors%multipleNMLSinFile = nml_iostats(1)
    nml_errors%badType1 = nml_iostats(2)
    nml_errors%badType2 = nml_iostats(3)
    nml_errors%missingVar = nml_iostats(4)

    do_nml_error_init = .FALSE.
  END SUBROUTINE nml_error_init

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
       ! only allow stdlog messages on root pe
         if ( mpp_pe() /= mpp_root_pe() ) return
     endif

     if (present(tag)) then
         write (logunit,'(/,80("="),/(a))') trim(version), trim(tag)
     else
         write (logunit,'(/,80("="),/(a))') trim(version)
     endif

 end subroutine write_version_number
! </SUBROUTINE>
! ========================================================================================================================



! ========================================================================================================================
! FROM FMS_IO_MOD
!use fms_io_mod,          only: save_restart, restart_file_type

! ========================================================================================================================




! ========================================================================================================================
! FROM DIAG_MANAGER_MOD
!use diag_manager_mod,    only: register_diag_field, send_data
! ========================================================================================================================





! ========================================================================================================================
! FROM MPP_DOMAINS_MOD
!use mpp_domains_mod,     only: mpp_define_domains, domain2d

! ========================================================================================================================




! ========================================================================================================================
! FROM MPP_MOD
!use mpp_mod,             only: input_nml_file, mpp_error, FATAL, WARNING, NOTE, stdout, stdlog
!use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE

! ========================================================================================================================




! ========================================================================================================================
! FROM OCEAN_DOMAINS_MOD
!use ocean_domains_mod,     only: get_local_indices

! ========================================================================================================================




! ========================================================================================================================
! FROM OCEAN_OBC_MOD
!use ocean_obc_mod,         only: store_ocean_obc_tracer_flux, ocean_obc_zero_boundary

! ========================================================================================================================




! ========================================================================================================================
! FROM OCEAN_PARAMETERS_MOD
!use ocean_parameters_mod,  only: ADVECT_UPWIND, ADVECT_QUICKER
!use ocean_parameters_mod,  only: TWO_LEVEL, missing_value

! ========================================================================================================================




! ========================================================================================================================
! FROM OCEAN_TYPES_MOD
!use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type, ocean_density_type, ocean_time_type
!use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_thickness_type, ocean_adv_vel_type

! ========================================================================================================================




! ========================================================================================================================
! FROM OCEAN_WORKSPACE_MOD
!use ocean_workspace_mod,   only: wrk1, wrk2, wrk3, wrk4

! ========================================================================================================================






! ========================================================================================================================
! FROM OCEAN_TRACER_UTIL_MOD
!use ocean_tracer_util_mod, only: tracer_psom_chksum

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


type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type)  , pointer :: Grd =>NULL()
type(ocean_domain_type), save    :: Dom_quicker
type(ocean_domain_type), save    :: Dom_flux

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

public  horz_advect_tracer
private horz_advect_tracer_upwind

public  vert_advect_tracer
private vert_advect_tracer_upwind

private advection_diag_init

private compute_adv_diss

public ocean_tracer_advect_init
public ocean_tracer_advect_end
public ocean_tracer_advect_restart

logical :: module_is_initialized   = .false.
logical :: have_obc                = .false.

logical :: limit_with_upwind       = .false.
logical :: debug_this_module       = .false.
logical :: advect_sweby_all        = .false.
logical :: zero_tracer_advect_horz = .false.
logical :: zero_tracer_advect_vert = .false.
logical :: write_a_restart         = .true.
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
subroutine ocean_tracer_advect_init (Grid, Domain, Time, Dens, T_prog, dtime, obc, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  real,                         intent(in)           :: dtime
  logical,                      intent(in)           :: obc
  logical,                      intent(in), optional :: debug

  integer :: n
  integer :: ioun, io_status, ierr

  integer :: stdoutunit,stdlogunit
  stdoutunit=stdout();stdlogunit=stdlog()

  write( stdlogunit,'(/a/)') trim(version)

  module_is_initialized = .true.

  have_obc = obc

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ocean_tracer_advect_nml, iostat=io_status)
  ierr = check_nml_error(io_status,'ocean_tracer_advect_nml')
#else
  ioun = open_namelist_file()
  read  (ioun, ocean_tracer_advect_nml,iostat=io_status)
  ierr = check_nml_error(io_status, 'ocean_tracer_advect_nml')
  call close_file(ioun)
#endif
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_tracer_advect_nml)
  write (stdlogunit, ocean_tracer_advect_nml)

#ifndef MOM_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk  = Grid%nk
  allocate(flux_x(isd:ied,jsd:jed,nk))
  allocate(flux_y(isd:ied,jsd:jed,nk))
  allocate(flux_z(isd:ied,jsd:jed,nk))
  allocate(advect_tendency(isd:ied,jsd:jed,nk))
#endif

  Dom => Domain
  Grd => Grid

  flux_x = 0.0
  flux_y = 0.0
  flux_z = 0.0

  advect_tendency(:,:,:)     = 0.0

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif
  if(debug_this_module) then
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_tracer_advect_mod with debug_this_module=.true.'
  endif

  num_prog_tracers = size(T_prog(:))
  do n=1,num_prog_tracers
     !if (trim(T_prog(n)%name) == 'temp') index_temp = n
     !if (trim(T_prog(n)%name) == 'passive_temp_sq') index_temp_sq = n
  enddo

  if(.not. write_a_restart) then
    write(stdoutunit,'(a)') '==>Note: running ocean_tracer_advect with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, so cannot restart if using PSOM advection.'
  endif

  if(zero_tracer_advect_horz) then
     call mpp_error(WARNING,'==>ocean_tracer_advect_mod: have turned OFF horizontal tracer advection')
     write(stdoutunit,*)'==>WARNING: have turned off horizontal tracer advection. Unrealistic simulation.'
  endif
  if(zero_tracer_advect_vert) then
     call mpp_error(WARNING,'==>ocean_tracer_advect_mod: have turned OFF vertical tracer advection')
     write(stdoutunit,*)'==>WARNING: have turned off vertical tracer advection. Unrealistic simulation.'
  endif

  if(advect_sweby_all) then
      if(async_domain_update) then
        write(stdoutunit,'(a)') &
        '==>Note: using asynchrnous domain update for MDFL SWEBY_all'
      endif

      write(stdoutunit,'(/a)')'==>ocean_tracer_advect_mod: advect_sweby_all=.true. so all tracers advected'
      write(stdoutunit,'(a)') '   with mdfl_sweby, regardless the settings in field_table.'
      write(stdoutunit,'(a/)')'   This method exploits mpp_update_domain capabilities and is faster in some cases.'

  else

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

  endif

  call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), Domain%layout, Dom_flux%domain2d, maskmap=Domain%maskmap&
             , xflags = Domain%xflags, yflags = Domain%yflags, xhalo=1, yhalo=1,name='flux'&
             , x_cyclic_offset = Domain%x_cyclic_offset, y_cyclic_offset = Domain%y_cyclic_offset)

   ! for binning some diagnostics to neutral density space
   neutralrho_nk = size(Dens%neutralrho_ref(:))
       allocate( nrho_work(isd:ied,jsd:jed,neutralrho_nk) )
       nrho_work(:,:,:) = 0.0

  call advection_diag_init(Time, Dens, T_prog(:))

  ! initialize clock ids
  id_clock_up_horz        = mpp_clock_id('(Ocean advect: horz up)          ',grain=CLOCK_ROUTINE)
  id_clock_up_vert        = mpp_clock_id('(Ocean advect: vert up)          ',grain=CLOCK_ROUTINE)


end subroutine ocean_tracer_advect_init
! </SUBROUTINE> NAME="ocean_tracer_advect_init"


!#######################################################################
! <SUBROUTINE NAME="advection_diag_init">
!
! <DESCRIPTION>
! Initialize the main tracer advection diagnostics.
! </DESCRIPTION>
!
subroutine advection_diag_init (Time, Dens, T_prog)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_density_type),     intent(in)  :: Dens
  type(ocean_prog_tracer_type), intent(in)  :: T_prog(:)

  integer :: n
  integer :: stdoutunit
  stdoutunit=stdout()

  allocate (id_tracer_advection(num_prog_tracers))
  allocate (id_tracer2_advection(num_prog_tracers))
  allocate (id_tracer_adv_diss(num_prog_tracers))
  allocate (id_horz_advect(num_prog_tracers))
  allocate (id_vert_advect(num_prog_tracers))
  allocate (id_advection_x(num_prog_tracers))
  allocate (id_advection_y(num_prog_tracers))
  allocate (id_advection_z(num_prog_tracers))
  allocate (id_xflux_adv(num_prog_tracers))
  allocate (id_yflux_adv(num_prog_tracers))
  allocate (id_zflux_adv(num_prog_tracers))
  allocate (id_xflux_adv_int_z(num_prog_tracers))
  allocate (id_yflux_adv_int_z(num_prog_tracers))

  id_tracer_advection =-1
  id_tracer2_advection=-1
  id_tracer_adv_diss  =-1
  id_horz_advect      =-1
  id_vert_advect      =-1
  id_advection_x      =-1
  id_advection_y      =-1
  id_advection_z      =-1
  id_xflux_adv        =-1
  id_yflux_adv        =-1
  id_zflux_adv        =-1
  id_xflux_adv_int_z  =-1
  id_yflux_adv_int_z  =-1

  do n=1,num_prog_tracers

    if(n==index_temp) then
      id_tracer_advection(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_advection', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*advection tendency',             &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_tracer2_advection(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sq_advection', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*advection tendency for (cp*temp)^2',    &
                   'kg/(m^2*sec)*(J/kg)^2', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_horz_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_horz_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*horz advect tendency',        &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_vert_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_vert_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*vert advect tendency',        &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_advection_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_x_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*x-advective heating',   &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_advection_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_y_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*y-advective heating',   &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_advection_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_z_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'cp*rho*dzt*z-advective heating',   &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_xflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_adv', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'cp*rho*dzt*dyt*u*temp',        &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_yflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_adv', &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time, 'cp*rho*dzt*dxt*v*temp',        &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_zflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_zflux_adv', &
                   Grd%tracer_axes_wt(1:3), Time%model_time, 'cp*rho*dxt*dyt*wt*temp',           &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_xflux_adv_int_z(n) = register_diag_field ('ocean_model', &
                   trim(T_prog(n)%name)//'_xflux_adv_int_z',      &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time,   &
                   'z-integral of cp*rho*dyt*u*temp',              &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_yflux_adv_int_z(n) = register_diag_field ('ocean_model', &
                   trim(T_prog(n)%name)//'_yflux_adv_int_z',      &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time,   &
                   'z-integral of cp*rho*dxt*v*temp',              &
                   'Watts', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_tracer_adv_diss(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_adv_diss', &
      Grd%tracer_axes(1:3), Time%model_time,                                                         &
      'dissipation of squared '//trim(T_prog(n)%name)//' from advection errors',                     &
      '(Watt/m^2)^2', missing_value=missing_value,                                                   &
       range=(/-1.e18,1.e18/))

    else

      id_tracer_advection(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_advection', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*advection tendency',         &
                   'kg/(sec*m^2)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_tracer2_advection(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sq_advection',    &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*advection tendency for tracer sqrd',&
                   'kg/(m^2*s)*(tracer squared)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_horz_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_horz_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*horz advect tendency',    &
                   'kg/(sec*m^2)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_vert_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_vert_advect', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*vert advect tendency',    &
                   'kg/(sec*m^2)', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_advection_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_x_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*x-advective tendency',     &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_advection_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_y_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*y-advective tendency',     &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_advection_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_z_adv', &
                   Grd%tracer_axes(1:3), Time%model_time, 'rho*dzt*z-advective tendency',     &
                   trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_xflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_adv', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho*dzt*dyt*u*tracer',         &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_yflux_adv(n)   = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_adv', &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time, 'rho*dzt*dxt*v*tracer',           &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_zflux_adv(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_zflux_adv', &
                   Grd%tracer_axes_wt(1:3), Time%model_time, 'rho*dxt*dyt*wt*tracer',            &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_xflux_adv_int_z(n) = register_diag_field ('ocean_model', &
                   trim(T_prog(n)%name)//'_xflux_adv_int_z',      &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time,   &
                   'z-integral of rho*dzt*dyt*u*tracer',           &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))
      id_yflux_adv_int_z(n) = register_diag_field ('ocean_model', &
                   trim(T_prog(n)%name)//'_yflux_adv_int_z',      &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time,   &
                   'z-integral of rho*dzt*dxt*v*tracer',           &
                   'kg/sec', missing_value=missing_value, range=(/-1.e18,1.e18/))

      id_tracer_adv_diss(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_adv_diss', &
      Grd%tracer_axes(1:3), Time%model_time,                                                         &
      'dissipation of squared '//trim(T_prog(n)%name)//' from advection errors',                     &
      '[kg/(m^2*sec)]^2', missing_value=missing_value,                                               &
       range=(/-1.e18,1.e18/))


    endif

  enddo


end subroutine advection_diag_init
! </SUBROUTINE> NAME="advection_diag_init"



!#######################################################################
! <SUBROUTINE NAME="horz_advect_tracer">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers
! </DESCRIPTION>
!
subroutine horz_advect_tracer(Time, Adv_vel, Thickness, Dens, T_prog, Tracer, ntracer, dtime, store_flux)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
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
          call mpp_error(FATAL,&
          '==>Error from ocean_tracer_advect_mod (horz_advect_tracer): chose invalid horz advection scheme')
      end select

      if (have_obc) call ocean_obc_zero_boundary(Tracer%wrk1, "T")

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


      ! send some diagnostics

      if (id_horz_advect(ntracer) > 0 .and. Tracer%horz_advect_scheme /= ADVECT_PSOM) then
           used = send_data(id_horz_advect(ntracer), Tracer%conversion*Tracer%wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                                        &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_psom_advect(ntracer) > 0 .and. Tracer%horz_advect_scheme == ADVECT_PSOM) then
           used = send_data(id_psom_advect(ntracer), Tracer%conversion*Tracer%wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                                        &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_xflux_adv(ntracer) > 0) then
           used = send_data(id_xflux_adv(ntracer), Tracer%conversion*flux_x(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_yflux_adv(ntracer) > 0) then
           used = send_data(id_yflux_adv(ntracer), Tracer%conversion*flux_y(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_xflux_adv_int_z(ntracer) > 0) then
          tmp_flux(:,:) = 0.0
          do k=1,nk
             tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) +  flux_x(isc:iec,jsc:jec,k)
          enddo
          used = send_data(id_xflux_adv_int_z(ntracer), Tracer%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif
      if (id_yflux_adv_int_z(ntracer) > 0) then
          tmp_flux(:,:) = 0.0
          do k=1,nk
             tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) +  flux_y(isc:iec,jsc:jec,k)
          enddo
          used = send_data(id_yflux_adv_int_z(ntracer), Tracer%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif

      if (have_obc) then
        call store_ocean_obc_tracer_flux(Time,Tracer,flux_x,ntracer,'z','adv')
        call store_ocean_obc_tracer_flux(Time,Tracer,flux_y,ntracer,'m','adv')
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
subroutine vert_advect_tracer(Time, Adv_vel, Dens, Thickness, T_prog, Tracer, ntracer, dtime)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_density_type),     intent(in)    :: Dens
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
        call mpp_error(FATAL,&
        '==>Error from ocean_tracer_advect_mod (vert_advect_tracer): invalid advection scheme chosen')
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

      if(id_tracer_advection(ntracer) > 0) then
          used = send_data(id_tracer_advection(ntracer), Tracer%conversion*advect_tendency(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:),                                            &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif

      ! send some more diagnostics
      if (id_vert_advect(ntracer) > 0) then
           used = send_data(id_vert_advect(ntracer), Tracer%conversion*Tracer%wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                                        &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_zflux_adv(ntracer) > 0)  then
          used = send_data(id_zflux_adv(ntracer), Tracer%conversion*flux_z(:,:,:),&
          Time%model_time, rmask=Grd%tmask(:,:,:),                                &
          is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif


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

  call mpp_clock_begin(id_clock_up_horz)

  if ( .not. module_is_initialized ) then
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_upwind): ocean_tracer_advect_mod not yet initialized')
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

  call mpp_clock_end(id_clock_up_horz)


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

  call mpp_clock_begin(id_clock_up_vert)

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

  call mpp_clock_end(id_clock_up_vert)

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

  call mpp_clock_begin(id_clock_adv_diss)

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

  call mpp_clock_end(id_clock_adv_diss)


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
! <SUBROUTINE NAME="ocean_tracer_advect_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_tracer_advect_restart(T_prog, time_stamp)
  type(ocean_prog_tracer_type), intent(in)           :: T_prog(:)
  character(len=*),             intent(in), optional :: time_stamp
   integer :: tau, taup1

  if(ANY(T_prog(1:num_prog_tracers)%horz_advect_scheme == ADVECT_PSOM) )then
     call save_restart(Adv_restart)
  end if

end subroutine ocean_tracer_advect_restart
! </SUBROUTINE> NAME="ocean_tracer_advect_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_advect_end">
!
! <DESCRIPTION>
! Write the PSOM moments for restarts.
! </DESCRIPTION>
subroutine ocean_tracer_advect_end(Time, T_prog)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  integer :: n
  character(len=128) :: filename

  integer :: stdoutunit
  stdoutunit=stdout()

  filename = 'RESTART/ocean_psom_moments.res'

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') '==>Warning from ocean_tracer_advect_end: NO restart written for PSOM moments.'
    call mpp_error(WARNING,'==>Warning from ocean_tracer_advect_end: NO restart written for PSOM moments.')
    return
  endif

  call ocean_tracer_advect_restart(T_prog)

  do n=1,num_prog_tracers
    if (T_prog(n)%horz_advect_scheme == ADVECT_PSOM) then
      call tracer_psom_chksum(Time, T_prog(n))
      write (stdoutunit,*) 'Completed write of psom restart for tracer ', trim(T_prog(n)%name)
    endif

  enddo

  return

end subroutine ocean_tracer_advect_end
! </SUBROUTINE> NAME="ocean_tracer_advect_end"



end module ocean_tracer_advect_mod
