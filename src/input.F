
!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************

! NAME:
!       INPUT

! AUTHOR:
!       Jose Lopez Miralles
!       Astronomy and Astrophysics Department, University of Valencia
!       46100, Burjassot (Valencia)

! SUPERVISOR:
!       Dr. Manel Perucho Pla
!       Astronomy and Astrophysics Department, University of Valencia
!       46100, Burjassot (Valencia)

!       Dr. Jose Maria Marti Puig
!       Astronomy and Astrophysics Department, University of Valencia
!       46100, Burjassot (Valencia)

!       Dr. Simone Migliari
!       European Space Agency / European Space Astronomy Center
!       28692, Villafranca del Castillo (Madrid)

! DESCRIPTION:
!       Load input parameters.

! ACADEMIC REFERENCE:


! ===============================================================================      
SUBROUTINE INPUT (path, suffix, suffixd, basenm, limiter,            &
                  nend, irstrt, nrstrt, nout, itstp, iars,           &
                  irk, iccv, niter, nsdim, nx, ny, nz,               &
                  igeomx, igeomy, igeomz, eqs,                       &
                  tmax, trstrt, tout, cfl, dtini,                    &
                  small, smlrho, smallp, smallu,                     &
                  gridlx, gridly, gridlz, stelar, cool)
! ===============================================================================
                  
        USE parallel
        USE math_constants
        USE physic_constants
        
        implicit none
# include "types.h"

        !------------------------------------------------------ Input  Variables

        character(256):: path

        !------------------------------------------------------ Output Variables

        character(2):: suffix
        character(3):: suffixd
        character(4):: basenm, limiter

        integer(is) :: nend, irstrt, nrstrt, nout, itstp, iars,        &
                       irk,  iccv, niter, nsdim, nx, ny, nz, trstrt,   &                           
                       igeomx, igeomy, igeomz, eqs, stelar, cool

        real(rs)    :: tmax, tout, cfl, dtini,                 &
                       small, smlrho, smallp, smallu,          &
                       gridlx, gridly, gridlz                          

        !------------------------------------------------------- Local Variables

        character(72) :: text

        character(8)  :: label, label1
        character(256):: fich

        integer(is)   :: ierr

	character(37) :: message1
        
  
        !-------------------------------------------------------------------------- BEGIN

        ierr  = 0
        
        ! Search input file in path.
        fich= TRIM(path)//'/'//'inpt.dat'
  
        open (1, file=fich, form='formatted', status='old', iostat=ierr)
        IF (ierr .NE. 0) THEN
           call ERROR ('INPUT: Error openning FILE', fich, __FILE__, __LINE__)
        ENDIF

        read (1,'(a72)') text
        read (1,'(a72)') text
        read (1,'(a72)') text

        ! ===========================
        ! INPUT STARTS
        ! ===========================

        ! BASENM is the root for the output, plot and restart file names. 
        read (1,'(a8,1x,a4)') label, basenm

        label1 = 'basenm'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------
        
        !NEND is the total number of timestps.
        read (1,'(a8,1x,i10)') label, nend

        label1 = 'nend'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF
        
        !------------------
        
        ! The program stops when time is .GE. than TMAX. 
        read (1,'(a8,1x,d23.2)') label, tmax

        label1 = 'tmax'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF


        !------------------

        !SUFFIX is the suffix for the output and restart file names.
        read (1,'(a8,1x,a2)') label, suffix

        label1 = 'suffix'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------        

        !SUFFIXD is the suffix for the 3D output file names.
        read (1,'(a8,1x,a3)') label,suffixd

        label1 = 'suffixd'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        ! If IRSTRT == 1, the program restarts from the file BASENM//'R'//SUFFIX 
        read (1,'(a8,1x,i10)') label, irstrt

        label1 = 'irstrt'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !A restart file is written every NRSTRT timesteps.
        read (1,'(a8,1x,i10)') label, nrstrt

        label1 = 'nrstrt'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !A restart file is written every TRSTRT time units.
        read (1,'(a8,1x,i10)') label, trstrt

        label1 = 'trstrt'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !An output file is written every NOUT timesteps. (Only for 1D problems).
        read (1,'(a8,1x,i10)') label, nout

        label1 = 'nout'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !An output file is written every TOUT time units. (Only for 1D problems).
        read (1,'(a8,1x,d23.2)') label, tout

        label1 = 'tout'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !DT is written on the screen every ITSTP timesteps.
        read (1,'(a8,1x,i10)') label, itstp

        label1 = 'itstp'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !CFL (Courant–Friedrichs–Lewy condition) is the 
        !time-limiting factor (.LT. 1).
        read (1,'(a8,1x,d23.2)') label, cfl

        label1 = 'cfl'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !IARS is the Riemann Solver (HLL=2,HLLC=3).
        read (1,'(a8,1x,i10)') label, iars
        
        label1 = 'iars'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF


        !------------------

        !IRK is the order of the Runge-Kutta used for time integration.
        read (1,'(a8,1x,i10)') label, irk

        label1 = 'irk'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !ICCV is the correction of conserved variables for
        !multidimensional simulations.
        read (1,'(a8,1x,i10)') label, iccv

        label1 = 'iccv'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !LIMITER for reconstruction ('VANL','VLMM','MC','MCMM','MINMOD').
        read (1,'(a8,1x,a4)') label, limiter

        label1 = 'limiter'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !DTINI is the initial DT.
        read (1,'(a8,1x,d23.2)') label, dtini

        label1 = 'dtini'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        ! NITER is the number of bisections to get the primitive variables and
        ! the number of Newton-Rhapson iterations to get the extreme wavespeeds.
        read (1,'(a8,1x,i10)') label, niter

        label1 = 'niter'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !SMALL is the threshold for nondimensional numbers (i.e. velocity).
        read (1,'(a8,1x,d23.2)') label, small

        label1 = 'small'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !SMLRHO is the threshold for densities (rho, r).
        read (1,'(a8,1x,d23.2)') label, smlrho

        label1 = 'smlrho'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !SMALLP is the threshold for pressure.
        read (1,'(a8,1x,d23.2)') label, smallp

        label1 = 'smallp'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !SMALLU is the threshold for internal energy.
        read (1,'(a8,1x,d23.2)') label, smallu

        label1 = 'smallu'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        ! STELAR is the flag for stelar mass load.
        read (1,'(a8,1x,i10)') label, stelar

        label1 = 'stelar'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        
        !------------------

        ! COOL is the flag for cooling.
        read (1,'(a8,1x,i10)') label, cool

        label1 = 'cool'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !NSDIM is the number of spatial dimensions (1, 2 or 3).
        read (1,'(a8,1x,i10)') label, nsdim

        label1 = 'nsdim'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF


        !------------------

        ! EQS determines the Equation of State (EoS)
        read (1,'(a8,1x,i10)') label, eqs

        label1 = 'eqs'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------

        !GRIDLX is the length of the grid in the x-direction.
        read (1,'(a8,1x,d23.2)') label, gridlx
        
        label1 = 'gridlx'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

     
        !GRIDLY is the length of the grid in the y-direction.
        read (1,'(a8,1x,d23.2)') label, gridly

        label1 = 'gridly'
        IF (label .NE. label1) then
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        
        !GRIDLZ is the length of the grid in the z-direction.
        read (1,'(a8,1x,d23.2)') label, gridlz

        label1 = 'gridlz'
        IF (label .NE. label1) then
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        !------------------
        
        !NX is the total number of grid points in the x-direction.
        read (1,'(a8,1x,i10)') label, nx

        label1 = 'nx'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        IF (nx .LT. 4) THEN
           ierr = ierr + 1
        ENDIF
  
        !NY is the total number of grid points in the y-direction.
        read (1,'(a8,1x,i10)') label, ny

        label1 = 'ny'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        IF (nsdim .GT. 1 .AND. ny .LT. 4) THEN
           ierr = ierr + 1
        ENDIF


        !NZ is the total number of grid points in the z-direction.
        read (1,'(a8,1x,i10)') label, nz

        label1 = 'nz'
        IF (label .NE. label1) THEN
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        IF (nsdim .EQ. 3 .AND. nz .LT. 4) THEN
           ierr = ierr + 1
        ENDIF

        !------------------

        !IGEOMX defines the type of the x coordinate. 
        read (1,'(a8,1x,i10)') label, igeomx

        label1 = 'igeomx'
        IF (label .NE. label1) then
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        ! LOSTREGO only works with planar coordinates.
        IF (igeomx .NE. 0) THEN
           ierr = ierr + 1
        ENDIF


        !IGEOMY defines the type of the y coordinate. 
        read (1,'(a8,1x,i10)') label, igeomy

        label1 = 'igeomy'
        IF (label .NE. label1) then
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        ! LOSTREGO only works with planar coordinates.        
        IF (igeomy .NE. 0) THEN
           ierr = ierr + 1
        ENDIF


        !IGEOMZ defines the type of the z coordinate. 
        read (1,'(a8,1x,i10)') label, igeomz

        label1 = 'igeomz'
        IF (label .NE. label1) then
           call INPUTERROR ('INPUT', label, label1, ierr)
        ENDIF

        ! LOSTREGO only works with planar coordinates.        
        IF (igeomz .NE. 0) THEN
           ierr = ierr + 1
        ENDIF

        !------------------

        close (1)

        ! Debug
        IF (ierr .NE. 0) THEN
           call ERROR ('INPUT', 'At least one error in input file',     &
                        __FILE__, __LINE__ )
        ENDIF

        RETURN
        
        !-------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE INPUT
! ===============================================================================

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


! NAME:
!       INPUTERROR

! AUTHOR:
!       Jose Lopez Miralles
!       Astronomy and Astrophysics Department, University of Valencia
!       46100, Burjassot (Valencia)

! SUPERVISOR:
!       Dr. Manel Perucho Pla
!       Astronomy and Astrophysics Department, University of Valencia
!       46100, Burjassot (Valencia)

!       Dr. Jose Maria Marti Puig
!       Astronomy and Astrophysics Department, University of Valencia
!       46100, Burjassot (Valencia)

!       Dr. Simone Migliari
!       European Space Agency / European Space Astronomy Center
!       28692, Villafranca del Castillo (Madrid)

! DESCRIPTION:
!       Debug: Write an error message for the input file.

! ACADEMIC REFERENCE:


! ===============================================================================       
SUBROUTINE INPUTERROR (rutina, label, label1, ierr)
! =============================================================================== 

        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(10)::  rutina
        character(8) ::  label, label1

        !------------------------------------------------------ Output Variables

        integer(is)  ::  ierr


        !------------------------------------------------------------------------- BEGIN

        write(6,*) ' '
        write(6,*) rutina, ': Incorrect input deck'
        write(6,*) '       label = ',label, '  expected label = ', label1
        ierr = ierr+1

        !------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE INPUTERROR   
! =============================================================================== 
