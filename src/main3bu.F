
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+                                                                                +
!+ NAME:                                                                          +
!+ -------- MAIN PROGRAM: LOSTREGO v1.2--------                                   +
!+                                                                                +
!+ AUTHOR:                                                                        +
!+       Jose Lopez Miralles                                                      +
!+       Astronomy and Astrophysics Department, University of Valencia            +
!+       46100, Burjassot (Valencia)                                              +
!+                                                                                +
!+ SUPERVISOR:                                                                    +
!+       Dr. Manel Perucho Pla                                                    +
!+       Astronomy and Astrophysics Department, University of Valencia            +
!+       46100, Burjassot (Valencia)                                              +
!+                                                                                +
!+       Dr. Jose Maria Marti Puig                                                +
!+       Astronomy and Astrophysics Department, University of Valencia            +
!+       46100, Burjassot (Valencia)                                              +
!+                                                                                +
!+ DESCRIPTION:                                                                   +
!+       MAIN PROGRAM: 1D, 2D, 3D SRMHD SIMULATIONS.                              +
!+                     HYBRID OMP+MPI PARALLELIZATION.                            +
!+                                                                                +
!+ ACADEMIC REFERENCE:                                                            +
!+       See academic references in subroutines.                                  +
!+                                                                                +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     


      
PROGRAM RATPENAT

        USE math_constants
        USE physic_constants
        USE tipos
        USE memoria
        USE parallel
        USE workarrays
        
        implicit none


# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

 
# include "types.h"

        !------------------------------------------------------- Local Variables
      
        logical(ls):: NoSalir

        character(8)  :: datfil 
        character(7)  :: rstfil, outfil, tdlfil
        character(4)  :: basenm, limiter
        character(3)  :: suffixd
        character(2)  :: suffix

        integer(is) :: longsize(3), mpi_dims1, mpi_dims2, mpi_dims3, max_dims
	integer (is) :: mpi_coordx, mpi_coordy, mpi_coordz              
	
	integer (is) :: H_openr, H_opend, rst_file,rst_rec, rst_filed, curR_R, r_filer, r_filed

        ! Error message

        character(40) :: message1, message3
        character(37) :: message2
        character(42) :: message4
        character(36) :: message5
        character(34) :: message6      
	character(32) :: message7

        integer(is) :: nx, ny, nz, nxtot, nytot, nztot
      
        ! Local & maximum sizes for workarrays (+1, +5, and +6)

        integer(is):: mnx1, mnx5, mnx6, &
                      mny1, mny5, mny6, &
                      mnz1, mnz5, mnz6, &
                      mn1,  mn5,  mn6

        ! Input
                      
        integer(is):: ierr, i, j, k,                                  &
                      nbegin, rkin, eqs,                              &
                      igeomx, igeomy, igeomz,                         &
                      bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz, &
                      nend, irstrt, nrstrt, nout, iars, irk, iccv,    &
                      itstp, niter, nsdim, stelar, cool,              &
                      nstep, nout1, nrst, trstrt, final
                   
        real(rs)   :: alpha_R, beta_R,                          &
                      gamma_ad, tmax, tout, cfl, dtini, &
                      small, smlrho, smallp, smallu,            &
                      gridlx, gridly, gridlz,                   &
                      delx, dely, delz,                         &
                      timej, dt, tout1, trst,                   &
                      me, mp, radius, R_b, rho_0a, m_0

        real(rs) :: pc, me_cgs, c, kb, mp_cgs, KB1, KB2, meff  
      

        ! Allocatable
                      
        real(rs), allocatable:: xznl(:), xzn(:), xznr(:), &
                                yznl(:), yzn(:), yznr(:), &
                                zznl(:), zzn(:), zznr(:), &
                                grav(:)

        type(Tphysic), allocatable:: physic (:,:,:)

	type(ConsTransp), allocatable:: CT (:,:,:)
      
        type(TboundPosX), allocatable:: boundPosX(:)
        type(TboundPosY), allocatable:: boundPosY(:)
        type(TboundPosZ), allocatable:: boundPosZ(:)

        integer :: numArgs, xCPU, yCPU, zCPU
        character(len=256) :: program_name, tmpStr

        real(rs) :: b2, bct2, vel2

        INTEGER:: iTime1, iTime2, rate

	real(rs) :: uwi, vwi, shift, Rorb, sum

        real(rs) :: RHOB, VELB, PB, TT1, GAMM1, RBM, BMJ, WB, BPHI1        
	real(rs) :: SCRCH1, SCRCH2, B40, B41, B42, B43, B44, DIVE, dist, ran1


        integer(is) :: idum, ncl, ic(2000), jc(2000), kc(2000), ii, jj, kk


       	integer(is):: flag_zero (3, 20)

	real(rs) :: vb

        !-------------------------------------------------------------------------- BEGIN


        ! CPU TIME
        ! ===================
        call SYSTEM_CLOCK(COUNT_RATE=RATE)
        call SYSTEM_CLOCK(iTime1)


        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !                                          INPUT
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        ! Retrieve n-th argument that is passed on the command line when the program is invoked.
        outRoot =''
        path_wkd =''
        numArgs = COMMAND_ARGUMENT_COUNT ()
        
        IF (numArgs .NE. 5) THEN        
           call ERROR ('MAIN3BU', 'Incorrect number of arguments...', __FILE__, __LINE__)
        ENDIF

        call GET_COMMAND_ARGUMENT (0, program_name)
        call GET_COMMAND_ARGUMENT (1, outRoot)
        call GET_COMMAND_ARGUMENT (2, path_wkd)
        
        call GET_COMMAND_ARGUMENT (3, tmpstr)
        read (tmpstr, *) xCPU
        
        call GET_COMMAND_ARGUMENT (4, tmpstr)
        read (tmpstr, *) yCPU
        
        call GET_COMMAND_ARGUMENT (5, tmpstr)
        read (tmpstr, *) zCPU

             

        ! Read initial parameters
        call INPUT (path_wkd, suffix, suffixd, basenm, limiter,        &
                    nend, irstrt, nrstrt, nout, itstp, iars,           &
                    irk, iccv, niter, nsdim, nx, ny, nz,               &
                    igeomx, igeomy, igeomz, eqs,                       &
                    tmax, trstrt, tout, cfl, dtini,                    &
                    small, smlrho, smallp, smallu,                     &
                    gridlx, gridly, gridlz, stelar, cool)


        
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !                                          PARALLELIZATION
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      
        !    ====================
        !    OMP+MPI INITIALIZATION
        !    ====================

        call OMP_INIT ()

        ! Save total nx,ny,nz before localization.
        nxtot = nx
        nytot = ny 
        nztot = nz

        call DECOMPOSE_MPI_DOMAIN (xCPU, yCPU, zCPU, nx, ny, nz)
	
        call MPI_INIT3D (xCPU, yCPU, zCPU, nsdim, max_dims)
       
        call DISTRIBUTE_WORK (max_dims, longsize, nxtot, nytot, nztot,    &       
                              mn1, mn5, mn6, mnx1, mnx5, mnx6,            &
                              mny1, mny5, mny6, mnz1, mnz5, mnz6)
       



        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !                                WRITING WARNINGS IN OUTPUT FILE 
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



        IF (cartrank .EQ. 0) THEN
           print*, 'MAIN:  xCPUS yCPUS zCPUS.....',  xCPU, yCPU, zCPU, outRoot
        ENDIF


        ! For tests, tmax is defined by the initial conditions.

        IF ((BASENM .EQ. 'BA1_' .OR. BASENM .EQ. 'BA2_' .OR. BASENM .EQ. 'BA3_' .OR. BASENM .EQ. 'BA3S' .OR.  &
             BASENM .EQ. 'BA4_') .AND. TMAX .NE. 0.402D0) THEN

           tmax = 0.402D0

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', 'tmax reset to 0.402d0', __FILE__, __LINE__)
           ENDIF

        ELSE IF (BASENM .EQ. 'BA5_' .AND. TMAX .NE. 0.552D0) THEN

           tmax = 0.552D0

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', 'tmax reset to 0.552d0', __FILE__, __LINE__)
           ENDIF

        ELSE IF (BASENM .EQ. 'GR2_' .AND. TMAX .NE. 1.502D0) THEN

           tmax = 1.502D0

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', 'tmax reset to 1.502d0', __FILE__, __LINE__)
           ENDIF

        ELSE IF (BASENM(1:2) .EQ. 'KO' .AND. TMAX .NE. 1.02D0) THEN

           tmax = 1.02D0

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', 'tmax reset to 1.02d0', __FILE__, __LINE__)
           ENDIF

        ELSE IF (BASENM .EQ. 'CMBW' .OR. BASENM .EQ. 'SMBW') THEN

           tmax = 4.1D0

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', 'tmax reset to 4.1D0', __FILE__, __LINE__)
           ENDIF

        ENDIF


	! Using Riemann Solver....

        IF (iars .EQ. 2 .AND. cartrank .EQ. 0) WRITE(6,*) 'USING RIEMANN SOLVER: HLLE'
        IF (iars .EQ. 3 .AND. cartrank .EQ. 0) WRITE(6,*) 'USING RIEMANN SOLVER: HLLC'
        IF (iars .EQ. 4 .AND. cartrank .EQ. 0) WRITE(6,*) 'USING RIEMANN SOLVER: HLLD'


        ! For tests, nsdim is defined by the initial conditions.

        IF ((BASENM .EQ. 'BA1_' .OR. BASENM .EQ. 'BA2_' .OR. BASENM .EQ. 'BA3_' .OR. BASENM .EQ. 'BA3S' .OR.     &
             BASENM .EQ. 'BA4_' .OR. BASENM .EQ. 'BA5_' .OR. BASENM .EQ. 'GR2_' .OR.      &
             BASENM .EQ. 'KO1_' .OR. BASENM .EQ. 'KO2_' .OR. BASENM .EQ. 'DZAW')          &
             .AND. NSDIM .NE. 1) THEN

           nsdim = 1

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', 'nsdim reset to 1', __FILE__, __LINE__)
           ENDIF


        ELSE IF ((BASENM .EQ. 'CMBW' .OR. BASENM .EQ. 'RROT' .OR. BASENM .EQ. 'SYBW')     & 
                  .AND. NSDIM .NE. 2) THEN

           nsdim = 2

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', 'nsdim reset to 2', __FILE__, __LINE__)
           ENDIF

        ELSE IF (BASENM .EQ. 'SMBW' .AND. NSDIM .NE. 3) THEN

           nsdim = 3

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', 'nsdim reset to 3', __FILE__, __LINE__)
           ENDIF

        ENDIF


        IF (iars .EQ. 3 .AND. nsdim .EQ. 2) THEN

           message1='HLLC in 2D is only valid for velz=0, Bz=0'

           IF (cartrank .EQ. 0) THEN
              call WARNING ('INPUT', message1, __FILE__, __LINE__)
           ENDIF

        ENDIF

        IF (nsdim .GT. 3) THEN
           call ERROR ('INPUT','Unsuitable number of spatial dimensions' , __FILE__, __LINE__)
        ENDIF


        ! Load and write constants (defined in Modulos.F).
        call SHOW_PHYS_CONST (pc, me_cgs, c, kb, mp_cgs, KB1, KB2, meff)

	IF (cartrank .EQ. 0) THEN

                   write(*,*) '------------ Physical constants ----------'
                   write(*,*) 'Physical variables in use....:'
                   write(*,*) 'pc=', pc
                   write(*,*) 'me_cgs=', me_cgs
                   write(*,*) 'c=', c
                   write(*,*) 'kb=', kb
                   write(*,*) 'mp_cgs=', mp_cgs
                   write(*,*) 'KB1=', KB1
                   write(*,*) 'KB2=', KB2
                   write(*,*) 'kb=', kb
                   write(*,*) 'meff=', meff
                   write(*,*) '-------------------------------------------'

	ENDIF


	! ----------------------------------------------------------------
        radius = 0.008D0*DTAN(5.D0*PI/180)
        R_b    = radius * pc
        rho_0a = 0.5D0*mp_cgs
        m_0    = rho_0a * R_b**3
        me     = me_cgs / m_0
        mp     = mp_cgs / m_0
	! ----------------------------------------------------------------

        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !                                          ALLOCATION
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                         
        allocate (xznl(mnx5), xzn(mnx5), xznr(mnx5),   &
                  yznl(mny5), yzn(mny5), yznr(mny5),   &
                  zznl(mnz5), zzn(mnz5), zznr(mnz5),   &
                  stat=ierr)

        IF (ierr .NE. 0) THEN

           message1='Allocation of main variables: Positions.'
           call ERROR ('MAIN3BU', message1, __FILE__, __LINE__)

        ENDIF

        allocate (physic(-4:mnx5,-4:mny5,-4:mnz5), stat=ierr)

        IF (ierr .NE. 0) THEN

           message2='Allocation of main variables: physic.'
           call ERROR ('MAIN3BU', message2, __FILE__, __LINE__)

        ENDIF

        allocate (boundPosX(-4:5), boundPosY(-4:5), boundPosZ(-4:5), stat=ierr)
      
        IF (ierr .NE. 0) THEN

           message3='Allocation of main variables: Positions.'
           call ERROR ('MAIN3BU', message3, __FILE__, __LINE__)

        ENDIF

	allocate (grav(-4:mny5), stat=ierr)

	IF (ierr .NE. 0) THEN

           message6='Allocation of main variable: grav.'
	   call ERROR ('MAIN3BU', message6, __FILE__, __LINE__)

	ENDIF

	allocate (CT(-4:mnx6,-4:mny6,-4:mnz6), stat=ierr)

        IF (ierr .NE. 0) THEN

           message7='Allocation of main variable: CT.'
           call ERROR ('MAIN3BU', message7, __FILE__, __LINE__)

        ENDIF

       
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !                                          START / RESTART
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! Initialize HDF5.

#ifdef HDF
      call rrhdf5_init(H_openr, H_opend)
#endif


        ! No Restart file :: construct a new model.

        IF (irstrt .EQ. 0) THEN

           timej  = 0.D0
           nstep  = 0
           dt     = dtini
           nbegin = 1

#ifdef HDF
           rst_file  = 0
           rst_rec   = 0
           rst_filed = 0
#endif
        
           IF (suffix(2:2) .NE. 'A') THEN
              call ERROR ('MAIN', 'Check input file suffix', __FILE__, __LINE__)
           ENDIF

           IF (cartrank .EQ. 0) THEN
              
              write(6,*) ' '
              write(6,*) ' MAIN: Constructing new initial model...'
              write(6,*) ' '
              
           ENDIF

           ! MPI Barrier.
           call PAR_BARRIER ()

           ! Allocate working arrays.
           call GETWAMEMORYSPACE (mn6, mn5, mn1, mnx1, mny1, mnz1, nx, ny, nz, &
                                  mpi_dims(1), mpi_dims(2), mpi_dims(3))
                      
          
           call GRID (basenm, nx, ny, nz, nsdim, gridlx, gridly, gridlz,                   &
                      xCPU, yCPU, zCPU, bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                      delx, dely, delz, xznl, xzn, xznr, yznl, yzn, yznr,                  &
                      zznl, zzn, zznr)


           call INIT (basenm, nx, ny, nz, mnx5, mny5, mnz5, mnx6, mny6, mnz6, eqs, iars,     &
                      xzn, xznl, xznr, yznr, delx, yzn, yznl, zzn, zznl, physic, gamma_ad, gridlx, gridly, gridlz, CT,   &
                      me, mp, m_0, kb, c, uwi, vwi, shift, Rorb, xCPU, yCPU, zCPU,ic,jc,kc,ncl,timej,pc)




           IF (cartrank .EQ. 0) THEN
              
              write(6,*) ' '
              write(6,*) 'MAIN: Finishing new initial model....'
              write(6,*) ' '
              
           ENDIF

           call INTERCAMBIABOUNDPOS (nx, xzn, xznl, xznr, ny, yzn, yznl, yznr,      &
                                     nz, zzn, zznl, zznr, boundposX, boundposY,     &
                                     boundposZ)


           !Pseudo-gravity for steady jets in pressure decreasing atmospheres.
           IF (basenm(1:1) .EQ. 'J') THEN
                                        
              message4 = 'MAIN: Jet simulations are not implemented.'
              call ERROR ('MAIN3BU', message1, __FILE__, __LINE__)

              message5 = 'Stelar mass load is not implemented.'                      
              IF (stelar .EQ. 1) call ERROR ('SOURCES', message2, __FILE__, __LINE__)

           ENDIF

           
           ! MPI Barrier.
           call PAR_BARRIER ()



           IF (nstep .NE. 0) THEN

              nbegin = nstep + 1
              nend   = nstep + nend

           ENDIF

           suffixd = 'A'//suffix
           outfil  = basenm//'O'//suffix
           rstfil  = basenm//'R'//suffix
           datfil  = basenm//'D'//suffixd


              !r_filed = rst_filed
              !r_filer = rst_file
              !curR_R = 0
              !H_openr = 0
              !H_opend = 0

              !mpi_coordx = mpi_coords(1)
              !mpi_coordy = mpi_coords(2)
              !mpi_coordz = mpi_coords(3)



                 !call RRHDF5_RECORD(0, datfil, rstfil, outfil,                                         &
                 !                   basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                        &
                 !                   H_openr, H_opend, curR_R, r_filer, r_filed,                        &
                 !                   gamma_ad, xzn, yzn, zzn, xznr, xznl, yznr, yznl, zznr, zznl,       &
                 !                   bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                    &
                 !                  gridlx, gridly, gridlz, nx, ny, nz, nsdim,                         &
                 !                   nend, irstrt, nrstrt, trstrt, nout, tout, itstp, iars, irk, niter, &
                 !                   igeomx, igeomy, igeomz, eqs, tmax, cfl, dtini, dt, timej, nstep,   &
                 !                   small, smlrho, smallp, smallu,                                     &
                 !                   mpi_coordx, mpi_coordy, mpi_coordz, physic,CT,mnx6,mny6,mnz6)

                 !stop

        ELSE



           IF (BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'H') THEN

           ! JET RADIUS :: RB = 6 x 10^9 cm

               vwi      = 6.67D-3
               uwi      = 9.05D-10
               shift    = 10.D0
               Rorb     = 500.D0

	   ENDIF

         
           ! Read restart file (R)
           rstfil  = basenm//'R'//suffix
           datfil  = basenm//'D'//suffixd

           IF (cartrank .EQ. 0) THEN
              
              write(6,*)
              write(6,*) 'MAIN: Restarting from....', rstfil
              write(6,*)
              
           ENDIF

           call PAR_BARRIER ()

           ! Allocate working arrays.
           call GETWAMEMORYSPACE (mn6, mn5, mn1, mnx1, mny1, mnz1, nx, ny, nz, &
                                  mpi_dims(1), mpi_dims(2), mpi_dims(3))
         
           call GRID (basenm, nx, ny, nz, nsdim, gridlx, gridly, gridlz,                   &
                      xCPU, yCPU, zCPU, bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                      delx, dely, delz, xznl, xzn, xznr, yznl, yzn, yznr,                  &
                      zznl, zzn, zznr)	  



                           
# ifdef HDF

               IF (nsdim .NE. 3) THEN
                  call ERROR ('MAIN', 'HDF works with 3D simulations', __FILE__, __LINE__)
               ENDIF



               call RST_HDF5(0, datfil, rstfil, outfil,                         &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                        &
               H_openr, H_opend, curR_R, r_filer, r_filed,                        &
               gamma_ad, xzn, yzn, zzn, xznr, xznl, yznr, yznl, zznr, zznl,       &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                    &
               gridlx, gridly, gridlz, nx, ny, nz, nsdim,                         &
               nend, irstrt, nrstrt, trstrt, nout, tout, itstp, iars, irk, niter, &
               igeomx, igeomy, igeomz, eqs, tmax, cfl, dtini, dt, timej, nstep,   &
               small, smlrho, smallp, smallu,                                     &
               mpi_coordx, mpi_coordy, mpi_coordz, physic)


# endif



           ! STAGGERED MAGNETIC FIELD AND CONSERVED VARIABLES
           ! ------------------------------------------------

                 call RESTRT_3D (outRoot, basenm, 0, nx, ny, nz, mnx1, mny1, mnz1,                &
                                 mnx5, mny5, mnz5, mnx6,mny6,mnz6,limiter, nend, irstrt, nrstrt,  &
                                 datfil, rstfil, outfil, tdlfil, suffixd, suffix,                 &
                                 nrst, nout, nout1, nsdim, igeomx, igeomy, igeomz, nstep,         &
                                 itstp, iars, irk, niter, tmax, trstrt, tout, cfl, dtini,         &
                                 small, smlrho, smallp, smallu,                                   &
                                 bndmnz, bndmxz, bndmny, bndmxy, bndmnx,                          &
                                 bndmxx, eqs, tout1, trst, timej, gamma_ad,                       &
                                 gridlx, gridly, gridlz, xznl, yznl, zznl,                        &
                                 xzn, yzn, zzn, xznr, yznr, zznr, physic, max_dims,CT,delx,me,mp)


	   call PAR_BARRIER()


           call INTERCAMBIABOUNDPOS (nx, xzn, xznl, xznr, ny, yzn, yznl, yznr,      &
                                     nz, zzn, zznl, zznr, boundposX, boundposY,     &
                                     boundposZ)


               
           call PAR_BARRIER()

           nbegin = nstep + 1


#ifdef HDF
           rst_file  = irstrt+1
           rst_rec   = 1
           rst_filed = irstrt+1
#endif

           
        ENDIF

             
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !                                             MAIN LOOP
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	
        ! Reset file counter initialization.
        nout1 = 0
        tout1 = 0.D0
        nrst  = 0
        trst  = 0.D0

        ! Runge-Kutta parameters.
        IF (irk .EQ. 2) THEN
           
           alpha_R = 2.D0
           beta_R  = 1.D0
           
        ELSE IF (irk .EQ. 3) THEN
           
           alpha_R = 4.D0
           beta_R  = 3.D0
           
        ELSE IF (irk .EQ. 4) THEN

	   irk = irk + 1

	ENDIF
        
        ! Start main loop
        nstep   = nbegin
        NoSalir = .TRUE.
	final   = 0
 
        call PAR_BARRIER ()

	
        DO WHILE (nstep .LE. nend .AND. NoSalir)


           timej = timej + dt
           nout1 = nout1 + 1
           tout1 = tout1 + dt
           nrst  = nrst  + 1
           trst  = trst  + dt

           ! Break the LOOP if t>tmax.
           IF (timej .GT. tmax) THEN              
              NoSalir = .FALSE.            
	      final   = 1
           ELSE

              IF (nsdim .EQ. 2) THEN

!$OMP PARALLEL DEFAULT(shared)    &
!$OMP PRIVATE(I,J,K)

                 DO k=1,nz
!$OMP DO SCHEDULE(static)
                    DO j=1,ny
                       DO i=1,nx
                          CT(I,J,K)%obx = CT(I,J,K)%bx
               		  CT(I,J,K)%oby = CT(I,J,K)%by
                       ENDDO
                    ENDDO
!$OMP END DO
                 ENDDO


                 DO k=1,nz
!$OMP DO SCHEDULE(static)
                    DO i=1,nx
                       CT(I,NY+1,K)%oby = CT(I,NY+1,K)%by
                    ENDDO
!$OMP END DO
                 ENDDO

                 DO k=1,nz
!$OMP DO SCHEDULE(static)
                    DO j=1,ny
                       CT(NX+1,J,K)%obx = CT(NX+1,J,K)%bx
                    ENDDO
!$OMP END DO
                 ENDDO
!$OMP END PARALLEL

              ELSE IF (nsdim .EQ. 3) THEN


!$OMP PARALLEL DEFAULT(shared)    &
!$OMP PRIVATE(I,J,K)

!$OMP DO SCHEDULE(static)
                 DO k=1,nz
                    DO j=1,ny
                       DO i=1,nx
                          CT(I,J,K)%obx = CT(I,J,K)%bx
                          CT(I,J,K)%oby = CT(I,J,K)%by
	                  CT(I,J,K)%obz = CT(I,J,K)%bz
                       ENDDO
                    ENDDO
                 ENDDO
!$OMP END DO


!$OMP DO SCHEDULE(static)
                 DO k=1,nz
                    DO i=1,nx
                       CT(I,NY+1,K)%oby = CT(I,NY+1,K)%by
                    ENDDO
                 ENDDO
!$OMP END DO

!$OMP DO SCHEDULE(static)
                 DO k=1,nz
                    DO j=1,ny
                       CT(NX+1,J,K)%obx = CT(NX+1,J,K)%bx
                    ENDDO
                 ENDDO
!$OMP END DO


!$OMP DO SCHEDULE(static)
                 DO j=1,ny
                    DO i=1,nx
                       CT(I,J,NZ+1)%obz = CT(I,J,NZ+1)%bz
                    ENDDO
                 ENDDO
!$OMP END DO


!$OMP END PARALLEL


	      ENDIF


              ! Save variables for Runge-Kutta evolution.

	      IF (nsdim .LE. 2) THEN

!$OMP PARALLEL DEFAULT (shared)     &
!$OMP PRIVATE (i,j,k)

                 DO k= 1, nz
!$OMP DO SCHEDULE(static)
                    DO j= 1, ny
                       DO i= 1, nx
                     
                          physic(i,j,k)%omomenx  = physic(i,j,k)%momenx
                       	  physic(i,j,k)%omomeny  = physic(i,j,k)%momeny
                       	  physic(i,j,k)%omomenz  = physic(i,j,k)%momenz
                          physic(i,j,k)%oenergy  = physic(i,j,k)%energy
                          physic(i,j,k)%odenstw  = physic(i,j,k)%denstw
                          physic(i,j,k)%odenstew = physic(i,j,k)%denstew
                          physic(i,j,k)%opres    = physic(i,j,k)%pres
                          physic(i,j,k)%odensty  = physic(i,j,k)%densty
                          physic(i,j,k)%otracer  = physic(i,j,k)%tracer
                       
                          physic(i,j,k)%obxc     = physic(i,j,k)%bxc
                          physic(i,j,k)%obyc     = physic(i,j,k)%byc
                          physic(i,j,k)%obzc     = physic(i,j,k)%bzc

                          physic(i,j,k)%ovelx     = physic(i,j,k)%velx
                          physic(i,j,k)%ovely     = physic(i,j,k)%vely
                          physic(i,j,k)%ovelz     = physic(i,j,k)%velz

                 
                       ENDDO
                    ENDDO
!$OMP END DO
                 ENDDO

!$OMP END PARALLEL

              ELSE IF (nsdim .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (shared)     &
!$OMP PRIVATE (i,j,k)

!$OMP DO SCHEDULE(static)
                 DO k= 1, nz
                    DO j= 1, ny
                       DO i= 1, nx

                          physic(i,j,k)%omomenx  = physic(i,j,k)%momenx
                          physic(i,j,k)%omomeny  = physic(i,j,k)%momeny
                          physic(i,j,k)%omomenz  = physic(i,j,k)%momenz
                          physic(i,j,k)%oenergy  = physic(i,j,k)%energy
                          physic(i,j,k)%odenstw  = physic(i,j,k)%denstw
                          physic(i,j,k)%odenstew = physic(i,j,k)%denstew
                          physic(i,j,k)%opres    = physic(i,j,k)%pres
                          physic(i,j,k)%odensty  = physic(i,j,k)%densty
                          physic(i,j,k)%otracer  = physic(i,j,k)%tracer

                          physic(i,j,k)%obxc     = physic(i,j,k)%bxc
                          physic(i,j,k)%obyc     = physic(i,j,k)%byc
                          physic(i,j,k)%obzc     = physic(i,j,k)%bzc

                          physic(i,j,k)%ovelx     = physic(i,j,k)%velx
                          physic(i,j,k)%ovely     = physic(i,j,k)%vely
                          physic(i,j,k)%ovelz     = physic(i,j,k)%velz


                       ENDDO
                    ENDDO
                 ENDDO
!$OMP END DO

!$OMP END PARALLEL

	      ENDIF


              ! Runge-Kutta loop

              DO rkin= 1, irk
 
		
                 call INTERCAMBIABOUNDPHYSIC (physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy)


	         IF (nsdim .LE. 2) THEN

!$OMP PARALLEL DEFAULT (SHARED)  &
!$OMP PRIVATE (I,J,K)

                    DO k = 1, nz
!$OMP DO SCHEDULE (STATIC)
                       DO j = 1, ny
                          DO i = 1, nx

                             physic(i,j,k)%bxc2 = physic(i,j,k)%bxc
                             physic(i,j,k)%byc2 = physic(i,j,k)%byc
                             physic(i,j,k)%bzc2 = physic(i,j,k)%bzc
                             physic(i,j,k)%tracer2 = physic(i,j,k)%tracer


                          ENDDO
                       ENDDO
!$OMP END DO
                    ENDDO
!$OMP END PARALLEL

                 ELSE IF (nsdim .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED)  &
!$OMP PRIVATE (I,J,K)

!$OMP DO SCHEDULE (STATIC)
                    DO k = 1, nz
                       DO j = 1, ny
                          DO i = 1, nx

                             physic(i,j,k)%bxc2 = physic(i,j,k)%bxc
                             physic(i,j,k)%byc2 = physic(i,j,k)%byc
                             physic(i,j,k)%bzc2 = physic(i,j,k)%bzc
		             physic(i,j,k)%tracer2 = physic(i,j,k)%tracer

                          ENDDO
                       ENDDO
                    ENDDO
!$OMP END DO

!$OMP END PARALLEL

                 ENDIF

		 IF (nsdim .LE. 2) THEN

                     call SWEEPS_2D (basenm, limiter, irk, rkin, iars, iccv,                            & 
                                     igeomx, igeomy, igeomz, bndmnx, bndmxx, bndmny, bndmxy,            &
                                     bndmnz, bndmxz, nx, mnx5, mnx6, nsdim, ny, mny5, mny6,             &
                                     eqs, nz, mnz5, mnz6, stelar, cool, gamma_ad, alpha_R, beta_R, dt,  &
                              	     me, mp, xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,         &
                                     boundPosX, boundPosY, boundPosZ, grav, physic, CT, uwi,vwi,shift,rorb)

                     call GETPRFQ2DNR (nx, ny, nz, mnx1, mnx5, mny5, mnz5, eqs, niter,                    &
                                     gamma_ad, small, smallp, me, mp, physic)

		 ELSE IF (nsdim .EQ. 3) THEN

                     call SWEEPS_3D (basenm, limiter, irk, rkin, iars, iccv,                            &
                                     igeomx, igeomy, igeomz, bndmnx, bndmxx, bndmny, bndmxy,            &
                                     bndmnz, bndmxz, nx, mnx5, mnx6, nsdim, ny, mny5, mny6,             &
                                     eqs, nz, mnz5, mnz6, stelar, cool, gamma_ad, alpha_R, beta_R, dt,  &
                                     me, mp, xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,         &
                                     boundPosX, boundPosY, boundPosZ, grav, physic, CT, uwi, vwi, shift,&
                                     Rorb, kb, rho_0a, R_b)


                     call GETPRFQ3DNR (nx, ny, nz, mnx1, mnx5, mny5, mnz5, eqs, niter,                    &
                                       gamma_ad, small, smallp, me, mp, physic)




		 ENDIF





                 IF (nsdim .EQ. 2 .AND. iccv .EQ. 2) THEN

!$OMP PARALLEL DEFAULT(shared)       & 
!$OMP PRIVATE(I,J,K,BCT2,B2,VEL2)  

                  DO k = 1, nz
!$OMP DO SCHEDULE (STATIC)
                     DO j = 1, ny
                        DO i = 1, nx

                           VEL2 = physic(I,J,K)%velx**2 + physic(I,J,K)%vely**2
                           B2   = physic(I,J,K)%bxzn**2 + physic(I,J,K)%byzn**2
                           BCT2 = physic(I,J,K)%bxct**2 + physic(I,J,K)%byct**2

                           physic(I,J,K)%momenx = physic(I,J,K)%momenx - physic(I,J,K)%velx*(B2-BCT2) +                                                & 
                                                  physic(I,J,K)%velx*(physic(I,J,K)%bxzn**2-physic(I,J,K)%bxct**2) +                                   & 
                                                  physic(I,J,K)%vely*(physic(I,J,K)%byzn*physic(I,J,K)%bxzn-physic(I,J,K)%byct*physic(I,J,K)%bxct) +   &
                                                  physic(I,J,K)%velz*physic(I,J,K)%bzc*(physic(I,J,K)%bxzn-physic(I,J,K)%bxct)

                           physic(I,J,K)%momeny = physic(I,J,K)%momeny - physic(I,J,K)%vely*(B2-BCT2) +                                                & 
                                                  physic(I,J,K)%velx*(physic(I,J,K)%byzn*physic(I,J,K)%bxzn-physic(I,J,K)%byct*physic(I,J,K)%bxct) +   &
                                                  physic(I,J,K)%vely*(physic(I,J,K)%byzn**2-physic(I,J,K)%byct**2) +                                   & 
                                                  physic(I,J,K)%velz*physic(I,J,K)%bzc*(physic(I,J,K)%byzn-physic(I,J,K)%byct)

                           physic(I,J,K)%momenz = physic(I,J,K)%momenz - physic(I,J,K)%velz*(B2-BCT2) + physic(I,J,K)%velx*physic(I,J,K)%bzc *         &
                                                 (physic(I,J,K)%bxzn-physic(I,J,K)%bxct) + physic(I,J,K)%vely*physic(I,J,K)%bzc*(physic(I,J,K)%byzn-physic(I,J,K)%byct)
 
                           physic(I,J,K)%energy = physic(I,J,K)%energy - 0.5D0*VEL2*(B2-BCT2) +       & 
                                                  0.5D0*(physic(I,J,K)%velx*physic(I,J,K)%bxzn +      & 
                                                         physic(I,J,K)%vely*physic(I,J,K)%byzn +      &
                                                         physic(I,J,K)%velz*physic(I,J,K)%bzc)**2 -   & 
                                                  0.5D0*(physic(I,J,K)%velx*physic(I,J,K)%bxct +      & 
                                                         physic(I,J,K)%vely*physic(I,J,K)%byct +      & 
                                                         physic(I,J,K)%velz*physic(I,J,K)%bzc)**2

                        ENDDO
                     ENDDO
!$OMP END DO
                  ENDDO

!$OMP END PARALLEL


                   call GETPRFQ2DNR (nx, ny, nz, mnx1, mnx5, mny5, mnz5, eqs, niter,      &
                                     gamma_ad, small, smallp, me, mp, physic)

	         
                ELSE IF (nsdim .EQ. 3 .AND. iccv .EQ. 2) THEN


!$OMP PARALLEL DEFAULT(shared)       &
!$OMP PRIVATE(I,J,K,BCT2,B2,VEL2)


!$OMP DO SCHEDULE (STATIC)
                  DO k = 1, nz
                     DO j = 1, ny
                        DO i = 1, nx

                           VEL2 = physic(I,J,K)%velx**2 + physic(I,J,K)%vely**2 + physic(I,J,K)%velz**2
                           B2   = physic(I,J,K)%bxzn**2 + physic(I,J,K)%byzn**2 + physic(I,J,K)%bzzn**2
                           BCT2 = physic(I,J,K)%bxct**2 + physic(I,J,K)%byct**2 + physic(I,J,K)%bzct**2

                           physic(I,J,K)%momenx = physic(I,J,K)%momenx - physic(I,J,K)%velx*(B2-BCT2) +                                                &
                                                  physic(I,J,K)%velx*(physic(I,J,K)%bxzn**2-physic(I,J,K)%bxct**2) +                                   &
                                                  physic(I,J,K)%vely*(physic(I,J,K)%byzn*physic(I,J,K)%bxzn-physic(I,J,K)%byct*physic(I,J,K)%bxct) +   &
                                                  physic(I,J,K)%velz*(physic(I,J,K)%bzzn*physic(I,J,K)%bxzn-physic(I,J,K)%bzct*physic(i,j,k)%bxct)

                           physic(I,J,K)%momeny = physic(I,J,K)%momeny - physic(I,J,K)%vely*(B2-BCT2) +                                                &
                                                  physic(I,J,K)%velx*(physic(I,J,K)%byzn*physic(I,J,K)%bxzn-physic(I,J,K)%byct*physic(I,J,K)%bxct) +   &
                                                  physic(I,J,K)%vely*(physic(I,J,K)%byzn**2-physic(I,J,K)%byct**2) +                                   &
                                                  physic(I,J,K)%velz*(physic(I,J,K)%byzn*physic(I,J,K)%bzzn-physic(I,J,K)%byct*physic(I,J,K)%bzct)

                           physic(I,J,K)%momenz = physic(I,J,K)%momenz - physic(I,J,K)%velz*(B2-BCT2) +                                                & 
                                                  physic(I,J,K)%velx*(physic(I,J,K)%bzzn*physic(I,J,K)%bxzn-physic(I,J,K)%bzct*physic(I,J,K)%bxct) +   &
                                                  physic(I,J,K)%vely*(physic(I,J,K)%bzzn*physic(I,J,K)%byzn-physic(I,J,K)%bzct*physic(I,J,K)%byct) +   &
                                                  physic(I,J,K)%velz*(physic(I,J,K)%bzzn**2-physic(I,J,K)%bzct**2)


                           physic(I,J,K)%energy = physic(I,J,K)%energy - 0.5D0*VEL2*(B2-BCT2) +         &
                                                  0.5D0*(physic(I,J,K)%velx*physic(I,J,K)%bxzn +        &
                                                         physic(I,J,K)%vely*physic(I,J,K)%byzn +        &
                                                         physic(I,J,K)%velz*physic(I,J,K)%bzzn)**2 -    &
                                                  0.5D0*(physic(I,J,K)%velx*physic(I,J,K)%bxct +        &
                                                         physic(I,J,K)%vely*physic(I,J,K)%byct +        &
                                                         physic(I,J,K)%velz*physic(I,J,K)%bzct)**2
			  
                        ENDDO
                     ENDDO
                  ENDDO
!$OMP END DO
!$OMP END PARALLEL


                   call GETPRFQ3DNR (nx, ny, nz, mnx1, mnx5, mny5, mnz5, eqs, niter,      &
                                     gamma_ad, small, smallp, me, mp, physic)



		ENDIF	

              ENDDO




           ENDIF


	              
           ! ==========
           ! Write Restart/Output files.
           ! ==========

           
           ! Write output file for 1D simulations.

           IF ((nsdim .EQ. 1) .AND. ((nout1 .GE. nout) .OR. (tout1 .GE. tout))) THEN
                      
                           
	      call PLTOUT (outRoot, basenm, nsdim, nx, ny, nz,     &
                           mnx5, mny5, mnz5, xzn, physic,          &
                   	   rstfil, outfil, tdlfil, datfil,         &
                   	   nout1, tout1, suffix, suffixd,          &
                   	   limiter, nend, irstrt, nrstrt,          &
                   	   nout, itstp, iars, irk, niter,          &
                   	   igeomx, igeomy, igeomz, eqs,            &
                   	   tmax, trstrt, tout, cfl, dtini,         &
                   	   small, smlrho, smallp, smallu,          &
                   	   gridlx, gridly, gridlz)
 
           ENDIF

           
           ! Write restart file (R) for 2D simulations.

           IF (nsdim .EQ. 2) THEN
           
              IF ((nrst .GE. nrstrt) .OR. (trst .GE. trstrt)) THEN
                
                 call RESTRT_2D (outRoot, basenm, 2, nx, ny, nz, mnx1, mny1, mnz1,                  &
                                 mnx5, mny5, mnz5, mnx6,mny6,mnz6,limiter, nend, irstrt, nrstrt,    &
                                 datfil, rstfil, outfil, tdlfil, suffixd, suffix,                   &
                                 nrst, nout, nout1, nsdim, igeomx, igeomy, igeomz, nstep,           &
                                 itstp, iars, irk, niter, tmax, trstrt, tout, cfl, dtini,           &
                                 small, smlrho, smallp, smallu,                                     &
                                 bndmnz, bndmxz, bndmny, bndmxy, bndmnx,                            &
                                 bndmxx, eqs, tout1, trst, timej, gamma_ad,                         &
                                 gridlx, gridly, gridlz, xznl, yznl, zznl,                          &
                                 xzn, yzn, zzn, xznr, yznr, zznr, physic, max_dims,CT)
              ENDIF

              IF ((nout1 .GE. nout) .OR. (tout1 .GE. tout)) THEN
                
                 call RESTRT_2D (outRoot, basenm, 2, nx, ny, nz, mnx1, mny1, mnz1,                  &
                                 mnx5, mny5, mnz5, mnx6,mny6,mnz6,limiter, nend, irstrt, nrstrt,    &
                      		 datfil, rstfil, outfil, tdlfil, suffixd, suffix,                   &
                                 nrst, nout, nout1, nsdim, igeomx, igeomy, igeomz, nstep,           &
                                 itstp, iars, irk, niter, tmax, trstrt, tout, cfl, dtini,           &
                                 small, smlrho, smallp, smallu,                                     &
                                 bndmnz, bndmxz, bndmny, bndmxy, bndmnx,                            &
                                 bndmxx, eqs, tout1, trst, timej, gamma_ad,                         &
                                 gridlx, gridly, gridlz, xznl, yznl, zznl,                          &
                                 xzn, yzn, zzn, xznr, yznr, zznr, physic, max_dims,CT)
                 
                 print*, 'Writing model...'           
              
              ENDIF

           ELSE IF (nsdim .EQ. 3) THEN

#ifdef HDF

              r_filed = rst_filed
      	      r_filer = rst_file
              curR_R = 0
              H_openr = 0
              H_opend = 0

              mpi_coordx = mpi_coords(1)
              mpi_coordy = mpi_coords(2)
              mpi_coordz = mpi_coords(3)


              IF ((nrst .GE. nrstrt) .OR. (trst .GE. trstrt)) THEN

                 call RRHDF5_RECORD(0, datfil, rstfil, outfil,                                         &
               			    basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                        &
               			    H_openr, H_opend, curR_R, r_filer, r_filed,                        &
               			    gamma_ad, xzn, yzn, zzn, xznr, xznl, yznr, yznl, zznr, zznl,       &
               			    bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                    &
               			    gridlx, gridly, gridlz, nx, ny, nz, nsdim,                         &
               			    nend, irstrt, nrstrt, trstrt, nout, tout, itstp, iars, irk, niter, &
               			    igeomx, igeomy, igeomz, eqs, tmax, cfl, dtini, dt, timej, nstep,   &
               			    small, smlrho, smallp, smallu,                                     &
               			    mpi_coordx, mpi_coordy, mpi_coordz, physic,CT,mnx6,mny6,mnz6)


                 call RESTRT_3D (outRoot, basenm, 1, nx, ny, nz, mnx1, mny1, mnz1,                &
                                 mnx5, mny5, mnz5, mnx6,mny6,mnz6,limiter, nend, irstrt, nrstrt,  &
                                 datfil, rstfil, outfil, tdlfil, suffixd, suffix,                 &
                                 nrst, nout, nout1, nsdim, igeomx, igeomy, igeomz, nstep,         &
                                 itstp, iars, irk, niter, tmax, trstrt, tout, cfl, dtini,         &
                                 small, smlrho, smallp, smallu,                                   &
                                 bndmnz, bndmxz, bndmny, bndmxy, bndmnx,                          &
                                 bndmxx, eqs, tout1, trst, timej, gamma_ad,                       &
                                 gridlx, gridly, gridlz, xznl, yznl, zznl,                        &
                                 xzn, yzn, zzn, xznr, yznr, zznr, physic, max_dims,CT,delx,me,mp)


              ENDIF

              IF ((nout1 .GE. nout) .OR. (tout1 .GE. tout)) THEN

                 call RRHDF5_RECORD(0, datfil, rstfil, outfil,                                         &
                                    basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                        &
                                    H_openr, H_opend, curR_R, r_filer, r_filed,                        &
                                    gamma_ad, xzn, yzn, zzn, xznr, xznl, yznr, yznl, zznr, zznl,       &
                                    bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                    &
                                    gridlx, gridly, gridlz, nx, ny, nz, nsdim,                         &
                                    nend, irstrt, nrstrt, trstrt, nout, tout, itstp, iars, irk, niter, &
                                    igeomx, igeomy, igeomz, eqs, tmax, cfl, dtini, dt, timej, nstep,   &
                                    small, smlrho, smallp, smallu,                                     &
                                    mpi_coordx, mpi_coordy, mpi_coordz, physic,CT,mnx6,mny6,mnz6)


                 call RESTRT_3D (outRoot, basenm, 1, nx, ny, nz, mnx1, mny1, mnz1,                &
                                 mnx5, mny5, mnz5, mnx6,mny6,mnz6,limiter, nend, irstrt, nrstrt,  &
                                 datfil, rstfil, outfil, tdlfil, suffixd, suffix,                 &
                                 nrst, nout, nout1, nsdim, igeomx, igeomy, igeomz, nstep,         &
                                 itstp, iars, irk, niter, tmax, trstrt, tout, cfl, dtini,         &
                                 small, smlrho, smallp, smallu,                                   &
                                 bndmnz, bndmxz, bndmny, bndmxy, bndmnx,                          &
                                 bndmxx, eqs, tout1, trst, timej, gamma_ad,                       &
                                 gridlx, gridly, gridlz, xznl, yznl, zznl,                        &
                                 xzn, yzn, zzn, xznr, yznr, zznr, physic, max_dims,CT,delx,me,mp)





                 rst_file = r_filer
        
                 nout1 = 0
         	 tout1 = 0.0

                 IF (cartrank .EQ. 0) THEN
              
                    write(6,*)
                    write(6,*) 'MAIN: Writing model....'
                    write(6,*)
              
                 ENDIF
              
              ENDIF

#else

              IF ((nrst .GE. nrstrt) .OR. (trst .GE. trstrt)) THEN

                 call RESTRT_3D (outRoot, basenm, 2, nx, ny, nz, mnx1, mny1, mnz1,                &
                                 mnx5, mny5, mnz5, mnx6,mny6,mnz6,limiter, nend, irstrt, nrstrt,  &
                                 datfil, rstfil, outfil, tdlfil, suffixd, suffix,                 &
                                 nrst, nout, nout1, nsdim, igeomx, igeomy, igeomz, nstep,         &
                                 itstp, iars, irk, niter, tmax, trstrt, tout, cfl, dtini,         &
                                 small, smlrho, smallp, smallu,                                   &
                                 bndmnz, bndmxz, bndmny, bndmxy, bndmnx,                          &
                                 bndmxx, eqs, tout1, trst, timej, gamma_ad,                       &
                                 gridlx, gridly, gridlz, xznl, yznl, zznl,                        &
                                 xzn, yzn, zzn, xznr, yznr, zznr, physic, max_dims,CT,delx,me,mp)
              ENDIF

              IF ((nout1 .GE. nout) .OR. (tout1 .GE. tout)) THEN

                 call RESTRT_3D (outRoot, basenm, 2, nx, ny, nz, mnx1, mny1, mnz1,                &
                                 mnx5, mny5, mnz5, mnx6,mny6,mnz6,limiter, nend, irstrt, nrstrt,  &
                                 datfil, rstfil, outfil, tdlfil, suffixd, suffix,                 &
                                 nrst, nout, nout1, nsdim, igeomx, igeomy, igeomz, nstep,         &
                                 itstp, iars, irk, niter, tmax, trstrt, tout, cfl, dtini,         &
                                 small, smlrho, smallp, smallu,                                   &
                                 bndmnz, bndmxz, bndmny, bndmxy, bndmnx,                          &
                                 bndmxx, eqs, tout1, trst, timej, gamma_ad,                       &
                                 gridlx, gridly, gridlz, xznl, yznl, zznl,                        &
                                 xzn, yzn, zzn, xznr, yznr, zznr, physic, max_dims,CT,delx,me,mp)

                 IF (cartrank .EQ. 0) THEN
              
                    write(6,*)
                    write(6,*) 'MAIN: Writing model....'
                    write(6,*)
              
                 ENDIF

              ENDIF

# endif


	   ENDIF

		
	    IF (nsdim .LE. 2) THEN

                call TSTEP_2D (limiter, nx, ny, nz, mnx5, mny5, mnz5,mnx6, mny6, mnz6, nsdim,       &
                               itstp, nstep, eqs, timej, cfl, me, mp, gamma_ad,                     &
                               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                      &
                               xznl, xznr, xzn, yznl, yznr, yzn, zznl, zznr, zzn,                   &
                               physic, boundPosX, boundPosY, boundPosZ, dt,CT)

            ELSE IF (nsdim .EQ. 3) THEN

                call TSTEP_3D (limiter, nx, ny, nz, mnx5, mny5, mnz5,mnx6, mny6, mnz6, nsdim,       &
                               itstp, nstep, eqs, timej, cfl, me, mp, gamma_ad,                     &
                               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                      &
                               xznl, xznr, xzn, yznl, yznr, yzn, zznl, zznr, zzn,                   &
                               physic, boundPosX, boundPosY, boundPosZ, dt,CT, basenm, m_0, c, kb, rho_0a, R_b)

	   ENDIF

           nstep = nstep + 1

           call PAR_BARRIER ()

	   ! SS433 SIMULATION
           ! RESET THE VALUES INSIDE THE INJECTION SPHERE.
	   ! --------------------------------------------

           IF (BASENM(1:2) .EQ. 'S4' .OR. BASENM .EQ. 'SYBW') THEN

               call SPHERE (basenm, nx, ny, nz, mnx5, mny5, mnz5, mnx6, mny6, mnz6, eqs, iars,                            &
                            xzn, xznl,xznr,yznr, delx, yzn, yznl, zzn,zznl, physic, gamma_ad, gridlx, gridly, gridlz, CT, &
                            me, mp, m_0, kb, c, uwi, vwi, shift, Rorb, xCPU, yCPU, zCPU,ic,jc,kc,ncl,timej,pc, dt)
           ENDIF

           !!call PAR_BARRIER ()

        ENDDO

        call PAR_BARRIER ()
        
        IF (cartrank .NE. 0) THEN
           
           call FREEWAMEMORYSPACE ()
         
           deallocate (xznl, xzn, xznr, &
                       yznl, yzn, yznr, &
                       zznl, zzn, zznr, &
                       grav )
           
           deallocate (physic)
	   deallocate (CT)
           deallocate (boundPosX, boundPosY, boundPosZ)
           
        ENDIF

	IF (cartrank .EQ. 0) THEN

           ! CPU TIME
           ! ===================

           call SYSTEM_CLOCK(iTime2)

           print*, '---------------------------------------'
           print*, 'Time =',real(iTime2-iTime1)/real(RATE)
           print*, '---------------------------------------'

           print*, 'MAIN: Test finished succesfully!.......'

	ENDIF
 
        call PAR_END ()



        !-------------------------------------------------------------------------- END

END PROGRAM RATPENAT
