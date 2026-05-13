!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************

! NAME:
!       RESTRT_2D

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
!       Read and Write data files.

! ACADEMIC REFERENCE:


! ===============================================================================      
SUBROUTINE RESTRT_2D (path, basenm, ird, nx, ny, nz, mnx1, mny1, mnz1,                &
                      mnx5, mny5, mnz5, mnx6,mny6,mnz6,limiter, nend, irstrt, nrstrt, &
                      datfil, rstfil, outfil, tdlfil, suffixd, suffix,                &
                      nrst, nout, nout1, nsdim, igeomx, igeomy, igeomz, nstep,        &
		      itstp, iars, irk, niter, tmax, trstrt, tout, cfl, dtini,        &
		      small, smlrho, smallp, smallu,                                  &
                      bndmnz, bndmxz, bndmny, bndmxy, bndmnx,                         &
                      bndmxx, eqs, tout1, trst, timej, gamma_ad,                      &
                      gridlx, gridly, gridlz, xznl, yznl, zznl,                       &
                      xzn, yzn, zzn, xznr, yznr, zznr, physic, max_dims,CT)
! ===============================================================================
                   
        USE parallel
        USE tipos
        USE workarrays
      
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(256):: path
        character(4)  :: basenm

        integer(is):: ird, nx, ny, nz, mnx1, mny1, mnz1, max_dims
	integer(is) :: mnx5, mny5, mnz5,mnx6,mny6,mnz6 	



        !------------------------------------------------------ Output Variables

        character(8):: datfil
        character(7):: rstfil, outfil, tdlfil
        character(3):: suffixd
        character(2):: suffix

        integer(is) :: nrst, nout1,  nout, nsdim, nend, irstrt, nrstrt,       &
                       igeomx, igeomy, igeomz, nstep, bndmnz, bndmxz,         &
                       bndmny, bndmxy, bndmxx, bndmnx, eqs,trstrt

        real(rs) :: tout1, trst, timej,                     &
                    gamma_ad, gridlx, gridly, gridlz
                         
        real(rs) :: xznl(mnx1), xzn(nx), xznr(nx), &
                    yznl(mny1), yzn(ny), yznr(ny), &
                     zznl(mnz1), zzn(mnz1), zznr(mnz1)

        type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)
	type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)

	character(4)::limiter

	integer(is) :: itstp, iars, irk, niter

	real(rs) :: tmax, tout, cfl, dtini, small, smlrho, smallp, smallu

        !------------------------------------------------------- Local Variables

        character(256) :: fich

        character(4), external :: int2ch4
        character(27) :: message1

        integer(is)    :: tid, i, j, k, ierr

        real(rs), pointer :: rho_s(:), xe_s(:), tem_s(:),   & 
                             eps_s(:), p_s(:), gam_s(:),    &
                             h_s(:), cs_s(:),               &
                             dpdrh_s(:), dpdeps_s(:)

        real(rs), pointer :: scrch1(:), scrch2(:)


        real(rs) :: b40, b41, b42, b43, b44

        !-------------------------------------------------------------------------- BEGIN

        tid = 1

        !  1-vectors
        rho_s     => wa(tid)%wa1rest100
        eps_s     => wa(tid)%wa1rest101
        p_s       => wa(tid)%wa1rest102
        dpdrh_s   => wa(tid)%wa1rest103
        dpdeps_s  => wa(tid)%wa1rest104
        h_s       => wa(tid)%wa1rest105
        cs_s      => wa(tid)%wa1rest106
      
        ! -4:mn6 vectors
        scrch1 => wa(tid)%wa6rest100
        scrch2 => wa(tid)%wa6rest101

        ! ----------------------------------------------
        ! ird .EQ. 0    :: Read  RESTART file.
        ! ird .EQ. 1    :: Write RESTART file.
        ! ird .EQ. 2    :: Write RESTART file (datfiles).
        ! ----------------------------------------------

        
        ! =======================
        !     READ FILE
        ! =======================

        IF (ird .EQ. 0) THEN   

           ! Read Restart file.
           fich = TRIM(path)//rstfil//'_'//int2ch4(cartrank)
           open (unit= 4, file= fich, form= 'unformatted', iostat= ierr)

           IF (ierr .NE. 0) THEN
              call ERROR ('RESTRT_2D', 'Error opening file.', __FILE__, __LINE__)
           ENDIF

           read(4) suffixd

           ! Read Input parameters in Restart file.
           read(4) xznl, yznl, zznl, timej, gamma_ad,                 &
                   bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                   gridlx, gridly, gridlz, nx, ny, nz,                &
                   nsdim, nstep, igeomx, igeomy, igeomz, eqs
               
           ! Read Physic variables in Restart file.
           read(4) (((physic(i,j,k)%velx, physic(i,j,k)%vely, physic(i,j,k)%velz,    &
                      physic(i,j,k)%densty, physic(i,j,k)%denstye,                   &
                      physic(i,j,k)%pres, physic(i,j,k)%tracer,                      &
                      physic(i,j,k)%bxc, physic(i,j,k)%byc, physic(i,j,k)%bzc,       &
                      i=1, nx+1), j=1, ny+1), k=1, nz+1)

                      
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    ! This subroutine will be used only with tests (1D,2D).
                    IF (eqs .EQ. 2) THEN

                       message1='Test works with ideal gas.'
                       call ERROR ('RESTRT_2D', message1, __FILE__, __LINE__)

                    ELSE

                       call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                       physic(i,j,k)%eps = eps_s(1)

                    ENDIF

                    ! Lorentz Factor W
                    scrch1(i)   = 1.0

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               
                    b44= (physic(i,j,k)%bxc*physic(i,j,k)%bxc +physic(i,j,k)%byc*physic(i,j,k)%byc +   &
                          physic(i,j,k)%bzc*physic(i,j,k)%bzc)/(scrch1(i))**2                      +   &
                         (physic(i,j,k)%velx*physic(i,j,k)%bxc                                     +   &
                          physic(i,j,k)%vely*physic(i,j,k)%byc                                     +   &
                          physic(i,j,k)%velz*physic(i,j,k)%bzc)**2

                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO

        ! ===================
        ! 1 :: Write (rstfil)
        ! ===================  

        ELSE IF (ird .EQ. 1) THEN
           
           ! Write Restart
           fich = TRIM(path)//rstfil//'_'//int2ch4(cartrank)
           
           open (unit = 4, file= fich, form= 'unformatted', iostat= ierr)
           IF (ierr .NE. 0) THEN
              call ERROR ('RESTRT_2D', 'Error opening file', __FILE__, __LINE__)
           ENDIF

           write(4) suffixd
           
           ! Input parameters.
           write(4) xznl, yznl, zznl, timej, gamma_ad,                 &
                    bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                    gridlx, gridly, gridlz, nx, ny, nz,                &
                    nsdim, nstep, igeomx, igeomy, igeomz, eqs
               
           ! Physic variables.
           write(4) (((physic(i,j,k)%velx, physic(i,j,k)%vely, physic(i,j,k)%velz,   &
                      physic(i,j,k)%densty, physic(i,j,k)%denstye,                   &
                      physic(i,j,k)%pres, physic(i,j,k)%tracer,                      &
                      physic(i,j,k)%bxc, physic(i,j,k)%byc, physic(i,j,k)%bzc,       &
                      i=1, nx+1), j=1, ny+1), k=1, nz+1)

           close(4)

           ! Reset counters.
           nrst = 0
           trst = 0.0

           ! This subroutine is overwritting the file, so we don't need to call FILNAM.

        ! ===================
        ! 2 :: Write (datfil)
        ! =================== 

        ELSE IF (ird .EQ. 2) THEN

           IF (cartrank .EQ. 0) THEN

              ! Write Restart
              fich = TRIM(path)//datfil//'_'//int2ch4(cartrank)

              open (unit = 12, file= fich, form= 'unformatted', iostat= ierr)
              IF (ierr .NE. 0) THEN
                 call ERROR ('RESTRT_2D', 'Error opening file', __FILE__, __LINE__)
              ENDIF

              write(12) nend, irstrt, nrstrt, nout,           &
                        itstp, iars, irk, niter, nsdim, nx, ny, nz, igeomx, igeomy, igeomz,     &
                        eqs, tmax, trstrt, tout, cfl, dtini, small, smlrho, smallp, smallu,     &
                        gridlx, gridly, gridlz, bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz, &
		        xznl, xzn, xznr, yznl, yzn, yznr, gamma_ad


              write (12) ((physic(i,j,1)%pres,                            &
                           physic(i,j,1)%densty, physic(i,j,1)%denstye,   &
                           physic(i,j,1)%velx, physic(i,j,1)%vely,        &
                           physic(i,j,1)%velz, physic(i,j,1)%bxc,         &
                           physic(i,j,1)%byc, physic(i,j,1)%bzc,          &
                           physic(i,j,1)%tracer, i=1, nx), j=1, ny)

	      write (12) ((CT(i,j,1)%bx, CT(i,j,1)%by, CT(i,j,1)%bz, i=1, nx+1), j=1, ny+1)

              close(12)

           ENDIF


           ! Reset counters.
           nout1 = 0
           tout1 = 0.0

           call FILNAM ('d', basenm, datfil, rstfil, outfil, tdlfil, suffixd, suffix)
	   print*, suffixd

        ELSE
           call ERROR ('RESTRT_2D', 'Wrong ird parameter...', __FILE__, __LINE__)
        ENDIF

        return

        !-------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE RESTRT_2D
! =============================================================================== 
