!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************

! NAME:
!       INIT

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
!       Define initial data for a collection of standard 1D,2D and 3D test 
!       problems and jet simulations.

! ACADEMIC REFERENCE:

      
! ===============================================================================
SUBROUTINE INIT (basenm, nx, ny, nz, mnx5, mny5, mnz5, mnx6, mny6, mnz6, eqs, iars, &
                 xzn, xznl,xznr,yznr, delx, yzn, yznl, zzn,zznl, physic, gamma_ad, gridlx, gridly, gridlz, CT, &
                 me, mp, m_0, kb, c, uwi, vwi, shift, Rorb, xCPU, yCPU, zCPU,ic,jc,kc,ncl,timej,pc)
! ===============================================================================

        USE math_constants
        USE parallel
        USE tipos
        USE workarrays
	USE physic_constants

        implicit none

# include "types.h"

        !------------------------------------------------------- Input Variables
        
        character(4):: basenm

        integer(is):: nx, ny, nz, mnx5, mny5, mnz5, eqs, iars, mnx6, mny6, mnz6, xCPU, yCPU, zCPU

        real(rs)   :: xzn(*), xznl(*), yzn(*), yznl(*), zzn(*), zznl(*), xznr(*), yznr(*), delx, gridlx, gridly, gridlz


	real(rs) :: me, mp, kb, m_0, c, timej, pc

        !------------------------------------------------------ Output Variables
                
        type(Tphysic):: physic (-4:mnx5,-4:mny5,-4:mnz5)
	type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)

	real(rs) :: gamma_ad, vwi, uwi, shift, Rorb

        !------------------------------------------------------- Local Variables
      
        integer(is):: i, j, k, tid, ICO, JCO, KCO, NBOXX, NBOXY, NBOXZ, IREF, JREF, KREF, ierr

        real(rs) :: b40, b41, b42, b43, b44, VA, RI, RA, R, OMEGA, tt1
  
        real(rs), pointer :: rho_s(:), p_s(:), eps_s(:), dpdrh_s(:),   &
                             dpdeps_s(:), h_s(:), cs_s(:), xe_s(:),    &
			     tem_s(:), gam_s(:)

        real(rs), pointer :: scrch1(:), scrch2(:)

	real(rs) :: theta1, theta2, beta, sigmap, sigmat, dist, distc, ac, cc, rhocl, ran1, dive

	integer(is) :: ncl, ic(2000), jc(2000), kc(2000)

	integer(is) :: ii, jj, iic, jjc, kkc, ci, cj, ck
   
	real(rs) :: delxh, delyh, delzh, xznCL, yznCL, zznCL


        real(rs) :: RHOB, VELB, PB, RBM, BMJ, WB, BPHI1, B1, B2, B3, X1, X2, X3


	!real(rs) :: AY(-4:NX+5,-4:0,-4:NZ+5)

        real(rs), allocatable :: AX(:,:,:),AY(:,:,:), BX(:,:,:), BY(:,:,:), BZ(:,:,:),      &
                                 XZNLLOC(:), XZNRLOC(:), XZNLOC(:), ZZNLLOC(:), ZZNRLOC(:), &
                                 ZZNLOC(:),  YZNLLOC(:), YZNRLOC(:), YZNLOC(:), AZ(:,:,:)


	real(rs), allocatable:: scr(:,:,:)

	integer(is) :: ini, ifi, jni, jfi, kni, kfi, kk

        character(256):: fich


        real(rs) :: RB, VB, PERIOD, PHI, H(3), R3(3), N(3), X(3), XN, ABSX, ABSN, ALPHA

        real(rs) :: UX, UY, UZ, BX_ANALYTIC, BY_ANALYTIC, BZ_ANALYTIC, BMZ,XROT,YROT,ZROT

        integer(is) :: fznpx_L, nznpx_L, fznpy_L, nznpy_L, fznpz_L, nznpz_L, &
                       fznpx_G, nznpx_G, fznpy_G, nznpy_G, fznpz_G, nznpz_G


	integer(is) :: idum
	!COMMON /random/ idum

        !------------------------------------------------------------------------- BEGIN

        allocate (AX(1:300+1,1:300+1,1:300+1), stat=ierr)
        allocate (AY(1:300+1,1:300+1,1:300+1), stat=ierr)
        allocate (AZ(1:300+1,1:300+1,1:300+1), stat=ierr)

        allocate (BX(1:300+1,1:300+1,1:300+1), stat=ierr)
        allocate (BY(1:300+1,1:300+1,1:300+1), stat=ierr)
        allocate (BZ(1:300+1,1:300+1,1:300+1), stat=ierr)

        allocate (XZNLLOC(1:300+1), stat=ierr)
        allocate (XZNRLOC(1:300+1), stat=ierr)
        allocate (XZNLOC(1:300+1), stat=ierr)

        allocate (YZNLLOC(1:300+1), stat=ierr)
        allocate (YZNRLOC(1:300+1), stat=ierr)
        allocate (YZNLOC(1:300+1), stat=ierr)

        allocate (ZZNLLOC(1:300+1), stat=ierr)
        allocate (ZZNRLOC(1:300+1), stat=ierr)
        allocate (ZZNLOC(1:300+1), stat=ierr)


        tid = 1 
   
        !  1-vectors
        rho_s     => wa(tid)%wa1init100
        eps_s     => wa(tid)%wa1init101
        p_s       => wa(tid)%wa1init102
        dpdrh_s   => wa(tid)%wa1init103
        dpdeps_s  => wa(tid)%wa1init104
        h_s       => wa(tid)%wa1init105
        cs_s      => wa(tid)%wa1init106

	xe_s      => wa(tid)%wa1init107
        tem_s     => wa(tid)%wa1init108
	gam_s     => wa(tid)%wa1init109
      
        !-4:mn6 vectors
        scrch1 => wa(tid)%wa6init100
        scrch2 => wa(tid)%wa6init101

        !   -------------
        !   ALL THE VELOCITIES, MAGNETIC FIELDS AND TRACER
        !   INITIALIZED TO ZERO BY DEFAULT
        !   -------------

        DO k=-4,nz+5
           DO j=-4,ny+5
              DO i=-4,nx+5

                 physic(i,j,k)%velx   = 0.D0
                 physic(i,j,k)%vely   = 0.D0
                 physic(i,j,k)%velz   = 0.D0

                 physic(i,j,k)%bxc    = 0.D0
                 physic(i,j,k)%byc    = 0.D0
                 physic(i,j,k)%bzc    = 0.D0

                 CT(i,j,k)%bx    = 0.D0
                 CT(i,j,k)%by    = 0.D0
                 CT(i,j,k)%bz    = 0.D0

                 physic(i,j,k)%tracer = 0.D0


              ENDDO
           ENDDO
        ENDDO

        

        !   -------------
        !   BALSARA'S TEST #1
        !   (BALSARA, 2001, ApJ, vol. 132, 83)
        !   -------------

        IF (basenm .EQ. 'BA1_') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'BA1_ works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 2.D0
         
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

		    IF (iars .EQ. 2) THEN

                       IF (xzn(i) .LE. 0.D0) THEN
                       
                          physic(i,j,k)%densty  = 1.D0
                          physic(i,j,k)%denstye = 1.D0
                          physic(i,j,k)%pres    = 1.D0
                          physic(i,j,k)%tracer  = 1.D0
                          physic(i,j,k)%bxc     = 0.5D0
                          physic(i,j,k)%byc     = 1.D0
                          physic(i,j,k)%bzc     = 0.D0
                       
                       ELSE
                       
                          physic(i,j,k)%densty  = 0.125D0
                          physic(i,j,k)%denstye = 0.125D0
                          physic(i,j,k)%pres    = 0.1D0
                          physic(i,j,k)%tracer  = 0.D0
                          physic(i,j,k)%bxc     = 0.5D0
                          physic(i,j,k)%byc     = -1.D0
                          physic(i,j,k)%bzc     = 0.D0
                       
                       ENDIF

		    ELSE IF (iars .GE. 3) THEN

                       IF (xzn(i) .LE. xznl(1)+gridlx/2.d0) THEN

                          physic(i,j,k)%densty  = 1.D0
                          physic(i,j,k)%denstye = 1.D0
                          physic(i,j,k)%pres    = 1.D0
                          physic(i,j,k)%tracer  = 1.D0
                          physic(i,j,k)%bxc     = 0.5D0
                          physic(i,j,k)%byc     = 1.D0
                          physic(i,j,k)%bzc     = 0.D0

                       ELSE IF (xzn(i) .LE. xznl(1)+gridlx/2.d0 + delx) THEN

                          physic(i,j,k)%densty  = 0.5625D0
                          physic(i,j,k)%denstye = 0.5626D0
                          physic(i,j,k)%pres    = 0.55D0
                          physic(i,j,k)%tracer  = 0.5D0
                          physic(i,j,k)%bxc     = 0.5D0
                          physic(i,j,k)%byc     = 0.D0
                          physic(i,j,k)%bzc     = 0.D0

                       ELSE

                          physic(i,j,k)%densty  = 0.125D0
                          physic(i,j,k)%denstye = 0.125D0
                          physic(i,j,k)%pres    = 0.1D0
                          physic(i,j,k)%tracer  = 0.D0
                          physic(i,j,k)%bxc     = 0.5D0
                          physic(i,j,k)%byc     = -1.D0
                          physic(i,j,k)%bzc     = 0.D0

                       ENDIF

		    ENDIF	

                    CT(i,j,k)%bx = physic(i,j,k)%bxc
                   
                    physic(i,j,k)%velx = 0.D0
                    physic(i,j,k)%vely = 0.D0
                    physic(i,j,k)%velz = 0.D0
                   
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)
	            
                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

		    ! Hydromagnetic specific enthalpy * density                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conserved variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                                              
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF

        
        !   -------------
        !   BALSARA'S TEST #2
        !   (BALSARA, 2001, ApJ, vol. 132, 83)
        !   -------------

        IF (basenm .EQ. 'BA2_') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'BA2_ works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0
         
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    IF (xzn(i) .LE. 0.D0) THEN
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 3.D1
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%bxc     = 5.D0
                       physic(i,j,k)%byc     = 6.D0
                       physic(i,j,k)%bzc     = 6.D0
                       
                       
                    ELSE
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%bxc     = 5.D0
                       physic(i,j,k)%byc     = 0.7D0
                       physic(i,j,k)%bzc     = 0.7D0
                       
                    ENDIF

                    CT(i,j,k)%bx = physic(i,j,k)%bxc

                    physic(i,j,k)%velx = 0.D0
                    physic(i,j,k)%vely = 0.D0
                    physic(i,j,k)%velz = 0.D0
                   
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43
                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF


        !   -------------
        !   BALSARA'S TEST #3
        !   (BALSARA, 2001, ApJ, vol. 132, 83)
        !   -------------

        IF (basenm .EQ. 'BA3_') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'BA3_ works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0

                   
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    IF (xzn(i) .LE. 0) THEN
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 1.D3
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%bxc     = 10.D0
                       physic(i,j,k)%byc     = 7.D0
                       physic(i,j,k)%bzc     = 7.D0
                       
                    ELSE
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 0.1D0
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%bxc     = 10.D0
                       physic(i,j,k)%byc     = 0.7D0
                       physic(i,j,k)%bzc     = 0.7D0
                       
                    END IF

                    CT(i,j,k)%bx = physic(i,j,k)%bxc

                    physic(i,j,k)%velx = 0.D0
                    physic(i,j,k)%vely = 0.D0
                    physic(i,j,k)%velz = 0.D0
                   
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43
                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF


        !   -------------
        !   BALSARA'S TEST #3 + MULTICOMPONENT GAS
        !   (BALSARA, 2001, ApJ, vol. 132, 83)
        !   -------------

        IF (basenm .EQ. 'BA3S') THEN

           !IF (eqs .NE. 2) THEN
           !   call ERROR ('INIT', 'BA3S works with multicomponent gas', __FILE__, __LINE__)
           !ENDIF

           gamma_ad = 5.D0/3.D0

                   
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    IF (xzn(i) .LE. 0) THEN
                       
                       physic(i,j,k)%densty  = 5.2D4
                       physic(i,j,k)%denstye = physic(i,j,k)%densty*me/(me+mp)
                       physic(i,j,k)%pres    = 5.5D-5
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%bxc     = 0.D0
                       physic(i,j,k)%byc     = 0.D0
                       physic(i,j,k)%bzc     = 0.D0

                       physic(i,j,k)%velx = 0.3D0
                       physic(i,j,k)%vely = 0.D0
                       physic(i,j,k)%velz = 0.D0

                       
                    ELSE
                       
                       physic(i,j,k)%densty  = 5.2D4
                       physic(i,j,k)%denstye = physic(i,j,k)%densty*me/(me+mp)
		       !physic(i,j,k)%denstye = 1.D0	
                       physic(i,j,k)%pres    = 5.5D-5
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%bxc     = 0.D0
                       physic(i,j,k)%byc     = 0.D0
                       physic(i,j,k)%bzc     = 0.D0
                       
                       physic(i,j,k)%velx = 0.2D0
                       physic(i,j,k)%vely = 0.D0
                       physic(i,j,k)%velz = 0.D0


                    END IF

                    CT(i,j,k)%bx = physic(i,j,k)%bxc

                                       
                    rho_s(1) = physic(i,j,k)%densty
		    xe_s(1)  = physic(i,j,k)%denstye/physic(i,j,k)%densty


                    tem_s(1)     = physic(i,j,k)%pres / (physic(i,j,k)%densty* &
                                  (xe_s(1)/me + (1.D0-xe_s(1))/mp))



                    call EOS_SYNGE (1, tid, 1, 1, me, mp, rho_s, xe_s, tem_s, eps_s,   &
                                    p_s, gam_s, h_s, cs_s, dpdrh_s, dpdeps_s)


                    physic(i,j,k)%eps = eps_s(1)
		    physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43
                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO
        ENDIF



        !   -------------
        !   BALSARA'S TEST #4
        !   (BALSARA, 2001, ApJ, vol. 132, 83)
        !   -------------

        IF (basenm .EQ. 'BA4_') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'BA4_ works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0
         
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


		    IF (iars .EQ. 2) THEN	

                       IF (xzn(i) .LE. 0) THEN
                       
                          physic(i,j,k)%densty  = 1.D0
                          physic(i,j,k)%denstye = 1.D0
                          physic(i,j,k)%pres    = 0.1D0
                          physic(i,j,k)%tracer  = 1.D0
                          physic(i,j,k)%velx    = 0.999D0
                          physic(i,j,k)%bxc     = 10.D0
                          physic(i,j,k)%byc     = 7.D0
                          physic(i,j,k)%bzc     = 7.D0
                       
                       ELSE
                       
                          physic(i,j,k)%densty  = 1.D0
                          physic(i,j,k)%denstye = 1.D0
                          physic(i,j,k)%pres    = 0.1D0
                          physic(i,j,k)%tracer  = 0.D0
                          physic(i,j,k)%velx    = -0.999D0
                          physic(i,j,k)%bxc     = 10.D0
                          physic(i,j,k)%byc     = -7.D0
                          physic(i,j,k)%bzc     = -7.D0
                       
                       END IF

		    ELSE IF (iars .GE. 3) THEN

                       IF (xzn(i) .LE. xzn(1) + gridlx/2.D0) THEN                       

                          physic(i,j,k)%densty  = 1.D0
                          physic(i,j,k)%denstye = 1.D0
                          physic(i,j,k)%pres    = 0.1D0
                          physic(i,j,k)%tracer  = 1.D0
                          physic(i,j,k)%velx    = 0.999D0
                          physic(i,j,k)%bxc     = 10.D0
                          physic(i,j,k)%byc     = 7.D0
                          physic(i,j,k)%bzc     = 7.D0

		       ELSE IF (xzn(i) .LE. xzn(1) + gridlx/2.0 + delx) THEN

                          physic(i,j,k)%densty  = 1.D0
                          physic(i,j,k)%denstye = 1.D0
                          physic(i,j,k)%pres    = 0.1D0
                          physic(i,j,k)%tracer  = 0.5D0
                          physic(i,j,k)%velx    = 0.D0
                          physic(i,j,k)%bxc     = 10.D0
                          physic(i,j,k)%byc     = 0.D0
                          physic(i,j,k)%bzc     = 0.D0

                       ELSE

                          physic(i,j,k)%densty  = 1.D0
                          physic(i,j,k)%denstye = 1.D0
                          physic(i,j,k)%pres    = 0.1D0
                          physic(i,j,k)%tracer  = 0.D0
                          physic(i,j,k)%velx    = -0.999D0
                          physic(i,j,k)%bxc     = 10.D0
                          physic(i,j,k)%byc     = -7.D0
                          physic(i,j,k)%bzc     = -7.D0

                       ENDIF

		    ENDIF

                    CT(i,j,k)%bx = physic(i,j,k)%bxc
                  
                    physic(i,j,k)%vely = 0.D0
                    physic(i,j,k)%velz = 0.D0
                   
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)
		
			
                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               
	            ! Assuming the Minkowski metric to be (-,+,+,+)
           	    b44 = -b40*b40+b41*b41+b42*b42+b43*b43                     
     
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF

        
        !   -------------
        !   BALSARA'S TEST #5
        !   (BALSARA, 2001, ApJ, vol. 132, 83)
        !   -------------

        IF (basenm .EQ. 'BA5_') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'BA5_ works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0
         
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    IF (xzn(i) .LE. 0) THEN
                       
                       physic(i,j,k)%densty  = 1.08D0
                       physic(i,j,k)%denstye = 1.08D0
                       physic(i,j,k)%pres    = 0.95D0
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%velx    = 0.4D0
                       physic(i,j,k)%vely    = 0.3D0
                       physic(i,j,k)%velz    = 0.2D0
                       physic(i,j,k)%bxc     = 2.D0
                       physic(i,j,k)%byc     = 0.3D0
                       physic(i,j,k)%bzc     = 0.3D0
                       
                    ELSE
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%velx    = -0.45D0
                       physic(i,j,k)%vely    = -0.2D0
                       physic(i,j,k)%velz    = 0.2D0
                       physic(i,j,k)%bxc     = 2.D0
                       physic(i,j,k)%byc     = -0.7D0
                       physic(i,j,k)%bzc     = 0.5D0
                       
                    END IF

                    CT(i,j,k)%bx = physic(i,j,k)%bxc
                   
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)
                    
                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               
                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43
                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF



        !   -------------
        !   RIEMANN PROBLEM WITH SINUSOIDAL PERTURBATION
        !   (DEL ZANNA & BUCCIANTINI, 2002, A&A, vol. 390, 3)
        !   -------------

        IF (basenm .EQ. 'DZBU') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'DZBU works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0

           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    IF (xzn(i) .LE. 0.5D0) THEN

                       physic(i,j,k)%densty  = 5.D0
                       physic(i,j,k)%denstye = 5.D0
                       physic(i,j,k)%pres    = 50.D0
                       physic(i,j,k)%tracer  = 1.D0

                    ELSE

                       physic(i,j,k)%densty  = 2.D0+0.3D0*DSIN(50.D0*xzn(i))
                       physic(i,j,k)%denstye = 2.D0+0.3D0*DSIN(50.D0*xzn(i))
                       physic(i,j,k)%pres    = 5.D0
                       physic(i,j,k)%tracer  = 0.D0

                    END IF


                    CT(i,j,k)%bx = physic(i,j,k)%bxc

                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)

                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40


                 ENDDO
              ENDDO
           ENDDO



           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF



        !   -------------
        !   ROTATIONAL WAVE
        !   MIGNONE ET AL., 2007, MNRAS, vol. 000, 1-15)
        !   -------------

        
        IF (basenm .EQ. 'CONT') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'CONT works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0
         
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    IF (xzn(i) .LE. 0) THEN
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 1.D0
		       physic(i,j,k)%velx    = 0.4D0
		       physic(i,j,k)%vely    = -0.3D0
		       physic(i,j,k)%velz    = 0.5D0
		       physic(i,j,k)%bxc     = 2.4D0
		       physic(i,j,k)%byc     = 1.D0
		       physic(i,j,k)%bzc     = -1.6D0
                       
                    ELSE
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%velx    = 0.377347D0
                       physic(i,j,k)%vely    = -0.482389D0
                       physic(i,j,k)%velz    = 0.424190D0
                       physic(i,j,k)%bxc     = 2.4D0
                       physic(i,j,k)%byc     = -0.1D0
                       physic(i,j,k)%bzc     = -2.178213
                       
                    END IF

                    CT(i,j,k)%bx = physic(i,j,k)%bxc

                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   =  1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               
                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43
                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF

        
        !   -------------
        !   KOMISSAROV'S SHOCK TUBE #2
        !   (KOMISSAROV, 1999, MNRAS, vol. 303, 343;
        !   INITIAL DATA CORRECTED AS IN KOMISSAROV, 2002, astro-ph/0209213v1)
        !   -------------

        IF (basenm .EQ. 'KO2_') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'KO2_ works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0
         
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    IF (xzn(i) .LE. 0) THEN
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 3.D1
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%bxc     = 0.D0
                       physic(i,j,k)%byc     = 2.D1
                       physic(i,j,k)%bzc     = 0.D0
                       
                    ELSE
                       
                       physic(i,j,k)%densty  = 0.1D0
                       physic(i,j,k)%denstye = 0.1D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%bxc     = 0.D0
                       physic(i,j,k)%byc     = 0.D0
                       physic(i,j,k)%bzc     = 0.D0
                       
                    END IF

                    physic(i,j,k)%velx = 0.D0
                    physic(i,j,k)%vely = 0.D0
                    physic(i,j,k)%velz = 0.D0

                    CT(i,j,k)%bx = physic(i,j,k)%bxc
                   
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               
                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43
                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF

        
        !   -------------
        !   GENERIC ALFVEN TEST
        !   (GIACOMAZZO & REZZOLLA, 2006, JFM, vol. 562, 223)
        !   -------------

        IF (basenm .EQ. 'GR2_') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'GR2_ works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0
         
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    IF (xzn(i) .LE. 0) THEN
                       
                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 5.D0
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%velx    = 0.D0
                       physic(i,j,k)%vely    = 0.3D0
                       physic(i,j,k)%velz    = 0.4D0
                       physic(i,j,k)%bxc     = 1.D0
                       physic(i,j,k)%byc     = 6.D0
                       physic(i,j,k)%bzc     = 2.D0
                       
                    ELSE
                       
                       physic(i,j,k)%densty  = 0.9D0
                       physic(i,j,k)%denstye = 0.9D0
                       physic(i,j,k)%pres    = 5.3D0
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%velx    = 0.D0
                       physic(i,j,k)%vely    = 0.D0
                       physic(i,j,k)%velz    = 0.D0
                       physic(i,j,k)%bxc     = 1.D0
                       physic(i,j,k)%byc     = 5.D0
                       physic(i,j,k)%bzc     = 2.D0
                       
                    END IF

                    CT(i,j,k)%bx = physic(i,j,k)%bxc
                   
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)

                   ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               
                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43
                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF

        
        !     -------------
        !     LARGE AMPLITUDE CP ALFVEN WAVE TEST
        !     (DEL ZANNA ET AL. 2007, A&A, vol. 473, 11)
        !     -------------
        
        IF (basenm .EQ. 'DCAW') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'DCAW works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0

           !  For this test we normalize our physical quantities 
           !  assuming: rho=p=Bo=eta=1.

           ! Since eps=3 --> h=5:
           VA    = DSQRT(2.D0/(7.D0 + DSQRT(45.D0)))
         
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    physic(i,j,k)%densty  = 1.D0
                    physic(i,j,k)%denstye = 1.D0
                    physic(i,j,k)%pres    = 1.D0
                    physic(i,j,k)%tracer  = 1.D0
                    physic(i,j,k)%bxc     = 1.D0
                    physic(i,j,k)%byc     = DCOS(xzn(i))
                    physic(i,j,k)%bzc     = DSIN(xzn(i))
                    physic(i,j,k)%velx    = 0.D0
                    physic(i,j,k)%vely    = -VA*physic(i,j,k)%byc
                    physic(i,j,k)%velz    = -VA*physic(i,j,k)%bzc

                    CT(i,j,k)%bx = physic(i,j,k)%bxc
                                         
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres
                  
                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    
                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)
                    
                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40
               
                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43
                          
                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                          
                          
                    
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

        ENDIF




        !   -------------
        !   RELATIVISTIC ORSZAG-TANG VORTEX
        !   (MARQUINA)
        !   -------------

        IF (basenm .EQ. 'ORTA') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'IMPL works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0

           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


                    physic(i,j,k)%densty = 1.D0
                    physic(i,j,k)%denstye= 1.D0
                    physic(i,j,k)%pres   = 1.D0

                    physic(i,j,k)%velx   = -0.75D0/DSQRT(2.D0)*DSIN(YZN(J)) 
                    physic(i,j,k)%vely   =  0.75D0/DSQRT(2.D0)*DSIN(XZN(I))                        
                    physic(i,j,k)%velz   =  0.D0                       


                    physic(i,j,k)%bxc    = -1.D0*DSIN(YZN(J))      
                    physic(i,j,k)%byc    =  1.D0*DSIN(2.D0*XZN(I))    
                    physic(i,j,k)%bzc    =  0.D0


                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)


                    physic(i,j,k)%eps = eps_s(1)
                    !physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40



                 ENDDO
              ENDDO
           ENDDO


           !	  ------------------
           !	  STAGGERED MAGNETIC FIELDS
           !	  (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1
           !	  AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !	  CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE
           !	  OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !	  ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO


           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO


        ENDIF



        
        !   -------------
        !   RELATIVISTIC SHOCK-CLOUD INTERACTION
        !   (MIGNONE & BODO, 2006, MNRAS, vol. 368, 1040-1054)
        !   -------------

        IF (basenm .EQ. 'IMPL') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'IMPL works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0
           
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


		    IF (xzn(i) .GT. 0.6D0) THEN

			physic(i,j,k)%densty = 1.D0 
			physic(i,j,k)%denstye= 1.D0
			physic(i,j,k)%pres   = 1.D-3
			physic(i,j,k)%velx   = -DSQRT(1.D0-1.D0/100.D0)
			physic(i,j,k)%bzc    = 0.5D0

		    ELSE

                        physic(i,j,k)%densty = 42.5942D0
                        physic(i,j,k)%denstye= 42.5942D0
                        physic(i,j,k)%pres   = 127.9483D0
                        physic(i,j,k)%bzc    = -2.12971D0

		    ENDIF

		    ! High-Density Clump

		    RA = 0.15D0

                    ! Radius
                    R = DSQRT((xzn(i)-0.8D0)**2 + yzn(j)**2)

		    IF (R .LE. RA) THEN

		        physic(i,j,k)%densty = 10.D0
			physic(i,j,k)%denstye = 10.D0
	
		    ENDIF

                 
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    

                    physic(i,j,k)%eps = eps_s(1)
		    !physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44      

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                    
                 ENDDO
              ENDDO
           ENDDO


           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1 
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE 
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

       
           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO

           
        ENDIF



        !   -------------
        !   KOMISSAROV'S CYLINDRICAL MAGNETIZED BLAST WAVE
        !   (KOMISSAROV, 1999, MNRAS, vol. 303, 343)
        !   -------------

        IF (basenm .EQ. 'CMBW') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'CMBW works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0

           ! Radius' limiters.

           RI = 0.8D0
           RA = 1.0D0
                    
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    ! Radius
                    R = DSQRT(xzn(i)**2 + yzn(j)**2)
		    
                    IF (R .LE. RI) THEN
                       
                       physic(i,j,k)%densty  = 1.D-2
                       physic(i,j,k)%denstye = 1.D-2
                       !physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 1.D0
		       physic(i,j,k)%eps = 300.D0
                       
                    ELSE IF (R .GT. RI .AND. R .LE. RA) THEN
                       
                       physic(i,j,k)%densty  = 1.D-4 + (1.D-2  - 1.D-4)*(RA-R)/(RA-RI)
                       physic(i,j,k)%denstye = 1.D-4 + (1.D-2  - 1.D-4)*(RA-R)/(RA-RI)
                       physic(i,j,k)%eps    = 150.D0 + (300.D0 - 150.D0)*(RA-R)/(RA-RI)
                       physic(i,j,k)%tracer  = 0.5D0
                       
                    ELSE

                       physic(i,j,k)%densty  = 1.D-4
                       physic(i,j,k)%denstye = 1.D-4
                       !physic(i,j,k)%pres    = 3.D-5
                       physic(i,j,k)%tracer  = 0.D0
		       physic(i,j,k)%eps = 150.D0

                    ENDIF
	         
                    physic(i,j,k)%velx    = 0.D0
                    physic(i,j,k)%vely    = 0.D0
                    physic(i,j,k)%velz    = 0.D0

                    ! Relatively strong magnetic field.
                    !physic(i,j,k)%bxc     = 1.D0

                    !Relatevely weak magnetic field.
                    !physic(i,j,k)%byc     = 1.d-2

                    !Very weak magnetic field.
                    physic(i,j,k)%bxc     = 1.d0                    
                    !physic(i,j,k)%byc     = 1.D-8
                                          
                    rho_s(1) = physic(i,j,k)%densty
                    !p_s(1)   = physic(i,j,k)%pres
		    eps_s(1)  = physic(i,j,k)%eps                  

                    call EOS (2, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    

                    !physic(i,j,k)%eps = eps_s(1)
		    physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44      

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                    
                 ENDDO
              ENDDO
           ENDDO


           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1 
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE 
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

       
           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO

           
        ENDIF


        !   -------------
        !   KOMISSAROV'S CYLINDRICAL MAGNETIZED BLAST WAVE
        !   WITH A GENERAL EQUATION OF STATE (ELECTRON-PROTON PLASMA).
        !   -------------

        IF (basenm .EQ. 'SYBW') THEN


           gamma_ad = 5.D0/3.D0
           RB = 0.008D0*DTAN(5.D0*PI/180) ! PARSEC
           VB = 0.26D0

           Period = 162.D0/2.D0*86400.D0/(RB*pc/c)

           PHI = PI*timej/Period + PI/4.D0
           !PHI = PI/4.D0

           !H = (/0.D0,0.0085D0/RB,0.D0/)
           !R3 = (/H(2)*DTAN(0.D0*PI/180)*DCOS(PHI),0.D0,H(2)*DTAN(20.D0*PI/180)*DSIN(PHI)/)
           N =  (/0.008D0/RB*DCOS(PHI), 0.008D0/RB*DSIN(PHI), 0.D0/)
           ABSN = DSQRT(N(1)**2+N(2)**2)

                    
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


                       X = (/XZN(I),YZN(J), 0.D0/)

                    ABSX = DSQRT(X(1)**2+X(2)**2)

                    XN = (X(1)*N(1)+X(2)*N(2))/ABSN

                    !ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                    !    ux = -DSIN(PHI)
                    !    uz =  DCOS(PHI)

                    ! XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+&
                    !        (-UZ*DSIN(20.D0*PI/180))*YZN(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                    ! YROT = (UZ*DSIN(20.D0*PI/180))*XZN(I)+(DCOS(20.D0*PI/180))*YZN(J)+(-UX*DSIN(20.D0*PI/180))*ZZN(K)

                    ! ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+(UX*DSIN(20.D0*PI/180))*YZN(J)+&
                    !        (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                    !    R = DSQRT(XROT**2+ZROT**2)
			 R = DABS(-XZN(I)*DSIN(PHI)+YZN(J)*DCOS(PHI)) 
			XN = DCOS(PHI)*XZN(I)+DSIN(PHI)*YZN(J)


                    IF (R .LE. 0.008D0/RB*DTAN(6.D0*PI/180) .AND. XN .LE. 0.008D0/RB) THEN


                        physic(i,j,k)%velx = VB*N(1)/ABSN
                        physic(i,j,k)%vely = VB*N(2)/ABSN
                        

                        !WB = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)
			WB = 1.D0/DSQRT(1.D0-0.26D0**2)

                        physic(i,j,k)%densty = 3.12D2*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)
                        physic(i,j,k)%denstye= 3.12D2*me/(me+mp)*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)

                        IF (physic(i,j,k)%densty .LT. 1.D0) physic(i,j,k)%densty = 1.D0
                        IF (physic(i,j,k)%denstye .LT. 1.D0*me/(me+mp)) physic(i,j,k)%denstye = 1.D0*me/(me+mp)

                        physic(i,j,k)%tracer = 1.D0*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)

                        RBM = 0.37D0
                        BMJ = 9.17D-4
                         PB = 2.66D-7

                        BPHI1 = (2.D0*BMJ*(1.D0/RBM))/(1.D0+(1.D0/RBM)**2)

                        physic(i,j,k)%pres   = 2.D0*(BMJ/(WB*(1.D0+(R/RBM)**2)))**2+PB-BPHI1**2/(2.D0*WB**2)*(1.D0+RBM**2)

                        rho_s(1) = physic(i,j,k)%densty
                        xe_s(1)  = physic(i,j,k)%denstye/physic(i,j,k)%densty

                        tem_s(1)     = physic(i,j,k)%pres / (physic(i,j,k)%densty* &
                                       (xe_s(1)/me + (1.D0-xe_s(1))/mp))

                        call EOS_SYNGE (1, tid, 1, 1, me, mp, rho_s, xe_s, tem_s, eps_s,   &
                                        p_s, gam_s, h_s, cs_s, dpdrh_s, dpdeps_s)


                        physic(i,j,k)%pres = p_s(1)
                        physic(i,j,k)%eps  = eps_s(1)


                        ! Lorentz Factor W
                        scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

  
                        ! Magnetic field four-vector.
                        b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                          physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                          physic(i,j,k)%velz*physic(i,j,k)%bzc)

                        b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                        b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                        b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                        ! Assuming the Minkowski metric to be (-,+,+,+)
                        b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                        scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                        ! Conservative variables.
                        physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                        physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                        physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                        physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                        physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                        physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                                physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40


                    ELSE


                            physic(i,j,k)%velx = 0.D0
                            physic(i,j,k)%vely = 0.D0
                            physic(i,j,k)%velz = 0.D0


                            physic(i,j,k)%densty  = 1.D0
                            physic(i,j,k)%denstye = 1.D0*me/(me+mp)
                            physic(i,j,k)%tracer  = 0.D0
                            physic(i,j,k)%pres    = 2.66D-7


                            rho_s(1) = physic(i,j,k)%densty
                            xe_s(1)  = physic(i,j,k)%denstye/physic(i,j,k)%densty

                            tem_s(1)     = physic(i,j,k)%pres / (physic(i,j,k)%densty* &
                                          (xe_s(1)/me + (1.D0-xe_s(1))/mp))

                            call EOS_SYNGE (1, tid, 1, 1, me, mp, rho_s, xe_s, tem_s, eps_s,   &
                                            p_s, gam_s, h_s, cs_s, dpdrh_s, dpdeps_s)


                            physic(i,j,k)%pres = p_s(1)
                            physic(i,j,k)%eps  = eps_s(1)


                            ! Lorentz Factor W
                            scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                            ! Magnetic field four-vector.
                            b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                              physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                              physic(i,j,k)%velz*physic(i,j,k)%bzc)

                            b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                            b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                            b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                            ! Assuming the Minkowski metric to be (-,+,+,+)
                            b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                            scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44


                            ! Conservative variables.
                            physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                            physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                            physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                            physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                            physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                            physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                                    physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40


                    ENDIF



                 ENDDO
              ENDDO
           ENDDO
		

           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1 
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE 
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

       
           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO

        ENDIF


        !     -------------
        !     RELATIVISTIC ROTOR
        !     (DEL ZANNA ET AL. 2003, A&A, vol. 400, 397)
        !     -------------
        
        IF (basenm .EQ. 'RROT') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'RROT works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0

           RI = 0.1D0
           OMEGA = 9.95D0
                    
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    ! Radius
                    R = DSQRT(xzn(i)**2 + yzn(j)**2)

                    IF (R .LE. RI) THEN
                       
                       physic(i,j,k)%densty  = 10.D0
                       physic(i,j,k)%denstye = 10.D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%velx    = -yzn(j)*OMEGA
                       physic(i,j,k)%vely    =  xzn(i)*OMEGA
                       physic(i,j,k)%velz    = 0.D0
                       
                    ELSE

                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%velx    = 0.D0
                       physic(i,j,k)%vely    = 0.D0
                       physic(i,j,k)%velz    = 0.D0

                    ENDIF

                    ! Background magnetic field.
                    physic(i,j,k)%bxc     = 1.D0                   
                    physic(i,j,k)%byc     = 0.D0
                    physic(i,j,k)%bzc     = 0.D0

                                          
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)

                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO

           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO


           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO


        ENDIF



        !     -------------
        !     SLAB JET IN CARTESIAN COORDINATES
        !     (S.S. KOMISSAROV 1999, MNRAS, vol. 303, 343-366)
        !     -------------

        IF (basenm .EQ. 'SLKO') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'SLKO works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0


           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


                    physic(i,j,k)%densty  = 10.D0
                    physic(i,j,k)%denstye = 10.D0
                    physic(i,j,k)%pres    = 0.01D0
                    physic(i,j,k)%tracer  = 0.D0
                    physic(i,j,k)%velx    = 0.D0
                    physic(i,j,k)%vely    = 0.D0
                    physic(i,j,k)%velz    = 0.D0


                    ! Background magnetic field.
                    physic(i,j,k)%bxc     = 1.D0
                    physic(i,j,k)%byc     = 0.D0
                    physic(i,j,k)%bzc     = 0.D0


                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)

                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO


           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO


           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO


        ENDIF



        !   -------------
        !   RELATIVISTIC 2D KELVIN-HELMOTZ INSTABILITY
        !   (BECKWITH&STONE, 2011, APJ, vol. 193:6)
        !   -------------

        IF (basenm .EQ. 'KH2D') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'KH2D works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0

           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


                    IF (yzn(j) .GT. 0.D0) THEN

                        physic(i,j,k)%densty = 0.505D0+0.495D0*DTANH((YZN(J)-0.5D0)/0.01D0)
                        physic(i,j,k)%denstye= 0.505D0+0.495D0*DTANH((YZN(J)-0.5D0)/0.01D0)
                        physic(i,j,k)%pres   = 1.D0
                        physic(i,j,k)%velx   = 0.5D0*DTANH((YZN(J)-0.5D0)/0.01D0)
                        physic(i,j,k)%vely   = 0.1D0*0.5D0*DSIN(2.D0*PI*XZN(I))*EXP(-(YZN(J)-0.5D0)**2/0.1D0)                        

                    ELSE

                        physic(i,j,k)%densty = 0.505D0-0.495D0*DTANH((YZN(J)+0.5D0)/0.01D0)
                        physic(i,j,k)%denstye= 0.505D0-0.495D0*DTANH((YZN(J)+0.5D0)/0.01D0)
                        physic(i,j,k)%pres   = 1.D0
                        physic(i,j,k)%velx   = -0.5D0*DTANH((YZN(J)+0.5D0)/0.01D0)
                        physic(i,j,k)%vely   = -0.1D0*0.5D0*DSIN(2.D0*PI*XZN(I))*EXP(-(YZN(J)+0.5D0)**2/0.1D0)                       

                    ENDIF

                    physic(i,j,k)%bxc = 1.D-3

		    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)


                    physic(i,j,k)%eps = eps_s(1)
                    !physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO


           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO


           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO


        ENDIF



        !   -------------
        !   KOMISSAROV'S SPHERICAL  MAGNETIZED BLAST WAVE
        !   (KOMISSAROV, 1999, MNRAS, vol. 303, 343)
        !   -------------

        IF (basenm .EQ. 'SMBW') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'SMBW works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0

           ! Radius' limiters.

           RI = 0.8D0
           RA = 1.0D0

           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    ! Radius
                    R = DSQRT(xzn(i)**2 + yzn(j)**2 + zzn(k)**2)

                    IF (R .LE. RI) THEN

                       physic(i,j,k)%densty  = 1.D-2
                       physic(i,j,k)%denstye = 1.D-2
                       !physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%eps = 300.D0

                    ELSE IF (R .GT. RI .AND. R .LE. RA) THEN

                       physic(i,j,k)%densty  = 1.D-4 + (1.D-2  - 1.D-4)*(RA-R)/(RA-RI)
                       physic(i,j,k)%denstye = 1.D-4 + (1.D-2  - 1.D-4)*(RA-R)/(RA-RI)
                       physic(i,j,k)%eps    = 150.D0 + (300.D0 - 150.D0)*(RA-R)/(RA-RI)
                       physic(i,j,k)%tracer  = 0.5D0


                    ELSE

                       physic(i,j,k)%densty  = 1.D-4
                       physic(i,j,k)%denstye = 1.D-4
                       !physic(i,j,k)%pres    = 3.D-5
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%eps = 150.D0

                    ENDIF

                    physic(i,j,k)%velx    = 0.D0
                    physic(i,j,k)%vely    = 0.D0
                    physic(i,j,k)%velz    = 0.D0

                    ! Relatively strong magnetic field.
                    !physic(i,j,k)%bxc     = 1.D0

                    !Relatevely weak magnetic field.
                    !physic(i,j,k)%bxc     = 0.D0

                    !Very weak magnetic field.
                    !physic(i,j,k)%bxc     = 0.01D0

                    !physic(i,j,k)%bxc     = 0.D0
                    

                    physic(i,j,k)%bxc     = 1.D0*DSIN(PI/2.D0)*DCOS(PI/4.D0)
                    physic(i,j,k)%byc     = 1.D0*DSIN(PI/2.D0)*DSIN(PI/4.D0)
                    physic(i,j,k)%bzc     = 1.D0*DCOS(PI/2.D0)


                    rho_s(1) = physic(i,j,k)%densty
                    !p_s(1)   = physic(i,j,k)%pres
                    eps_s(1)  = physic(i,j,k)%eps

                    call EOS (2, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)

                    !physic(i,j,k)%eps = eps_s(1)
                    physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO


           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=2, nz
              DO j=1, ny
                 DO i=1, nx
                    CT(i,j,k)%bz = (physic(i,j,k-1)%bzc + physic(i,j,k)%bzc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO


           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO


           ! z-boundary
           DO j=1, ny
              DO i=1, nx

                 CT(i,j,1)%bz = physic(i,j,1)%bzc
                 CT(i,j,nz+1)%bz = physic(i,j,nz)%bzc

              ENDDO
           ENDDO


	ENDIF


        !     -------------
        !     3D RELATIVISTIC ROTOR
        !     (MIGNONE ET AL. 2009, MNRAS, 393, 1141–1156)
        !     -------------


        IF (basenm .EQ. '3ROT') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'RROT works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0

           RI = 0.1D0
           OMEGA = 9.95D0
                    
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx

                    ! Radius
                    R = DSQRT(xzn(i)**2 + yzn(j)**2+zzn(k)**2)

                    IF (R .LE. RI) THEN
                       
                       physic(i,j,k)%densty  = 10.D0
                       physic(i,j,k)%denstye = 10.D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 1.D0
                       physic(i,j,k)%velx    = -yzn(j)*OMEGA
                       physic(i,j,k)%vely    =  xzn(i)*OMEGA
                       physic(i,j,k)%velz    = 0.D0
                       
                    ELSE

                       physic(i,j,k)%densty  = 1.D0
                       physic(i,j,k)%denstye = 1.D0
                       physic(i,j,k)%pres    = 1.D0
                       physic(i,j,k)%tracer  = 0.D0
                       physic(i,j,k)%velx    = 0.D0
                       physic(i,j,k)%vely    = 0.D0
                       physic(i,j,k)%velz    = 0.D0

                    ENDIF

                    ! Background magnetic field.
                    physic(i,j,k)%bxc     = 1.D0                   
                    physic(i,j,k)%byc     = 0.D0
                    physic(i,j,k)%bzc     = 0.D0

                                          
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)

                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO

           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=2, nz
              DO j=1, ny
                 DO i=1, nx
                    CT(i,j,k)%bz = (physic(i,j,k-1)%bzc + physic(i,j,k)%bzc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO


           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO

           ! z-boundary
           DO j=1, ny
              DO i=1, nx

                 CT(i,j,1)%bz = physic(i,j,1)%bzc
                 CT(i,j,nz+1)%bz = physic(i,j,nz)%bzc

              ENDDO
           ENDDO


        ENDIF


        !   -------------
        !   RELATIVISTIC 3D SHOCK-CLOUD INTERACTION
        !   (MIGNONE ET AL., 2012, APJ, vol. 198:7)
        !   -------------

        IF (basenm .EQ. 'SC3D') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'SC3D works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0
           
           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


		    physic(i,j,k)%tracer = 0.D0

                    IF (xzn(i) .GT. 0.6D0) THEN

                        physic(i,j,k)%densty = 1.D0
                        physic(i,j,k)%denstye= 1.D0
                        physic(i,j,k)%pres   = 1.D-3
                        physic(i,j,k)%velx   = -DSQRT(1.D0-1.D0/100.D0)
                        physic(i,j,k)%bzc    = 0.5D0

                    ELSE

                        physic(i,j,k)%densty = 42.5942D0
                        physic(i,j,k)%denstye= 42.5942D0
                        physic(i,j,k)%pres   = 127.9483D0
                        physic(i,j,k)%bzc    = -2.12971D0

                    ENDIF

		    physic(i,j,k)%bzc = physic(i,j,k)%bzc/DSQRT(2.D0)
		    physic(i,j,k)%byc = physic(i,j,k)%bzc


		    ! High-Density Clump

		    RI = 0.15D0

                    ! Radius
                    R = DSQRT((xzn(i)-0.8D0)**2 + yzn(j)**2 + zzn(k)**2)


                    IF (R .LE. RI) THEN

                       physic(i,j,k)%densty  = 1.D1
                       physic(i,j,k)%denstye = 1.D1
		       physic(i,j,k)%tracer  = 1.D0

                    ENDIF


                 
                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    

                    physic(i,j,k)%eps = eps_s(1)
		    !physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44      

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40
                    
                 ENDDO
              ENDDO
           ENDDO


           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1 
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE 
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=2, nz
              DO j=1, ny
                 DO i=1, nx
                    CT(i,j,k)%bz = (physic(i,j,k-1)%bzc + physic(i,j,k)%bzc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO

       
           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO


           ! z-boundary
           DO j=1, ny
              DO i=1, nx

                 CT(i,j,1)%bz = physic(i,j,1)%bzc
                 CT(i,j,nz+1)%bz = physic(i,j,nz)%bzc

              ENDDO
           ENDDO
           
        ENDIF


        !   -------------
        !   RELATIVISTIC 3D KELVIN-HELMOTZ INSTABILITY
        !   (BECKWITH&STONE, 2011, APJ, vol. 193:6)
        !   -------------

        IF (basenm .EQ. 'KH3D') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'KH3D works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 4.D0/3.D0

           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx



                    IF (yzn(j) .GT. 0.D0) THEN

                        physic(i,j,k)%densty = 0.505D0+0.495D0*DTANH((YZN(J)-0.5D0)/0.01D0)
                        physic(i,j,k)%denstye= 0.505D0+0.495D0*DTANH((YZN(J)-0.5D0)/0.01D0)
                        physic(i,j,k)%pres   = 1.D0
                        physic(i,j,k)%velx   = 0.5D0*DTANH((YZN(J)-0.5D0)/0.01D0)
                        physic(i,j,k)%vely   = 0.1D0*0.5D0*DSIN(2.D0*PI*XZN(I))*EXP(-(YZN(J)-0.5D0)**2/0.1D0)

                    ELSE

                        physic(i,j,k)%densty = 0.505D0-0.495D0*DTANH((YZN(J)+0.5D0)/0.01D0)
                        physic(i,j,k)%denstye= 0.505D0-0.495D0*DTANH((YZN(J)+0.5D0)/0.01D0)
                        physic(i,j,k)%pres   = 1.D0
                        physic(i,j,k)%velx   = -0.5D0*DTANH((YZN(J)+0.5D0)/0.01D0)
                        physic(i,j,k)%vely   = -0.1D0*0.5D0*DSIN(2.D0*PI*XZN(I))*EXP(-(YZN(J)+0.5D0)**2/0.1D0)

                    ENDIF


                    physic(i,j,k)%bxc = 1.D-3

		    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)


                    physic(i,j,k)%eps = eps_s(1)
                    !physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO


           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=2, nz
              DO j=1, ny
                 DO i=1, nx
                    CT(i,j,k)%bz = (physic(i,j,k-1)%bzc + physic(i,j,k)%bzc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO


           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO


           ! z-boundary
           DO j=1, ny
              DO i=1, nx

                 CT(i,j,1)%bz = physic(i,j,1)%bzc
                 CT(i,j,nz+1)%bz = physic(i,j,nz)%bzc

              ENDDO
           ENDDO

        ENDIF



        !   -----------------------
        !   ADVECTION TOROIDAL FIELD
        !   -----------------------

        IF (basenm .EQ. 'LOOP') THEN


           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'KH3D works with ideal gas', __FILE__, __LINE__)
           ENDIF


           gamma_ad = 5.D0/3.D0
           !vwi      = 6.67D-3
           !uwi      = 9.05D-10
           !shift    = 10.D0
           !Rorb     = 500.D0


               RHOB = 10.D0
               VELB = 0.5D0
                PB  = 0.5D0
               RBM  = 0.37D0
                BMJ = 1.D0
                WB  = 1.D0/DSQRT(1.D0-VELB**2)


           DO k = 1, nz+1
              DO j = 1, ny
                 DO i = 1, nx+1

                     R = DSQRT(XZNL(I)**2+ZZNL(K)**2)

                     IF (R .LE. 1.D0) THEN

                         AY(I,J,K) = -0.37D0*LOG(0.1369D0+R**2)+0.37D0*LOG(0.1369D0+1.D0**2)

                     ELSE

                         AY(I,J,K) = 0.D0

                     ENDIF


                 ENDDO
              ENDDO
           ENDDO


           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx+1

                    CT(i,j,k)%bx = (AY(I,J,K+1) - AY(I,J,K))/DELX

                  ENDDO
                ENDDO
           ENDDO



           DO k = 1, nz+1
              DO j = 1, ny
                 DO i = 1, nx

                    CT(i,j,k)%bz = -(AY(I+1,J,K) - AY(I,J,K))/DELX

                  ENDDO
                ENDDO
           ENDDO



           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


                     physic(i,j,k)%bxc = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                     physic(i,j,k)%bzc = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)/2.D0


                  ENDDO
                ENDDO
           ENDDO



        DO k=1,nz-1
           DO j=-4,0
              DO i=1, nx-1


                DIVE = (physic(i+1,j,k)%bxc-physic(i,j,k)%bxc)/delx+&
                      (physic(i,j+1,k)%byc-physic(i,j,k)%byc)/delx+&
                      (physic(i,j,k+1)%bzc-physic(i,j,k)%bzc)/delx
                IF (DIVE .GT. 1.D-6) print*, 'divergence', I, j, K
             ENDDO
           ENDDO
        ENDDO




           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx




                    tt1 = DSQRT(ZZN(K)*ZZN(K) + XZN(I)*XZN(I))

                    IF (tt1 .LE. 1.D0) THEN

                        	BPHI1 = (2.D0*BMJ*(1.D0/RBM))/(1.D0+(1.D0/RBM)**2)

                               physic(i,j,k)%densty  = 10.D0
                               physic(i,j,k)%denstye = 10.D0
                               physic(i,j,k)%pres    = 2.D0*(BMJ/(WB*(1.D0+(TT1/RBM)**2)))**2+PB-BPHI1**2/(2.D0*WB**2)*(1.D0+RBM**2)
                               physic(i,j,k)%tracer    = 1.D0

                    physic(i,j,k)%velx = 0.D0
                    physic(i,j,k)%vely = 0.5D0
                    physic(i,j,k)%velz = 0.D0


                    ELSE


                               physic(i,j,k)%densty  = 1.D0
                               physic(i,j,k)%denstye = 1.D0
                               physic(i,j,k)%pres    = pb
                               physic(i,j,k)%tracer    = 0.D0

                    physic(i,j,k)%velx = 0.D0
                    physic(i,j,k)%vely = 0.D0
                    physic(i,j,k)%velz = 0.D0

		    ENDIF



		    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)


                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO



        ENDIF



       
        !   -------------
        !   RELATIVISTIC MHD JET SIMULATION
        !   (MIGNONE ET AL., 2010, MNRAS, 402,7-12)
        !   -------------

        IF (basenm .EQ. 'MI10') THEN

           IF (eqs .NE. 1) THEN
              call ERROR ('INIT', 'MI10 works with ideal gas', __FILE__, __LINE__)
           ENDIF

           gamma_ad = 5.D0/3.D0

           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


                    physic(i,j,k)%densty = 1.D0
                    physic(i,j,k)%denstye= 1.D0
                    physic(i,j,k)%pres   = 0.0000161428D0

                    physic(i,j,k)%velx   = 0.D0
                    physic(i,j,k)%vely   = 0.D0
		    physic(i,j,k)%velz   = 0.D0


                    physic(i,j,k)%bxc   = 0.D0                              
                    physic(i,j,k)%byc   = 0.D0                                                          
                    physic(i,j,k)%bzc   = 0.D0

                    rho_s(1) = physic(i,j,k)%densty
                    p_s(1)   = physic(i,j,k)%pres

                    call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)


                    physic(i,j,k)%eps = eps_s(1)
                    !physic(i,j,k)%pres = p_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

                 ENDDO
              ENDDO
           ENDDO


           !      ------------------
           !      STAGGERED MAGNETIC FIELDS
           !      (BXS AT INTERFACES I=1 AND NX+1, AND BYS AT INTERFACES J=1
           !      AND NY+1 ARE FIXED AS INITIAL CONDITIONS BUT THEY MUST BE FIXED
           !      CONSISTENTLY WITH THE BOUNDARY CONDITIONS (AS, E.G., IN THE CASE
           !      OF PERIODIC B.C., WHERE BXS(1,J,K)=BXS(NX+1,J,K), BYS(I,1,K)=BYS(I,NY+1,K))
           !      ------------------


           DO k=1, nz
              DO j=1, ny
                 DO i=2, nx
                    CT(i,j,k)%bx = (physic(i-1,j,k)%bxc + physic(i,j,k)%bxc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           DO k=1, nz
              DO j=2, ny
                 DO i=1, nx
                    CT(i,j,k)%by = (physic(i,j-1,k)%byc + physic(i,j,k)%byc)/2.D0
                 ENDDO
              ENDDO
           ENDDO


           DO k=2, nz
              DO j=1, ny
                 DO i=1, nx
                    CT(i,j,k)%bz = (physic(i,j,k-1)%bzc + physic(i,j,k)%bzc)/2.D0
                 ENDDO
              ENDDO
           ENDDO

           ! x-boundary
           DO k=1,nz
              DO j=1,ny

                 CT(1,j,k)%bx = physic(1,j,k)%bxc
                 CT(nx+1,j,k)%bx = physic(nx,j,k)%bxc

              ENDDO
           ENDDO


           ! y-boundary
           DO k=1, nz
              DO i=1, nx

                 CT(i,1,k)%by = physic(i,1,k)%byc
                 CT(i,ny+1,k)%by = physic(i,ny,k)%byc

              ENDDO
           ENDDO



           ! z-boundary
           DO j=1, ny
              DO i=1, nx

                 CT(i,j,1)%bz = physic(i,j,1)%bzc
                 CT(i,j,nz+1)%bz = physic(i,j,nz)%bzc

              ENDDO
           ENDDO

        ENDIF



	! ==================================================
        ! MICROQUASAR JET SIMULATIONS WITH HOMOGENEOUS WIND.
        !    The units of the code are chosen so that
        !       - RB     (beam radius)
        !       - C      (light speed)
        !       - RHOEXT (density of the external medium)
        !    are equal to 1.
        !
        !    This means that the Length Unit is RB,
        !                        Time   Unit is RB/C,
        !                        Mass   Unit is RHOEXT*RB**3.
        ! ==================================================


        IF (BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'H') THEN


           allocate (scr(120,120,120), stat= ierr)

           ini=mpi_coords(1)*nx
           ifi=ini+nx
           jni=mpi_coords(2)*ny
           jfi=jni+ny
           kni=mpi_coords(3)*nz
           kfi=kni+nz

           fich= '/storage/scratch/lv72/lv72764/RMHD-3D/MH00/dnsty'

           open( 1, file=fich, form='unformatted', action='read', status='old', access='stream', iostat=ierr )
           if (ierr/=0) then
               call ERROR( 'initbu: Error openning ', fich, __FILE__, __LINE__ )
           endif

           read(1) (((scr(i,j,k), i=1,120), j=1,120), k=1,120)

           close(1)


           do kk=1,nz
              do jj=1,ny
                 do ii=1,nx
                     physic(ii,jj,kk)%densty = scr(ii+ini,jj+jni,kk+kni)
                     physic(ii,jj,kk)%denstye = physic(ii,jj,kk)%densty
                 enddo
              enddo
           enddo

           deallocate(scr)



	   ! JET RADIUS :: RB = 6 x 10^9 cm

           gamma_ad = 5.D0/3.D0
           vwi      = 6.67D-3
           uwi      = 9.05D-10
           shift    = 120.D0
           Rorb     = 500.D0


           DO k= 1, nz
              DO j= 1, ny
                 DO i= 1, nx
			

	 	    
	            ! INCLINATION
                    theta1 = DATAN(zzn(k)/       &
                             DSQRT((yzn(j)+shift)**2+(Rorb-xzn(i))**2))

	            ! AZIMUTHAL ANGLE.
                    theta2 = DATAN((yzn(j)+shift)/(Rorb-xzn(i)))
	
	
          	    physic(i,j,k)%velx = -vwi*DCOS(theta1)*DCOS(theta2)
          	    physic(i,j,k)%vely =  vwi*DCOS(theta1)*DSIN(theta2)
          	    physic(i,j,k)%velz =  vwi*DSIN(theta1)


                    physic(i,j,k)%tracer = 0.D0

                    !physic(i,j,k)%eps    = uwi

                    !physic(i,j,k)%densty=1.D0
                    !physic(i,j,k)%denstye=1.D0


                    physic(i,j,k)%pres = (gamma_ad-1.D0)*uwi*1.D0

                    rho_s(1) = physic(i,j,k)%densty
                    !eps_s(1) = physic(i,j,k)%eps
		    p_s(1)   = physic(i,j,k)%pres

                    IF (eqs .EQ. 1) THEN
                       call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
                    ENDIF

                    !physic(i,j,k)%pres = p_s(1)
                    physic(i,j,k)%eps    = eps_s(1)


                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40


	 	 ENDDO
	      ENDDO
	   ENDDO




	  IF (mpi_coords(2) .EQ. 0) THEN


           !DO k = 1, nz+1
           !   DO j = -4, 0
           !      DO i = 1, nx+1

           !          R = DSQRT(XZNL(I)**2+ZZNL(K)**2)

           !          IF (R .LE. 2.D0) THEN

           !              AY(I,J,K) = -0.0395206D0*LOG(0.1369D0+R**2)/DCOSH(0.5D0*(R+0.1D0)**2)!+&
           !                           !0.0395206D0*LOG(0.1369D0+2.D0**2)/DCOSH(0.5D0*(2.D0+0.1D0)**2)
	   !	     ELSE

		         !AY(i,j,k) = 0.D0


	!	     ENDIF


         !        ENDDO
         !     ENDDO
         !  ENDDO




	   DELX = GRIDLZ/(NZ*ZCPU)

	   ZZNLLOC(1) = -10.D0

           DO  k=2,NZ*zCPU+1
              ZZNLLOC(K) = ZZNLLOC(K-1) + DELX
           ENDDO

           DO  k=1,NZ*zCPU+1
              ZZNRLOC(k) = ZZNLLOC(K) + DELX
           ENDDO

           DO  K=1,NZ*zCPU
              ZZNLOC(K)  = 0.5D0*(ZZNLLOC(K) + ZZNRLOC(K))
           ENDDO


           DELX	= GRIDLX/(NX*XCPU)


           XZNLLOC(1) = -10.D0

           DO  I=2,NX*XCPU+1
              XZNLLOC(I) = XZNLLOC(I-1) + DELX
           ENDDO

           DO  I=1,NX*XCPU+1
              XZNRLOC(I) = XZNLLOC(I) + DELX
           ENDDO

           DO  I=1,NX*XCPU
              XZNLOC(I)  = 0.5D0*(XZNLLOC(I) + XZNRLOC(I))
           ENDDO



	   DO k = 1, nz*zCPU
              DO j = -4, 0
                 DO i = 1, nx*xCPU+1

                     R = DSQRT(XZNLLOC(I)**2+ZZNLOC(K)**2)

		     IF (R .LE. 2.D0) THEN

                        BX(i,j,k) = -(0.577364D0*ZZNLOC(K)/DCOSH(R**4))/(1.D0+7.3046D0*R**2)

		     ELSE

                        BX(i,j,k) = 0.D0 

		     ENDIF

                 ENDDO
              ENDDO
           ENDDO



              DO j = -4, 0
                 DO i = 1, nx*xCPU

                     BZ(I,J,1) = 0.D0

                  ENDDO
                ENDDO



           DO k = 1, nz*zCPU
              DO j = -4, 0
                 DO i = 1, nx*xCPU

                     R = DSQRT(XZNLOC(I)**2+ZZNLLOC(K)**2)

		     BZ(i,j,k+1) =  BZ(i,j,k)+BX(i,j,k)-BX(i+1,j,k)


                  ENDDO
                ENDDO
           ENDDO



           DO k = 1, nz
              DO j = -4, 0
                 DO i = 1, nx+1

                    !R = DSQRT(XZNL(I)**2+ZZN(K)**2)


                    !CT(i,j,k)%bx = (AY(I,J,K+1) - AY(I,J,K))/DELX  !/DCOSH(R**4)
                    ! CT(i,j,k)%bx = -(0.577364D0*ZZN(K)/DCOSH((-0.05D0+R)**2))/(1.D0+7.3046D0*R**2)
 		    CT(i,j,k)%bx = BX(mpi_coords(1)*NX+i,j,mpi_coords(3)*NZ+k)



                  ENDDO
                ENDDO
           ENDDO


              !DO j = -4, 0
              !   DO i = 1, nx

                     !R = DSQRT(XZN(I)**2+ZZNL(K)**2)


                    !CT(i,j,k)%bz = -(AY(I+1,J,K) - AY(I,J,K))/DELX   !/DCOSH(R**4)
                    !CT(i,j,1)%bz = 0.D0




                  !ENDDO
                !ENDDO



           DO k = 1, nz+1
              DO j = -4, 0
                 DO i = 1, nx

                     !R = DSQRT(XZN(I)**2+ZZNL(K)**2)


                    !CT(i,j,k)%bz = -(AY(I+1,J,K) - AY(I,J,K))/DELX   !/DCOSH(R**4)
                     !CT(i,j,k+1)%bz = CT(i,j,k)%bz+CT(i,j,k)%bx-CT(i+1,j,k)%bx

 	            CT(i,j,k)%bz = BZ(mpi_coords(1)*NX+i,j,mpi_coords(3)*NZ+k)



                  ENDDO
                ENDDO
           ENDDO



           DO k = 1, nz
              DO j = -4, 0
                 DO i = 1, nx


                     physic(i,j,k)%bxc = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                     physic(i,j,k)%bzc = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)/2.D0


                  ENDDO
                ENDDO
           ENDDO



	  ENDIF



        DO k=1,nz-1
           DO j=-4,0
              DO i=1, nx-1


                DIVE = (CT(i+1,j,k)%bx-CT(i,j,k)%bx)/delx+&
                      (CT(i,j+1,k)%by-CT(i,j,k)%by)/delx+&
                      (CT(i,j,k+1)%bz-CT(i,j,k)%bz)/delx
                IF (DIVE .GT. 1.D-6) print*, 'divergence', I, j, K

             ENDDO
           ENDDO
        ENDDO


	ENDIF        



	! ==================================================
        ! MICROQUASAR JET SIMULATIONS WITH CLUMPY STELLAR WIND.
        !    The units of the code are chosen so that
        !       - RB     (beam radius)
        !       - C      (light speed)
        !       - RHOEXT (mean density of the external medium)
        !    are equal to 1.
        !
        !    This means that the Length Unit is RB,
        !                        Time   Unit is RB/C,
        !                        Mass   Unit is RHOEXT*RB**3.
        ! ==================================================


        IF (BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'C') THEN


           ! GENERATE THE POSITIONS OF THE CLOUDS.

           !idum=-3
           !ncl=240
           
           !IF (ncl .GT. 2000) THEN
           !   call ERROR ('Main', 'You tend to have too many clouds', __FILE__, __LINE__)
           !ENDIF
           
           !DO jj = 1, ncl

           !   ic(jj) = ran1(idum)*nx*xCPU
           !   jc(jj) = ran1(idum)*ny*yCPU
           !   kc(jj) = ran1(idum)*nz*zCPU

              
           !   IF (jj .GE. 2) THEN                 
           !      DO ii = 1, jj-1
                    

	   !         dist = sqrt( (real(ic(jj)-ic(ii)))**2 + (real(jc(jj)-jc(ii)))**2 +   &
           !                       real((kc(jj)-kc(ii)))**2)

           !         IF (dist .LT. 55) THEN

           !             ic(jj) = ran1(idum)*nx*xCPU
           !   		jc(jj) = ran1(idum)*ny*yCPU
           !   		kc(jj) = ran1(idum)*nz*zCPU

           !         ENDIF


           !      ENDDO                 
           !   ENDIF



                    ! =========================
                    ! MOVE BOUNDARY CLUMPS
                    ! =========================


            !        IF (mpi_coords(1) .EQ. 0 .AND. ic(jj) .LE. 40) ic(jj) = 40
            !        IF (mpi_coords(2) .EQ. 0 .AND. jc(jj) .LE. 40) jc(jj) = 40
            !        IF (mpi_coords(3) .EQ. 0 .AND. kc(jj) .LE. 40) kc(jj) = 40


            !        IF (mpi_coords(1) .EQ. (mpi_dims(1)-1)  .AND. ic(jj) .GE. (nx*xCPU - 40)) ic(jj) = nx*xCPU - 40
            !        IF (mpi_coords(2) .EQ. (mpi_dims(2)-1)  .AND. jc(jj) .GE. (ny*yCPU - 40)) jc(jj) = ny*yCPU - 40
            !        IF (mpi_coords(3) .EQ. (mpi_dims(3)-1)  .AND. kc(jj) .GE. (nz*zCPU - 40)) kc(jj) = nz*zCPU - 40




           !ENDDO 

           ac = 20.D0
           cc = 2.D0

           ! Intialization of density.

           DO k=1,nz
              DO j=1,ny
                 DO i=1,nx
                    physic(i,j,k)%densty = 0.D0
		    physic(i,j,k)%denstye = 0.D0
                 ENDDO
              ENDDO
           ENDDO


           DO jj = 1, ncl
              
              ci=ic(jj)
              cj=jc(jj)
              ck=kc(jj)


	      delxh = gridlx /(nx*xCPU)
              delyh = gridly /(ny*yCPU)
              delzh = gridlz /(nz*zCPU)

	      xznCL = ic(jj)*delxh + 0.5D0*delxh - gridlx/2.D0
              yznCL = jc(jj)*delyh + 0.5D0*delyh 
              zznCL = kc(jj)*delzh + 0.5D0*delzh - gridlz/2.D0
 
	 

              DO kkc = 1, nz
                 DO jjc = 1, ny
                    DO iic = 1, nx

                       distc = DSQRT ((xzn(iic)-xznCL)**2+(yzn(jjc)-yznCL)**2 + &
                                      (zzn(kkc)-zznCL)**2)
                              
                       IF (distc .NE. 0.D0) THEN
                          rhocl = ac/exp(distc**2/(2.D0*cc**2))
                       ELSE 
                          rhocl = ac
                       ENDIF
                              
                       !IF (rhocl .LE. 0.5D0) rhocl = 0.D0

                       physic(iic,jjc,kkc)%densty = physic(iic,jjc,kkc)%densty + rhocl
                       physic(iic,jjc,kkc)%denstye = physic(iic,jjc,kkc)%denstye + rhocl
 
                    ENDDO
                 ENDDO
              ENDDO

           ENDDO

           
	   ! JET RADIUS :: RB = 6 x 10^9 cm

           gamma_ad = 5.D0/3.D0
           vwi      = 6.67D-3
           uwi      = 9.05D-10
           shift    = 10.D0
           Rorb     = 500.D0

           
           DO k= 1, nz
              DO j= 1, ny
                 DO i= 1, nx
			
	 	    
	            ! INCLINATION
                    theta1 = DATAN(zzn(k)/       &
                             DSQRT((yzn(j)+shift)**2+(Rorb-xzn(i))**2))

	            ! AZIMUTHAL ANGLE.
                    theta2 = DATAN((yzn(j)+shift)/(Rorb-xzn(i)))
	
	
          	    !physic(i,j,k)%velx = -vwi*DCOS(theta1)*DCOS(theta2)
          	    !physic(i,j,k)%vely =  vwi*DCOS(theta1)*DSIN(theta2)
          	    !physic(i,j,k)%velz =  vwi*DSIN(theta1)


                    physic(i,j,k)%tracer = 0.D0

                    
                    IF (physic(i,j,k)%densty .LT. 0.1D0) THEN
                       physic(i,j,k)%densty  = 0.1D0
                       physic(i,j,k)%denstye = 0.1D0
                    ENDIF
                    

                    physic(i,j,k)%pres = (gamma_ad-1.D0)*uwi


                    rho_s(1) = physic(i,j,k)%densty
		    p_s(1) = physic(i,j,k)%pres

		    IF (eqs .EQ. 1) THEN
                       call EOS (1, 1, 1, gamma_ad, eps_s, rho_s, p_s, dpdrh_s, dpdeps_s, h_s, cs_s)
		    ENDIF

                    physic(i,j,k)%eps = eps_s(1)

                    ! Lorentz Factor W
                    scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                    ! Magnetic field four-vector.
                    b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                      physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                      physic(i,j,k)%velz*physic(i,j,k)%bzc)

                    b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                    b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                    b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                    ! Conservative variables.
                    physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                    physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                    physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                    physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                    physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                    physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                            physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40


	 	 ENDDO
	      ENDDO
	   ENDDO



          IF (mpi_coords(2) .EQ. 0) THEN


           DO k = 1, nz+1
              DO j = -4, 0
                 DO i = 1, nx+1

                     R = DSQRT(XZNL(I)**2+ZZNL(K)**2)

                     IF (R .LE. 1.5D0) THEN

                         AY(I,J,K) = -0.0395206D0*LOG(0.1369D0+R**2)+0.0395206D0*LOG(0.1369D0+1.5D0**2)

                     ELSE

                         AY(I,J,K) = 0.D0

                     ENDIF


                 ENDDO
              ENDDO
           ENDDO


          DO k = 1, nz
              DO j = -4, 0
                 DO i = 1, nx+1

                    R = DSQRT(XZNL(I)**2+ZZN(K)**2)


                    CT(i,j,k)%bx = (AY(I,J,K+1) - AY(I,J,K))/DELX/DCOSH(R**4)



                  ENDDO
                ENDDO
           ENDDO



           DO k = 1, nz+1
              DO j = -4, 0
                 DO i = 1, nx

                     R = DSQRT(XZN(I)**2+ZZNL(K)**2)


                    CT(i,j,k)%bz = -(AY(I+1,J,K) - AY(I,J,K))/DELX/DCOSH(R**4)

                  ENDDO
                ENDDO
           ENDDO



           DO k = 1, nz
              DO j = -4, 0
                 DO i = 1, nx


                     physic(i,j,k)%bxc = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                     physic(i,j,k)%bzc = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)/2.D0


                  ENDDO
                ENDDO
           ENDDO



          ENDIF


	ENDIF        


        ! ==================================================
        ! SIMULATIONS OF SS433 SUBPARSEC JETS
        !    The units of the code are chosen so that
        !	- RB     (beam radius)
        !	- C	 (light speed)
        !	- RHOEXT (density of the external medium)
        !    are equal to 1.
        !
	!    This means that the Length Unit is RB,
        !                        Time   Unit is RB/C,
        !                        Mass   Unit is RHOEXT*RB**3.
        ! ==================================================


        IF (BASENM(1:2) .EQ. 'S4') THEN


           gamma_ad = 5.D0/3.D0
           RB = 0.008D0*DTAN(5.D0*PI/180) ! PARSEC
           VB = 0.26D0

	   !idum = -3

           Period = 162.D0*86400.D0/(RB*pc/c)

           PHI = 2.D0*PI*timej/Period + PI/4.D0 
           !PHI = 0.D0

           H = (/0.D0,0.0085D0/RB,0.D0/)
           R3 = (/H(2)*DTAN(20.D0*PI/180)*DCOS(PHI),0.D0,H(2)*DTAN(20.D0*PI/180)*DSIN(PHI)/)
           N =  H + R3 !AXIS OF PRECESSION CONE
           ABSN = DSQRT(N(1)**2+N(2)**2+N(3)**2)


           DO k= 1, nz
              DO j= 1, ny
                 DO i= 1, nx


                       X = (/XZN(I),YZN(J),ZZN(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)

                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)

                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZN(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZN(I)+(DCOS(20.D0*PI/180))*YZN(J)+(-UX*DSIN(20.D0*PI/180))*ZZN(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+(UX*DSIN(20.D0*PI/180))*YZN(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                        R = DSQRT(XROT**2+ZROT**2)


                    IF (R .LE. 0.008D0/RB*DTAN(6.D0*PI/180) .AND. XN .LE. 0.008D0/RB) THEN

			IF (BASENM .EQ. 'S430') THEN

                            physic(i,j,k)%velx = VB*N(1)/ABSN
                            physic(i,j,k)%vely = VB*N(2)/ABSN
                            physic(i,j,k)%velz = VB*N(3)/ABSN

			ELSE IF (BASENM .EQ. 'S433') THEN

                            physic(i,j,k)%velx = 0.D0
                            physic(i,j,k)%vely = 0.D0
                            physic(i,j,k)%velz = 0.D0

			ENDIF


                        !WB = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)
			!WB = 1.D0/DSQRT(1.D0-0.26D0**2)

                        !physic(i,j,k)%densty = 3.12D2*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)
                        !physic(i,j,k)%denstye= 3.12D2*me/(me+mp)*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)

                        !IF (physic(i,j,k)%densty .LT. 1.D0) physic(i,j,k)%densty = 1.D0
                        !IF (physic(i,j,k)%denstye .LT. 1.D0*me/(me+mp)) physic(i,j,k)%denstye = 1.D0*me/(me+mp)

                        !physic(i,j,k)%tracer = 1.D0*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)

                        !RBM = 0.37D0
                        !BMJ = 9.17D-4
                        ! PB = 2.66D-7

                        !BPHI1 = (2.D0*BMJ*(1.D0/RBM))/(1.D0+(1.D0/RBM)**2)

                        !physic(i,j,k)%pres   = 2.D0*(BMJ/(WB*(1.D0+(R/RBM)**2)))**2+PB-BPHI1**2/(2.D0*WB**2)*(1.D0+RBM**2)


                        physic(i,j,k)%densty  = 1.D0
                        physic(i,j,k)%denstye = 1.D0*me/(me+mp)
                        physic(i,j,k)%tracer  = 0.D0
                        physic(i,j,k)%pres    = 2.66D-7


                        rho_s(1) = physic(i,j,k)%densty
                        xe_s(1)  = physic(i,j,k)%denstye/physic(i,j,k)%densty

                        tem_s(1)     = physic(i,j,k)%pres / (physic(i,j,k)%densty* &
                                       (xe_s(1)/me + (1.D0-xe_s(1))/mp))

                        call EOS_SYNGE (1, tid, 1, 1, me, mp, rho_s, xe_s, tem_s, eps_s,   &
                                        p_s, gam_s, h_s, cs_s, dpdrh_s, dpdeps_s)


                        physic(i,j,k)%pres = p_s(1)
                        physic(i,j,k)%eps  = eps_s(1)


                        ! Lorentz Factor W
                        scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)

  
                        ! Magnetic field four-vector.
                        b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                          physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                          physic(i,j,k)%velz*physic(i,j,k)%bzc)

                        b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                        b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                        b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                        ! Assuming the Minkowski metric to be (-,+,+,+)
                        b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                        scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44

                        ! Conservative variables.
                        physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                        physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                        physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                        physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                        physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                        physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                                physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40


                    ELSE


                            physic(i,j,k)%velx = 0.D0
                            physic(i,j,k)%vely = 0.D0
                            physic(i,j,k)%velz = 0.D0


                            physic(i,j,k)%densty  = 1.D0
                            physic(i,j,k)%denstye = 1.D0*me/(me+mp)
                            physic(i,j,k)%tracer  = 0.D0
                            physic(i,j,k)%pres    = 2.66D-7


                            rho_s(1) = physic(i,j,k)%densty
                            xe_s(1)  = physic(i,j,k)%denstye/physic(i,j,k)%densty

                            tem_s(1)     = physic(i,j,k)%pres / (physic(i,j,k)%densty* &
                                          (xe_s(1)/me + (1.D0-xe_s(1))/mp))

                            call EOS_SYNGE (1, tid, 1, 1, me, mp, rho_s, xe_s, tem_s, eps_s,   &
                                            p_s, gam_s, h_s, cs_s, dpdrh_s, dpdeps_s)


                            physic(i,j,k)%pres = p_s(1)
                            physic(i,j,k)%eps  = eps_s(1)


                            ! Lorentz Factor W
                            scrch1(i)   = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)


                            ! Magnetic field four-vector.
                            b40 = scrch1(i)*( physic(i,j,k)%velx*physic(i,j,k)%bxc +   &
                                              physic(i,j,k)%vely*physic(i,j,k)%byc +   &
                                              physic(i,j,k)%velz*physic(i,j,k)%bzc)

                            b41 = physic(i,j,k)%bxc/scrch1(i) + physic(i,j,k)%velx*b40
                            b42 = physic(i,j,k)%byc/scrch1(i) + physic(i,j,k)%vely*b40
                            b43 = physic(i,j,k)%bzc/scrch1(i) + physic(i,j,k)%velz*b40

                            ! Assuming the Minkowski metric to be (-,+,+,+)
                            b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                            scrch2(i)  = physic(i,j,k)%densty*h_s(1) + b44


                            ! Conservative variables.
                            physic(i,j,k)%denstw  = physic(i,j,k)%densty*scrch1(i)
                            physic(i,j,k)%denstew = physic(i,j,k)%denstye*scrch1(i)
                            physic(i,j,k)%momenx  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velx - b40*b41
                            physic(i,j,k)%momeny  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%vely - b40*b42
                            physic(i,j,k)%momenz  = scrch2(i)*scrch1(i)*scrch1(i)*physic(i,j,k)%velz - b40*b43
                            physic(i,j,k)%energy  = scrch2(i)*scrch1(i)*scrch1(i)-physic(i,j,k)%pres -           &
                                                    physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40


                    ENDIF


                 ENDDO
              ENDDO
           ENDDO

           DELX = 0.1D0

           ZZNLLOC(1) = -15.D0

           DO  k=2,300+1
               ZZNLLOC(K) = ZZNLLOC(K-1) + DELX
           ENDDO

           DO  k=1,300
              ZZNRLOC(k) = ZZNLLOC(K) + DELX
           ENDDO

           DO  K=1,300
              ZZNLOC(K)  = 0.5D0*(ZZNLLOC(K) + ZZNRLOC(K))
           ENDDO


           DELX = 0.1D0


           XZNLLOC(1) = -15.D0

           DO  I=2,300+1
              XZNLLOC(I) = XZNLLOC(I-1) + DELX
           ENDDO

           DO  I=1,300
              XZNRLOC(I) = XZNLLOC(I) + DELX
           ENDDO

           DO  I=1,300
              XZNLOC(I)  = 0.5D0*(XZNLLOC(I) + XZNRLOC(I))
           ENDDO

           DELX = 0.1D0


           YZNLLOC(1) = 0.D0

           DO  J=2,150+1
              YZNLLOC(J) = YZNLLOC(J-1) + DELX
           ENDDO

           DO  J=1,150
              YZNRLOC(J) = YZNLLOC(J) + DELX
           ENDDO

           DO  J=1,150
              YZNLOC(J)  = 0.5D0*(YZNLLOC(J) + YZNRLOC(J))
           ENDDO


           DO k = 1, 300+1
              DO j = 1, 150+1
                  DO i = 1, 300+1

                       X = (/XZNLLOC(I),YZNLLOC(J),ZZNLLOC(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)

                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)

                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZNLLOC(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZNLLOC(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZNLLOC(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZNLLOC(I)+(DCOS(20.D0*PI/180))*YZNLLOC(J)+(-UX*DSIN(20.D0*PI/180))*ZZNLLOC(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZNLLOC(I)+(UX*DSIN(20.D0*PI/180))*YZNLLOC(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZNLLOC(K)

                        R = DSQRT(XROT**2+ZROT**2)


                    IF (timej .eq. 0.d0) THEN


                       AZ(I,J,K) = -(-0.00034D0*LOG(0.1369D0+R**2)*UX*DSIN(20.D0*PI/180))  &
                                   *((DTANH(-20.D0*(XN-0.008D0/RB))+1.D0)/2.D0)

                    ELSE

                       AZ(I,J,K) = -(-0.00034D0*LOG(0.1369D0+R**2)*UX*DSIN(20.D0*PI/180))  &
                                   *((DTANH(-20.D0*(XN-1.1D0*0.008D0/RB))+1.D0)/2.D0)

                    ENDIF

                 ENDDO
              ENDDO
           ENDDO


           DO k = 1, 300+1
              DO j = 1, 150+1
                  DO i = 1, 300+1

                    X = (/XZNLLOC(I),YZNLLOC(J),ZZNLLOC(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)

                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)

                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZNLLOC(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZNLLOC(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZNLLOC(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZNLLOC(I)+(DCOS(20.D0*PI/180))*YZNLLOC(J)+(-UX*DSIN(20.D0*PI/180))*ZZNLLOC(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZNLLOC(I)+(UX*DSIN(20.D0*PI/180))*YZNLLOC(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZNLLOC(K)

                        R = DSQRT(XROT**2+ZROT**2)


                    IF (timej .eq. 0.d0) THEN


                        AY(I,J,K) = (-0.00034D0*LOG(0.1369D0+R**2)*DCOS(20.D0*PI/180))  &
                                   *((DTANH(-20.D0*(XN-0.008D0/RB))+1.D0)/2.D0)

                    ELSE

                        AY(I,J,K) = (-0.00034D0*LOG(0.1369D0+R**2)*DCOS(20.D0*PI/180))  &
                                   *((DTANH(-20.D0*(XN-1.1D0*0.008D0/RB))+1.D0)/2.D0)

                    ENDIF


                 ENDDO
              ENDDO
           ENDDO


           DO k = 1, 300+1
              DO j = 1, 150+1
                  DO i = 1, 300+1

                       X = (/XZNLLOC(I),YZNLLOC(J),ZZNLLOC(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)

                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)

                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZNLLOC(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZNLLOC(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZNLLOC(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZNLLOC(I)+(DCOS(20.D0*PI/180))*YZNLLOC(J)+(-UX*DSIN(20.D0*PI/180))*ZZNLLOC(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZNLLOC(I)+(UX*DSIN(20.D0*PI/180))*YZNLLOC(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZNLLOC(K)

                        R = DSQRT(XROT**2+ZROT**2)


                     IF (timej .eq. 0.d0) THEN


                         AX(I,J,K) = (-0.00034D0*LOG(0.1369D0+R**2)*UZ*DSIN(20.D0*PI/180)) &
                                      *((DTANH(-20.D0*(XN-0.008D0/RB))+1.D0)/2.D0)

                     ELSE

                         AX(I,J,K) = (-0.00034D0*LOG(0.1369D0+R**2)*UZ*DSIN(20.D0*PI/180)) &
                                      *((DTANH(-20.D0*(XN-1.1D0*0.008D0/RB))+1.D0)/2.D0)


                     ENDIF


                 ENDDO
              ENDDO
           ENDDO


           DO k = 1, 300
              DO j = 1, 150
                  DO i = 1, 300+1

                       X = (/XZNLLOC(I),YZNLOC(J),ZZNLOC(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)


                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)


                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZNLLOC(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZNLOC(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZNLOC(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZNLLOC(I)+(DCOS(20.D0*PI/180))*YZNLOC(J)+(-UX*DSIN(20.D0*PI/180))*ZZNLOC(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZNLLOC(I)+(UX*DSIN(20.D0*PI/180))*YZNLOC(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZNLOC(K)

                        R = DSQRT(XROT**2+ZROT**2)


                     BX(i,j,k) = ((AY(I,J,K+1)+AY(I,J+1,K+1))/2.D0 - (AY(I,J,K)+AY(I,J+1,K))/2.D0)/DELX-&
                                 ((AZ(I,J+1,K)+AZ(I,J+1,K+1))/2.D0 - (AZ(I,J,K)+AZ(I,J,K+1))/2.D0)/DELX

                     BX(i,j,k) = BX(i,j,k)/DCOSH(R**8)



                  ENDDO
                ENDDO
           ENDDO


           DO k = 1, 300
              DO j = 1, 150+1
                  DO i = 1, 300


                       X = (/XZNLOC(I),YZNLLOC(J),ZZNLOC(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)


                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)


                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZNLOC(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZNLLOC(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZNLOC(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZNLOC(I)+(DCOS(20.D0*PI/180))*YZNLLOC(J)+(-UX*DSIN(20.D0*PI/180))*ZZNLOC(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZNLOC(I)+(UX*DSIN(20.D0*PI/180))*YZNLLOC(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZNLOC(K)

                        R = DSQRT(XROT**2+ZROT**2)


                     BY(i,j,k) = -((AX(I,J,K+1)+AX(I+1,J,K+1))/2.D0 - (AX(I,J,K)+AX(I+1,J,K))/2.D0)/DELX + &
                                  ((AZ(I+1,J,K)+AZ(I+1,J,K+1))/2.D0 - (AZ(I,J,K)+AZ(I,J,K+1))/2.D0)/DELX

                        BY(i,j,k) = BY(i,j,k)/DCOSH(R**8)


                  ENDDO
                ENDDO
           ENDDO


              DO j = 1, 150
                 DO i = 1, 300

                     BZ(I,J,1) = 0.D0

                  ENDDO
                ENDDO



           DO j = 1, 150
              DO k = 1, 300
                 DO i = 1, 300

                     BZ(i,j,k+1) =  BZ(i,j,k)+BX(i,j,k)+BY(i,j,k)-BX(i+1,j,k)-BY(i,j+1,k)

                  ENDDO
                ENDDO
           ENDDO


           IF ((XZNL(1) .LE. XZNLLOC(300+1) .AND. XZNLLOC(1) .LE. XZNL(NX+1)) .AND. &
                (YZN(1)  .LE. YZNLOC(150)  .AND. YZNLOC(1)  .LE. YZN(NY))  .AND. &
                (ZZN(1)  .LE. ZZNLOC(300)  .AND. ZZNLOC(1)  .LE. ZZN(NZ))) THEN

                IF (XZNL(1) .LT. XZNLLOC(1)) THEN

			 IF (mpi_coords(1) .EQ. 0) THEN
                             FZNPX_L = FLOOR((XZNLLOC(1)-XZNL(11))/DELX+0.5D0)+11
			 ELSE
                             FZNPX_L = FLOOR((XZNLLOC(1)-XZNL(1))/DELX+0.5D0)+1
			 ENDIF

                         FZNPX_G = 1
                ELSE
                         FZNPX_L = 1
                         FZNPX_G = FLOOR((XZNL(1)-XZNLLOC(1))/DELX+0.5D0)+1
                ENDIF

                IF (YZN(1) .LT. YZNLOC(1)) THEN
                       FZNPY_L = FLOOR((YZNLOC(1)-YZN(1))/DELX+0.5D0)+1
                       FZNPY_G = 1
                ELSE
                        FZNPY_L = 1
                        FZNPY_G = FLOOR((YZN(1)-YZNLOC(1))/DELX+0.5D0)+1
                ENDIF

                IF (ZZN(1) .LT. ZZNLOC(1)) THEN

                         IF (mpi_coords(3) .EQ. 0) THEN
                             FZNPZ_L = FLOOR((ZZNLOC(1)-ZZN(11))/DELX+0.5D0)+11
                         ELSE
                             FZNPZ_L = FLOOR((ZZNLOC(1)-ZZN(1))/DELX+0.5D0)+1
                         ENDIF

                       FZNPZ_G = 1
                ELSE
                        FZNPZ_L = 1
                        FZNPZ_G = FLOOR((ZZN(1)-ZZNLOC(1))/DELX+0.5D0)+1
                ENDIF


                IF (XZNL(NX+1) .GT. XZNLLOC(300+1)) THEN
		     
		       IF (mpi_coords(1) .EQ. mpi_dims(1)-1) THEN
                           NZNPX_L = (NX)-FLOOR((XZNL(NX+1-10)-XZNLLOC(300+1))/DELX+0.5D0)+1-10
		       ELSE
                           NZNPX_L = (NX)-FLOOR((XZNL(NX+1)-XZNLLOC(300+1))/DELX+0.5D0)+1
		       ENDIF

                ELSE
                        NZNPX_L = NX+1
                ENDIF

                IF (YZN(NY) .GT. YZNLOC(150)) THEN
                        NZNPY_L = NY-FLOOR((YZN(NY)-YZNLOC(150))/DELX+0.5D0)
                ELSE
                        NZNPY_L = NY
                ENDIF

                IF (ZZN(NZ) .GT. ZZNLOC(300)) THEN

                       IF (mpi_coords(3) .EQ. mpi_dims(3)-1) THEN
                           NZNPZ_L = (NZ)-FLOOR((ZZN(NZ-10)-ZZNLOC(300))/DELX+0.5D0)-10
                       ELSE
                           NZNPZ_L = (NZ)-FLOOR((ZZN(NZ)-ZZNLOC(300))/DELX+0.5D0)
                       ENDIF

                ELSE
                        NZNPZ_L = NZ
                ENDIF


            ELSE

              fznpx_L = 1
              nznpx_L = -1

              fznpy_L = 1
              nznpy_L = -1

              fznpz_L = 1
              nznpz_L = -1

              fznpx_G = 1
              nznpx_G = -1

              fznpy_G = 1
              nznpy_G = -1

              fznpz_G = 1
              nznpz_G = -1


           ENDIF


           DO k = fznpz_l, nznpz_l
              DO j = fznpy_l, nznpy_l
                 DO i = fznpx_l, nznpx_l


                       X = (/XZNL(I),YZN(J),ZZN(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)


                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)


                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZNL(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZN(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZNL(I)+(DCOS(20.D0*PI/180))*YZN(J)+(-UX*DSIN(20.D0*PI/180))*ZZN(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZNL(I)+(UX*DSIN(20.D0*PI/180))*YZN(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                        R = DSQRT(XROT**2+ZROT**2)


                    IF (timej .EQ. 0.D0) THEN

                        CT(i,j,k)%bx = BX(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

                    ELSE


                        IF (R .Le. 0.008D0/RB*DTAN(6.D0*PI/180) .and. xn .le. 0.008d0/RB) THEN

                            CT(i,j,k)%bx = BX(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

                        ELSE IF (R .gt. 0.008D0/RB*DTAN(6.D0*PI/180) .and. absx .le. 0.008d0/RB) THEN

                            CT(i,j,k)%bx = 0.D0

                        ENDIF


                    ENDIF

                  ENDDO
                ENDDO
           ENDDO


           IF ((XZN(1) .LE. XZNLLOC(300) .AND. XZNLOC(1) .LE. XZN(NX)) .AND. &
                (YZNL(1)  .LE. YZNLLOC(150+1)  .AND. YZNLLOC(1)  .LE. YZNL(NY+1))  .AND. &
                (ZZN(1)  .LE. ZZNLOC(300)  .AND. ZZNLOC(1)  .LE. ZZN(NZ))) THEN

                IF (XZN(1) .LT. XZNLOC(1)) THEN

                         IF (mpi_coords(1) .EQ. 0) THEN
                             FZNPX_L = FLOOR((XZNLOC(1)-XZN(11))/DELX+0.5D0)+11
                         ELSE
                             FZNPX_L = FLOOR((XZNLOC(1)-XZN(1))/DELX+0.5D0)+1
                         ENDIF

                         FZNPX_G = 1
                ELSE
                         FZNPX_L = 1
                         FZNPX_G = FLOOR((XZN(1)-XZNLOC(1))/DELX+0.5D0)+1
                ENDIF

                IF (YZNL(1) .LT. YZNLLOC(1)) THEN
                       FZNPY_L = FLOOR((YZNLLOC(1)-YZNL(1))/DELX+0.5D0)+1
                       FZNPY_G = 1
                ELSE
                        FZNPY_L = 1
                        FZNPY_G = FLOOR((YZNL(1)-YZNLLOC(1))/DELX+0.5D0)+1
                ENDIF

                IF (ZZN(1) .LT. ZZNLOC(1)) THEN

                         IF (mpi_coords(3) .EQ. 0) THEN
                             FZNPZ_L = FLOOR((ZZNLOC(1)-ZZN(11))/DELX+0.5D0)+11
                         ELSE
                             FZNPZ_L = FLOOR((ZZNLOC(1)-ZZN(1))/DELX+0.5D0)+1
                         ENDIF

                       FZNPZ_G = 1
                ELSE
                        FZNPZ_L = 1
                        FZNPZ_G = FLOOR((ZZN(1)-ZZNLOC(1))/DELX+0.5D0)+1
                ENDIF


                IF (XZN(NX) .GT. XZNLOC(300)) THEN

                       IF (mpi_coords(1) .EQ. mpi_dims(1)-1) THEN
                           NZNPX_L = (NX)-FLOOR((XZN(NX-10)-XZNLOC(300))/DELX+0.5D0)-10
                       ELSE
                           NZNPX_L = (NX)-FLOOR((XZN(NX)-XZNLOC(300))/DELX+0.5D0)
                       ENDIF

                ELSE
                        NZNPX_L = NX
                ENDIF

                IF (YZNL(NY+1) .GT. YZNLLOC(150+1)) THEN
                        NZNPY_L = NY-FLOOR((YZNL(NY+1)-YZNLLOC(150+1))/DELX+0.5D0)+1
                ELSE
                        NZNPY_L = NY+1
                ENDIF

                IF (ZZN(NZ) .GT. ZZNLOC(300)) THEN

                       IF (mpi_coords(3) .EQ. mpi_dims(3)-1) THEN
                           NZNPZ_L = (NZ)-FLOOR((ZZN(NZ-10)-ZZNLOC(300))/DELX+0.5D0)-10
                       ELSE
                           NZNPZ_L = (NZ)-FLOOR((ZZN(NZ)-ZZNLOC(300))/DELX+0.5D0)
                       ENDIF

                ELSE
                        NZNPZ_L = NZ
                ENDIF

            ELSE

              fznpx_L = 1
              nznpx_L = -1

              fznpy_L = 1
              nznpy_L = -1

              fznpz_L = 1
              nznpz_L = -1

              fznpx_G = 1
              nznpx_G = -1

              fznpy_G = 1
              nznpy_G = -1

              fznpz_G = 1
              nznpz_G = -1


            ENDIF

           DO k = fznpz_l, nznpz_l
              DO j = fznpy_l, nznpy_l
                 DO i = fznpx_l, nznpx_l


                       X = (/XZN(I),YZNL(J),ZZN(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)


                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)


                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZNL(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZN(I)+(DCOS(20.D0*PI/180))*YZNL(J)+(-UX*DSIN(20.D0*PI/180))*ZZN(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+(UX*DSIN(20.D0*PI/180))*YZNL(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                        R = DSQRT(XROT**2+ZROT**2)


                     IF (timej .EQ. 0.D0) THEN

                         CT(i,j,k)%by = BY(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

                     ELSE


                        IF (R .Le. 0.008D0/RB*DTAN(6.D0*PI/180) .and. xn .le. 0.008d0/RB) THEN

                            CT(i,j,k)%by = BY(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

                        ELSE IF (R .gt. 0.008D0/RB*DTAN(6.D0*PI/180) .and. absx .le. 0.008d0/RB) THEN

                            CT(i,j,k)%by = 0.d0

                        ENDIF

                    ENDIF

                  ENDDO
                ENDDO
           ENDDO

           IF ((XZN(1) .LE. XZNLOC(300) .AND. XZNLOC(1) .LE. XZN(NX)) .AND. &
                (YZN(1)  .LE. YZNLOC(150)  .AND. YZNLOC(1)  .LE. YZN(NY))  .AND. &
                (ZZNL(1)  .LE. ZZNLLOC(300+1)  .AND. ZZNLLOC(1)  .LE. ZZNL(NZ+1))) THEN

                IF (XZN(1) .LT. XZNLOC(1)) THEN

                         IF (mpi_coords(1) .EQ. 0) THEN
                             FZNPX_L = FLOOR((XZNLOC(1)-XZN(11))/DELX+0.5D0)+11
                         ELSE
                             FZNPX_L = FLOOR((XZNLOC(1)-XZN(1))/DELX+0.5D0)+1
                         ENDIF

                         FZNPX_G = 1
                ELSE
                         FZNPX_L = 1
                         FZNPX_G = FLOOR((XZN(1)-XZNLOC(1))/DELX+0.5D0)+1
                ENDIF

                IF (YZN(1) .LT. YZNLOC(1)) THEN
                       FZNPY_L = FLOOR((YZNLOC(1)-YZN(1))/DELX+0.5D0)+1
                       FZNPY_G = 1
                ELSE
                        FZNPY_L = 1
                        FZNPY_G = FLOOR((YZN(1)-YZNLOC(1))/DELX+0.5D0)+1
                ENDIF

                IF (ZZNL(1) .LT. ZZNLLOC(1)) THEN

                         IF (mpi_coords(3) .EQ. 0) THEN
                             FZNPZ_L = FLOOR((ZZNLLOC(1)-ZZNL(11))/DELX+0.5D0)+11
                         ELSE
                             FZNPZ_L = FLOOR((ZZNLLOC(1)-ZZNL(1))/DELX+0.5D0)+1
                         ENDIF

                       FZNPZ_G = 1
                ELSE
                        FZNPZ_L = 1
                        FZNPZ_G = FLOOR((ZZNL(1)-ZZNLLOC(1))/DELX+0.5D0)+1
                ENDIF


                IF (XZN(NX) .GT. XZNLOC(300)) THEN

                       IF (mpi_coords(1) .EQ. mpi_dims(1)-1) THEN
                           NZNPX_L = (NX)-FLOOR((XZN(NX-10)-XZNLOC(300))/DELX+0.5D0)-10
                       ELSE
                           NZNPX_L = (NX)-FLOOR((XZN(NX)-XZNLOC(300))/DELX+0.5D0)
                       ENDIF

                ELSE
                        NZNPX_L = NX
                ENDIF

                IF (YZN(NY) .GT. YZNLOC(150)) THEN
                        NZNPY_L = NY-FLOOR((YZN(NY)-YZNLOC(150))/DELX+0.5D0)
                ELSE
                        NZNPY_L = NY
                ENDIF

                IF (ZZNL(NZ+1) .GT. ZZNLLOC(300+1)) THEN

                       IF (mpi_coords(3) .EQ. mpi_dims(3)-1) THEN
                           NZNPZ_L = (NZ)-FLOOR((ZZNL(NZ+1-10)-ZZNLLOC(300+1))/DELX+0.5D0)+1-10
                       ELSE
                           NZNPZ_L = (NZ)-FLOOR((ZZNL(NZ+1)-ZZNLLOC(300+1))/DELX+0.5D0)+1
                       ENDIF

                ELSE
                        NZNPZ_L = NZ+1
                ENDIF

            ELSE

              fznpx_L = 1
              nznpx_L = -1

              fznpy_L = 1
              nznpy_L = -1

              fznpz_L = 1
              nznpz_L = -1

              fznpx_G = 1
              nznpx_G = -1

              fznpy_G = 1
              nznpy_G = -1

              fznpz_G = 1
              nznpz_G = -1

          ENDIF

           DO k = fznpz_l, nznpz_l
              DO j = fznpy_l, nznpy_l
                 DO i = fznpx_l, nznpx_l

                       X = (/XZN(I),YZN(J),ZZNL(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)


                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)


                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZN(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZNL(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZN(I)+(DCOS(20.D0*PI/180))*YZN(J)+(-UX*DSIN(20.D0*PI/180))*ZZNL(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+(UX*DSIN(20.D0*PI/180))*YZN(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZNL(K)

                        R = DSQRT(XROT**2+ZROT**2)



                     IF (timej .EQ. 0.D0) THEN

                         CT(i,j,k)%bz = BZ(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

                     ELSE


                        IF (R .Le. 0.008D0/RB*DTAN(6.D0*PI/180) .and. xn .le. 0.008d0/RB) THEN

                            CT(i,j,k)%bz = BZ(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

                        ELSE IF (R .gt. 0.008D0/RB*DTAN(6.D0*PI/180) .and. absx .le. 0.008d0/RB) THEN

                            CT(i,j,k)%bz = 0.d0

                        ENDIF

                    ENDIF


                  ENDDO
                ENDDO
           ENDDO

           DO k = 1, nz
              DO j = 1, ny
                 DO i = 1, nx


                       X = (/XZN(I),YZN(J),ZZN(K)/)

                    ABSX = DSQRT(X(1)**2+X(2)**2+X(3)**2)


                    XN = (X(1)*N(1)+X(2)*N(2)+X(3)*N(3))/ABSN

                    ALPHA = DACOS(XN*ABSN/(ABSX*ABSN))

                        ux = -DSIN(PHI)
                        uz =  DCOS(PHI)

                     XROT = (DCOS(20.D0*PI/180)+UX**2*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+&
                            (-UZ*DSIN(20.D0*PI/180))*YZN(J)+(UX*UZ*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                     YROT = (UZ*DSIN(20.D0*PI/180))*XZN(I)+(DCOS(20.D0*PI/180))*YZN(J)+(-UX*DSIN(20.D0*PI/180))*ZZN(K)

                     ZROT = (UZ*UX*(1.D0-DCOS(20.D0*PI/180)))*XZN(I)+(UX*DSIN(20.D0*PI/180))*YZN(J)+&
                            (DCOS(20.D0*PI/180)+UZ**2*(1.D0-DCOS(20.D0*PI/180)))*ZZN(K)

                        R = DSQRT(XROT**2+ZROT**2)


                     IF (timej .EQ. 0.D0) THEN

                         physic(i,j,k)%bxc = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                         physic(i,j,k)%byc = (CT(i,j,k)%by + CT(i,j+1,k)%by)/2.D0
                         physic(i,j,k)%bzc = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)/2.D0


                     ELSE


                        IF (R .Le. 0.008D0/RB*DTAN(6.D0*PI/180) .and. xn .le. 0.008d0/RB) THEN

                            physic(i,j,k)%bxc = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                            physic(i,j,k)%byc = (CT(i,j,k)%by + CT(i,j+1,k)%by)/2.D0
                            physic(i,j,k)%bzc = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)/2.D0


                        ELSE IF (R .ge. 0.008D0/RB*DTAN(6.D0*PI/180) .and. absx .le. 0.008d0/RB) THEN

                            physic(i,j,k)%bxc = 0.D0
                            physic(i,j,k)%byc = 0.D0
                            physic(i,j,k)%bzc = 0.D0


                        ENDIF


                     ENDIF


                  ENDDO
                ENDDO
           ENDDO


        ENDIF


        deallocate(BX,BY,BZ, AX, AY, AZ)        

        !------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE INIT
! ===============================================================================


      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


! NAME:
!       RAN1

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
!       �\200\234Minimal�\200\235 random number generator of Park and Miller with Bays-Durham shuffle and
!        added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
!        the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
!        alter idum between successive deviates in a sequence. RNMX should approximate the largest
!        floating value that is less than 1.

! ACADEMIC REFERENCE:


     
! ===============================================================================      
FUNCTION RAN1(idum)
! ===============================================================================

      
# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

# include "types.h"

        
        !------------------------------------------------------- Local Variables

       
        INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
        REAL(rs) ran1,AM,EPS,RNMX
        PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
        NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)


        INTEGER j,k,iv(NTAB),iy
        SAVE iv,iy
        DATA iv /NTAB*0/, iy /0/


        !------------------------------------------------------------------------- BEGIN

        
        if (idum.le.0.or.iy.eq.0) then

           idum=max(-idum,1)

           do j=NTAB+8,1,-1
              k=idum/IQ
              idum=IA*(idum-k*IQ)-IR*k
              if (idum.lt.0) idum=idum+IM
              if (j.le.NTAB) iv(j)=idum
           enddo
           
           iy=iv(1)
           
        endif

        k=idum/IQ

        idum=IA*(idum-k*IQ)-IR*k

        if (idum.lt.0) idum=idum+IM
        
        j=1+iy/NDIV
        iy=iv(j)
        iv(j)=idum
        ran1=min(AM*iy,RNMX)
        
        RETURN


        !------------------------------------------------------------------------- END

! ===============================================================================            
END      
! ===============================================================================      
