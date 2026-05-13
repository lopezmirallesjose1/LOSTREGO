
!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************

! NAME:
!       SPHERE

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
SUBROUTINE SPHERE (basenm, nx, ny, nz, mnx5, mny5, mnz5, mnx6, mny6, mnz6, eqs, iars,                         &
                 xzn, xznl,xznr,yznr, delx, yzn, yznl,zzn,zznl, physic, gamma_ad, gridlx, gridly, gridlz, CT, &
                 me, mp, m_0, kb, c, uwi, vwi, shift, Rorb, xCPU, yCPU, zCPU,ic,jc,kc,ncl,timej,pc,dt)
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


	real(rs) :: me, mp, kb, m_0, c, TIMEJ,pc, dt

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

	real(rs) :: B0, XI, KPARA, KPERP, SIGMA, SIGMAHOT, EX, EY, EZ, RHO0, H0

	real(rs) :: RB, PERIOD, PHI, H(3), R3(3), N(3), X(3), XN, ABSX, ABSN, ALPHA

	real(rs) :: UX, UY, UZ, BX_ANALYTIC, BY_ANALYTIC, BZ_ANALYTIC, BMZ,XROT,YROT,ZROT, reset, VB, LB

	integer(is) :: fznpx_L, nznpx_L, fznpy_L, nznpy_L, fznpz_L, nznpz_L, &
                       fznpx_G, nznpx_G, fznpy_G, nznpy_G, fznpz_G, nznpz_G

# ifdef OPENMP
	integer :: omp_get_thread_num, omp_get_num_threads,                 &
                   omp_get_max_threads,curth
# endif


	logical(ls):: inject
	integer(is) :: NN, nblobs
	real(rs) :: phi2, OPHI

	!COMMON /random/ idum

        logical, save :: first_time=.FALSE.
	integer(is), save :: idum=-3
	real(rs), save :: OVB=0.3D0

!$omp threadprivate(first_time,idum,OVB)

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



        ! ==================================================
        ! SIMULATIONS OF SS433 SUBPARSEC JETS
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


        IF (BASENM(1:2) .EQ. 'S4') THEN


           IF (BASENM .EQ. 'S433' .OR. BASENM .EQ. 'S431' .OR. BASENM .EQ. 'S432') THEN
               !idum = -3
               !VB = ran1(idum)*0.08D0 + 0.22D0
	        VB = 0.26D0
                LB = OVB*5.D0 - 0.4D0	
           ENDIF

!$OMP PARALLEL DEFAULT(shared)                                      & 
!$OMP PRIVATE(tid,i,j,k,rho_s,xe_s,tem_s,eps_s)                     &
!$OMP PRIVATE(p_s,gam_s,h_s,cs_s,dpdrh_s,dpdeps_s)                  &
!$OMP PRIVATE(SCRCH1, SCRCH2, X, ABSX, XN, ALPHA, R, DIVE, PHI)     &
!$OMP PRIVATE(WB, B40, B41, B42, B43, B44,UX,UZ,XROT,YROT,ZROT,NN,IDUM)           
           
# ifdef OPENMP
            tid = omp_get_thread_num()+1
# else
            tid = 1
# endif           
   
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
           
           gamma_ad = 5.D0/3.D0
           RB = 0.008D0*DTAN(5.D0*PI/180) ! PARSEC
           !VB = 0.26D0

	   !IF (BASENM .EQ. 'S433') THEN
           !    idum = -3
           !    VB = ran1(idum)*0.08D0 + 0.22D0
	   !    print*, VB, cartrank

	   !ENDIF

           Period = 162.D0*86400.D0/(RB*pc/c)

           PHI = 2.D0*PI*timej/Period + PI/4.D0
           OPHI = 2.D0*PI*(timej-dt)/Period + PI/4.D0


	   IF (BASENM .EQ. 'S433' .OR. BASENM .EQ. 'S431' .OR. BASENM .EQ. 'S432') THEN

		inject = .FALSE.

		IF (first_time) THEN

		    first_time = .FALSE.
		    inject = .TRUE.

		    OVB = VB

!$OMP MASTER
                        print*, 'injecting first blob with velocity:', VB
!$OMP END MASTER

		ENDIF

		nblobs = 25

		DO NN = 1, nblobs*10

		    PHI2 = NN*2.D0*PI/nblobs + PI/4.D0

		    !IF (PHI/PHI2 .GE. 1.D0 .AND. OPHI/PHI2 .LT. 1.D0) THEN
                    !IF (PHI/PHI2 .GE. 1.D0 .AND. OPHI/PHI2 .LT. 1.D0 .AND. &
		    !	(timej .LT. 1010.D0 .OR. (timej .GT. 1050.D0 .AND. &
                    !     timej .LT. 1070.D0) .OR. timej .GT. 1110.D0)) THEN


                    IF (PHI/PHI2 .GE. 1.D0 .AND. OPHI/PHI2 .LT. 1.D0 .AND. &
                        (timej .LT. 1010.D0 .OR. timej .GT. 1050.D0)) THEN


		        inject = .TRUE.
			!PHI = PHI - (2.D0*PI/nblobs)

			OVB = VB

!$OMP MASTER
		        print*, 'injecting blob with velocity:', VB
!$OMP END MASTER 
                        EXIT

		    ENDIF
		
		ENDDO
	

	   ENDIF


           H = (/0.D0,0.0085D0/RB,0.D0/)
           R3 = (/H(2)*DTAN(20.D0*PI/180)*DCOS(PHI),0.D0,H(2)*DTAN(20.D0*PI/180)*DSIN(PHI)/)
           N =  H + R3 !AXIS OF PRECESSION CONE
           ABSN = DSQRT(N(1)**2+N(2)**2+N(3)**2)


!$OMP DO SCHEDULE (STATIC)
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



		    IF (R .le. 0.008D0/RB*DTAN(6.D0*PI/180) .AND. XN .LE. (0.008D0/RB + LB) .AND. XN .GT. 0.008D0/RB) THEN

			IF (inject) THEN

                	    physic(i,j,k)%velx = VB*N(1)/ABSN
                            physic(i,j,k)%vely = VB*N(2)/ABSN
                            physic(i,j,k)%velz = VB*N(3)/ABSN

			    !if (alpha .le. 5.D0*PI/180.D0) THEN

                            !physic(i,j,k)%velx = VB*X(1)/ABSX
                            !physic(i,j,k)%vely = VB*X(2)/ABSX
                            !physic(i,j,k)%velz = VB*X(3)/ABSX

			    !endif



                        !WB = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)
                        wb = 1.d0/dsqrt(1.d0-VB**2)

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

			!ELSE

                        !    physic(i,j,k)%velx = 0.D0
                        !    physic(i,j,k)%vely = 0.D0
                        !    physic(i,j,k)%velz = 0.D0

                        !    physic(i,j,k)%densty  = 1.D0
                        !    physic(i,j,k)%denstye = 1.D0*me/(me+mp)
                        !    physic(i,j,k)%tracer  = 0.D0
                        !    physic(i,j,k)%pres    = 2.66D-7

			!ENDIF


			!WB = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2) 
			!wb = 1.d0/dsqrt(1.d0-0.26d0**2)

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


		    ELSE IF ((R .le. 0.008D0/RB*DTAN(6.D0*PI/180) .AND. XN .LE. 0.008D0/RB) .AND. &
                            (timej .GE. 1020.D0 .AND. timej .LE. 1040.D0)) THEN


                        physic(i,j,k)%velx = 0.75D0*N(1)/ABSN
                        physic(i,j,k)%vely = 0.75D0*N(2)/ABSN
                        physic(i,j,k)%velz = 0.75D0*N(3)/ABSN


                        !WB = 1.D0/DSQRT(1.D0-physic(i,j,k)%velx**2-physic(i,j,k)%vely**2-physic(i,j,k)%velz**2)
                        wb = 1.d0/dsqrt(1.d0-0.75D0**2)

                        physic(i,j,k)%densty = 2.16536D0*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)
                        physic(i,j,k)%denstye= 2.16536D0*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)

                        IF (physic(i,j,k)%densty .LT. 1.D-6) physic(i,j,k)%densty = 1.D-6
                        IF (physic(i,j,k)%denstye .LT. 1.D-6) physic(i,j,k)%denstye = 1.D-6

                        physic(i,j,k)%tracer = 1.D0*((DTANH(-10.D0*(R-0.008D0/RB*DTAN(5.D0*PI/180)))+1.D0)/2.D0)

                        RBM = 0.37D0
                        BMJ = 1.33901D0
                         PB = 0.265721D0

                        BPHI1 = (2.D0*BMJ*(1.D0/RBM))/(1.D0+(1.D0/RBM)**2)

			IF (timej .LE. 1040.D0) THEN

			    physic(i,j,k)%pres = PB

			ELSE

                            physic(i,j,k)%pres   = 2.D0*(BMJ/(WB*(1.D0+(R/RBM)**2)))**2+PB-BPHI1**2/(2.D0*WB**2)*(1.D0+RBM**2)

		        ENDIF

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

		
			!IF (R .gt. 0.008D0/RB*DTAN(6.D0*PI/180) .AND. ABSX .LE. 0.008D0/RB) THEN
                        IF (ABSX .LT. 0.008D0/RB .OR. (ABSX .LT. 0.008D0/RB+0.7D0 .AND. inject)) THEN

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


		    ENDIF


                 ENDDO
              ENDDO
           ENDDO
!$OMP END DO

!$OMP BARRIER

           
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

           

!$OMP DO SCHEDULE (STATIC)
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


                    !ELSE IF ((timej .GE. 1020.D0 .AND. timej .LE. 1040.D0) .OR. (timej .GE. 1080.D0 .AND. timej .LE. 1100.D0)) THEN

                    !   AZ(I,J,K) = -(-0.495432D0*LOG(0.1369D0+R**2)*UX*DSIN(20.D0*PI/180))  &
                    !               *((DTANH(-20.D0*(XN-0.008D0/RB))+1.D0)/2.D0)

                    ELSE

                       AZ(I,J,K) = -(-0.00034D0*LOG(0.1369D0+R**2)*UX*DSIN(20.D0*PI/180))  &
                                   *((DTANH(-20.D0*(XN-0.008D0/RB-LB))+1.D0)/2.D0)

                    ENDIF

                 ENDDO
              ENDDO
           ENDDO
!$OMP END DO

!$OMP BARRIER
           
!$OMP DO SCHEDULE (STATIC)
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


                    !ELSE IF ((timej .GE. 1020.D0 .AND. timej .LE. 1040.D0) .OR. (timej .GE. 1080.D0 .AND. timej .LE. 1100.D0)) THEN

                    !    AY(I,J,K) = (-0.495432D0*LOG(0.1369D0+R**2)*DCOS(20.D0*PI/180))  &
                    !               *((DTANH(-20.D0*(XN-0.008D0/RB))+1.D0)/2.D0)

                    ELSE

                        AY(I,J,K) = (-0.00034D0*LOG(0.1369D0+R**2)*DCOS(20.D0*PI/180))  &
                                   *((DTANH(-20.D0*(XN-0.008D0/RB-LB))+1.D0)/2.D0)

                    ENDIF


                 ENDDO
              ENDDO
           ENDDO
!$OMP END DO

!$OMP BARRIER
                      
!$OMP DO SCHEDULE (STATIC)
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


                    !ELSE IF ((timej .GE. 1020.D0 .AND. timej .LE. 1040.D0) .OR. (timej .GE. 1080.D0 .AND. timej .LE. 1100.D0)) THEN

                    !     AX(I,J,K) = (-0.495432D0*LOG(0.1369D0+R**2)*UZ*DSIN(20.D0*PI/180)) &
                    !                  *((DTANH(-20.D0*(XN-0.008D0/RB))+1.D0)/2.D0)


                     ELSE

                         AX(I,J,K) = (-0.00034D0*LOG(0.1369D0+R**2)*UZ*DSIN(20.D0*PI/180)) &
                                      *((DTANH(-20.D0*(XN-0.008D0/RB-LB))+1.D0)/2.D0)


                     ENDIF


                 ENDDO
              ENDDO
           ENDDO
!$OMP END DO

!$OMP BARRIER
           
!$OMP DO SCHEDULE (STATIC)
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
!$OMP END DO

!$OMP BARRIER
           
!$OMP DO SCHEDULE (STATIC)

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
!$OMP END DO

!$OMP BARRIER
           
!$OMP DO SCHEDULE (STATIC)
              DO j = 1, 150
                 DO i = 1, 300

                     BZ(I,J,1) = 0.D0

                  ENDDO
                ENDDO
!$OMP END DO

!$OMP BARRIER
           
!$OMP DO SCHEDULE (STATIC)                
           DO j = 1, 150
              DO k = 1, 300
                 DO i = 1, 300

                     BZ(i,j,k+1) =  BZ(i,j,k)+BX(i,j,k)+BY(i,j,k)-BX(i+1,j,k)-BY(i,j+1,k)

                  ENDDO
                ENDDO
           ENDDO
!$OMP END DO

!$OMP BARRIER

           
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

              
!$OMP DO SCHEDULE (STATIC)                 
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

                        CT(i,j,k)%bx = BX(mpi_coords(1)*NX+i,mpi_coords(2)*NY+j,mpi_coords(3)*NZ+k)

                    ELSE


                        IF (R .Le. 0.008D0/RB*DTAN(6.D0*PI/180) .and. XN .LE. (0.008D0/RB + LB) .AND. XN .GE. 0.008D0/RB) THEN

		            IF (inject) THEN

                            CT(i,j,k)%bx = BX(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

			    !ELSE

			    !ct(i,j,k)%bx = 0.D0

			    ENDIF

                        !ELSE IF (R .LE. 0.008D0/RB*DTAN(6.D0*PI/180) .and. XN .le. 0.008d0/RB) THEN

                        !    IF (timej .GE. 1080.D0 .AND. timej .LE. 1100.D0) THEN

                        !    CT(i,j,k)%bx = BX(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

			!    ELSE

                        !    CT(i,j,k)%bx = 0.D0

			!    ENDIF

                        ELSE 

                            IF (ABSX .LT. 0.008D0/RB) THEN

                            CT(i,j,k)%bx = 0.D0

			    ENDIF

                        ENDIF


                    ENDIF

                  ENDDO
                ENDDO
           ENDDO
!$OMP END DO
!$OMP BARRIER



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


!$OMP DO SCHEDULE (STATIC)                 
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

                         CT(i,j,k)%by = BY(mpi_coords(1)*NX+i,mpi_coords(2)*NY+j,mpi_coords(3)*NZ+k)

                     ELSE


                        IF (R .Le. 0.008D0/RB*DTAN(6.D0*PI/180) .and. XN .LE. (0.008D0/RB + LB) .AND. XN .GE. 0.008D0/RB) THEN

			    IF (inject) THEN

                            CT(i,j,k)%by = BY(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

			    !ELSE

			    !CT(i,j,k)%by = 0.D0

			    ENDIF

                        !ELSE IF (R .LE. 0.008D0/RB*DTAN(6.D0*PI/180) .and. XN .le. 0.008d0/RB) THEN

                        !    IF (timej .GE. 1080.D0 .AND. timej .LE. 1100.D0) THEN

                        !    CT(i,j,k)%by = BY(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

                        !    ELSE

                        !    CT(i,j,k)%by = 0.D0

                        !    ENDIF


                        !ELSE IF (R .gt. 0.008D0/RB*DTAN(6.D0*PI/180) .and. absx .le. 0.008d0/RB) THEN
                        ELSE 

			    IF (absx .le. 0.008d0/RB) THEN

                            CT(i,j,k)%by = 0.d0

			    ENDIF

                        ENDIF

                    ENDIF

                  ENDDO
                ENDDO
           ENDDO
!$OMP END DO

!$OMP BARRIER


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


!$OMP DO SCHEDULE (STATIC)                 

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

                         CT(i,j,k)%bz = BZ(mpi_coords(1)*NX+i,mpi_coords(2)*NY+j,mpi_coords(3)*NZ+k)

                     ELSE


                        IF (R .Le. 0.008D0/RB*DTAN(6.D0*PI/180) .and. XN .LE. (0.008D0/RB + LB) .AND. XN .GE. 0.008D0/RB) THEN

			    IF (inject) THEN

                            CT(i,j,k)%bz = BZ(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

			    !ELSE

			    !CT(i,j,k)%bz = 0.D0

			    ENDIF

                        !ELSE IF (R .LE. 0.008D0/RB*DTAN(6.D0*PI/180) .and. XN .le. 0.008d0/RB) THEN

                        !    IF (timej .GE. 1080.D0 .AND. timej .LE. 1100.D0) THEN

                        !    CT(i,j,k)%bz = BZ(i-fznpx_l+fznpx_g,j-fznpy_l+fznpy_g,k-fznpz_l+fznpz_g)

                        !    ELSE

                        !    CT(i,j,k)%bz = 0.D0

                        !    ENDIF


                        !ELSE IF (R .gt. 0.008D0/RB*DTAN(6.D0*PI/180) .and. absx .le. 0.008d0/RB) THEN
                        ELSE 

			    IF (absx .le. 0.008d0/RB) THEN

                            CT(i,j,k)%bz = 0.d0

			    ENDIF

                        ENDIF

                    ENDIF


                  ENDDO
                ENDDO
           ENDDO
!$OMP END DO

!$OMP BARRIER

!$OMP DO SCHEDULE (STATIC)
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


                        IF (R .Le. 0.008D0/RB*DTAN(6.D0*PI/180) .and. XN .LE. (0.008D0/RB + LB) .AND. XN .GE. 0.008D0/RB) THEN

			    IF (inject) THEN

                            physic(i,j,k)%bxc = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                            physic(i,j,k)%byc = (CT(i,j,k)%by + CT(i,j+1,k)%by)/2.D0
                            physic(i,j,k)%bzc = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)/2.D0

			    !else

                            !physic(i,j,k)%bxc = 0.D0
                            !physic(i,j,k)%byc = 0.D0
                            !physic(i,j,k)%bzc = 0.D0

			    endif
			  
                        ELSE IF (R .LE. 0.008D0/RB*DTAN(6.D0*PI/180) .and. XN .le. 0.008d0/RB) THEN

                            IF (timej .GE. 1080.D0 .AND. timej .LE. 1100.D0) THEN

                            physic(i,j,k)%bxc = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                            physic(i,j,k)%byc = (CT(i,j,k)%by + CT(i,j+1,k)%by)/2.D0
                            physic(i,j,k)%bzc = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)/2.D0


                            ELSE

                            physic(i,j,k)%bxc = 0.D0
                            physic(i,j,k)%byc = 0.D0
                            physic(i,j,k)%bzc = 0.D0


                            ENDIF



                        !ELSE IF (R .ge. 0.008D0/RB*DTAN(6.D0*PI/180) .and. absx .le. 0.008d0/RB) THEN
			ELSE

                            IF (absx .le. 0.008d0/RB) THEN

                            physic(i,j,k)%bxc = 0.D0
                            physic(i,j,k)%byc = 0.D0
                            physic(i,j,k)%bzc = 0.D0

			    ENDIF


                        ENDIF


                     ENDIF


                  ENDDO
                ENDDO
           ENDDO
!$OMP END DO

!$OMP END PARALLEL

	ENDIF


	deallocate(BX,BY,BZ, AX, AY, AZ)

        !------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE SPHERE
! ===============================================================================

