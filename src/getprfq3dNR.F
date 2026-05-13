! NAME:
!       GETPRFQ2DNR

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
!       Recover primitive variables in 3D simulations.

! ACADEMIC REFERENCE:
!       " Relativistic MHD simulations of extragalactic jets."
!       T. Leismann, L.Anton, M.A.Aloy et al. 2005. A&A 436, 503-526 (2005).
!       See also: "An HLLC Riemann Solver for Relativistic Flows-II, Magnetohydrodynamics."
!       A. Mignone & G. Bodo. 2006. MNRAS. 368, 1040-1054


! ===============================================================================
SUBROUTINE GETPRFQ2DNR (nx, ny, nz, mnx1, mnx5, mny5, mnz5, eqs, niter,      &
                      gamma_ad, small, smallp, me, mp, physic)
! ===============================================================================  
        
        USE tipos
        USE parallel
        USE workarrays
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        integer(is):: nx, ny, nz, mnx1, mnx5, mny5, mnz5, eqs, niter

        real(rs)   :: small, smallp, me, mp, gamma_ad
        
        !------------------------------------------------------ Output Variables

        type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

        !------------------------------------------------------- Local Variables

        integer(is):: tid, i, j, k, ll, nl, nelem, iter, ierr

        character(26) :: message1
      
# ifdef OPENMP
        integer :: omp_get_thread_num, omp_get_num_threads,                 &
                   omp_get_max_threads,curth
# endif

       
        integer(is):: flag_zero (3, 20)
                             
                
        real(rs), pointer   :: xe_so(:), rho_so(:), tem_so(:), eps_so(:), p_so(:),   &
                               gam_so(:), h_so(:), cs_so(:), dpdrho_so(:), dpdeps_so(:)

        real(rs)   :: nve, HT, WT, B22, M2, BM, S2, ZETA, ERRP, E, R, G1G, &
		      SCRCH1, SCRCH2, FUN, DFUN, OZETA, B, C, CHI, ZETAP,  &
		      ZETA2, S2_ZETA2, VEL2, W2, DV2_DZETA, DCHI_DZETA, DRHO_DZETA, DP_DCHI, DP_DRHO, DP 		      


        real(rs), pointer :: rho_s(:),xe_s(:),tem_s(:),eps_s(:),p_s(:),         &
                             gam_s(:), h_s(:),cs_s(:),dpdrho_s(:),dpdeps_s(:)
                           

	real(rs), pointer :: w(:), denstw_s(:), denstew_s(:), energy_s(:),   &
			     momenx_s(:), momeny_s(:), momenz_s(:), pmin(:), &
			     velx_s(:), vely_s(:), velz_s(:), bx_s(:), by_s(:), bz_s(:)



        real(rs),pointer    :: rho_b(:), eps_b(:), p_b(:), h_b(:), cs_b(:),  &
                               dpdrho_b(:), dpdeps_b(:), gam_b(:), xe_b(:), tem_bb(:)


	real(rs) :: vx, vy, vz, LOR

        !------------------------------------------------------------------------- BEGIN

        tid = 1

        niter = 10000
        nelem = 0
        maxele = 100
	

        ! SMALLP prevents the magnetic field to take unphysical small values.
        
!$OMP PARALLEL DEFAULT(shared)   &
!$OMP PRIVATE(i,j,k)             


        DO k = 1, nz
!$OMP DO SCHEDULE(STATIC)
           DO j = 1, ny
              DO i = 1, nx

                 IF ((physic(i,j,k)%bxc*physic(i,j,k)%bxc) .LE. smallp) physic(i,j,k)%bxc = 0.D0
                 IF ((physic(i,j,k)%byc*physic(i,j,k)%byc) .LE. smallp) physic(i,j,k)%byc = 0.D0
                 IF ((physic(i,j,k)%bzc*physic(i,j,k)%bzc) .LE. smallp) physic(i,j,k)%bzc = 0.D0

              ENDDO
           ENDDO
!$OMP END DO
        ENDDO

!$OMP END PARALLEL



!$OMP PARALLEL DEFAULT(shared)                                               & 
!$OMP PRIVATE(tid,i,j,k,rho_s,rho_b,xe_s,xe_b,tem_s,tem_bb,eps_s,eps_b,iter) &
!$OMP PRIVATE(p_s,gam_s,h_s,cs_s,dpdrho_s,dpdeps_s)                          &
!$OMP PRIVATE(p_b,gam_b,h_b,cs_b,dpdrho_b,dpdeps_b)                          &
!$OMP PRIVATE(w,denstw_s,denstew_s, energy_s,momenx_s)                       &
!$OMP PRIVATE(momeny_s,momenz_s, velx_s, vely_s, velz_s)                     &
!$OMP PRIVATE(bx_s, by_s, bz_s, B22, M2, BM, S2, ZETA, OZETA)                &
!$OMP PRIVATE(ERRP,E,R,G1G,SCRCH1,SCRCH2,B,C,CHI,ZETAP,ZETA2)                &
!$OMP PRIVATE(FUN, DFUN,S2_ZETA2,VEL2,W2,DV2_DZETA,DCHI_DZETA)               &
!$OMP PRIVATE(pmin,nl, message1,DRHO_DZETA,DP_DCHI,DP_DRHO,DP)               &
!$OMP PRIVATE(vx,vy,vz,lor)                                                  &
!$OMP REDUCTION(+:nelem)


# ifdef OPENMP
            tid = omp_get_thread_num()+1
# else
            tid = 1
# endif


      
        ! -4:mn5 vectors
        p_s      => wa(tid)%wa5get100
        rho_s    => wa(tid)%wa5get101
        xe_s     => wa(tid)%wa5get102
        tem_s    => wa(tid)%wa5get103
        eps_s    => wa(tid)%wa5get104
        gam_s    => wa(tid)%wa5get105
        h_s      => wa(tid)%wa5get106
        cs_s     => wa(tid)%wa5get107
        dpdrho_s => wa(tid)%wa5get108
        dpdeps_s => wa(tid)%wa5get109


        ! -4:mn5 vectors
        p_b      => wa(tid)%wa5bis100
        rho_b    => wa(tid)%wa5bis101
        xe_b     => wa(tid)%wa5bis102
        tem_bb   => wa(tid)%wa5bis103
        eps_b    => wa(tid)%wa5bis104
        gam_b    => wa(tid)%wa5bis105
        h_b	 => wa(tid)%wa5bis106
        cs_b     => wa(tid)%wa5bis107
        dpdrho_b => wa(tid)%wa5bis108
        dpdeps_b => wa(tid)%wa5bis109


        w          => wa(tid)%wa5get110
        denstw_s   => wa(tid)%wa5get111
        denstew_s  => wa(tid)%wa5get112
        energy_s   => wa(tid)%wa5get113
        momenx_s   => wa(tid)%wa5get114
        momeny_s   => wa(tid)%wa5get115
        momenz_s   => wa(tid)%wa5get116
        pmin       => wa(tid)%wa5get117
        velx_s     => wa(tid)%wa5get118
        vely_s     => wa(tid)%wa5get119
        velz_s     => wa(tid)%wa5get120
        bx_s       => wa(tid)%wa5get121
        by_s       => wa(tid)%wa5get122
        bz_s       => wa(tid)%wa5get123

        ! Extract one-dimensional vectors.

        DO k=1, nz
!$OMP DO SCHEDULE (STATIC)
           DO j= 1, ny
              DO i= 1, nx

                 p_s(i)       = physic(i,j,k)%pres
                 denstw_s(i)  = physic(i,j,k)%denstw
                 denstew_s(i) = physic(i,j,k)%denstew
                 momenx_s(i)  = physic(i,j,k)%momenx
                 momeny_s(i)  = physic(i,j,k)%momeny
                 momenz_s(i)  = physic(i,j,k)%momenz
	         energy_s(i)  = physic(i,j,k)%energy
                 bx_s(i)      = physic(i,j,k)%bxc
                 by_s(i)      = physic(i,j,k)%byc
                 bz_s(i)      = physic(i,j,k)%bzc
                 xe_s(i)      = physic(i,j,k)%denstew/physic(i,j,k)%denstw

		 ! FIRST GUESS FOR THE TEMPERATURE BASED ON PREVIOUS PRIMITIVES (SYNGE).
                 tem_s(i)     = physic(i,j,k)%pres / (physic(i,j,k)%densty*                        &
                                (xe_s(i)/me + (1.D0-xe_s(i))/mp))    

    
                 pmin(i)     = DMAX1 (DSQRT(physic(i,j,k)%momenx**2+physic(i,j,k)%momeny**2 +         &
                                            physic(i,j,k)%momenz**2) - physic(i,j,k)%energy -         &
                                            physic(i,j,k)%denstw, smallp)        


              ENDDO

	      
              ! EOS_SYNGE :: OBTAIN FIRST GUESS FOR GAMMA FACTOR.
              IF (eqs .EQ. 2) THEN                 
                 call EOS_SYNGE (3, tid, 1, nx, me, mp, rho_s, xe_s, tem_s, eps_s,                 &
                                 p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)
              ENDIF


              DO i= 1, nx

                 B22 = bx_s(i)*bx_s(i)+by_s(i)*by_s(i) +bz_s(i)*bz_s(i)

                 M2 = momenx_s(i)*momenx_s(i) + momeny_s(i)*momeny_s(i) +                          &
                      momenz_s(i)*momenz_s(i)
                 BM = bx_s(i)*momenx_s(i) +  by_s(i)*momeny_s(i) +                                 &
                      bz_s(i)*momenz_s(i)
		 S2 = BM*BM    


                 iter = 0
		 ERRP = 1.D0

		 E = energy_s(i) 
		 R = denstw_s(i)


		 IF (eqs .EQ. 2) THEN
                     G1G = (gam_s(i) - 1.D0)/gam_s(i)
		 ELSE
                     G1G = (GAMMA_AD - 1.D0)/GAMMA_AD
		 ENDIF


		 ! Provide initial guess for ZETA.

		 
		 B = -4.D0*(E + R - B22)
		 C =  M2 + B22*B22 - 2.D0*B22*(E + R)

		 CHI = B*B - 12.D0*C
		 CHI = DMAX1(CHI, 0.D0)

	         ZETA = (-B + DSQRT(CHI))/6.D0
		 ZETA = DMAX1(R, ZETA)

		 ZETAP = ZETA -  R

		 GAM_B(I) = 0.D0
		
                 DO WHILE (ERRP .GT. 1.D-8  .AND. iter .LT. niter)

                    IF (iter + 1 .GT. niter) THEN

                        ierr = ierr+1

                        IF (ierr .NE. 0) THEN
                            call ERROR ('GETPRFQ3D', 'MAXIMUM NUMBER OF ITERATIONS',   &
                                         __FILE__, __LINE__)
                        ENDIF

                    ENDIF

		    ZETA = ZETAP + R

		    ZETA2 = ZETA*ZETA
		    S2_ZETA2 = S2/ZETA2

		    SCRCH1 = 1.D0/(ZETA + B22)
		    SCRCH2 = SCRCH1 * SCRCH1

		    VEL2 = S2_ZETA2*SCRCH1*(SCRCH1*ZETA + 1.D0) + M2*SCRCH2
		      W2 = 1.D0/(1.D0-VEL2)
		     W(I)= DSQRT(W2)

		    IF (VEL2 .GT. 1.D0) THEN

                        ierr = ierr+1

                        IF (ierr .NE. 0) THEN
                            call ERROR ('GETPRFQ3D', 'VEL2 > 0',   &
                                         __FILE__, __LINE__)
                        ENDIF

		    ENDIF

		    CHI = ZETAP/W2 - R*VEL2/(W(I)+1.D0)

		    DV2_DZETA = -2.D0*SCRCH2*(3.D0*S2_ZETA2 + SCRCH1*(S2_ZETA2*B22*B22/ZETA + M2))

		    RHO_S(I) = R/W(I)


		    ! KINEMATICAL TERMS

		    dchi_dZETA =  1.D0 - VEL2 - 0.5D0*W(I)*(R + 2.D0*chi*W(I))*dv2_dZETA
                    drho_dZETA = -0.5D0*R*W(I)*dv2_dZETA


		    ! ASSUMING IDEAL GAS
                    dp_dchi = G1G
                    dp_drho = 0.0D0
                    p_S(I)  = chi*dp_dchi
		    DP = dp_dchi*dchi_dZETA + dp_drho*drho_dZETA

		    FUN = ZETAP - (E+P_S(I)) + 0.5D0*(B22 + (B22*M2 - S2)*SCRCH2)
		    DFUN= 1.D0 - DP - (B22*M2 - S2)*SCRCH1*SCRCH2

		    OZETA = ZETAP			

		    ZETAP = ZETAP - FUN/DFUN

		    ERRP = ABS(1.D0-ZETAP/OZETA)

                    ITER = ITER + 1


                    ! =====================================================================
		    IF (eqs .EQ. 2) THEN

 	                gam_b(I)  = gam_s(I)
           	    	xe_b(I)   = xe_s(I)
            	    	tem_bb(I) = tem_s(I)

                        ! Since ZETAP is updated in each interation, we need to update
                        ! the gamma factor using the Synge equation.

           	    	vx = (momenx_s(i) + bm*bx_s(i)/(zetap+R))/((zetap+R)+b22)
           	    	vy = (momeny_s(i) + bm*by_s(i)/(zetap+R))/((zetap+R)+b22)
           	    	vz = (momenz_s(i) + bm*bz_s(i)/(zetap+R))/((zetap+R)+b22)

           	    	LOR	= 1.D0/DSQRT(1.D0-vx**2-vy**2-vz**2)
           	    	rho_b(I) = R/LOR

           	    	eps_b(I) = ((zetap+R)/rho_b(I)/lor/lor - 1.D0 ) / gam_b(I)


                        ! ------------------------------------------------------------


           	        call EOS_SYNGE (2, tid, I, I, me, mp, rho_b, xe_b, tem_bb, eps_b,   &
                                        p_b, gam_b, h_b, cs_b, dpdrho_b, dpdeps_b)

		        GAM_S(I) = gam_b(I)

                    	G1G = (GAM_S(I) - 1.D0)/GAM_S(I)
		
		    ENDIF
                    ! =======================================================================


                 ENDDO

		 ZETA = ZETAP + R
	        

                    ! Recover primitive velocity.
                    velx_s(i) = (momenx_s(i) + bm*bx_s(i)/zeta)/(zeta+b22)
                    vely_s(i) = (momeny_s(i) + bm*by_s(i)/zeta)/(zeta+b22)
                    velz_s(i) = (momenz_s(i) + bm*bz_s(i)/zeta)/(zeta+b22)

	         
                    ! SMALL prevents the velocity to take unphysical small values.
                    IF (ABS(velx_s(i)) .LE. small) velx_s(i) = 0.D0
                    IF (ABS(vely_s(i)) .LE. small) vely_s(i) = 0.D0
                    IF (ABS(velz_s(i)) .LE. small) velz_s(i) = 0.D0
                
                    ! Recover Lorentz factor.
                    w(i) = 1.D0/DSQRT(1.D0 - velx_s(i)*velx_s(i) - vely_s(i)*vely_s(i) - velz_s(i)*velz_s(i))
        
                    ! Recover primitive density.
                    rho_s(i)  = denstw_s(i)/w(i)
                 
                    ! Calculate the internal energy eps.

                    IF (eqs .EQ. 2) THEN
                       eps_s(i) = (zeta/rho_s(i)/w(i)/w(i) - 1.D0) / gam_s(i)
                    ELSE
                       eps_s(i) = (zeta/rho_s(i)/w(i)/w(i) - 1.D0) / gamma_ad
                    ENDIF

              ENDDO

              ! Calculate the pressure p.

              IF (eqs .EQ. 2) THEN
                 call EOS_SYNGE (2, tid, 1, nx, me, mp, rho_s, xe_s, tem_s, eps_s,        &
                                 p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)

              ELSE
                 call EOS (2, 1, nx, gamma_ad, eps_s, rho_s, p_s, dpdrho_s, dpdeps_s, h_s, cs_s)
              ENDIF


	      DO i = 1, nx
	
	         IF (p_s(i) .LE. SMALLP) THEN

                    !WRITE(6,*) 'GETPRFQ_MHD_2D: MINIMUM PRESSURE REACHED AT POINT:'
                    !print*, i,j,k
                    !WRITE(6,*) '             PRESSURE RESET TO SMALLP'
                    !WRITE(6,*) '             REMAINING THERMODYNAMICAL VARIABLES MADE CONSISTENT'

                    p_s(i) = SMALLP

                    call EOS (1, i, i, gamma_ad, eps_s, rho_s, p_s, dpdrho_s, dpdeps_s, h_s, cs_s)

                 END IF

	      ENDDO


              DO i = 1, nx      
              
                 physic(i,j,k)%pres = p_s(i)
                 physic(i,j,k)%densty = rho_s(i)
                 physic(i,j,k)%denstye = xe_s(i)*rho_s(i)
                 physic(i,j,k)%eps = eps_s(i)   
                 physic(i,j,k)%velx = velx_s(i)
                 physic(i,j,k)%vely = vely_s(i)
                 physic(i,j,k)%velz = velz_s(i)

              ENDDO
           ENDDO
!$OMP END DO
        ENDDO

!$OMP END PARALLEL



        !------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE GETPRFQ2DNR
! ===============================================================================



! NAME:
!       GETPRFQ3D

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
!       Recover primitive variables in 3D simulations.

! ACADEMIC REFERENCE:
!       " Relativistic MHD simulations of extragalactic jets."
!       T. Leismann, L.Anton, M.A.Aloy et al. 2005. A&A 436, 503-526 (2005).
!       See also: "An HLLC Riemann Solver for Relativistic Flows-II, Magnetohydrodynamics."
!       A. Mignone & G. Bodo. 2006. MNRAS. 368, 1040-1054


! ===============================================================================
SUBROUTINE GETPRFQ3DNR (nx, ny, nz, mnx1, mnx5, mny5, mnz5, eqs, niter,      &
                      gamma_ad, small, smallp, me, mp, physic)
! ===============================================================================  
        
        USE tipos
        USE parallel
        USE workarrays
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        integer(is):: nx, ny, nz, mnx1, mnx5, mny5, mnz5, eqs, niter

        real(rs)   :: small, smallp, me, mp, gamma_ad
        
        !------------------------------------------------------ Output Variables

        type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

        !------------------------------------------------------- Local Variables

        integer(is):: tid, i, j, k, ll, nl, nelem, iter, ierr

        character(26) :: message1
      
# ifdef OPENMP
        integer :: omp_get_thread_num, omp_get_num_threads,                 &
                   omp_get_max_threads,curth
# endif

       
        integer(is):: flag_HLL (3, 20)
                             
                
        real(rs), pointer   :: xe_so(:), rho_so(:), tem_so(:), eps_so(:), p_so(:),   &
                               gam_so(:), h_so(:), cs_so(:), dpdrho_so(:), dpdeps_so(:)

        real(rs)   :: nve, HT, WT, B22, M2, BM, S2, ZETA, ERRP, E, R, G1G, &
		      SCRCH1, SCRCH2, FUN, DFUN, OZETA, B, C, CHI, ZETAP,  &
		      ZETA2, S2_ZETA2, VEL2, W2, DV2_DZETA, DCHI_DZETA, DRHO_DZETA, DP_DCHI, DP_DRHO, DP 		      


        real(rs), pointer :: rho_s(:),xe_s(:),tem_s(:),eps_s(:),p_s(:),         &
                             gam_s(:), h_s(:),cs_s(:),dpdrho_s(:),dpdeps_s(:)
                           

        real(rs),pointer    :: rho_b(:), eps_b(:), p_b(:), h_b(:), cs_b(:),  &
                               dpdrho_b(:), dpdeps_b(:), gam_b(:), xe_b(:), tem_bb(:)


	real(rs), pointer :: w(:), denstw_s(:), denstew_s(:), energy_s(:),   &
			     momenx_s(:), momeny_s(:), momenz_s(:), pmin(:), &
			     velx_s(:), vely_s(:), velz_s(:), bx_s(:), by_s(:), bz_s(:)


        real(rs), pointer :: scrch11(:), scrch22(:)

        real(rs) :: b40, b41, b42, b43, b44

	real(rs) :: v2min, v2max, wmin, wmax, fmin, fmax, lor, vx, vy, vz


        !------------------------------------------------------------------------- BEGIN

        tid = 1

        ! 1D-vector
        p_so      => wa(tid)%wa1get100
        rho_so    => wa(tid)%wa1get101
        xe_so     => wa(tid)%wa1get102
        tem_so    => wa(tid)%wa1get103
        eps_so    => wa(tid)%wa1get104
        gam_so    => wa(tid)%wa1get105
        h_so      => wa(tid)%wa1get106
        cs_so     => wa(tid)%wa1get107
        dpdrho_so => wa(tid)%wa1get108
        dpdeps_so => wa(tid)%wa1get109
        

        niter = 10000
        nelem = 0
        maxele = 100

        ! SMALLP prevents the magnetic field to take unphysical small values.
        
!$OMP PARALLEL DEFAULT(shared)   &
!$OMP PRIVATE(i,j,k)             

!$OMP DO SCHEDULE(STATIC)
        DO k = 1, nz
           DO j = 1, ny
              DO i = 1, nx


                  IF (physic(i,j,k)%bxc*physic(i,j,k)%bxc .LE. SMALL) physic(i,j,k)%bxc = 0.D0
                  IF (physic(i,j,k)%byc*physic(i,j,k)%byc .LE. SMALL) physic(i,j,k)%byc = 0.D0
                  IF (physic(i,j,k)%bzc*physic(i,j,k)%bzc .LE. SMALL) physic(i,j,k)%bzc = 0.D0

              ENDDO
           ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL DEFAULT(shared)                                                           & 
!$OMP PRIVATE(tid,i,j,k,rho_s,rho_b,xe_s,xe_b,tem_s,tem_bb,eps_s,eps_b,iter)             &
!$OMP PRIVATE(p_s,p_b,gam_s,gam_b,h_s,h_b,cs_s,cs_b,dpdrho_s,dpdrho_b,dpdeps_s,dpdeps_b) &
!$OMP PRIVATE(w,denstw_s,denstew_s, energy_s,momenx_s)              &
!$OMP PRIVATE(momeny_s,momenz_s, velx_s, vely_s, velz_s)            &
!$OMP PRIVATE(bx_s, by_s, bz_s, B22, M2, BM, S2, ZETA, OZETA)       &
!$OMP PRIVATE(ERRP,E,R,G1G,SCRCH1,SCRCH2,B,C,CHI,ZETAP,ZETA2)       &
!$OMP PRIVATE(FUN, DFUN,S2_ZETA2,VEL2,W2,DV2_DZETA,DCHI_DZETA )     &
!$OMP PRIVATE(pmin,nl, message1,DRHO_DZETA,DP_DCHI,DP_DRHO,DP)      &
!$OMP PRIVATE(scrch11, scrch22, b40, b41, b42, b43, b44, flag_HLL)  &
!$OMP PRIVATE(v2min,v2max,wmin,wmax,fmin,fmax,lor,ll,vx,vy,vz)      &
!$OMP REDUCTION(+:nelem)


# ifdef OPENMP
            tid = omp_get_thread_num()+1
# else
            tid = 1
# endif


        ! -4:mn6 vectors
        scrch11 => wa(tid)%wa6rest100
        scrch22 => wa(tid)%wa6rest101

      
        ! -4:mn5 vectors
        p_s      => wa(tid)%wa5get100
        rho_s    => wa(tid)%wa5get101
        xe_s     => wa(tid)%wa5get102
        tem_s    => wa(tid)%wa5get103
        eps_s    => wa(tid)%wa5get104
        gam_s    => wa(tid)%wa5get105
        h_s      => wa(tid)%wa5get106
        cs_s     => wa(tid)%wa5get107
        dpdrho_s => wa(tid)%wa5get108
        dpdeps_s => wa(tid)%wa5get109


        ! -4:mn5 vectors
        p_b      => wa(tid)%wa5bis100
        rho_b    => wa(tid)%wa5bis101
        xe_b     => wa(tid)%wa5bis102
        tem_bb   => wa(tid)%wa5bis103
        eps_b    => wa(tid)%wa5bis104
        gam_b    => wa(tid)%wa5bis105
        h_b      => wa(tid)%wa5bis106
        cs_b     => wa(tid)%wa5bis107
        dpdrho_b => wa(tid)%wa5bis108
        dpdeps_b => wa(tid)%wa5bis109

        w          => wa(tid)%wa5get110
        denstw_s   => wa(tid)%wa5get111
        denstew_s  => wa(tid)%wa5get112
        energy_s   => wa(tid)%wa5get113
        momenx_s   => wa(tid)%wa5get114
        momeny_s   => wa(tid)%wa5get115
        momenz_s   => wa(tid)%wa5get116
        pmin       => wa(tid)%wa5get117
        velx_s     => wa(tid)%wa5get118
        vely_s     => wa(tid)%wa5get119
        velz_s     => wa(tid)%wa5get120
        bx_s       => wa(tid)%wa5get121
        by_s       => wa(tid)%wa5get122
        bz_s       => wa(tid)%wa5get123

        ! Extract one-dimensional vectors.


!$OMP DO SCHEDULE (STATIC)
        DO k=1, nz
           DO j= 1, ny
              DO i= 1, nx

                 p_s(i)       = physic(i,j,k)%pres
                 denstw_s(i)  = physic(i,j,k)%denstw
                 denstew_s(i) = physic(i,j,k)%denstew
                 momenx_s(i)  = physic(i,j,k)%momenx
                 momeny_s(i)  = physic(i,j,k)%momeny
                 momenz_s(i)  = physic(i,j,k)%momenz
	         energy_s(i)  = physic(i,j,k)%energy
                 bx_s(i)      = physic(i,j,k)%bxc
                 by_s(i)      = physic(i,j,k)%byc
                 bz_s(i)      = physic(i,j,k)%bzc
                 xe_s(i)      = physic(i,j,k)%denstew/physic(i,j,k)%denstw
                 tem_s(i)     = physic(i,j,k)%pres / (physic(i,j,k)%densty*                        &
                                (xe_s(i)/me + (1.D0-xe_s(i))/mp))    

    
                 pmin(i)     = DMAX1 (DSQRT(physic(i,j,k)%momenx**2+physic(i,j,k)%momeny**2 +      &
                                            physic(i,j,k)%momenz**2) - physic(i,j,k)%energy -      &
                                            physic(i,j,k)%denstw, smallp)        


              ENDDO

	      
              ! EOS_SYNGE :: Obtain gamma factor for BISECT.
              IF (eqs .EQ. 2) THEN                 
                 call EOS_SYNGE (3, tid, 1, nx, me, mp, rho_s, xe_s, tem_s, eps_s,                 &
                                 p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)
              ENDIF

              DO i= 1, nx

                 B22 = bx_s(i)*bx_s(i)+by_s(i)*by_s(i) +bz_s(i)*bz_s(i)

                 M2 = momenx_s(i)*momenx_s(i) + momeny_s(i)*momeny_s(i) +                          &
                      momenz_s(i)*momenz_s(i)
                 BM = bx_s(i)*momenx_s(i) +  by_s(i)*momeny_s(i) +                                 &
                      bz_s(i)*momenz_s(i)
		 S2 = BM*BM    


                 iter = 0
		 ERRP = 1.D0

		 E = energy_s(i) 
		 R = denstw_s(i)

                 IF (eqs .EQ. 2) THEN
                     G1G = (gam_s(i) - 1.D0)/gam_s(i)
                 ELSE
                     G1G = (GAMMA_AD - 1.D0)/GAMMA_AD
                 ENDIF

		 ! Provide initial guess.

		 
		 B = -4.D0*(E + R - B22)
		 C =  M2 + B22*B22 - 2.D0*B22*(E + R)

		 CHI = B*B - 12.D0*C
		 CHI = DMAX1(CHI, 0.D0)

	         ZETA = (-B + DSQRT(CHI))/6.D0
		 ZETA = DMAX1(R, ZETA)

		 ZETAP = ZETA -  R

                 DO WHILE (ERRP .GT. 1.D-10 .AND. iter .LT. niter)

                    IF (iter + 1 .GT. niter) THEN

                        ierr = ierr+1

                        IF (ierr .NE. 0) THEN
                            call ERROR ('GETPRFQ3D', 'MAXIMUM NUMBER OF ITERATIONS',   &
                                         __FILE__, __LINE__)
                        ENDIF

                    ENDIF

		    ZETA = ZETAP + R

		    ZETA2 = ZETA*ZETA
		    S2_ZETA2 = S2/ZETA2

		    SCRCH1 = 1.D0/(ZETA + B22)
		    SCRCH2 = SCRCH1 * SCRCH1

		    VEL2 = S2_ZETA2*SCRCH1*(SCRCH1*ZETA + 1.D0) + M2*SCRCH2
		      W2 = 1.D0/(1.D0-VEL2)
		     W(I)= DSQRT(W2)

		    IF (VEL2 .GT. 1.D0) THEN


			print*, i,j,k
			print*, mpi_coords(1), mpi_coords(2), mpi_coords(3)

                        ierr = ierr+1

                        IF (ierr .NE. 0) THEN
                            call ERROR ('GETPRFQ3D', 'VEL2 > 0',   &
                                         __FILE__, __LINE__)

                        ENDIF

		    ENDIF

		    CHI = ZETAP/W2 - R*VEL2/(W(I)+1.D0)

		    DV2_DZETA = -2.D0*SCRCH2*(3.D0*S2_ZETA2 + SCRCH1*(S2_ZETA2*B22*B22/ZETA + M2))

		    RHO_S(I) = R/W(I)


		    ! KINEMATICAL TERMS

		    dchi_dZETA =  1.D0 - VEL2 - 0.5D0*W(I)*(R + 2.D0*chi*W(I))*dv2_dZETA
                    drho_dZETA = -0.5D0*R*W(I)*dv2_dZETA


		    ! ASSUMING IDEAL GAS
                    dp_dchi = G1G
                    dp_drho = 0.0D0
                    p_S(I)  = chi*dp_dchi
		    DP = dp_dchi*dchi_dZETA + dp_drho*drho_dZETA

		    FUN = ZETAP - (E+P_S(I)) + 0.5D0*(B22 + (B22*M2 - S2)*SCRCH2)
		    DFUN= 1.D0 - DP - (B22*M2 - S2)*SCRCH1*SCRCH2

		    OZETA = ZETAP			

		    ZETAP = ZETAP - FUN/DFUN

		    ERRP = ABS(1.D0-ZETAP/OZETA)

                    ITER = ITER + 1

                    ! =====================================================================
                    IF (eqs .EQ. 2) THEN

                        gam_b(I)  = gam_s(I)
                        xe_b(I)   = xe_s(I)
                        tem_bb(I) = tem_s(I)

                        ! Since ZETAP is updated in each interation, we need to update
                        ! the gamma factor using the Synge equation.

                        vx = (momenx_s(i) + bm*bx_s(i)/(zetap+R))/((zetap+R)+b22)
                        vy = (momeny_s(i) + bm*by_s(i)/(zetap+R))/((zetap+R)+b22)
                        vz = (momenz_s(i) + bm*bz_s(i)/(zetap+R))/((zetap+R)+b22)

                        LOR     = 1.D0/DSQRT(1.D0-vx**2-vy**2-vz**2)
                        rho_b(I) = R/LOR

                        eps_b(I) = ((zetap+R)/rho_b(I)/lor/lor - 1.D0 ) / gam_b(I)


                        ! ------------------------------------------------------------


                        call EOS_SYNGE (2, tid, I, I, me, mp, rho_b, xe_b, tem_bb, eps_b,   &
                                        p_b, gam_b, h_b, cs_b, dpdrho_b, dpdeps_b)

                        GAM_S(I) = gam_b(I)

                        G1G = (GAM_S(I) - 1.D0)/GAM_S(I)

                    ENDIF
                    ! =======================================================================


                 ENDDO

		 ZETA = ZETAP + R
	        

                    ! Recover primitive velocity.
                    velx_s(i) = (momenx_s(i) + bm*bx_s(i)/zeta)/(zeta+b22)
                    vely_s(i) = (momeny_s(i) + bm*by_s(i)/zeta)/(zeta+b22)
                    velz_s(i) = (momenz_s(i) + bm*bz_s(i)/zeta)/(zeta+b22)

	         
                    ! SMALL prevents the velocity to take unphysical small values.
                    IF (ABS(velx_s(i)) .LE. small) velx_s(i) = 0.D0
                    IF (ABS(vely_s(i)) .LE. small) vely_s(i) = 0.D0
                    IF (ABS(velz_s(i)) .LE. small) velz_s(i) = 0.D0
                
                    ! Recover Lorentz factor.
                    w(i) = 1.D0/DSQRT(1.D0 - velx_s(i)*velx_s(i) - vely_s(i)*vely_s(i) - velz_s(i)*velz_s(i))
        
                    ! Recover primitive density.
                    rho_s(i)  = denstw_s(i)/w(i)
                 
                    ! Calculate the internal energy eps.

                    IF (eqs .EQ. 2) THEN
                       eps_s(i) = (zeta/rho_s(i)/w(i)/w(i) - 1.D0) / gam_s(i)
                    ELSE
                       eps_s(i) = (zeta/rho_s(i)/w(i)/w(i) - 1.D0) / gamma_ad
                    ENDIF

              ENDDO

              ! Calculate the pressure p.

              IF (eqs .EQ. 2) THEN
                 call EOS_SYNGE (2, tid, 1, nx, me, mp, rho_s, xe_s, tem_s, eps_s,        &
                                 p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)

              ELSE
                 call EOS (2, 1, nx, gamma_ad, eps_s, rho_s, p_s, dpdrho_s, dpdeps_s, h_s, cs_s)
              ENDIF


	      DO i = 1, nx


	
	         IF (p_s(i) .LE. SMALLP) THEN

		    !IF (mpi_coords(2) .NE. 0) THEN
                    !WRITE(6,*) 'GETPRFQ_MHD_2D: MINIMUM PRESSURE REACHED AT POINT:'
                    !print*, i,j,k
                    !WRITE(6,*) '             PRESSURE RESET TO SMALLP'
                    !WRITE(6,*) '             SOLVING FOR V2'
		    !ENDIF                
		    !stop

                    p_s(i) = SMALLP


		    ! USE BISECTION ALGORITHM TO SOLVE FOR V2
		    ! ---------------------------------------

                    IF (eqs .EQ. 2) THEN
                       G1G = (gam_s(i) - 1.D0)/gam_s(i)
                    ELSE
                       G1G = (GAMMA_AD - 1.D0)/GAMMA_AD
                    ENDIF

                    B22 = bx_s(i)*bx_s(i)+by_s(i)*by_s(i) +bz_s(i)*bz_s(i)

                    M2 = momenx_s(i)*momenx_s(i) + momeny_s(i)*momeny_s(i) +                          &
                         momenz_s(i)*momenz_s(i)
                    BM = bx_s(i)*momenx_s(i) +  by_s(i)*momeny_s(i) +                                 &
                         bz_s(i)*momenz_s(i)
                    S2 = BM*BM


		    V2MIN = 0.D0
		    V2MAX = 1.D0-1.D-6

		    WMIN  = 1.D0/DSQRT(1.D0-V2MIN)
		    WMAX  = 1.D0/DSQRT(1.D0-V2MAX)

		    ZETA  = (DENSTW_S(I)+P_S(I)*WMIN/G1G)*WMIN
                    ZETA2 = ZETA*ZETA
                    S2_ZETA2 = S2/ZETA2

                    SCRCH1 = 1.D0/(ZETA + B22)
                    SCRCH2 = SCRCH1 * SCRCH1

		    FMIN  = S2_ZETA2*SCRCH1*(SCRCH1*ZETA + 1.D0) + M2*SCRCH2 - V2MIN


                    ZETA  = (DENSTW_S(I)+P_S(I)*WMAX/G1G)*WMAX
                    ZETA2 = ZETA*ZETA
                    S2_ZETA2 = S2/ZETA2

                    SCRCH1 = 1.D0/(ZETA + B22)
                    SCRCH2 = SCRCH1 * SCRCH1

                    FMAX  = S2_ZETA2*SCRCH1*(SCRCH1*ZETA + 1.D0) + M2*SCRCH2 - V2MAX

		
		    DO ll = 1, NITER

			VEL2 = 0.5D0*(V2MIN+V2MAX)
                        LOR    = 1.D0/DSQRT(1.D0-VEL2)

                        ZETA  = (DENSTW_S(I)+P_S(I)*LOR/G1G)*LOR

                    	ZETA2 = ZETA*ZETA
                    	S2_ZETA2 = S2/ZETA2

                    	SCRCH1 = 1.D0/(ZETA + B22)
                    	SCRCH2 = SCRCH1 * SCRCH1

                    	FUN = S2_ZETA2*SCRCH1*(SCRCH1*ZETA + 1.D0) + M2*SCRCH2 - VEL2

			IF (FUN*FMIN .GT. 0.D0) THEN
			
			    V2MIN = VEL2 
			    FMIN  = FUN

			ELSE

                            V2MAX = VEL2   
                            FMAX  = FUN
			
			ENDIF


			IF (ABS(FUN) .LT. 1.D-8) THEN

                        ! REDEFINE ENERGY

                        PHYSIC(I,J,K)%ENERGY = ZETA-DENSTW_S(I)-P_S(I)+0.5D0*(1.D0+VEL2)*B22-0.5D0*S2_ZETA2

                        ! RECOMPUTE DENSITY

                        RHO_S(I) = DENSTW_S(I)/LOR

                        tem_s(i) = p_s(i) / (rho_s(i)*   &
                                  (xe_s(i)/me + (1.D0-xe_s(i))/mp))  

		       
			IF (eqs .EQ. 1) THEN

                            call EOS (1, i, i, gamma_ad, eps_s, rho_s, p_s, dpdrho_s, dpdeps_s, h_s, cs_s)

		        ELSE

                            call EOS_SYNGE (1, tid, i, i, me, mp, rho_s, xe_s, tem_s, eps_s,        &
                                            p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)

                            G1G = (gam_s(i) - 1.D0)/gam_s(i)


			ENDIF


                        velx_s(i) = (momenx_s(i) + bm*bx_s(i)/zeta)/(zeta+b22)
                        vely_s(i) = (momeny_s(i) + bm*by_s(i)/zeta)/(zeta+b22)
                        velz_s(i) = (momenz_s(i) + bm*bz_s(i)/zeta)/(zeta+b22)

			EXIT

			ENDIF


			IF (VEL2 .GE. 1.D0 .OR. (ll+1) .GT. NITER) THEN
			
			    
                            tem_s(i) = p_s(i) / (rho_s(i)*   &
                                      (xe_s(i)/me + (1.D0-xe_s(i))/mp))


                            IF (eqs .EQ. 1) THEN

                                call EOS (1, i, i, gamma_ad, eps_s, rho_s, p_s, dpdrho_s, dpdeps_s, h_s, cs_s)

                            ELSE

                                call EOS_SYNGE (1, tid, i, i, me, mp, rho_s, xe_s, tem_s, eps_s,        &
                                                p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)


                            ENDIF


                            ! Lorentz Factor W
                    	    scrch11(i)   = 1.D0/DSQRT(1.D0-velx_s(i)**2-vely_s(i)**2-velz_s(i)**2)


                            ! Magnetic field four-vector.
                            b40 = scrch11(i)*(velx_s(i)*bx_s(i) +   &
                                              vely_s(i)*by_s(i) +   &
                                              velz_s(i)*bz_s(i))

                            b41 = bx_s(i)/scrch11(i) + velx_s(i)*b40
                            b42 = by_s(i)/scrch11(i) + vely_s(i)*b40
                            b43 = bz_s(i)/scrch11(i) + velz_s(i)*b40

                            ! Assuming the Minkowski metric to be (-,+,+,+)
                            b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                            scrch22(i)  = rho_s(i)*h_s(i) + b44

                            ! Conservative variables.
                            physic(i,j,k)%denstw  = rho_s(i)*scrch11(i)
                            physic(i,j,k)%denstew = rho_s(i)*scrch11(i)
                            physic(i,j,k)%momenx  = scrch22(i)*scrch11(i)*scrch11(i)*velx_s(i) - b40*b41
                            physic(i,j,k)%momeny  = scrch22(i)*scrch11(i)*scrch11(i)*vely_s(i) - b40*b42
                            physic(i,j,k)%momenz  = scrch22(i)*scrch11(i)*scrch11(i)*velz_s(i) - b40*b43
                            physic(i,j,k)%energy  = scrch22(i)*scrch11(i)*scrch11(i)-p_s(i) -           &
                                                    physic(i,j,k)%denstw - 0.5D0*b44 - b40*b40

			    EXIT


			ENDIF


			! END..............


		    ENDDO


                 END IF

	      ENDDO


              DO i = 1, nx      
              
                 physic(i,j,k)%pres = p_s(i)
                 physic(i,j,k)%densty = rho_s(i)
                 physic(i,j,k)%denstye = xe_s(i)*rho_s(i)
                 physic(i,j,k)%eps = eps_s(i)   
                 physic(i,j,k)%velx = velx_s(i)
                 physic(i,j,k)%vely = vely_s(i)
                 physic(i,j,k)%velz = velz_s(i)

              ENDDO
           ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL




        !------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE GETPRFQ3DNR
! ===============================================================================



