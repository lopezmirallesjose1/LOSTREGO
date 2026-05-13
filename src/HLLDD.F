!*********************************
!     3D-RMHD CODE: LOSTREGO
!*********************************

! NAME:
!       HLLD (Approximate Riemann Solver)

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
!       Calculate the fluxes in the direction of the sweep using an approximate
!       Riemann solver (HLL/HLLC/HLLD).

! ACADEMIC REFERENCE:
!       "An HLLC Riemann Solver for Relativistic Flows:I. Hydrodynamics."
!        A. Mignone & G. Bodo. 2005. MNRAS. 000, 1-13.
!       "An HLLC Riemann Solver for Relativistic Flows-II, Magnetohydrodynamics."
!        A. Mignone & G. Bodo. 2006. MNRAS. 368, 1040-1054
!       "The HLLD Riemann solver based on magnetic field decomposition method for
!        the numerical simulation of MHD". X. Guo, V. Florinski, C. Wang.
!        2016 . JCP 327 543-552.
!        See also: "Riemann Solvers and Numerical Methods for Fluid Dynamics". E.F.Toro. 
!        See also: "Magnetohidrodinamica relativista numerica:
!        Aplicaciones en Relatividad Especial y General". L.A. Ruiz. 2008


! ===================================================================================================
SUBROUTINE HLLDD (ii,tid, xyzswp,basenm,iars, fznp, nznp, rhol, rhor,                      &
                vell, velr, veltl, veltr, velttl, velttr, pl, pr, bl, br, btl, btr,        &
                bttl, bttr, b40l, b40r, b41l, b41r, b42l, fl, fr,                          &
                b42r, b43l, b43r, b44l, b44r, rl, rr, wl, wr, hl, hr, csl, csr,            &
                rflx, mflx, mtflx, mttflx, eflx, btflx, bttflx,gamma_ad,ishock)
! ===================================================================================================


        USE workarrays
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        integer(is):: ii,tid, iars, fznp, nznp, xyzswp

        character(4) :: basenm

        ! size -4:mn5
        real(rs)    :: pl(-4:*), pr(-4:*), rhol(-4:*), rhor(-4:*),     &
                       vell(-4:*), velr(-4:*), veltl(-4:*),            &
                       veltr(-4:*), velttl(-4:*), velttr(-4:*),        &
                       wl(-4:*), wr(-4:*), hl(-4:*), hr(-4:*),         &
                       csl(-4:*), csr(-4:*),                           &
                       rl(-4:*), rr(-4:*), bl(-4:*), br(-4:*),         &
                       btl(-4:*), btr(-4:*), bttl(-4:*),               &
                       bttr(-4:*), b40l(-4:*), b40r(-4:*),             &
                       b41l(-4:*), b41r(-4:*), b42l(-4:*),             &
                       b42r(-4:*), b43l(-4:*), b43r(-4:*),             &
                       b44l(-4:*), b44r(-4:*), fl(-4:*), fr(-4:*)


        real(rs) :: gamma_ad               

	integer(is) :: ishock (-4:*)
                       
        !------------------------------------------------------- Output Variables

        ! size -4:mn6
        real(rs)   :: rflx(-4:*), mflx(-4:*), mtflx(-4:*),               &
                      mttflx(-4:*), eflx(-4:*), btflx(-4:*),             &
                      bttflx(-4:*)


       !------------------------------------------------------- Local Variables

        integer(is):: i, niter, ITER, IERR, iter2, EXITTRUE, total

	integer(is), pointer :: HLLswitch(:)
	!integer(is) :: HLLswitch(NZNP+1)


        real(rs)   :: tt2, tt4,                           &
                      ca0, ca1, ca2, ca3, ca4,            &
                      lamb1l, lamb1r, lamb4l, lamb4r, calf2, omega2, aux1, aux2

        real(rs)   :: AH, BH, CH, bvt, bvtt

        real(rs), pointer :: p(:), vel(:), velt(:),                      &
                             veltt(:), w(:), r(:), m(:),                 &
                             mt(:), mtt(:), e(:), b(:), bt(:), btt(:),   &
                             rls(:), rrs(:),                             &
                             btls(:), btrs(:), bttls(:), bttrs(:),       &
                             ml(:), mtl(:), mttl(:), mr(:), mtr(:),      &
                             mttr(:), el(:), er(:), frl(:), fml(:),      &
                             fmtl(:), fmttl(:), fel(:),                  &
                             fbtl(:), fbttl(:), frr(:), fmr(:),          &
                             fmtr(:), fmttr(:), fer(:),                  &
                             fbtr(:), fbttr(:), lbm(:), lbn(:)

        real(rs) :: PB0, B22, BM, M2, ZETA, velh, velth, veltth, rho, eps, phll,            &
                    DELTA0, RELERR, RELERR0, OLDP, F1, F0, A_, G_, C_, Q_, X_,              &
                    VEL_ALF_L, VEL_ALF_TL, VEL_ALF_TTL, VEL_ALF_R, VEL_ALF_TR, VEL_ALF_TTR, &
                    B_ALF_TL, B_ALF_TTL, B_ALF_TR, B_ALF_TTR, WWL, WWR, ETA_L, ETA_R,       &
                    KL, KTL, KTTL, KR, KTR, KTTR, LBN_ALF, LBM_ALF, scrch1,                 &                    
                    b40, b41, b42, b43, b44, B_C, B_C_T, B_C_TT, YL, YR, b40h, b41h,        &
                    b42h, b43h, b44h, VEL_ALF, LBA, dK

                    
        real(rs), pointer :: RD_L(:), RM_L(:), RMT_L(:), RMTT_L(:), RE_L(:), &
                             RB_L(:), RBT_L(:), RBTT_L(:)

        real(rs), pointer :: RD_R(:), RM_R(:), RMT_R(:), RMTT_R(:), RE_R(:), &
                             RB_R(:), RBT_R(:), RBTT_R(:)

        real(rs), pointer :: lbna(:), lbma(:), lbc(:)                     


        real(rs), pointer :: frl2(:), fml2(:),                             &
                             fmtl2(:), fmttl2(:), fel2(:),                 &
                             fbtl2(:), fbttl2(:), frr2(:), fmr2(:),        &
                             fmtr2(:), fmttr2(:), fer2(:),                 &
                             fbtr2(:), fbttr2(:)


        real(rs) :: RL_A, RR_A, EL_A, ER_A, ML_A, MR_A, MTL_A, MTR_A, MTTL_A, MTTR_A
        real(rs) :: RL_C, RR_C, EL_C, ER_C, ML_C, MR_C, MTL_C, MTR_C, MTTL_C, MTTR_C    
	real(rs) :: R_A, E_A, M_A, MT_A, MTT_A, R_C, E_C, M_C, MT_C, MTT_C    

	real(rs) :: VEL_C_L, VEL_C_R, VEL_C_TL, VEL_C_TR, VEL_C_TTL, VEL_C_TTR                             
	real(rs) :: KLBC, KRBC, VEL_C1, dP, vB, vBC, SCRCHL, SCRCHR, VEL_C, VEL_C_T, VEL_C_TT        

        real(rs)   :: S2, ERRP, G1G, P_S, signB,                     &
                      SCRCH2, FUN, DFUN, OZETA, B1, C1, CHI, ZETAP,  &
                      ZETA2, S2_ZETA2, VEL2, W2, DV2_DZETA, DCHI_DZETA, DRHO_DZETA, DP_DCHI, DP_DRHO

        !-------------------------------------------------------------------------- BEGIN

        ! size -4:mn5
        p    => wa(tid)%wa5hlle100(:)
        vel  => wa(tid)%wa5hlle101(:)
        velt => wa(tid)%wa5hlle102(:)
        veltt=> wa(tid)%wa5hlle103(:)
        w    => wa(tid)%wa5hlle104(:)

        r   => wa(tid)%wa5hlle105(:)
        m   => wa(tid)%wa5hlle106(:)
        mt  => wa(tid)%wa5hlle107(:)
        mtt => wa(tid)%wa5hlle108(:)
        e   => wa(tid)%wa5hlle109(:)
        b   => wa(tid)%wa5hlle110(:)
        bt  => wa(tid)%wa5hlle111(:)
        btt => wa(tid)%wa5hlle112(:)

        ml  => wa(tid)%wa5hlle113(:)
        mtl => wa(tid)%wa5hlle114(:)
        mttl=> wa(tid)%wa5hlle115(:)
        mr  => wa(tid)%wa5hlle116(:)
        mtr => wa(tid)%wa5hlle117(:)
        mttr=> wa(tid)%wa5hlle118(:)
        el  => wa(tid)%wa5hlle119(:)
        er  => wa(tid)%wa5hlle120(:)

        rls => wa(tid)%wa5hlle121(:)
        rrs => wa(tid)%wa5hlle122(:)

        btls  =>  wa(tid)%wa5hlle123(:)
        btrs  =>  wa(tid)%wa5hlle124(:)
        bttls =>  wa(tid)%wa5hlle125(:)
        bttrs =>  wa(tid)%wa5hlle126(:)

        frl  => wa(tid)%wa5hlle127(:)
        fml  => wa(tid)%wa5hlle128(:)
        fmtl => wa(tid)%wa5hlle129(:)
        fmttl=> wa(tid)%wa5hlle130(:)
        fel  => wa(tid)%wa5hlle131(:)
        fbtl => wa(tid)%wa5hlle132(:)
        fbttl=> wa(tid)%wa5hlle133(:)

        frr  => wa(tid)%wa5hlle134(:)
        fmr  => wa(tid)%wa5hlle135(:)
        fmtr => wa(tid)%wa5hlle136(:)
        fmttr=> wa(tid)%wa5hlle137(:)
        fer  => wa(tid)%wa5hlle138(:)
        fbtr => wa(tid)%wa5hlle139(:)
        fbttr=> wa(tid)%wa5hlle140(:)

        lbm => wa(tid)%wa5hlle141(:)
        lbn => wa(tid)%wa5hlle142(:)


        RD_L      =>wa(tid)%wa5hlld100(:)
        RM_L      =>wa(tid)%wa5hlld101(:)
        RMT_L     =>wa(tid)%wa5hlld102(:)
        RMTT_L    =>wa(tid)%wa5hlld103(:)
        RE_L      =>wa(tid)%wa5hlld104(:)
        RB_L      =>wa(tid)%wa5hlld105(:)
        RBT_L     =>wa(tid)%wa5hlld106(:)
        RBTT_L    =>wa(tid)%wa5hlld107(:)

        RD_R      =>wa(tid)%wa5hlld108(:)
        RM_R      =>wa(tid)%wa5hlld109(:)
        RMT_R     =>wa(tid)%wa5hlld110(:)
        RMTT_R    =>wa(tid)%wa5hlld111(:)
        RE_R      =>wa(tid)%wa5hlld112(:)
        RB_R      =>wa(tid)%wa5hlld113(:)
        RBT_R     =>wa(tid)%wa5hlld114(:)
        RBTT_R    =>wa(tid)%wa5hlld115(:)

        lbna      =>wa(tid)%wa5hlld116(:)
        lbma      =>wa(tid)%wa5hlld117(:)
        lbc       =>wa(tid)%wa5hlld118(:)

        
        frl2  => wa(tid)%wa5hlld119(:)
        fml2  => wa(tid)%wa5hlld120(:)
        fmtl2 => wa(tid)%wa5hlld121(:)
        fmttl2=> wa(tid)%wa5hlld122(:)
        fel2  => wa(tid)%wa5hlld123(:)
        fbtl2 => wa(tid)%wa5hlld124(:)
        fbttl2=> wa(tid)%wa5hlld125(:)

        frr2  => wa(tid)%wa5hlld126(:)
        fmr2  => wa(tid)%wa5hlld127(:)
        fmtr2 => wa(tid)%wa5hlld128(:)
        fmttr2=> wa(tid)%wa5hlld129(:)
        fer2  => wa(tid)%wa5hlld130(:)
        fbtr2 => wa(tid)%wa5hlld131(:)
        fbttr2=> wa(tid)%wa5hlld132(:)

	HLLswitch =>wa(tid)%wa5hlld133(:)
        
        DO I=FZNP,NZNP+1
	    				
	    CALF2  = B44L(I)/(RHOL(I)*HL(I) + B44L(I))
            OMEGA2 = CSL(I)**2*(1.D0 - CALF2) + CALF2 
            AUX1   = OMEGA2/wL(I)**2
            AUX2   = CSL(I)**2/(RHOL(I)*HL(I) + B44L(I))/wL(I)**4
		
            CA0    = VELL(I)**4*(1.D0 - OMEGA2) - VELL(I)**2*AUX1 + B41L(I)**2*AUX2
            CA1    = -4.D0*VELL(I)**3*(1.D0 - OMEGA2) + 2.D0*VELL(I)*AUX1 -          &
                     2.D0*B40L(I)*B41L(I)*AUX2		
            CA2    = 6.D0*VELL(I)**2*(1.D0 - OMEGA2) - (1.D0 - VELL(I)**2)*AUX1 +    &
                    (B40L(I)**2 - B41L(I)**2)*AUX2	
	    CA3    = -4.D0*VELL(I)*(1.D0 - OMEGA2) - 2*VELL(I)*AUX1 +                & 
      	              2.D0*B40L(I)*B41L(I)*AUX2
	    CA4    = 1.D0 - OMEGA2*(VELL(I)**2 + VELTL(I)**2 + VELTTL(I)**2) -       &
                     B40L(I)**2*AUX2	

            call FAST (ca0,ca1,ca2,ca3,ca4,lamb1l,lamb4l)
		 
				
	    CALF2  = B44R(I)/(RHOR(I)*HR(I) + B44R(I))
	    OMEGA2 = CSR(I)**2*(1.D0 - CALF2) + CALF2
	    AUX1   = OMEGA2/WR(I)**2
	    AUX2   = CSR(I)**2/(RHOR(I)*HR(I) + B44R(I))/WR(I)**4
		
            CA0    = VELR(I)**4*(1.D0 - OMEGA2) - VELR(I)**2*AUX1 + B41R(I)**2*AUX2
            CA1    = -4.D0*VELR(I)**3*(1.D0 - OMEGA2) + 2.D0*VELR(I)*AUX1 -          &
                     2.D0*B40R(I)*B41R(I)*AUX2		
	    CA2    = 6.D0*VELR(I)**2*(1.D0 - OMEGA2) - (1.D0 - VELR(I)**2)*AUX1 +    &
                    (B40R(I)**2 - B41R(I)**2)*AUX2	
	    CA3    = -4.D0*VELR(I)*(1.D0 - OMEGA2) - 2*VELR(I)*AUX1 +                & 
      	             2.D0*B40R(I)*B41R(I)*AUX2
	    CA4    = 1.D0 - OMEGA2*(VELR(I)**2 + VELTR(I)**2 + VELTTR(I)**2) -       &
                     B40R(I)**2*AUX2	

            call FAST (ca0,ca1,ca2,ca3,ca4,lamb1r,lamb4r)


            lbm(i) = MAX (lamb4l, lamb4r)
            lbn(i) = MIN (lamb1l, lamb1r)
            
            
            !!! Approximate Riemann Solvers employed:

            !   iars = 2 ----> HLLE approximate Riemann Solver.
            !   iars = 3 ----> HLLC approximate Riemann Solver.
            !   iars = 4 ----> HLLD approximate Riemann Solver.


            !  -----------------------------------------------------------------
            !  APPROXIMATE RIEMANN SOLVER: HLL (iars = 2).
            !  -----------------------------------------------------------------

            !                      t
            !     lbn*T              ^               lbm*T
            !-------\----------------|----------------/--------
            !         \              |              /
            !           \            |            /
            !             \          |          /
            !     UL        \       Uhll      /        UR
            !                 \      |      /
            !                   \    |    /
            !                     \  |  /
            !------------------------*------------------------->
            !                        O                         x

            ! The fluxes are calculated in the direction of the sweep.
            ! The corresponding interface numerical flux at t=T is defined as:


            ! FL = F (UL)    if   lbn .GE. 0             (1)
            ! Fhll           if   lbn .LE. 0 .LE. lbm    (2)
            ! FR = F (UR)    if   lbm .LE. 0             (3)


            ! IMPORTANT TO NOTE: 
            ! The magnetic field quadrivector (b40, b41, b42, b43) has been previously
            ! defined in terms of the field components (b, bt, btt).
            

            !	 -----------
            !	 FLUXES IN THE HLL INTERMEDIATE STATE
            !	 -----------

	    RFLX(I)  = (lbm(I)*RHOL(I)*velL(I)*wl(i) - lbn(I)*RHOR(I)*velr(i)*wr(I) +   &
                        lbm(I)*lbn(I)*(RHOR(I)*wR(I) - RHOL(I)*wL(I)))/                 &
     	               (lbm(I) - lbn(I))
	  
       	    MFLX(I) = (lbm(I)*((RHOL(I)*HL(I) + B44L(I))*vell(i)*wl(i)*vell(i)*wl(i) +   &
                       PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I)) -                        &
                       lbn(I)*((RHOR(I)*HR(I) + B44R(I))*wr(i)*velr(i)*wr(i)*velr(i) +   &
     	               PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I)) +                        &
                       lbm(I)*lbn(I)*                                                    &
                     ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*velr(i) -                    & 
                       B40R(I)*B41R(I) -                                                 &
                      (RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*vell(i) +                    &
                       B40L(I)*B41L(I)))/(lbm(I) - lbn(I))

	    MTFLX(I) = (lbm(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*vell(i)*wl(I)*veltl(i) -   &
     	                B41L(I)*B42L(I)) -                                                 &
                        lbn(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*velr(i)*WR(I)*veltr(i)  -  &
     		        B41R(I)*B42R(I)) +                                                 &
    		        lbm(I)*lbn(I)*                                                     &
                      ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*veltr(i) -                    & 
                        B40R(I)*B42R(I) -                                                  &
                       (RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*veltl(i) +                    &
                        B40L(I)*B42L(I)))/(lbm(I) - lbn(I))

            MTTFLX(I) = (lbm(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*vell(i)*wl(i)*velttl(I) - &
     	                 B41L(I)*B43L(I)) -                                                &
                         lbn(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*velr(i)*wr(i)*velttr(I) - &
     	          	 B41R(I)*B43R(I)) +                                                &
      		  	 lbm(I)*lbn(I)*                                                    &
                 	 ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(i)*velttr(I) -                & 
                         B40R(I)*B43R(I) -                                                 &
                         (RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(i)*velttL(I) +                 &
                 	 B40L(I)*B43L(I)))/(lbm(I) - lbn(I))
 
	    EFLX(I)  = (lbm(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*vell(i) -            &
     	                B40L(I)*B41L(I)) -                         			   &
                        lbn(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*velr(i) -  	   &
     		  	B40R(I)*B41R(I)) + 						   &
     		  	lbm(I)*lbn(I)* 							   &
                 	((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I) - 			   & 
                 	PR(I) - 0.5D0*B44R(I) - B40R(I)*B40R(I) - 			   &
                 	(RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I) +  			   &
                 	PL(I) + 0.5D0*B44L(I) + B40L(I)*B40L(I)))/ 			   &
                 	(lbm(I) - lbn(I)) - RFLX(I)


            !FBX(I) = 0.D0
	  
	    BTFLX(I) = (lbm(I)*(BtL(I)*VELL(I) - BL(I)*VELtL(I)) -   &
                        lbn(I)*(BtR(I)*VELR(I) - BR(I)*VELtR(I)) +   &
                        lbm(I)*lbn(I)*(BtR(I) - BtL(I)))/            &
                       (lbm(I) - lbn(I)) 	 

	    BTTFLX(I) = (lbm(I)*(BttL(I)*VELL(I) - BL(I)*VELttL(I)) -   &
                         lbn(I)*(BttR(I)*VELR(I) - BR(I)*VELttR(I)) +   &
                         lbm(I)*lbn(I)*(BttR(I) - BttL(I)))/            &
                        (lbm(I) - lbn(I)) 	 


        ENDDO


        !IF (iars .EQ. 2 .OR. (iars .EQ. 4 .AND. (BASENM .EQ. 'CMBW'    .OR.   & 
        !                                         BASENM(1:2) .EQ. 'MH' .OR.   &
	!					 BASENM .EQ. 'LOOP' ))) THEN

        IF (iars .EQ. 2) THEN


           DO I = FZNP, NZNP+1

              IF (lbn(I).GT.0.D0) THEN 

                 RFLX(I)  = RHOL(I)*WL(I)*VELL(I)

                 mflx(I) = (RHOL(I)*HL(I) + B44L(I))*wL(I)*vell(i)*WL(I)*VELL(I) +     &
                            PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I)
                 mtflx(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) -   &
     	                     B41L(I)*B42L(I)
                 mttflx(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) - &
     	                      B41L(I)*B43L(I)

                 eflx(I)  = (RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) -            &
     	                     B40L(I)*B41L(I) - RFLX(I)

                 !FBX(I) = 0.D0
                 btflx(I) = BTL(I)*VELL(I) - BL(I)*VELTL(I)
                 bttflx(I) = BTTL(I)*VELL(I) - BL(I)*VELTTL(I)

              END IF


              IF (LBM(I).LE.0.D0) THEN

                 RFLX(I)  = RHOR(I)*WR(I)*VELR(I)

                 MFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELR(I) +     &
                            PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I)
                 MTFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) -   &
     	                     B41R(I)*B42R(I)
                 MTTFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) - &
     	                      B41R(I)*B43R(I)

                 EFLX(I)  = (RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) -            &
     	                     B40R(I)*B41R(I) - RFLX(I)

                 !FBX(I) = 0.D0
                 BTFLX(I) = BTR(I)*VELR(I) - BR(I)*VELTR(I)
                 BTTFLX(I) = BTTR(I)*VELR(I) - BR(I)*VELTTR(I)

              END IF

           ENDDO

           IF (IARS .EQ. 2) RETURN

        END IF
      

        !     ------
        !     HLLC CONTINUES
        !     ------


        !  -----------------------------------------------------------------
        !  APPROXIMATE RIEMANN SOLVER: HLLC (iars = 3).
        !  -----------------------------------------------------------------

        !                      t
        !     lbn*T              ^     lb'*T     lbm*T
        !-------\----------------|-------.--------/--------
        !         \              |      .       /
        !           \            |     .      /
        !             \          |    .     /
        !     UL        \      UL*    UR* /        UR
        !                 \      |  .   /
        !                   \    | .  /
        !                     \  |. /
        !------------------------*------------------------->
        !                        O                         x


        ! The solution is approximated to the initial problem with two constant states
        ! bounded by two fast shocks and a contact discontinuity in the middle (with velocity lb').

        ! The corresponding interface numerical flux at t=T is defined as:

        ! FL = F (UL)    if   lbn .GE. 0             (1)
        ! FL'            if   lbn .LE. 0 .LE. lb'    (2)
        ! FR'            if   lb' .LE. 0 .LE. lbm    (3)
        ! FR = F (UR)    if   lbm .LE. 0

        ! IMPORTANT TO NOTE: 
        ! The magnetic field quadrivector (b40, b41, b42, b43) has been previously
        ! defined in terms of the field components (b, bt, btt).

      
        IF (IARS .EQ. 3) THEN

      
           DO I = FZNP, NZNP+1



              R(I)  = (lbm(I)*RHOR(I)*WR(I) - lbn(I)*RHOL(I)*WL(I) +                             &
              	        RHOL(I)*WL(I)*VELL(I) - RHOR(I)*WR(I)*VELR(I))/(lbm(I) - lbn(I))

              M(I) = (lbm(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*velr(i) - B40R(I)*B41R(I)) - & 
               	      lbn(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*vell(i) - B40L(I)*B41L(I)) + &
                     ((RHOL(I)*HL(I) + B44L(I))*WL(I)*vell(i)*WL(I)*vell(i) +                    &
                       PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I)) -                                &
                     ((RHOR(I)*HR(I) + B44R(I))*WR(I)*velr(i)*WR(I)*velr(i) +                    &
                       PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I)))/(lbm(I) - lbn(I))

              MT(I) = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTR(I) -                     & 
                       B40R(I)*B42R(I)) - LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTL(I) -  & 
                       B40L(I)*B42L(I)) + ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) - &
     	               B41L(I)*B42L(I)) - ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) - &
     	               B41R(I)*B42R(I)))/(LBM(I) - LBN(I))

              MTT(I) = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTTR(I) -                     & 
                        B40R(I)*B43R(I)) - LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTTL(I) -  & 
                        B40L(I)*B43L(I)) + ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) - &
     	                B41L(I)*B43L(I)) - ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) - &
     	                B41R(I)*B43R(I)))/(LBM(I) - LBN(I))

              E(I)  = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I) -                    & 
                       PR(I) - 0.5D0*B44R(I) - B40R(I)*B40R(I)) -                         &
                       LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I) -                    & 
                       PL(I) - 0.5D0*B44L(I) - B40L(I)*B40L(I)) +                         &
                      ((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) - B40L(I)*B41L(I)) - & 
                      ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) - B40R(I)*B41R(I)))/ &
                      (LBM(I) - LBN(I)) - R(I)

              B(I) = (BL(I)+BR(I))/2.D0

              BT(I) = (LBM(I)*BTR(I) - LBN(I)*BTL(I) + (BTL(I)*VELL(I) - BL(I)*VELTL(I)) -       &
                      (BtR(I)*VELR(I) - BR(I)*VELTR(I)))/(LBM(I) - LBN(I))

              BTT(I) = (LBM(I)*BTTR(I) - LBN(I)*BTTL(I) + (BTTL(I)*VELL(I) - BL(I)*VELTTL(I)) -  &
                       (BTTR(I)*VELR(I) - BR(I)*VELTTR(I)))/(LBM(I) - LBN(I))


              !	    ------
	      !       GETTING VELX FROM THE QUADRATIC EQUATION 
	      !       AND THE REST OF VARIABLES AT THE HLLC INTERMEDIATE STATE
	      !     ------
       

              IF (B(I).EQ.0.D0) THEN

                 AH = EFLX(I) + RFLX(I)

           	 BH = - MFLX(I) - (E(I) + R(I))

           	 CH = M(I)

           	 IF (AH .NE. 0.D0) VEL(I) = (- BH - DSQRT(BH**2 - 4.D0*AH*CH))/2.D0/AH
	   	 IF (AH .EQ. 0.D0) VEL(I) = -CH/BH


           	 P(I)    = MFLX(I) - VEL(I)*(EFLX(I) + RFLX(I))

    	         RLs(I)   = (lbn(I) - VELL(I))*RHOL(I)*WL(I)/(LBN(I) - VEL(I))

                 MTL(I)  = (LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTL(I) - B40L(I)*B42L(I)) - &
                           ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) - B41L(I)*B42L(I)))/ &
                           (LBN(I) - VEL(I))

                 MTTL(I)  = (LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTTL(I) - B40L(I)*B43L(I)) - &
                            ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) - B41L(I)*B43L(I)))/ &
                            (LBN(I) - VEL(I))

                 EL(I)   = (LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I) -                     &
                            PL(I) - 0.5D0*B44L(I) - B40L(I)*B40L(I)) -                          &
                            ((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) - B40L(I)*B41L(I)) + & 
                            P(I)*VEL(I))/(LBN(I) - VEL(I))

                 ML(I)  = (EL(I) + P(I))*VEL(I)

           	 BTLS(I) = (LBN(I)-VELL(I))*BTL(I)/(LBN(I) - VEL(I))

           	 BTTLS(I) = (LBN(I)-VELL(I))*BTTL(I)/(LBN(I) - VEL(I))

           	 RRS(I)   = (LBM(I) - VELR(I))*RHOR(I)*WR(I)/(LBM(I) - VEL(I))

           	 MTR(I)  = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTR(I) - B40R(I)*B42R(I)) - &
                           ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) - B41R(I)*B42R(I)))/ &
                 	    (LBM(I) - VEL(I))

           	 MTTR(I)  = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTTR(I) - B40R(I)*B43R(I)) -  &
                             ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) - B41R(I)*B43R(I)))/ &
                  	     (LBM(I) - VEL(I))

           	 ER(I)   = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I) -                    & 
                 	    PR(I) - 0.5D0*B44R(I) - B40R(I)*B40R(I)) -                         & 
                    	   ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) - B40R(I)*B41R(I)) + & 
                    	    P(I)*VEL(I))/(LBM(I) - VEL(I))

           	 MR(I)  = (ER(I) + P(I))*VEL(I)

           	 BTRS(I) = (LBM(I)-VELR(I))*BTR(I)/(LBM(I) - VEL(I))

           	 BTTRS(I) = (LBM(I)-VELR(I))*BTTR(I)/(LBM(I) - VEL(I))


	         !          -------------
	         !          FLUXES AT THE HLLC INTERMEDIATE STATES
	         !          -------------

	         FRL(I)  = RLS(I)*VEL(I)

           	 FML(I) = ML(I)*VEL(I) + P(I)

           	 FMTL(I) = MTL(I)*VEL(I) 

           	 FMTTL(I) = MTTL(I)*VEL(I) 

           	 FEL(I)  = ML(I)

           	 FBTL(I) = BTLS(I)*VEL(I)

           	 FBTTL(I) = BTTLS(I)*VEL(I) 

            	 FRR(I)  = RRS(I)*VEL(I)

           	 FMR(I) = MR(I)*VEL(I) + P(I)

           	 FMTR(I) = MTR(I)*VEL(I) 

           	 FMTTR(I) = MTTR(I)*VEL(I)

           	 FER(I)  = MR(I)

           	 FBTR(I) = BTRS(I)*VEL(I) 

           	 FBTTR(I) = BTTRS(I)*VEL(I) 

              ELSE

                 AH = EFLX(I) + RFLX(I) - BT(I)*BTFLX(I) - BTT(I)*BTTFLX(I)

                 BH = - MFLX(I) - (E(I) + R(I)) + (BT(I)*BT(I) + BTT(I)*BTT(I)) + & 
                      (BTFLX(I)*BTFLX(I) + BTTFLX(I)*BTTFLX(I))

                 CH = M(I) - BT(I)*BTFLX(I) - BTT(I)*BTTFLX(I)

	   	 IF (AH .NE. 0.D0)  VEL(I) = (- BH - DSQRT(BH**2 - 4.D0*AH*CH))/2.D0/AH
	   	 IF (AH .EQ. 0.D0)  VEL(I) = -CH/BH

           	 BVT = BT(I)*VEL(I) - BTFLX(I)
	   	 BVTT = BTT(I)*VEL(I) - BTTFLX(I)


                 IF (BASENM.EQ.'CMBW' .AND. XYZSWP.EQ.1) THEN 

                    VELT(I) = 0.D0           !NOT NEEDED (ALWAYS APPEARS MULTIPLIED BY ZERO)
              	    VELTT(I) = 0.D0
              	    BTT(I)   = 0.D0
              	    BVTT = 0.D0

                 ELSE IF (BASENM.EQ.'CMBW' .AND. XYZSWP.EQ.2) THEN

                    VELTT(I) = 0.D0           !NOT NEEDED (ALWAYS APPEARS MULTIPLIED BY ZERO) 
                    VELT(I) = 0.D0
                    BT(I)   = 0.D0
                    BVT = 0.D0

	         ELSE
 
                    VELT(I) = (BT(I)*VEL(I) - BTFLX(I))/B(I)
                    VELTT(I) = (BTT(I)*VEL(I) - BTTFLX(I))/B(I)

                 END IF

                 P(I)    = MFLX(I) + B(I)**2 - BVT**2 - BVTT**2 -             &
                           VEL(I)*(EFLX(I) + RFLX(I) - BVT*BT(I) - BVTT*BTT(I))

                 RLS(I)   = (LBN(I) - VELL(I))*RHOL(I)*WL(I)/(LBN(I) - VEL(I))

                 MTL(I)  = (- B(I)*BT(I)*(1.D0 - VEL(I)**2) + BVTT*BT(I)*VELTT(I) -                     &
                            BVT*(B(I)*VEL(I) + BTT(I)*VELTT(I)) +                                       &
                    	    LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTL(I) - B40L(I)*B42L(I)) - &
                           ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) - B41L(I)*B42L(I)))/ &
                    	   (LBN(I) - VEL(I))

           	 MTTL(I)  = (- B(I)*BTT(I)*(1.D0 - VEL(I)**2) + BVT*BTT(I)*VELT(I) -                       &
                             BVTT*(B(I)*VEL(I) + BT(I)*VELT(I)) +                                          &
                             LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTTL(I) - B40L(I)*B43L(I)) -  &
                             ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) - B41L(I)*B43L(I)))/ &
                             (LBN(I) - VEL(I))

                 EL(I)   = (LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I) -                    & 
                            PL(I) - 0.5D0*B44L(I) - B40L(I)*B40L(I)) -                         &
                           ((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) - B40L(I)*B41L(I)) + &
                            P(I)*VEL(I) - B(I)**2*VEL(I) - BVT*BT(I) - BVTT*BTT(I))/           &
                           (LBN(I) - VEL(I))

                 ML(I)  = (EL(I) + P(I))*VEL(I) - B(I)**2*VEL(I) - BVT*BT(I) - BVTT*BTT(I)

                 RRS(I)   = (LBM(I) - VELR(I))*RHOR(I)*WR(I)/(LBM(I) - VEL(I))

                 MTR(I)  = (- B(I)*BT(I)*(1.D0 - VEL(I)**2) + BVTT*BT(I)*VELTT(I) -                     &
                           BVT*(B(I)*VEL(I) + BTT(I)*VELTT(I)) +                                        &
                           LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTR(I) - B40R(I)*B42R(I)) -  &
                           ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) - B41R(I)*B42R(I)))/ &
                           (LBM(I) - VEL(I))

                 MTTR(I)  = (- B(I)*BTT(I)*(1.D0 - VEL(I)**2) + BVT*BTT(I)*VELT(I) -                      &
                            BVTT*(B(I)*VEL(I) + BT(I)*VELT(I)) +                                          &
                    	    LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTTR(I) - B40R(I)*B43R(I)) -  &
                    	    ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) - B41R(I)*B43R(I)))/ &
                    	    (LBM(I) - VEL(I))

           	 ER(I)   = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I) -                               &
                           PR(I) - 0.5D0*B44R(I) - B40R(I)*B40R(I)) -                                     &
                           ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) - B40R(I)*B41R(I)) +            &  
                           P(I)*VEL(I) - B(I)**2*VEL(I) - BVT*BT(I) - BVTT*BTT(I))/                       &
                           (LBM(I) - VEL(I))

           	 MR(I)  = (ER(I) + P(I))*VEL(I) - B(I)**2*VEL(I) - BVT*BT(I) - BVTT*BTT(I)


                 !          -------------
                 !          FLUXES AT THE HLLC INTERMEDIATE STATES
                 !          -------------

           	 FRL(I)  = RLS(I)*VEL(I)

           	 FML(I) = ML(I)*VEL(I) - B(I)**2 + BVT**2 + BVTT**2 - & 
                          BVT*BT(I)*VEL(I) - BVTT*BTT(I)*VEL(I) + P(I)

                 FMTL(I) = MTL(I)*VEL(I) - B(I)*BT(I)*(1.D0-VEL(I)**2) + &
                           BVTT*BT(I)*VELTT(I) -                         &
                           BVT*(B(I)*VEL(I) + BTT(I)*VELTT(I))

                 FMTTL(I) = MTTL(I)*VEL(I) - B(I)*BTT(I)*(1.D0-VEL(I)**2) + &
                            BVT*BTT(I)*VELT(I) -                            &
                            BVTT*(B(I)*VEL(I) + BT(I)*VELT(I))

           	 FEL(I)  = ML(I)

           	 FBTL(I) = BT(I)*VEL(I) - BVT

           	 FBTTL(I) = BTT(I)*VEL(I) - BVTT

           	 FRR(I)  = RRS(I)*VEL(I)

           	 FMR(I) = MR(I)*VEL(I) - B(I)**2 + BVT**2 + BVTT**2 - & 
                          BVT*BT(I)*VEL(I) - BVTT*BTT(I)*VEL(I) + P(I)

           	 FMTR(I) = MTR(I)*VEL(I) - B(I)*BT(I)*(1.D0-VEL(I)**2) + &
                           BVTT*BT(I)*VELTT(I) -                         &
                    	   BVT*(B(I)*VEL(I) + BTT(I)*VELTT(I))

           	 FMTTR(I) = MTTR(I)*VEL(I) - B(I)*BTT(I)*(1.D0-VEL(I)**2) + &
                            BVT*BTT(I)*VELT(I) -                            &
                    	    BVTT*(B(I)*VEL(I) + BT(I)*VELT(I))

                 FER(I)  = MR(I)

                 FBTR(I) = BT(I)*VEL(I) - BVT

                 FBTTR(I) = BTT(I)*VEL(I) - BVTT

              END IF

           ENDDO

           !     --------
	   !     HLLC INTERCELL NUMERICAL FLUXES 
	   !     --------

           DO I=FZNP,NZNP+1


              IF (LBN(I).GT.0.D0) THEN 

                 RFLX(I)  = RHOL(I)*WL(I)*VELL(I)

                 MFLX(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELL(I) + &
                            PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I)
                 MTFLX(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) - &
     	                     B41L(I)*B42L(I)
            	 MTTFLX(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) - &
     	                      B41L(I)*B43L(I)

            	 EFLX(I)  = (RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) - &
     	                     B40L(I)*B41L(I) - RFLX(I)

            	 BTFLX(I) = BTL(I)*VELL(I) - BL(I)*VELTL(I)
            	 BTTFLX(I) = BTTL(I)*VELL(I) - BL(I)*VELTTL(I)

              ELSE IF (LBN(I).LE.0.D0.AND.VEL(I).GT.0.D0) THEN

                 RFLX(I)  = FRL(I)

            	 MFLX(I) = FML(I) 
            	 MTFLX(I) = FMTL(I) 
            	 MTTFLX(I) = FMTTL(I) 

            	 EFLX(I)  = FEL(I) - RFLX(I)

            	 BTFLX(I) = FBTL(I)
            	 BTTFLX(I) = FBTTL(I)

              ELSE IF (VEL(I).LE.0.D0.AND.LBM(I).GT.0.D0) THEN

            	 RFLX(I)  = FRR(I)

            	 MFLX(I) = FMR(I)
            	 MTFLX(I) = FMTR(I) 
            	 MTTFLX(I) = FMTTR(I) 

            	 EFLX(I)  = FER(I) - RFLX(I)

            	 BTFLX(I) = FBTR(I)
            	 BTTFLX(I) = FBTTR(I)

              ELSE 

                 RFLX(I)  = RHOR(I)*WR(I)*VELR(I)

            	 MFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELR(I) + &
                            PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I)
                 MTFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) - &
     	                     B41R(I)*B42R(I)
                 MTTFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) - &
     	                      B41R(I)*B43R(I)

            	 EFLX(I)  = (RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) - &
     	                     B40R(I)*B41R(I) - RFLX(I)

                 BTFLX(I) = BTR(I)*VELR(I) - BR(I)*VELTR(I)
                 BTTFLX(I) = BTTR(I)*VELR(I) - BR(I)*VELTTR(I)

              END IF

           ENDDO

	   RETURN

        ENDIF


        !     ------
        !     HLLD CONTINUES
        !     ------

        !  -----------------------------------------------------------------
        !  APPROXIMATE RIEMANN SOLVER: HLLD (iars = 4).
        !  -----------------------------------------------------------------

        !                      t
        !lbn*T   lbna*T          ^     lbc*T   lbma*T  lbm*T
        !---.----\---------------|-------.--------/--.------
        !     .   \              |      .       /  .
        !       .   \            |     .      /  .
        !         .   \          |    .     /  .
        !     UL    .   \      UL*    UR* /  .       UR
        !              .  \      |  .   /  .
        !                 . \    | .  / .
        !                   . \  |. /.
        !------------------------*------------------------->
        !                        O                         x


        ! The corresponding interface numerical flux at t=T is defined as:

        ! FL = F (UL)         if   lbn .GE. 0             
	! Fal                 if   lbn .LT. 0 .LE. lbna    
        ! Fal+lbc*(ucl-ual)   if   lbna .LT. 0 .LE. lbc    
        ! Far+lbc*(ucr-uar)   if   lbc .LT. 0 .LE. lbma
        ! Far                 if   lbma .LT. 0 .LE. lbm    
        ! FR = F (UR)         if   lbm .LE. 0

      
        IF (IARS .EQ. 4) THEN



        DO I = FZNP, NZNP+1


	   HLLswitch(i) = 0

           EXITTRUE = 0

		
	   ! switch to HLL in proximity of strong shocks.

           IF ((basenm(1:2) .EQ. 'MH' .OR. basenm(1:2) .EQ. 'MC' .or. basenm .eq. 'LOOP' .OR. BASENM(1:2) .EQ.  &
                'S4') .AND. ishock(i) .EQ. 1) THEN
               HLLswitch(i) = 1 ; CYCLE
           ENDIF

           !IF (ishock(i) .EQ. 1) THEN
           !    HLLswitch(i) = 1 ; CYCLE
           !ENDIF

	
           ! ============
           !HLL STATES
           ! ============          

           R(I)  = (lbm(I)*RHOR(I)*WR(I) - lbn(I)*RHOL(I)*WL(I) +                             &
                    RHOL(I)*WL(I)*VELL(I) - RHOR(I)*WR(I)*VELR(I))/(lbm(I) - lbn(I))

           M(I) = (lbm(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*velr(i) - B40R(I)*B41R(I)) - &
                   lbn(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*vell(i) - B40L(I)*B41L(I)) + &
                  ((RHOL(I)*HL(I) + B44L(I))*WL(I)*vell(i)*WL(I)*vell(i) +                    &
                    PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I)) -                                &
                  ((RHOR(I)*HR(I) + B44R(I))*WR(I)*velr(i)*WR(I)*velr(i) +                    &
                    PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I)))/(lbm(I) - lbn(I))

           MT(I) = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTR(I) -                     &
                    B40R(I)*B42R(I)) - LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTL(I) -  &
                    B40L(I)*B42L(I)) + ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) - &
                    B41L(I)*B42L(I)) - ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) - &
                    B41R(I)*B42R(I)))/(LBM(I) - LBN(I))

           MTT(I) = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTTR(I) -                      &
                     B40R(I)*B43R(I)) - LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTTL(I) -   &
                     B40L(I)*B43L(I)) + ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) -  &
                     B41L(I)*B43L(I)) - ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) -  &
                     B41R(I)*B43R(I)))/(LBM(I) - LBN(I))

           E(I)  = (LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I) -                   &
                    PR(I) - 0.5D0*B44R(I) - B40R(I)*B40R(I)) -                        &
                    LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I) -                   &
                    PL(I) - 0.5D0*B44L(I) - B40L(I)*B40L(I)) +                        &
                  ((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) - B40L(I)*B41L(I)) - &
                  ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) - B40R(I)*B41R(I)))/ &
                   (LBM(I) - LBN(I)) - R(I)

           B(I) = (BL(I)+BR(I))/2.D0

           BT(I) = (LBM(I)*BTR(I) - LBN(I)*BTL(I) + (BTL(I)*VELL(I) - BL(I)*VELTL(I)) -       &
                   (BtR(I)*VELR(I) - BR(I)*VELTR(I)))/(LBM(I) - LBN(I))

           BTT(I) = (LBM(I)*BTTR(I) - LBN(I)*BTTL(I) + (BTTL(I)*VELL(I) - BL(I)*VELTTL(I)) -  &
                    (BTTR(I)*VELR(I) - BR(I)*VELTTR(I)))/(LBM(I) - LBN(I))


	
           IF (B(I)*B(I) .LT. 1.D-15) B(I)=0.D0
           IF (BT(I)*BT(I) .LT. 1.D-15) BT(I)=0.D0
           IF (BTT(I)*BTT(I) .LT. 1.D-15) BTT(I)=0.D0


           B22 = b(i)*b(i)+bt(i)*bt(i) +btt(i)*btt(i)
           M2 = m(i)*m(i) + mt(i)*mt(i) +                          &
                mtt(i)*mtt(i)
           BM = b(i)*m(i) +  bt(i)*mt(i) +                         &
                btt(i)*mtt(i)
	   S2 = BM*BM

           niter = 500
           iter = 0
           ERRP = 1.D0

           G1G = (GAMMA_AD - 1.D0)/GAMMA_AD


           ! Provide initial guess.


           B1 = -4.D0*(E(I) + R(I) - B22)
           C1 =  M2 + B22*B22 - 2.D0*B22*(E(I) + R(I))

           CHI = B1*B1 - 12.D0*C1
           CHI = DMAX1(CHI, 0.D0)

           ZETA = (-B1 + DSQRT(CHI))/6.D0
           ZETA = DMAX1(R(I), ZETA)

           ZETAP = ZETA -  R(I)

           DO WHILE (ERRP .GT. 1.D-11 .AND. iter .LT. niter)

                    IF (iter + 1 .GT. niter) THEN

                        ierr = ierr+1

                        IF (ierr .NE. 0) THEN
                            call ERROR ('GETPRFQ3D', 'MAXIMUM NUMBER OF ITERATIONS',   &
                                         __FILE__, __LINE__)
                        ENDIF

                    ENDIF

                    ZETA = ZETAP + R(I)

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
                            call ERROR ('HLLD', 'VEL2 > 0',   &
                                         __FILE__, __LINE__)
                        ENDIF

                    ENDIF

                    CHI = ZETAP/W2 - R(I)*VEL2/(W(I)+1.D0)

                    DV2_DZETA = -2.D0*SCRCH2*(3.D0*S2_ZETA2 + SCRCH1*(S2_ZETA2*B22*B22/ZETA + M2))

                    RHO = R(I)/W(I)

                    ! KINEMATICAL TERMS

                    dchi_dZETA =  1.D0 - VEL2 - 0.5D0*W(I)*(R(I) + 2.D0*chi*W(I))*dv2_dZETA
                    drho_dZETA = -0.5D0*R(I)*W(I)*dv2_dZETA


                    ! ASSUMING IDEAL GAS
                    dp_dchi = G1G
                    dp_drho = 0.0D0
                    p_S  = chi*dp_dchi
                    DP = dp_dchi*dchi_dZETA + dp_drho*drho_dZETA

                    FUN = ZETAP - (E(I)+P_S) + 0.5D0*(B22 + (B22*M2 - S2)*SCRCH2)
                    DFUN= 1.D0 - DP - (B22*M2 - S2)*SCRCH1*SCRCH2

                    OZETA = ZETAP

                    ZETAP = ZETAP - FUN/DFUN

                    ERRP = ABS(1.D0-ZETAP/OZETA)

                    ITER = ITER + 1

           ENDDO

           ZETA = ZETAP + R(I)

           velH = (m(i) + bm*b(i)/zeta)/(zeta+b22)
           veltH = (mt(i) + bm*bt(i)/zeta)/(zeta+b22)
           velttH = (mtt(i) + bm*btt(i)/zeta)/(zeta+b22)


           ! SMALL prevents the velocity to take unphysical small values.
           IF (ABS(velH) .LE. 1.D-15) velH = 0.D0
           IF (ABS(veltH) .LE. 1.D-15) veltH = 0.D0
           IF (ABS(velttH) .LE. 1.D-15) velttH = 0.D0
                
           ! Recover the Lorentz factor in the HLLE state.
           w(i) = 1.D0/DSQRT(1.D0 - velH*velH - veltH*veltH - velttH*velttH)
        
           ! Recover primitive density.
           rho  = r(i)/w(i)
                 
           ! Calculate the internal energy eps.
           eps = (zeta/rho/w(i)/w(i) - 1.D0) / gamma_ad
                

           ! Magnetic field four-vector.
           b40H = w(i)*( velH*b(i) +     &
                         veltH*bt(i) +   &
                         velttH*btt(i))

           b41H = b(i)/w(i) + velH*b40H
           b42H = bt(i)/w(i) + veltH*b40H
           b43H = btt(i)/w(i) + velttH*b40H

           ! Assuming the Minkowski metric to be (-,+,+,+)
           b44H = -b40H*b40H+b41H*b41H+b42H*b42H+b43H*b43H


           phll = ((gamma_ad-1.D0)*rho*eps)+(0.5D0*b44h)
	   

           ! ============
           ! RL AND RR VECTORS
           ! ============   


           RD_L(I)   =  lbn(I)*RHOL(I)*WL(I) - &
                        (RHOL(I)*WL(I)*VELL(I))
           
           RM_L(I)   =  lbn(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*vell(i) - B40L(I)*B41L(I)) -   &
                        ((RHOL(I)*HL(I) + B44L(I))*wL(I)*vell(i)*WL(I)*VELL(I) +                     &
                        PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I))

           RMT_L(I)  =  LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTL(I) - B40L(I)*B42L(I)) -  &
                        ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) -                    &
     	                B41L(I)*B42L(I))
                        
           RMTT_L(I) =  LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTTL(I) - B40L(I)*B43L(I)) - &
                        ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) -                   &
     	                B41L(I)*B43L(I))

           
           RE_L(I)   =  LBN(I)*((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I) - PL(I) - 0.5D0*B44L(I) - B40L(I)*B40L(I)) - &
                        ((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) -                                           &
     	                B40L(I)*B41L(I))

           
           RB_L(I)   =  LBN(I)*B(I)
           RBT_L(I)  =  LBN(I)*BTL(I)-(BTL(I)*VELL(I) - BL(I)*VELTL(I))
           RBTT_L(I) =  LBN(I)*BTTL(I)-(BTTL(I)*VELL(I) - BL(I)*VELTTL(I))


           RD_R(I)   =  lbM(I)*RHOR(I)*WR(I) - &
                        (RHOR(I)*WR(I)*VELR(I))
           
           RM_R(I)   =  lbM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*velR(i) - B40R(I)*B41R(I)) - &
                        ((RHOR(I)*HR(I) + B44R(I))*wR(I)*velR(i)*WR(I)*VELR(I) +                   &
                        PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I))

           RMT_R(I)  =  LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTR(I) - B40R(I)*B42R(I)) - &
                        ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) -                   &
     	                B41R(I)*B42R(I))
                        
           RMTT_R(I) =  LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTTR(I) - B40R(I)*B43R(I)) - &
                        ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) -                   &
     	                B41R(I)*B43R(I))

           
           RE_R(I)   =  LBM(I)*((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I) - PR(I) - 0.5D0*B44R(I) - B40R(I)*B40R(I)) - &
                        ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) -                                           &
     	                B40R(I)*B41R(I))

           
           RB_R(I)   =  LBM(I)*B(I)
           RBT_R(I)  =  LBM(I)*BTR(I)-(BTR(I)*VELR(I) - BR(I)*VELTR(I))
           RBTT_R(I) =  LBM(I)*BTTR(I)-(BTTR(I)*VELR(I) - BR(I)*VELTTR(I))


           ! *****************************************************
           ! We start calculating the pressure.

           AH = LBM(I) - LBN(I)
           BH = RE_R(I)-RE_L(I)+LBM(I)*RM_L(I)-LBN(I)*RM_R(I)
           CH = RM_L(I)*RE_R(I)-RM_R(I)*RE_L(I)

	   scrch1 = BH**2-4.D0*AH*CH
	   scrch1 = DMAX1(scrch1, 0.D0)

           ! Since AH ne 0 always, then:
	   PB0 = 0.5D0*(- BH + DSQRT(scrch1))*(1.D0/(LBM(I) - LBN(I)))
	

	   IF (B(I) .EQ. 0.D0) THEN

                 P(I) = PB0

	   ELSE
           
              ! ====================
              ! ITERATIVE PROCESS
              ! ====================
              
              ! We need to start with an initial guess for P0, according to MB.
	                    
	      scrch1 = MAX(PL(I),PR(I))

              IF ((B(I)**2/scrch1) .LT. 1.D-2) THEN
                 P(I) = PB0
              ELSE   
                 P(I) = PHLL
              ENDIF   

              !-------------------------------------

	      ! ROOT FINDER: SECANT METHOD

              ITER   = 0
              IERR   = 0
              relerr = 1.D0
	      relerr0 =1.D0
	      EXITTRUE = 0             
	      niter = 20	

              oldP = P(I)
              dP = 100000000000.D0
              
              ! =========      
              ! STEP 1      
              ! =========  

              A_ = RM_L(I)-LBN(I)*RE_L(I)+P(I)*(1.D0-LBN(I)**2)
              G_ = RBT_L(I)*RBT_L(I) + RBTT_L(I)*RBTT_L(I)
              C_ = RMT_L(I)*RBT_L(I) + RMTT_L(I)*RBTT_L(I)
              Q_ = -A_ - G_ + B(I)**2*(1.D0-LBN(I)**2)
              X_ = B(I)*(A_*LBN(I)*B(I)+C_) - (A_+G_)*(LBN(I)*P(I)+RE_L(I))

              VEL_ALF_L   = B(I)*(A_*B(I)+LBN(I)*C_)-(A_+G_)*(P(I)+RM_L(I))                 
              VEL_ALF_TL  = Q_*RMT_L(I)+RBT_L(I)*(C_+B(I)*(LBN(I)*RM_L(I)-RE_L(I)))
              VEL_ALF_TTL = Q_*RMTT_L(I)+RBTT_L(I)*(C_+B(I)*(LBN(I)*RM_L(I)-RE_L(I)))

              
              scrch1 = VEL_ALF_L*RM_L(I)+VEL_ALF_TL*RMT_L(I)+VEL_ALF_TTL*RMTT_L(I)
              scrch1 = X_*RE_L(I) - scrch1

              WWL = P(I) + scrch1/(X_*LBN(I)-VEL_ALF_L)

              VEL_ALF_L = VEL_ALF_L/X_
              VEL_ALF_TL=VEL_ALF_TL/X_
              VEL_ALF_TTL=VEL_ALF_TTL/X_


              B_ALF_TL = -(RBT_L(I)*(LBN(I)*P(I)+RE_L(I))-B(I)*RMT_L(I))/A_
              B_ALF_TTL = -(RBTT_L(I)*(LBN(I)*P(I)+RE_L(I))-B(I)*RMTT_L(I))/A_


              A_ = RM_R(I)-LBM(I)*RE_R(I)+P(I)*(1.D0-LBM(I)**2)
              G_ = RBT_R(I)*RBT_R(I) + RBTT_R(I)*RBTT_R(I)
              C_ = RMT_R(I)*RBT_R(I) + RMTT_R(I)*RBTT_R(I)
              Q_ = -A_ - G_ + B(I)**2*(1.D0-LBM(I)**2)
              X_ = B(I)*(A_*LBM(I)*B(I)+C_) - (A_+G_)*(LBM(I)*P(I)+RE_R(I))

              VEL_ALF_R   = B(I)*(A_*B(I)+LBM(I)*C_)-(A_+G_)*(P(I)+RM_R(I))
              VEL_ALF_TR  = Q_*RMT_R(I)+RBT_R(I)*(C_+B(I)*(LBM(I)*RM_R(I)-RE_R(I)))
              VEL_ALF_TTR = Q_*RMTT_R(I)+RBTT_R(I)*(C_+B(I)*(LBM(I)*RM_R(I)-RE_R(I)))

              
              scrch1 = VEL_ALF_R*RM_R(I)+VEL_ALF_TR*RMT_R(I)+VEL_ALF_TTR*RMTT_R(I)
              scrch1 = X_*RE_R(I) - scrch1

              WWR = P(I) + scrch1/(X_*LBM(I)-VEL_ALF_R)


              VEL_ALF_R = VEL_ALF_R/X_
              VEL_ALF_TR=VEL_ALF_TR/X_
              VEL_ALF_TTR=VEL_ALF_TTR/X_

              B_ALF_TR = -(RBT_R(I)*(LBM(I)*P(I)+RE_R(I))-B(I)*RMT_R(I))/A_
              B_ALF_TTR = -(RBTT_R(I)*(LBM(I)*P(I)+RE_R(I))-B(I)*RMTT_R(I))/A_
                 

              IF (WWL .LT. 0.D0 .OR. WWR .LT. 0.D0) EXITTRUE = 1 


	      ! if W>0, we can go ahead...
              ! ......................................................................
	      
              IF (EXITTRUE .EQ. 0) THEN

              ! =========      
              ! STEP 2      
              ! =========                   

	      IF (B(I) .GT. 0.D0) THEN
	          signB = 1.D0
	      ELSE
		  signB = -1.D0
	      ENDIF

              ETA_L = - signB*DSQRT(WWL)      
              ETA_R =   signB*DSQRT(WWR)  


              scrch1 = 1.D0/(LBN(I)*P(I)+RE_L(I)+ETA_L*B(I))

              KL    = scrch1*(RM_L(I)+P(I)+ETA_L*RB_L(I))
              KTL   = scrch1*(RMT_L(I)+ETA_L*RBT_L(I))
              KTTL  = scrch1*(RMTT_L(I)+ETA_L*RBTT_L(I))

              scrch1 = 1.D0/(LBM(I)*P(I)+RE_R(I)+ETA_R*B(I))

              KR    = scrch1*(RM_R(I)+P(I)+ETA_R*RB_R(I))
              KTR   = scrch1*(RMT_R(I)+ETA_R*RBT_R(I))
              KTTR  = scrch1*(RMTT_R(I)+ETA_R*RBTT_R(I))

              LBN_ALF = KL
              LBM_ALF = KR


	      ! ========================================

              dK = KR-KL+1.D-12

              B_C    = B(I)*dK

              B_C_T  = (B_ALF_TR*(KR-VEL_ALF_R)+B(I)*VEL_ALF_TR)-  &
                       (B_ALF_TL*(KL-VEL_ALF_L)+B(I)*VEL_ALF_TL)
                           

              B_C_TT = (B_ALF_TTR*(KR-VEL_ALF_R)+B(I)*VEL_ALF_TTR)-  &
                       (B_ALF_TTL*(KL-VEL_ALF_L)+B(I)*VEL_ALF_TTL)
                 

              ! =========      
              ! STEP 3      
              ! =========                   


              KLBC = KL*B_C + KTL*B_C_T + KTTL*B_C_TT
              KRBC = KR*B_C + KTR*B_C_T + KTTR*B_C_TT

              
              VEL_C_L = KL - ((dK*B(I)*(1.D0-KL**2-KTL**2-KTTL**2))/(ETA_L*dK-KLBC))
              VEL_C_R = KR - ((dK*B(I)*(1.D0-KR**2-KTR**2-KTTR**2))/(ETA_R*dK-KRBC))	

              VEL_C1 = (VEL_C_L+VEL_C_R)*0.5D0

	     
              ! Check conditions to switch HLL.
              IF ((VEL_C_L-KL) .GT. -1.D-6 .AND. (KR-VEL_C_R) .GT. -1.D-6 .AND.             &
                  (LBN(I)-VEL_ALF_L) .LT. 1.D-6 .AND. (VEL_ALF_R-LBM(I)) .LT. 1.D-6 .AND.   &
                  (WWR-P(I)) .GT. 0.D0 .AND. (WWL-P(I)) .GT. 0.D0) THEN

                  EXITTRUE = 0	 
	
              ELSE

                  EXITTRUE = 1


	      ENDIF


              F0 = VEL_C_L-VEL_C_R; relerr0 = ABS(F0)

              !IF (relerr0 .GE. 1.D-12) THEN
              !        P(I) = 1.025D0*P(I) 
              !ENDIF   

	      ENDIF
	      ! ......................................................................              

	      IF (relerr0 .GE. 1.D-12 .and. EXITTRUE .EQ. 0) THEN


                      P(I) = 1.025D0*P(I)

                      F1 = F0


	      DO iter = 0, niter

		 !F1 = F0

                 ! =========      
                 ! STEP 1      
                 ! =========  

                 A_ = RM_L(I)-LBN(I)*RE_L(I)+P(I)*(1.D0-LBN(I)**2)
                 G_ = RBT_L(I)*RBT_L(I) + RBTT_L(I)*RBTT_L(I)
                 C_ = RMT_L(I)*RBT_L(I) + RMTT_L(I)*RBTT_L(I)
                 Q_ = -A_ - G_ + B(I)**2*(1.D0-LBN(I)**2)
                 X_ = B(I)*(A_*LBN(I)*B(I)+C_) - (A_+G_)*(LBN(I)*P(I)+RE_L(I))

                 VEL_ALF_L   = B(I)*(A_*B(I)+LBN(I)*C_)-(A_+G_)*(P(I)+RM_L(I))                 
                 VEL_ALF_TL  = Q_*RMT_L(I)+RBT_L(I)*(C_+B(I)*(LBN(I)*RM_L(I)-RE_L(I)))
                 VEL_ALF_TTL = Q_*RMTT_L(I)+RBTT_L(I)*(C_+B(I)*(LBN(I)*RM_L(I)-RE_L(I)))


		 scrch1 = VEL_ALF_L*RM_L(I)+VEL_ALF_TL*RMT_L(I)+VEL_ALF_TTL*RMTT_L(I)
	         scrch1 = X_*RE_L(I) - scrch1

	         WWL = P(I) + scrch1/(X_*LBN(I)-VEL_ALF_L)


	         VEL_ALF_L = VEL_ALF_L/X_
	         VEL_ALF_TL=VEL_ALF_TL/X_
		 VEL_ALF_TTL=VEL_ALF_TTL/X_

                 B_ALF_TL = -(RBT_L(I)*(LBN(I)*P(I)+RE_L(I))-B(I)*RMT_L(I))/A_
                 B_ALF_TTL = -(RBTT_L(I)*(LBN(I)*P(I)+RE_L(I))-B(I)*RMTT_L(I))/A_


                 A_ = RM_R(I)-LBM(I)*RE_R(I)+P(I)*(1.D0-LBM(I)**2)
                 G_ = RBT_R(I)*RBT_R(I) + RBTT_R(I)*RBTT_R(I)
                 C_ = RMT_R(I)*RBT_R(I) + RMTT_R(I)*RBTT_R(I)
                 Q_ = -A_ - G_ + B(I)**2*(1.D0-LBM(I)**2)
                 X_ = B(I)*(A_*LBM(I)*B(I)+C_) - (A_+G_)*(LBM(I)*P(I)+RE_R(I))

                 VEL_ALF_R   = B(I)*(A_*B(I)+LBM(I)*C_)-(A_+G_)*(P(I)+RM_R(I))
                 VEL_ALF_TR  = Q_*RMT_R(I)+RBT_R(I)*(C_+B(I)*(LBM(I)*RM_R(I)-RE_R(I)))
                 VEL_ALF_TTR = Q_*RMTT_R(I)+RBTT_R(I)*(C_+B(I)*(LBM(I)*RM_R(I)-RE_R(I)))


                 scrch1 = VEL_ALF_R*RM_R(I)+VEL_ALF_TR*RMT_R(I)+VEL_ALF_TTR*RMTT_R(I)
                 scrch1 = X_*RE_R(I) - scrch1

                 WWR = P(I) + scrch1/(X_*LBM(I)-VEL_ALF_R)

                 IF (WWL .LT. 0.D0 .OR. WWR .LT. 0.D0) THEN
                        EXITTRUE = 1 ; EXIT
                 ENDIF

                 VEL_ALF_R = VEL_ALF_R/X_
                 VEL_ALF_TR=VEL_ALF_TR/X_
                 VEL_ALF_TTR=VEL_ALF_TTR/X_

                 B_ALF_TR = -(RBT_R(I)*(LBM(I)*P(I)+RE_R(I))-B(I)*RMT_R(I))/A_
                 B_ALF_TTR = -(RBTT_R(I)*(LBM(I)*P(I)+RE_R(I))-B(I)*RMTT_R(I))/A_
                 

                 ! =========      
                 ! STEP 2      
                 ! =========                   


                 IF (B(I) .GT. 0.D0) THEN
                     signB = 1.D0
                 ELSE
                     signB = -1.D0
                 ENDIF

                 ETA_L = - signB*DSQRT(WWL)
                 ETA_R =   signB*DSQRT(WWR)

		 scrch1 = 1.D0/(LBN(I)*P(I)+RE_L(I)+ETA_L*B(I))

                 KL    = scrch1*(RM_L(I)+P(I)+ETA_L*RB_L(I))
                 KTL   = scrch1*(RMT_L(I)+ETA_L*RBT_L(I))
                 KTTL  = scrch1*(RMTT_L(I)+ETA_L*RBTT_L(I))

                 scrch1 = 1.D0/(LBM(I)*P(I)+RE_R(I)+ETA_R*B(I))

                 KR    = scrch1*(RM_R(I)+P(I)+ETA_R*RB_R(I))
                 KTR   = scrch1*(RMT_R(I)+ETA_R*RBT_R(I))
                 KTTR  = scrch1*(RMTT_R(I)+ETA_R*RBTT_R(I))

                 LBN_ALF = KL
                 LBM_ALF = KR


	         ! ========================================

	 	 dK = KR-KL+1.D-12

		 B_C    = B(I)*dK

                 B_C_T  = (B_ALF_TR*(KR-VEL_ALF_R)+B(I)*VEL_ALF_TR)-  &
                          (B_ALF_TL*(KL-VEL_ALF_L)+B(I)*VEL_ALF_TL)
                           

                 B_C_TT = (B_ALF_TTR*(KR-VEL_ALF_R)+B(I)*VEL_ALF_TTR)-  &
                          (B_ALF_TTL*(KL-VEL_ALF_L)+B(I)*VEL_ALF_TTL)
                 

                 ! =========      
                 ! STEP 3      
                 ! =========                   


		 KLBC = KL*B_C + KTL*B_C_T + KTTL*B_C_TT
	 	 KRBC = KR*B_C + KTR*B_C_T + KTTR*B_C_TT


		 VEL_C_L = KL - ((dK*B(I)*(1.D0-KL**2-KTL**2-KTTL**2))/(ETA_L*dK-KLBC))
                 VEL_C_R = KR - ((dK*B(I)*(1.D0-KR**2-KTR**2-KTTR**2))/(ETA_R*dK-KRBC))	

		 VEL_C1 = (VEL_C_L+VEL_C_R)*0.5D0



              ! Check conditions to switch HLL.
              IF ((VEL_C_L-KL) .GT. -1.D-6 .AND. (KR-VEL_C_R) .GT. -1.D-6 .AND.             &
                  (LBN(I)-VEL_ALF_L) .LT. 1.D-6 .AND. (VEL_ALF_R-LBM(I)) .LT. 1.D-6 .AND.   &
                  (WWR-P(I)) .GT. 0.D0 .AND. (WWL-P(I)) .GT. 0.D0) THEN

                  EXITTRUE = 0

              ELSE

                  EXITTRUE = 1 ; EXIT


              ENDIF


                
                 !=================================================


                 F1 = VEL_C_L-VEL_C_R
		 
                 IF (abs(f1) .gt. abs(f0) .and. iter .gt. 4) THEN
			EXITTRUE=1 ; EXIT
		 ENDIF

                 dP = (P(I) - oldP)/(F1 - F0)*F1

                 oldP = P(I) ; F0 = F1

                 P(I) = P(I) - dP
		

                 IF (P(I) .LT. 0.D0) P(I) = 1.D-6*ABS(P(I))

              
                 IF (abs(dP) .LT. 1.d-5*p(I) .OR. abs(f1) .LT. 1.d-6) EXIT

		 ! The algorithm is not converging.....
		 IF (iter + 1 .GT. niter) THEN
		    EXITTRUE = 1 ;EXIT
		 ENDIF

              ENDDO

	   ELSE

		P(I) = OLDP

	   ENDIF
              
	   ENDIF


	   ! Check failures in the root finder.
	   IF (EXITTRUE .NE. 0) THEN
	      HLLswitch(I) = 1 ; CYCLE
	   ENDIF


           ! Once P is calculated, we need to solve the system...
	           
           ! ================
           ! ALFVEN STATES
           ! ================


           A_ = RM_L(I)-LBN(I)*RE_L(I)+P(I)*(1.D0-LBN(I)**2)
           G_ = RBT_L(I)*RBT_L(I) + RBTT_L(I)*RBTT_L(I)
           C_ = RMT_L(I)*RBT_L(I) + RMTT_L(I)*RBTT_L(I)
           Q_ = -A_ - G_ + B(I)**2*(1.D0-LBN(I)**2)
           X_ = B(I)*(A_*LBN(I)*B(I)+C_) - (A_+G_)*(LBN(I)*P(I)+RE_L(I))


           VEL_ALF_L   = B(I)*(A_*B(I)+LBN(I)*C_)-(A_+G_)*(P(I)+RM_L(I))                 
           VEL_ALF_TL  = Q_*RMT_L(I)+RBT_L(I)*(C_+B(I)*(LBN(I)*RM_L(I)-RE_L(I)))
           VEL_ALF_TTL = Q_*RMTT_L(I)+RBTT_L(I)*(C_+B(I)*(LBN(I)*RM_L(I)-RE_L(I)))


	   scrch1 = VEL_ALF_L*RM_L(I)+VEL_ALF_TL*RMT_L(I)+VEL_ALF_TTL*RMTT_L(I)
	   scrch1 = X_*RE_L(I)-scrch1


	   WWL = P(I) + scrch1/(X_*LBN(I)-VEL_ALF_L)	


	   VEL_ALF_L   = VEL_ALF_L/X_
	   VEL_ALF_TL  = VEL_ALF_TL/X_
	   VEL_ALF_TTL = VEL_ALF_TTL/X_

           
           A_ = RM_R(I)-LBM(I)*RE_R(I)+P(I)*(1.D0-LBM(I)**2)
           G_ = RBT_R(I)*RBT_R(I) + RBTT_R(I)*RBTT_R(I)
           C_ = RMT_R(I)*RBT_R(I) + RMTT_R(I)*RBTT_R(I)
           Q_ = -A_ - G_ + B(I)**2*(1.D0-LBM(I)**2)
           X_ = B(I)*(A_*LBM(I)*B(I)+C_) - (A_+G_)*(LBM(I)*P(I)+RE_R(I))


           VEL_ALF_R   = B(I)*(A_*B(I)+LBM(I)*C_)-(A_+G_)*(P(I)+RM_R(I))    
           VEL_ALF_TR  = Q_*RMT_R(I)+RBT_R(I)*(C_+B(I)*(LBM(I)*RM_R(I)-RE_R(I)))    
           VEL_ALF_TTR = Q_*RMTT_R(I)+RBTT_R(I)*(C_+B(I)*(LBM(I)*RM_R(I)-RE_R(I)))    


           scrch1 = VEL_ALF_R*RM_R(I)+VEL_ALF_TR*RMT_R(I)+VEL_ALF_TTR*RMTT_R(I)
           scrch1 = X_*RE_R(I)-scrch1


           WWR = P(I) + scrch1/(X_*LBM(I)-VEL_ALF_R)


           VEL_ALF_R   = VEL_ALF_R/X_
           VEL_ALF_TR  = VEL_ALF_TR/X_
           VEL_ALF_TTR = VEL_ALF_TTR/X_

           IF (WWL .LT. 0.D0 .OR. WWR .LT. 0.D0) THEN
               HLLswitch(i) =1 ; cycle
	   ENDIF


              ! Check conditions to switch HLL.
              IF ((LBN(I)-VEL_ALF_L) .LT. 1.D-6 .AND. (VEL_ALF_R-LBM(I)) .LT. 1.D-6 .AND.   &
                  (WWR-P(I)) .GT. 0.D0 .AND. (WWL-P(I)) .GT. 0.D0) THEN

                  CONTINUE

              ELSE

                  HLLswitch(i) =1 ; cycle


              ENDIF


              ! Check conditions to switch HLL.
              !IF ((LBN(I)-VEL_ALF_L) .GT. 0.D0 .OR. (LBM(I)-VEL_ALF_R) .LT. 0.D0 .OR. &
              !    (WWR-P(I)) .LT. 0.D0 .OR. (WWL-P(I)) .LT. 0.D0) THEN

              ! HLLswitch(i) =1 ; cycle


              !ENDIF


           IF (B(I) .GT. 0.D0) THEN
                  signB = 1.D0
           ELSE
                  signB = -1.D0
           ENDIF

           ETA_L = - signB*DSQRT(WWL)
           ETA_R =   signB*DSQRT(WWR)

           scrch1 = 1.D0/(LBN(I)*P(I)+RE_L(I)+B(I)*ETA_L)

           KL    = scrch1*(RM_L(I)+P(I)+ETA_L*RB_L(I))
           KTL   = scrch1*(RMT_L(I)+ETA_L*RBT_L(I))
           KTTL  = scrch1*(RMTT_L(I)+ETA_L*RBTT_L(I))

           scrch1 = 1.D0/(LBM(I)*P(I)+RE_R(I)+B(I)*ETA_R)

           KR    = scrch1*(RM_R(I)+P(I)+ETA_R*RB_R(I))
           KTR   = scrch1*(RMT_R(I)+ETA_R*RBT_R(I))
           KTTR  = scrch1*(RMTT_R(I)+ETA_R*RBTT_R(I))


           ! Saving velocities.....
           LBNA(I) = KL
           LBMA(I) = KR



           ! ------------------      
           ! Left state      
           ! ------------------

	   SCRCH1 = 1.D0/(LBN(I)-VEL_ALF_L)                 

           RL_A   = RD_L(I)*SCRCH1

           B_ALF_TL  = (RBT_L(I)-B(I)*VEL_ALF_TL)*scrch1            
           B_ALF_TTL = (RBTT_L(I)-B(I)*VEL_ALF_TTL)*scrch1 

	   VB = VEL_ALF_L*B(I)+VEL_ALF_TL*B_ALF_TL+VEL_ALF_TTL*B_ALF_TTL

           EL_A   = (RE_L(I)+P(I)*VEL_ALF_L-B(I)*VB)*scrch1
           
           ML_A   = (EL_A+P(I))*VEL_ALF_L-B(I)*vB
           MTL_A  = (EL_A+P(I))*VEL_ALF_TL-B_ALF_TL*vB
           MTTL_A = (EL_A+P(I))*VEL_ALF_TTL-B_ALF_TTL*vB
           


           ! ------------------      
           ! Right state      
           ! ------------------
                 
           SCRCH1 = 1.D0/(LBM(I)-VEL_ALF_R)

           RR_A   = RD_R(I)*SCRCH1

           B_ALF_TR  = (RBT_R(I)-B(I)*VEL_ALF_TR)*scrch1
           B_ALF_TTR = (RBTT_R(I)-B(I)*VEL_ALF_TTR)*scrch1

           VB = VEL_ALF_R*B(I)+VEL_ALF_TR*B_ALF_TR+VEL_ALF_TTR*B_ALF_TTR

           ER_A   = (RE_R(I)+P(I)*VEL_ALF_R-B(I)*VB)*scrch1

           MR_A   = (ER_A+P(I))*VEL_ALF_R-B(I)*vB                                                        
           MTR_A  = (ER_A+P(I))*VEL_ALF_TR-B_ALF_TR*vB                                                   
           MTTR_A = (ER_A+P(I))*VEL_ALF_TTR-B_ALF_TTR*vB 


           ! ================
           ! CONTACT STATES
           ! ================


	   IF (LBNA(I) .LT. -1.D-6 .OR. lbma(I) .GT. 1.d-6) THEN


	   dK = (KR-KL)+1.D-12
           
	   B_C     = B(I)*dK

           B_C_T   = (B_ALF_TR*(KR-VEL_ALF_R)+B(I)*VEL_ALF_TR)-  &
                     (B_ALF_TL*(KL-VEL_ALF_L)+B(I)*VEL_ALF_TL)
                      

           B_C_TT  = (B_ALF_TTR*(KR-VEL_ALF_R)+B(I)*VEL_ALF_TTR)-  &
                     (B_ALF_TTL*(KL-VEL_ALF_L)+B(I)*VEL_ALF_TTL)


	   B_C   = B(I)                      
	   B_C_T = B_C_T/dK
	   B_C_TT = B_C_TT/dK


	   KLBC = KL*B_C + KTL*B_C_T + KTTL*B_C_TT
	   KRBC = KR*B_C + KTR*B_C_T + KTTR*B_C_TT

	   SCRCHL = (1.D0-KL**2-KTL**2-KTTL**2)/(ETA_L-KLBC)
	   SCRCHR = (1.D0-KR**2-KTR**2-KTTR**2)/(ETA_R-KRBC)


           VEL_C_L =   KL-B_C*SCRCHL 
           VEL_C_TL =  KTL-B_C_T*SCRCHL
           VEL_C_TTL = KTTL-B_C_TT*SCRCHL

           VEL_C_R =   KR-B_C*SCRCHR                                                                                 
           VEL_C_TR =  KTR-B_C_T*SCRCHR                                                                              
           VEL_C_TTR = KTTR-B_C_TT*SCRCHR

	   VEL_C = 0.5D0*(VEL_C_L+VEL_C_R) 
	   VEL_C_T = 0.5D0*(VEL_C_TL+VEL_C_TR)
	   VEL_C_TT = 0.5D0*(VEL_C_TTL+VEL_C_TTR)

	   LBC(I) = VEL_C


           !IF ((VEL_C_L-KL) .LT. 0.D0 .OR. (KR-VEL_C_R) .LT. 0.D0 .OR.             &
           !    (LBN(I)-VEL_ALF_L) .GT. 0.D0 .OR. (LBM(I)-VEL_ALF_R) .LT. 0.D0 .OR. &
           !    (WWR-P(I)) .LT. 0.D0 .OR. (WWL-P(I)) .LT. 0.D0) THEN

           !     HLLswitch(i) =1 ; cycle

           !ENDIF



              IF ((VEL_C_L-KL) .GT. -1.D-6 .AND. (KR-VEL_C_R) .GT. -1.D-6 .AND.             &
                  (LBN(I)-VEL_ALF_L) .LT. 1.D-6 .AND. (VEL_ALF_R-LBM(I)) .LT. 1.D-6 .AND.   &
                  (WWR-P(I)) .GT. 0.D0 .AND. (WWL-P(I)) .GT. 0.D0) THEN

                  CONTINUE

              ELSE

                  HLLswitch(i) =1 ; cycle


              ENDIF


	   IF (VEL_C .GT. 0.D0) THEN

	      R_A = RL_A
	      E_A = EL_A

	      M_A = ML_A
	      MT_A = MTL_A
	      MTT_A = MTTL_A

	      LBA = KL

	      VEL_ALF = VEL_ALF_L 

	   ELSE

              R_A = RR_A
              E_A = ER_A

              M_A = MR_A
              MT_A = MTR_A
              MTT_A = MTTR_A

              LBA = KR     

              VEL_ALF = VEL_ALF_R

	   ENDIF	


	   !IF (LBA .EQ. VEL_C) GO TO 10
           IF (ABS(LBA-VEL_C).LE.1.D-15) GO TO 10


           vBc = VEL_C*B_C + VEL_C_T*B_C_T + VEL_C_TT*B_C_TT



           R_C   = R_A*(LBA-VEL_ALF)/(LBA-VEL_C)

           E_C   = (LBA*E_A-M_A+P(I)*VEL_C-B(I)*vBc)&
                     /(LBA-VEL_C)


           M_C   = (E_C+P(I))*VEL_C-B(I)*vBc
           MT_C  = (E_C+P(I))*VEL_C_T-B_C_T*vBc
           MTT_C = (E_C+P(I))*VEL_C_TT-B_C_TT*vBc


           !=============
           ! COMPUTE FLUXES
           !=============

	   10 continue

	   ENDIF           
           ! ------------------
           ! Left state
           ! ------------------

           
           FRL(I)    = RHOL(I)*WL(I)*VELL(I) + LBN(I)*(RL_A-RHOL(I)*WL(I))

           FML(I)    =  ((RHOL(I)*HL(I) + B44L(I))*wL(I)*vell(i)*WL(I)*VELL(I) +      & 
                        PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I)) +                    & 
	 		LBN(I)*(ML_A-((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*vell(i) - &
                        B40L(I)*B41L(I)))


           FMTL(I)   =  ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) -       &
                        B41L(I)*B42L(I)) +                                              &
	 		LBN(I)*(MTL_A-((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTL(I) - &
                        B40L(I)*B42L(I)))


           FMTTL(I)   =  ((RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) -       &
                        B41L(I)*B43L(I)) +                                                &
                        LBN(I)*(MTTL_A-((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELTTL(I) - &
                        B40L(I)*B43L(I)))

 
           FEL(I)    =  ((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) -                  &
                        B40L(I)*B41L(I)) +                                                &
		        LBN(I)*(EL_A-((RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I) - PL(I) -     &
                        0.5D0*B44L(I) - B40L(I)*B40L(I)))    


           FBTL(I)   = (BTL(I)*VELL(I) - BL(I)*VELTL(I)) + &
		        LBN(I)*(B_ALF_TL-BTL(I))	

           FBTTL(I)  = (BTTL(I)*VELL(I) - BL(I)*VELTTL(I)) + &
                        LBN(I)*(B_ALF_TTL-BTTL(I))

           
           ! ------------------      
           ! Right state      
           ! ------------------           
       

           FRR(I)    = RHOR(I)*WR(I)*VELR(I) + LBM(I)*(RR_A-RHOR(I)*WR(I))

           FMR(I)    =  ((RHOR(I)*HR(I) + B44R(I))*wR(I)*velR(i)*WR(I)*VELR(I) +       &
                        PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I)) +                     &
                        LBM(I)*(MR_A-((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*velR(i) -  &
                        B40R(I)*B41R(I)))


           FMTR(I)   =  ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) -       &
                        B41R(I)*B42R(I)) +                                              &
                        LBM(I)*(MTR_A-((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTR(I) - &
                        B40R(I)*B42R(I)))


           FMTTR(I)   =  ((RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) -       &
                        B41R(I)*B43R(I)) +                                                &
                        LBM(I)*(MTTR_A-((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELTTR(I) - &
                        B40R(I)*B43R(I)))


           FER(I)    =  ((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) -                 &
                        B40R(I)*B41R(I)) +                                               &
                        LBM(I)*(ER_A-((RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I) - PR(I) -    &
                        0.5D0*B44R(I) - B40R(I)*B40R(I)))     


           FBTR(I)   = (BTR(I)*VELR(I) - BR(I)*VELTR(I)) + &
                        LBM(I)*(B_ALF_TR-BTR(I))

           FBTTR(I)  = (BTTR(I)*VELR(I) - BR(I)*VELTTR(I)) + &
                        LBM(I)*(B_ALF_TTR-BTTR(I))

	   
           IF (LBNA(I) .LT. -1.D-6 .OR. lbma(I) .GT. 1.d-6) THEN

	   !if (lbna(i) .EQ. vel_c) go to 11 
           IF (ABS(LBNA(I)-VEL_C).LE.1.D-15) GO TO 11

           
	   IF (VEL_C .GT. 0.D0) THEN

              FRL2(I)    = FRL(I)+LBNA(I)*(R_C-RL_A)                            
              FML2(I)    = FML(I)+LBNA(I)*(M_C-ML_A)           
              FMTL2(I)   = FMTL(I)+LBNA(I)*(MT_C-MTL_A)           
              FMTTL2(I)  = FMTTL(I)+LBNA(I)*(MTT_C-MTTL_A)           
              FEL2(I)    = FEL(I)+LBNA(I)*(E_C-EL_A)           
              FBTL2(I)   = FBTL(I)+LBNA(I)*(B_C_T-B_ALF_TL)           
              FBTTL2(I)  = FBTTL(I)+LBNA(I)*(B_C_TT-B_ALF_TTL)

	   ELSE

              FRR2(I)    = FRR(I)+LBMA(I)*(R_C-RR_A)                             
              FMR2(I)    = FMR(I)+LBMA(I)*(M_C-MR_A)           
              FMTR2(I)   = FMTR(I)+LBMA(I)*(MT_C-MTR_A)           
              FMTTR2(I)  = FMTTR(I)+LBMA(I)*(MTT_C-MTTR_A)           
              FER2(I)    = FER(I)+LBMA(I)*(E_C-ER_A)           
              FBTR2(I)   = FBTR(I)+LBMA(I)*(B_C_T-B_ALF_TR)          
              FBTTR2(I)  = FBTTR(I)+LBMA(I)*(B_C_TT-B_ALF_TTR)

	   ENDIF

	11 continue


	   ENDIF
	
        ENDDO

	! Debug errors in the HLL switch.

	!IF (SUM(HLLswitch) .GT. 0.1*NZNP) THEN
        !    call ERROR ('HLLD',             &
	!                'HLLSWITCH > 0.01 NZNP', __FILE__, __LINE__)
	!ENDIF
	!print*, SUM(HLLswitch)

        !     --------
        !     HLLD INTERCELL NUMERICAL FLUXES 
        !     --------

	
        DO I=FZNP,NZNP+1

           !IF ((FL(I)+FR(I))/2.D0 .LT. 1.D-10 .AND. (BASENM .EQ. 'CMBW' .OR. BASENM(1:2) .EQ. 'MH' .OR. BASENM .EQ. 'LOOP')) CYCLE

	 IF (HLLswitch(i) .EQ. 1) THEN


            IF (lbn(I).GT.0.D0) THEN

               RFLX(I)  = RHOL(I)*WL(I)*VELL(I)

               mflx(I) = (RHOL(I)*HL(I) + B44L(I))*wL(I)*vell(i)*WL(I)*VELL(I) +     &
                       PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I)
               mtflx(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) -   &
                        B41L(I)*B42L(I)
               mttflx(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) - &
                        B41L(I)*B43L(I)

               eflx(I)  = (RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) -            &
                        B40L(I)*B41L(I) - RFLX(I)


               !FBX(I) = 0.D0
               btflx(I) = BTL(I)*VELL(I) - BL(I)*VELTL(I)
               bttflx(I) = BTTL(I)*VELL(I) - BL(I)*VELTTL(I)

            END IF

            IF (LBM(I).LE.0.D0) THEN

               RFLX(I)  = RHOR(I)*WR(I)*VELR(I)

               MFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELR(I) +     &
                       PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I)
               MTFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) -   &
                        B41R(I)*B42R(I)
               MTTFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) - &
                        B41R(I)*B43R(I)

               EFLX(I)  = (RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) -            &
                        B40R(I)*B41R(I) - RFLX(I)


               !FBX(I) = 0.D0
               BTFLX(I) = BTR(I)*VELR(I) - BR(I)*VELTR(I)
               BTTFLX(I) = BTTR(I)*VELR(I) - BR(I)*VELTTR(I)

            END IF

	 ELSE

               IF (LBN(I).GE. 0.d0) THEN 

                   RFLX(I)  = RHOL(I)*WL(I)*VELL(I)

                   MFLX(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELL(I) +     &
                              PL(I) + 0.5D0*B44L(I) - B41L(I)*B41L(I)
                   MTFLX(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTL(I) -   &
     	                       B41L(I)*B42L(I)
                   MTTFLX(I) = (RHOL(I)*HL(I) + B44L(I))*WL(I)*VELL(I)*WL(I)*VELTTL(I) - &
     	                        B41L(I)*B43L(I)

                   EFLX(I)  = (RHOL(I)*HL(I) + B44L(I))*WL(I)*WL(I)*VELL(I) -            &
     	                       B40L(I)*B41L(I) - RFLX(I)

                   BTFLX(I) = BTL(I)*VELL(I) - BL(I)*VELTL(I)
                   BTTFLX(I) = BTTL(I)*VELL(I) - BL(I)*VELTTL(I)


               ELSE IF (lbm(I) .LE. 0.d0) THEN

                   RFLX(I)  = RHOR(I)*WR(I)*VELR(I)

                   MFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELR(I) +     &
                              PR(I) + 0.5D0*B44R(I) - B41R(I)*B41R(I)
                   MTFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTR(I) -   &
                               B41R(I)*B42R(I)
                   MTTFLX(I) = (RHOR(I)*HR(I) + B44R(I))*WR(I)*VELR(I)*WR(I)*VELTTR(I) - &
                                B41R(I)*B43R(I)

                   EFLX(I)  = (RHOR(I)*HR(I) + B44R(I))*WR(I)*WR(I)*VELR(I) -            &
                               B40R(I)*B41R(I) - RFLX(I)

                   BTFLX(I) = BTR(I)*VELR(I) - BR(I)*VELTR(I)
                   BTTFLX(I) = BTTR(I)*VELR(I) - BR(I)*VELTTR(I)


	       ELSE


                 IF (LBNA(I).GE.-1.d-6) THEN


                   RFLX(I)  = FRL(I)

                   MFLX(I) = FML(I) 
            	   MTFLX(I) = FMTL(I) 
            	   MTTFLX(I) = FMTTL(I) 

            	   EFLX(I)  = FEL(I) - RFLX(I)

            	   BTFLX(I) = FBTL(I)
            	   BTTFLX(I) = FBTTL(I)

                 ELSE IF (LBMA(I).LE.1.d-6) THEN


                   RFLX(I)  = FRR(I)

                   MFLX(I) = FMR(I)
                   MTFLX(I) = FMTR(I)
                   MTTFLX(I) = FMTTR(I)

                   EFLX(I)  = FER(I) - RFLX(I)

                   BTFLX(I) = FBTR(I)
                   BTTFLX(I) = FBTTR(I)

		 ELSE


                    IF (LBC(I).GT.0.D0) THEN

                      RFLX(I)  = FRL2(I)

            	      MFLX(I) = FML2(I)
            	      MTFLX(I) = FMTL2(I) 
            	      MTTFLX(I) = FMTTL2(I) 

            	      EFLX(I)  = FEL2(I) - RFLX(I)

            	      BTFLX(I) = FBTL2(I)
            	      BTTFLX(I) = FBTTL2(I)
            
	            ELSE

                      RFLX(I)  = FRR2(I)

                      MFLX(I) = FMR2(I) 
            	      MTFLX(I) = FMTR2(I) 
            	      MTTFLX(I) = FMTTR2(I) 

            	      EFLX(I)  = FER2(I) - RFLX(I)

            	      BTFLX(I) = FBTR2(I)
            	      BTTFLX(I) = FBTTR2(I)

		ENDIF
              ENDIF
	   ENDIF


	 ENDIF

         ENDDO

         RETURN
         
      ENDIF

END SUBROUTINE HLLDD


      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


! NAME:
!       FAST

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
!       Obtain the extreme fast magnetosonic eigenvalues using the method
!       of Newton-Rhapson.
      
! ACADEMIC REFERENCE:
!       "Magnetohidrodinamica relativista numerica:
!       Aplicaciones en Relatividad Especial y General". L.A.Ruiz. 2008
!       See also: "Relativistic fluids and magneto-fluids". A.M. Anile. 1989.


! ===============================================================================       
SUBROUTINE FAST (ca0,ca1,ca2,ca3,ca4,lamb1,lamb4)
! ===============================================================================
      
      
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables
        real(rs)   :: ca0, ca1, ca2, ca3, ca4

        !------------------------------------------------------- Output Variables
        ! size -4:mn6
        real(rs)   :: lamb1, lamb4

        !------------------------------------------------------- Local Variables
        integer(is):: niter, iter

        real(rs)   :: deriv, quartic

        !-------------------------------------------------------------------------- BEGIN

        iter = 0

        ! Initialize the Newton-Raphson process.
        lamb1 = -1.D0     ! speed of light
        quartic = 1.D0    !--> 0 (err. limit 1.e-10).

        DO WHILE (iter .LT. 100 .AND. ABS(quartic) .GE. 1.D-10)

           iter = iter + 1

           ! General form of the quartic equation.
           quartic = (((ca4*lamb1 + ca3)*lamb1 + ca2)* lamb1 + ca1)* lamb1 + ca0

           ! General form of the derivative of the quartic equation.
           deriv = ((4.D0*ca4*lamb1 + 3.D0*ca3)*lamb1 + 2.D0*ca2)*lamb1 + ca1

           lamb1 = lamb1 - quartic/deriv

        ENDDO

        ! Stop the Newton Raphson if the max. number of iterations is reached
        ! before the error limit.
        IF (iter .EQ. 500) THEN

           print*,'fast: reached iteration limit in lambda1',lamb1,quartic
           stop

        ENDIF

        iter = 0

        ! Initialize the Newton-Raphson process.
        lamb4 = 1.D0     ! speed of light
        quartic = 1.D0   ! --> (err. limit 1.D-10).

        DO WHILE (iter .LT. 100 .AND. ABS(quartic) .GE. 1.D-10)

           iter = iter + 1

           ! General form of the quartic equation.
           quartic = (((ca4*lamb4 + ca3)*lamb4 + ca2)* lamb4 + ca1)* lamb4 + ca0

           ! General form of the derivative of the quartic equation.
           deriv = ((4.D0*ca4*lamb4 + 3.D0*ca3)*lamb4 + 2.D0*ca2)*lamb4 + ca1

           lamb4 = lamb4 - quartic/deriv

        ENDDO


        ! Stop the Newton Raphson if the max. number of iterations is reached
        ! before the error limit.
        IF (iter .EQ. 500) THEN

           print*,'fast: reached iteration limit in lambda4',lamb4,quartic
           stop

        ENDIF

        !-------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE FAST
! =============================================================================== 
