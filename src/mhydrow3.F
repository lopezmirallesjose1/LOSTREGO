!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************

! NAME:
!       MHYDROW3

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
!       Calculate the numerical fluxes in the direction of the sweep.

! ACADEMIC REFERENCE:
!       "An HLLC Riemann Solver for Relativistic Flows:I. Hydrodynamics."
!        A. Mignone & G. Bodo. 2005. MNRAS. 000, 1-13.
!       "An HLLC Riemann Solver for Relativistic Flows-II, Magnetohydrodynamics."
!        A. Mignone & G. Bodo. 2006. MNRAS. 368, 1040-1054
!       "The HLLD Riemann solver based on magnetic field decomposition method for
!        the numerical simulation of MHD".
!        X. Guo, V. Florinski, C. Wang. 2016. JCP 327 543-552.
!        See also: Riemann Solvers and Numerical Methods for Fluid Dynamics. E.F.Toro. 
!        See also: Magnetohidrodinamica relativista numerica: Aplicaciones en Relatividad
!                  Especial y General. L.A. Ruiz. 2008


! =============================================================================== 
SUBROUTINE MHYDROW3 (ii,tid,xyzswp,basenm, limiter, fznp, nznp, iars,           &
                     eqs, gamma_ad, me, mp, x, xl, xr, areal,                   &
                     p, rho, rhoe, rhopr, vel, velt, veltt, f, b, bs, bt, btt,  &
                     rflx, reflx, mflx, mtflx, mttflx,                          &
                     eflx, btflx, bttflx, fflx, yt)
! =============================================================================== 

        USE workarrays
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(4) :: limiter, basenm

        integer(is) :: tid, fznp, nznp, iars, eqs, xyzswp,ii

        real(rs)   ::  me, mp, gamma_ad, yt


        ! size -4:mn5
        real(rs)   :: x(-4:*), xl(-4:*), xr(-4:*), areal(-4:*)

        ! size -4:mn5
        real(rs)   :: p(-4:*), rho(-4:*), rhoe(-4:*), rhopr(-4:*), vel(-4:*),  &
                      velt(-4:*), veltt(-4:*), f(-4:*), b(-4:*), bt(-4:*),     &
                      btt(-4:*), bs(-4:*) 

        !------------------------------------------------------- Output Variables
 
        ! size -4:mn6
        real(rs)   :: rflx(-4:*), reflx(-4:*),  mflx(-4:*), mtflx(-4:*), &
                      mttflx(-4:*), eflx(-4:*), fflx(-4:*), btflx(-4:*), & 
                      bttflx(-4:*)

        !------------------------------------------------------- Local Variables

        integer(is):: i

        real(rs)   :: tt1, tt2, tt3, norm, b40, b41, b42, b43, b44, VEL2

        integer(is), pointer :: ishock(:)

        real(rs), pointer :: rpflx(:)

        real(rs), pointer :: beta(:)

        real(rs), pointer :: wm(:), wp(:)

        real(rs), pointer :: pl(:), pr(:), rhol(:), rhor(:),       &
                             rhoel(:), rhoer(:),                   &
                             rhoprl(:), rhoprr(:),                 &
                             xel(:), xer(:), teml(:),              &
                             temr(:), vell(:), velr(:),            &
                             veltl(:), veltr(:), velttl(:),        &  
                             velttr(:), wl(:),                     &
                             wr(:), ul(:), ur(:), csl(:),          &
                             csr(:), hl(:), hr(:),                 &
                             dpdrhl(:), dpdrhr(:),                 &
                             dpdul(:), dpdur(:),                   &
                             gaml(:), gamr(:),                     &
                             rl(:), rr(:), fl(:), fr(:),           &
                             btl(:), btr(:), bttl(:),              &
                             bttr(:), b40l(:), b40r(:),            &
                             b41l(:), b41r(:), b42l(:),            &
                             b42r(:), b43l(:), b43r(:),            &
                             b44l(:), b44r(:), bl(:), br(:),       &
                             rel(:), rpl(:), rer(:), rpr(:)

        real(rs), pointer :: w(:), wvel(:), wvelt(:), wveltt(:)

        real(rs), pointer :: pm(:), pp(:), rhom(:), rhop(:), &
                             rhoem(:), rhoep(:),             &
                             rhopm(:), rhopp(:),             &
                             velm(:), velp(:), veltm(:),     &
                             veltp(:), velttm(:), velttp(:), &
                             fm(:), fp(:), btm(:),           & 
                             btp(:), bttm(:), bttp(:)
 
        !-------------------------------------------------------------------------- BEGIN

        ! -4:mn6 vectors

        rpflx  => wa(tid)%wa6mhydrow100

        wm =>     wa(tid)%wa6mhydrow101
        wp =>     wa(tid)%wa6mhydrow102

        pl     => wa(tid)%wa6mhydrow103
        pr     => wa(tid)%wa6mhydrow104
        rhol   => wa(tid)%wa6mhydrow105
        rhor   => wa(tid)%wa6mhydrow106
        rhoel  => wa(tid)%wa6mhydrow107
        rhoer  => wa(tid)%wa6mhydrow108
        rhoprl => wa(tid)%wa6mhydrow109
        rhoprr => wa(tid)%wa6mhydrow110
        xel    => wa(tid)%wa6mhydrow111
        xer    => wa(tid)%wa6mhydrow112
        teml   => wa(tid)%wa6mhydrow113
        temr   => wa(tid)%wa6mhydrow114
        vell   => wa(tid)%wa6mhydrow115
        velr   => wa(tid)%wa6mhydrow116
        veltl  => wa(tid)%wa6mhydrow117
        veltr  => wa(tid)%wa6mhydrow118
        velttl => wa(tid)%wa6mhydrow119
        velttr => wa(tid)%wa6mhydrow120
        wl     => wa(tid)%wa6mhydrow121
        wr     => wa(tid)%wa6mhydrow122
        ul     => wa(tid)%wa6mhydrow123
        ur     => wa(tid)%wa6mhydrow124
        csl    => wa(tid)%wa6mhydrow125
        csr    => wa(tid)%wa6mhydrow126
        hl     => wa(tid)%wa6mhydrow127
        hr     => wa(tid)%wa6mhydrow128
        dpdrhl => wa(tid)%wa6mhydrow129
        dpdrhr => wa(tid)%wa6mhydrow130
        dpdul  => wa(tid)%wa6mhydrow131
        dpdur  => wa(tid)%wa6mhydrow132
        gaml   => wa(tid)%wa6mhydrow133
        gamr   => wa(tid)%wa6mhydrow134 
        rl     => wa(tid)%wa6mhydrow135
        rr     => wa(tid)%wa6mhydrow136
        fl     => wa(tid)%wa6mhydrow137
        fr     => wa(tid)%wa6mhydrow138
        btl    => wa(tid)%wa6mhydrow139
        btr    => wa(tid)%wa6mhydrow140
        bttl   => wa(tid)%wa6mhydrow141
        bttr   => wa(tid)%wa6mhydrow142
        b40l   => wa(tid)%wa6mhydrow143
        b40r   => wa(tid)%wa6mhydrow144
        b41l   => wa(tid)%wa6mhydrow145
        b41r   => wa(tid)%wa6mhydrow146
        b42l   => wa(tid)%wa6mhydrow147
        b42r   => wa(tid)%wa6mhydrow148
        b43l   => wa(tid)%wa6mhydrow149
        b43r   => wa(tid)%wa6mhydrow150
        b44l   => wa(tid)%wa6mhydrow151
        b44r   => wa(tid)%wa6mhydrow152
        bl     => wa(tid)%wa6mhydrow153
        br     => wa(tid)%wa6mhydrow154
        rel    => wa(tid)%wa6mhydrow155
        rpl    => wa(tid)%wa6mhydrow156
        rer    => wa(tid)%wa6mhydrow157
        rpr    => wa(tid)%wa6mhydrow158

        
        ! -4:mn5 vectors
        w      => wa(tid)%wa5mhydrow100
        wvel   => wa(tid)%wa5mhydrow101
        wvelt  => wa(tid)%wa5mhydrow102
        wveltt => wa(tid)%wa5mhydrow103


        pm     => wa(tid)%wa5mhydrow104
        pp     => wa(tid)%wa5mhydrow105
        rhom   => wa(tid)%wa5mhydrow106
        rhop   => wa(tid)%wa5mhydrow107
        rhoem  => wa(tid)%wa5mhydrow108
        rhoep  => wa(tid)%wa5mhydrow109
        rhopm  => wa(tid)%wa5mhydrow110
        rhopp  => wa(tid)%wa5mhydrow111    
        velm   => wa(tid)%wa5mhydrow112
        velp   => wa(tid)%wa5mhydrow113
        veltm  => wa(tid)%wa5mhydrow114
        veltp  => wa(tid)%wa5mhydrow115
        velttm => wa(tid)%wa5mhydrow116
        velttp => wa(tid)%wa5mhydrow117
        fm     => wa(tid)%wa5mhydrow118
        fp     => wa(tid)%wa5mhydrow119
        btm    => wa(tid)%wa5mhydrow120
        btp    => wa(tid)%wa5mhydrow121
        bttm   => wa(tid)%wa5mhydrow122
        bttp   => wa(tid)%wa5mhydrow123
        
        ishock => wa(tid)%wa5mhydrow124
	beta => wa(tid)%wa5mhydrow125


        ! Relativistic velocity and magnetic fields.
        DO i= fznp-5, nznp+5

           w(i)      = 1.0/DSQRT(1.D0-vel(i)*vel(i)-velt(i)*velt(i)-veltt(i)*veltt(i))
           wvel(i)   = w(i)*vel(i)
           wvelt(i)  = w(i)*velt(i)
           wveltt(i) = w(i)*veltt(i)

           b40  = w(i)*(vel(i)*b(i) + velt(i)*bt(i) + veltt(i)*btt(i))  
           b41  = b(i)/w(i) + b40*vel(i)
           b42  = bt(i)/w(i) + b40*velt(i)
           b43  = btt(i)/w(i) + b40*veltt(i)

           ! Assuming the Minkowski metric to be (-,+,+,+)
           b44 = -b40*b40+b41*b41+b42*b42+b43*b43

           ! Beta factor (magnetic pressure fraction) for shock detection.
           beta(i) =  b44/2.D0/p(i)
           ishock(i) = 0
	   
        ENDDO


        ! 'Detect' only works with MCMM/VLMM limiter, which requires
        ! the detection of a shock.

        IF (limiter .EQ. 'MCMM' .OR. limiter .EQ. 'VLMM' .OR. iars .EQ. 4) THEN 
            call DETECT (nznp, tid, vel, beta, p, ishock)
        ENDIF


        DO i= fznp, nznp+1

            IF (p(i) .EQ. 1.D-15) THEN
                ishock(i) = 1
                ishock(i+1) = 1
            ENDIF

            IF ((BASENM(1:2) .EQ. 'S4' .OR. BASENM .EQ. 'SYBW') .AND. ((XYZSWP .EQ. 1 .AND. YT   .LT. 12.D0) .OR.  &
                                                                       (XYZSWP .EQ. 2 .AND. X(I) .LT. 12.D0) .OR.  &
                                                                       (XYZSWP .EQ. 3 .AND. YT   .LT. 12.D0))) THEN
                ISHOCK(I) = 1


            ENDIF

        ENDDO

	IF (limiter .EQ. 'MMOD' .OR. limiter .EQ. 'MC00' .OR. limiter .EQ. 'VANL' .OR. &
	    limiter .EQ. 'MCMM' .OR. limiter .EQ. 'VLMM') THEN	

            ! PRIMITIVE variables reconstruction.
            call PLM (limiter, ishock, nznp, x, xl, xr, rhom, rho, rhop)
            call PLM (limiter, ishock, nznp, x, xl, xr, rhoem, rhoe, rhoep)
            call PLM (limiter, ishock, nznp, x, xl, xr, rhopm, rhopr, rhopp)
            call PLM (limiter, ishock, nznp, x, xl, xr, pm, p, pp)
            call PLM (limiter, ishock, nznp, x, xl, xr, btm, bt, btp)
            call PLM (limiter, ishock, nznp, x, xl, xr, bttm, btt, bttp)
            call PLM (limiter, ishock, nznp, x, xl, xr, fm, f, fp)

            call PLM (limiter, ishock, nznp, x, xl, xr, velm, vel, velp)
            call PLM (limiter, ishock, nznp, x, xl, xr, veltm, velt, veltp)
            call PLM (limiter, ishock, nznp, x, xl, xr, velttm, veltt, velttp)

	ELSE IF (limiter .EQ. 'MP3_' .OR. limiter .EQ. 'MP5_' .OR. &
		 limiter .EQ. 'MP7_' .OR. limiter .EQ. 'MP9_') THEN

            ! PRIMITIVE variables reconstruction.
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, rhom, rho, rhop)
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, rhoem, rhoe, rhoep)
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, rhopm, rhopr, rhopp)
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, pm, p, pp)
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, btm, bt, btp)
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, bttm, btt, bttp)
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, fm, f, fp)

            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, velm, vel, velp)
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, veltm, velt, veltp)
            call MPREC (tid, limiter, ishock, nznp, x, xl, xr, velttm, veltt, velttp)

	ENDIF

        DO i= fznp-1, nznp+1

	   VEL2 = VELM(I)*VELM(I) + VELTM(I)*VELTM(I) + VELTTM(I)*VELTTM(I)

           IF (VEL2.GE.1.D0) THEN

              RHOM(I) = RHO(I)
	      RHOEM(I)= RHOE(I)
	      RHOPM(I)= RHOPR(I)
              PM(I)   = P(I)
              FM(I)   = F(I)

              BTM(I)  = BT(I)
              BTTM(I) = BTT(I)

              VELM(I)   = VEL(I)
              VELTM(I)  = VELT(I)
              VELTTM(I) = VELTT(I)

           ENDIF 
	   
           WM(I)  = 1.D0/DSQRT(1.D0 - VELM(I)*VELM(I) - &
                    VELTM(I)*VELTM(I) - VELTTM(I)*VELTTM(I))


	   VEL2 = VELP(I)*VELP(I) + VELTP(I)*VELTP(I) + VELTTP(I)*VELTTP(I)


           IF (VEL2.GE.1.D0) THEN

              RHOP(I) = RHO(I)
              RHOEP(I)= RHOE(I)
              RHOPP(I)= RHOPR(I)
              PP(I)   = P(I)
              FP(I)   = F(I)

              BTP(I)  = BT(I)
              BTTP(I) = BTT(I)

              VELP(I)   = VEL(I)
              VELTP(I)  = VELT(I)
              VELTTP(I) = VELTT(I)

           ENDIF
	 
	   WP(I)  = 1.D0/DSQRT(1.D0 - VELP(I)*VELP(I) - &
                    VELTP(I)*VELTP(I) - VELTTP(I)*VELTTP(I))

        ENDDO

   
        DO i= fznp, nznp+1
 

           !!------------------!------------------!-------
           !!                  !                  !
           !!                  !                  !
           !!U1M    U 1     U1P!U2M     U 2    U2P!U3M
           !!        *         !         *        !
           !!                  !                  !
           !!                  !                  !
           !!------------------!------------------!--------
           !1                  2                  3
           !                                     
           
           ! In the interface:
           ! UL (i) = UP (i-1)
           ! UR (i) = UM (i)

           ! Left states:
           rhol(i)  = rhop(i-1)
           rhoel(i) = rhoep(i-1)
           rhoprl(i)= rhopp(i-1)
	   
           pl(i)    = pp(i-1)
           tt1      = pl(i) / rhol(i)

           vell(i)  = velp(i-1)
           veltl(i) = veltp(i-1)
           velttl(i)= velttp(i-1)
           wl(i)    = wp(i-1)
           xel(i)   = rhoel(i)/rhol(i)
           teml(i)  = pl(i) / (rhol(i)*( xel(i)/me + (1.D0-xel(i))/mp ))    

	   bl(i)    = bs(i)	   
           btl(i)   = btp(i-1)
           bttl(i)  = bttp(i-1)

           b40l(i)  = wl(i)*(vell(i)*bl(i) + veltl(i)*btl(i) + velttl(i)*bttl(i))  
           b41l(i)  = bl(i)/wl(i) + b40l(i)*vell(i)
           b42l(i)  = btl(i)/wl(i) + b40l(i)*veltl(i)
           b43l(i)  = bttl(i)/wl(i) + b40l(i)*velttl(i)
           
           ! Assuming the metric to be (-,+,+,+).
           b44l(i) = -b40l(i)*b40l(i)+b41l(i)*b41l(i)+b42l(i)*b42l(i)+b43l(i)*b43l(i)

           fl(i) = fp(i-1)
           

           ! Right states:
           rhor(i)  = rhom(i)
           rhoer(i) = rhoem(i)
           rhoprr(i)= rhopm(i) 
           pr(i)    = pm(i)
           tt3      = pr(i) / rhor(i)

           velr(i)  = velm(i)
           veltr(i) = veltm(i)
           velttr(i) = velttm(i)
           wr(i)    = wm(i)
           xer(i)   = rhoer(i)/rhor(i)
           temr(i)  = pr(i) / (rhor(i)*( xer(i)/me + (1.D0-xer(i))/mp ))

           br(i)    = bl (i)
           btr(i)   = btm(i)
           bttr(i)  = bttm(i)
           
           b40r(i)  = wr(i)*(velr(i)*br(i) + veltr(i)*btr(i) + velttr(i)*bttr(i))  
           b41r(i)  = br(i)/wr(i) + b40r(i)*velr(i)
           b42r(i)  = btr(i)/wr(i) + b40r(i)*veltr(i)
           b43r(i)  = bttr(i)/wr(i) + b40r(i)*velttr(i)
           
           ! Assuming the metric to be (-,+,+,+).
           b44r(i) = -b40r(i)*b40r(i)+b41r(i)*b41r(i)+b42r(i)*b42r(i)+b43r(i)*b43r(i)

           fr(i) = fm(i)
		           
        ENDDO


        ! Equation of state for the unknown thermodynamical quantities.
        IF (eqs .EQ. 2) THEN

           call EOS_SYNGE (1, tid, fznp, nznp+1, me, mp, rhol, xel, teml, ul,   &
                           pl, gaml, hl, csl, dpdrhl, dpdul)	
           call EOS_SYNGE (1, tid, fznp, nznp+1, me, mp, rhor, xer, temr, ur,   &
                           pr, gamr, hr, csr, dpdrhr, dpdur)

        ELSE

           ! Test
           call EOS (1, fznp, nznp+1, gamma_ad, ul, rhol, pl, dpdrhl, dpdul, hl, csl)   
           call EOS (1, fznp, nznp+1, gamma_ad, ur, rhor, pr, dpdrhr, dpdur, hr, csr) 

        ENDIF

	
        ! Definition of the conserved densities D;
        DO i= fznp, nznp+1  

           rl(i)    = rhol(i)*wl(i)
           rel(i)   = rhoel(i)*wl(i)
           rpl(i)   = rhoprl(i)*wl(i)
           rr(i)    = rhor(i)*wr(i)
           rer(i)   = rhoer(i)*wr(i)
           rpr(i)   = rhoprr(i)*wr(i)
	   	   	   	   
        ENDDO



        ! HLL/HLLC/HLLD Approximate Riemann Solver.
        call HLLDD (ii,tid,xyzswp,basenm, iars, fznp, nznp, rhol, rhor,                       &
                   vell, velr, veltl, veltr, velttl, velttr, pl, pr, bl, br, btl, btr,        &
                   bttl, bttr, b40l, b40r, b41l, b41r, b42l,fl,fr,                            &
                   b42r, b43l, b43r, b44l, b44r, rl, rr, wl, wr, hl, hr, csl, csr,            &
                   rflx, mflx, mtflx, mttflx, eflx, btflx, bttflx, gamma_ad, ishock)


        DO i=fznp,nznp+1
	    
            IF (rflx(i) .GT. 0.D0) THEN
  
               reflx(i) = rel(i)*rflx(i)/rl(i)
               rpflx(i) = rpl(i)*rflx(i)/rl(i)
               	      
            ELSE

               reflx(i) = rer(i)*rflx(i)/rr(i)
               rpflx(i) = rpr(i)*rflx(i)/rr(i)
               
            ENDIF

         ENDDO


         ! Renormalize fluxes (error correction):
         DO i=fznp,nznp+1
	    
            IF (rflx(i) .NE. 0.D0) THEN
	       	
               norm = (reflx(i)+rpflx(i))/rflx(i)
               reflx(i) = reflx(i)/norm

            ENDIF

         ENDDO
	 

         DO i= fznp, nznp+1
		
            rflx(i)   = rflx(i)   *areal(i)
            reflx(i)  = reflx(i)  *areal(i)
            mflx(i)   = mflx(i)   *areal(i)
            mtflx(i)  = mtflx(i)  *areal(i)
	    mttflx(i) = mttflx(i) *areal(i)
            eflx(i)   = eflx(i)   *areal(i)
	    btflx(i)  = btflx(i)  *areal(i)
	    bttflx(i) = bttflx(i) *areal(i)

            IF (rflx(i) .GT. 0.D0) THEN
                fflx(i) = fl(i)*rflx(i)
            ELSE
      	  	fflx(i) = fr(i)*rflx(i)
    	    ENDIF
	    	    
         ENDDO

         !---------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE MHYDROW3
! =============================================================================== 


