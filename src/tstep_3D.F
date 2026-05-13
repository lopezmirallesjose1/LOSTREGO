!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************

! NAME:
!       TSTEP (3D VERSION)

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
!       Compute the minimum DT using the Courant condition.

! ACADEMIC REFERENCE:
!       " Relativistic MHD simulations of extragalactic jets."
!       T. Leismann, L.Anton, M.A.Aloy et al. 2005. A&A 436, 503-526 (2005).
!       See also: "An HLLC Riemann Solver for Relativistic Flows-II, Magnetohydrodynamics."
!       A. Mignone & G. Bodo. 2006. MNRAS. 368, 1040-1054


! =============================================================================== 
SUBROUTINE TSTEP_3D (limiter, nx, ny, nz, mnx5, mny5, mnz5, mnx6, mny6, mnz6, nsdim,       &
                  itstp, nstep, eqs, timej, cfl, me, mp, gamma_ad,                         &
                  bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                          &
                  xznl, xznr, xzn, yznl, yznr, yzn, zznl, zznr, zzn,                       &
                  physic, boundPosX, boundPosY, boundPosZ, dt, CT, basenm, m_0, cl, kb, rho_0a, R_b)
! =============================================================================== 


        USE tipos
        USE parallel
        USE workarrays
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(4):: limiter, basenm

        integer(is):: nx, ny, nz, mnx5, mny5, mnz5, mnx6, mny6, mnz6,     &
                      nsdim, itstp, nstep, eqs,                           &
                      bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz  

        real(rs):: timej, cfl, me, mp, gamma_ad                        

        real(rs):: xznl(*), xznr(*), xzn(*),       &
                   yznl(*), yznr(*), yzn(*),       &
                   zznl(*), zznr(*), zzn(*)

        type(Tphysic) :: physic(-4:mnx5,-4:mny5,-4:mnz5)
	type(ConsTransp) :: CT(-4:mnx6,-4:mny6,-4:mnz6)

        type(TboundPosX) :: boundPosX(-4:5)
        type(TboundPosY) :: boundPosY(-4:5)
        type(TboundPosZ) :: boundPosZ(-4:5)

	real(rs) :: m_0, cl, kb, rho_0a, R_b

        !------------------------------------------------------- Output Variables 

        real(rs) :: dt

        !------------------------------------------------------- Local Variables

        integer(is) :: tid, i, j, k, ic, jc, kc, pmin, reducemindt, kk, xyzswp, ierr

        integer(is), pointer:: ishock(:)  

        integer(is):: ict(128), jct(128), kct(128)

        integer(is) :: fznp, nznp, fzn, nzn, nxp, nyp, nzp, bndmin, bndmax
      
        real(rs):: dtt(128), soundt(128), dtcc_omp(128), dtestmax, dtc_omp(128)

        real(rs):: lamb1xl, lamb1xr, lamb4xl, lamb4xr,           &
                   lamb1yl, lamb1yr, lamb4yl, lamb4yr,           &
                   lamb1zl, lamb1zr, lamb4zl, lamb4zr,           &
                   lambda, sound, v, olddt, tt1, tt2, tt3 

        real(rs):: dtcc, b40, b41, b42, b43, b44

        real(rs) :: y, yl, yr, dy,        &
                    z, zl, zr, dz




        real(rs):: ca0, ca1, ca2, ca3, ca4, calf2, omega2, aux1, aux2		 

        real(rs), pointer :: scrch1(:)

# ifdef OPENMP
        integer :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
# endif
        integer :: np

        real(rs), pointer :: vel(:), velt(:), veltt(:),           &
                             rho(:), rhoe(:), rhopr(:), p(:),     &
                             b(:), bt(:), btt(:), bs(:),          &
                             x(:), xl(:), xr(:), dx(:), f(:)       

        real(rs), pointer :: wm(:), wp(:)

        real(rs), pointer :: pl(:), pr(:), rhol(:), rhor(:),     &
                             rhoel(:), rhoer(:),                 &
                             rhoprl(:), rhoprr(:),               &
                             xel(:), xer(:), teml(:),            &
                             temr(:), vell(:), velr(:),          &
                             veltl(:), veltr(:), velttl(:),      &  
                             velttr(:), wl(:),                   &
                             wr(:), ul(:), ur(:), csl(:),        &
                             csr(:), hl(:), hr(:),               &
                             dpdrhl(:), dpdrhr(:),               &
                             dpdul(:), dpdur(:),                 &
                             gaml(:), gamr(:),                   &
                             btl(:), btr(:), bttl(:),            &
                             bttr(:), b40l(:), b40r(:),          &
                             b41l(:), b41r(:), b42l(:),          &
                             b42r(:), b43l(:), b43r(:),          &
                             b44l(:), b44r(:), fl(:), fr(:),     &
                             bl(:), br(:)

        real(rs), pointer :: w(:), wvel(:), wvelt(:), wveltt(:)

        real(rs), pointer :: pm(:), pp(:), rhom(:), rhop(:),     &
                             rhoem(:), rhoep(:), fm(:), fp(:),   &
                             rhopm(:), rhopp(:),                 &
                             velm(:), velp(:), veltm(:),         &
                             veltp(:), velttm(:), velttp(:),     &
                             btm(:), btp(:), bttm(:), bttp(:)


	real(rs), pointer :: trac(:,:,:)

	real(rs) :: vwi, uwi, shift, rorb

	real(rs) :: xe, tem, dedt, dtco, dtc, dtcmin

        !--------------------------------------------------------------------------- BEGIN

        dtcc    = 0.D0
	 dtc    = 1.D0
        olddt   = dt

        vwi      = 6.67D-3
        uwi      = 9.05D-10
        shift    = 10.D0
        Rorb     = 500.D0


!$OMP PARALLEL DEFAULT(shared)                                                        &  
!$OMP PRIVATE(tid,i,j,k,ishock)                                                       &
!$OMP PRIVATE(scrch1,bndmin,bndmax,fzn,nzn,fznp,nznp,nyp,nzp,xyzswp)                  &
!$OMP PRIVATE(lamb1xl,lamb1xr,lamb4xl,lamb4xr,lamb1yl,lamb1yr,lamb4yl,lamb4yr)        &
!$OMP PRIVATE(lamb1zl, lamb1zr, lamb4zl, lamb4zr, b40, b41, b42, b43, b44)            &
!$OMP PRIVATE(vel,velt,veltt,rho,rhoe,rhopr,p,f,b,bt,btt,bs,x,xl,xr,dx)               &
!$OMP PRIVATE(y,yl,yr,dy,z,zl,zr,dz)                                                  &
!$OMP PRIVATE(wm,wp,pm,pp,rhom,rhop,rhoem,rhoep,rhopm,rhopp,velm,velp,veltm,veltp)    &
!$OMP PRIVATE(fm,fp,velttm,velttp,btm,btp,bttm,bttp,w,wvel,wvelt,wveltt)              &
!$OMP PRIVATE(pl,pr,rhol,rhor,rhoel,rhoer,rhoprl,rhoprr,xel,xer,teml,temr)            &
!$OMP PRIVATE(vell,velr,veltl,veltr,velttl,velttr,wl,wr,ul,ur,csl,csr,fl,fr,trac)     &
!$OMP PRIVATE(hl,hr,dpdrhl,dpdrhr,dpdul,dpdur,gaml,gamr,btl,btr,bttl,bttr)            &
!$OMP PRIVATE(b40l,b40r,b41l,b41r,b42l,b42r,b43l,b43r,b44l,b44r,bl,br)                &
!$OMP PRIVATE(dtestmax, calf2,omega2,aux1,aux2,ca0, ca1, ca2, ca3, ca4, tt1, tt2,tt3) &
!$OMP PRIVATE(lambda, xe, tem, dedt, dtco, dtcmin)                                    &
!$OMP REDUCTION(max:dtcc)                                                             &                                                            
!$OMP REDUCTION(min:dtc)



# ifdef OPENMP
            tid = omp_get_thread_num()+1
            np = omp_get_num_threads()
# else
            tid = 1
            np  = 1
# endif


                  
        ! -4:mn5 vectors

        vel   => wa(tid)%wa5tstep100
        velt  => wa(tid)%wa5tstep101
        veltt => wa(tid)%wa5tstep102 
        rho   => wa(tid)%wa5tstep103
        rhoe  => wa(tid)%wa5tstep104
        rhopr => wa(tid)%wa5tstep105
        p     => wa(tid)%wa5tstep106
        x     => wa(tid)%wa5tstep107
        xl    => wa(tid)%wa5tstep108
        xr    => wa(tid)%wa5tstep109
        dx    => wa(tid)%wa5tstep110

        f    => wa(tid)%wa5tstep111

        b    => wa(tid)%wa5tstep112
        bt   => wa(tid)%wa5tstep113
        btt  => wa(tid)%wa5tstep114

        w      => wa(tid)%wa5tstep115
        wvel   => wa(tid)%wa5tstep116
        wvelt  => wa(tid)%wa5tstep117
        wveltt => wa(tid)%wa5tstep118
 
        pm     => wa(tid)%wa5tstep119
        pp     => wa(tid)%wa5tstep120
        rhom   => wa(tid)%wa5tstep121
        rhop   => wa(tid)%wa5tstep122
        rhoem  => wa(tid)%wa5tstep123
        rhoep  => wa(tid)%wa5tstep124
        rhopm  => wa(tid)%wa5tstep125
        rhopp  => wa(tid)%wa5tstep126    
        velm   => wa(tid)%wa5tstep127
        velp   => wa(tid)%wa5tstep128
        veltm  => wa(tid)%wa5tstep129
        veltp  => wa(tid)%wa5tstep130
        velttm => wa(tid)%wa5tstep131
        velttp => wa(tid)%wa5tstep132
        fm     => wa(tid)%wa5tstep133
        fp     => wa(tid)%wa5tstep134
        btm    => wa(tid)%wa5tstep135
        btp    => wa(tid)%wa5tstep136
        bttm   => wa(tid)%wa5tstep137
        bttp   => wa(tid)%wa5tstep138

        ishock => wa(tid)%wa5tstep139
	scrch1 => wa(tid)%wa5tstep140

	bs     => wa(tid)%wa5tstep141
        
        ! -4:mn6
        
        wm =>     wa(tid)%wa6tstep101
        wp =>     wa(tid)%wa6tstep102

        pl     => wa(tid)%wa6tstep103
        pr     => wa(tid)%wa6tstep104
        rhol   => wa(tid)%wa6tstep105
        rhor   => wa(tid)%wa6tstep106
        rhoel  => wa(tid)%wa6tstep107
        rhoer  => wa(tid)%wa6tstep108
        rhoprl => wa(tid)%wa6tstep109
        rhoprr => wa(tid)%wa6tstep110
        xel    => wa(tid)%wa6tstep111
        xer    => wa(tid)%wa6tstep112
        teml   => wa(tid)%wa6tstep113
        temr   => wa(tid)%wa6tstep114
        vell   => wa(tid)%wa6tstep115
        velr   => wa(tid)%wa6tstep116
        veltl  => wa(tid)%wa6tstep117
        veltr  => wa(tid)%wa6tstep118
        velttl => wa(tid)%wa6tstep119
        velttr => wa(tid)%wa6tstep120
        wl     => wa(tid)%wa6tstep121
        wr     => wa(tid)%wa6tstep122
        ul     => wa(tid)%wa6tstep123
        ur     => wa(tid)%wa6tstep124
        csl    => wa(tid)%wa6tstep125
        csr    => wa(tid)%wa6tstep126
        hl     => wa(tid)%wa6tstep127
        hr     => wa(tid)%wa6tstep128
        dpdrhl => wa(tid)%wa6tstep129
        dpdrhr => wa(tid)%wa6tstep130
        dpdul  => wa(tid)%wa6tstep131
        dpdur  => wa(tid)%wa6tstep132
        gaml   => wa(tid)%wa6tstep133
        gamr   => wa(tid)%wa6tstep134 
        btl    => wa(tid)%wa6tstep135
        btr    => wa(tid)%wa6tstep136
        bttl   => wa(tid)%wa6tstep137
        bttr   => wa(tid)%wa6tstep138
        b40l   => wa(tid)%wa6tstep139
        b40r   => wa(tid)%wa6tstep140
        b41l   => wa(tid)%wa6tstep141
        b41r   => wa(tid)%wa6tstep142
        b42l   => wa(tid)%wa6tstep143
        b42r   => wa(tid)%wa6tstep144
        b43l   => wa(tid)%wa6tstep145
        b43r   => wa(tid)%wa6tstep146
        b44l   => wa(tid)%wa6tstep147
        b44r   => wa(tid)%wa6tstep148
        bl     => wa(tid)%wa6tstep149
        br     => wa(tid)%wa6tstep150
        fl     => wa(tid)%wa6tstep151
        fr     => wa(tid)%wa6tstep152


	! 3D Working Array
	trac   => wa(tid)%wa3tstep100
  
        ict(tid) = 0
        jct(tid) = 0 
        kct(tid) = 0
      
        dtcc_omp(tid) = 0.D0
	 dtc_omp(tid) = 1.D0

	dtcmin = 10.D0

        ! =============================
        !           X-SWEEP
        ! =============================

        ! XSweep own parameters.
        xyzswp  = 1
        bndmin = bndmnx
        bndmax = bndmxx
        fzn    = 1
        nzn    = nx
        fznp   = fzn
        nznp   = nzn
        nyp    = ny
        nzp    = nz

!$OMP DO SCHEDULE (static)
        DO k=1, nzp
           DO j= 1, nyp

             DO i=fznp, nznp+1
                bs(i) = CT(i,j,k)%bx
             ENDDO

              ! Take the primitive variables at the cell centers.

              DO i= fznp, nznp

                 vel(i)   = physic(i,j,k)%velx
                 velt(i)  = physic(i,j,k)%vely
                 veltt(i) = physic(i,j,k)%velz 

                 rho(i)  = physic(i,j,k)%densty
                 rhoe(i) = physic(i,j,k)%denstye
                 rhopr(i)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                 p(i)    = physic(i,j,k)%pres

                 b(i)   = physic(i,j,k)%bxc
                 bt(i)  = physic(i,j,k)%byc
                 btt(i) = physic(i,j,k)%bzc

                 f(i)   = physic(i,j,k)%tracer
		 trac(i,j,k) = physic(i,j,k)%tracer
        
                 xl(i)   = xznl(i)
                 x(i)    = xzn(i)
                 xr(i)   = xznr(i)
                 dx(i)   = xznr(i) - xznl(i)

              ENDDO


              ! We need to provide boundary conditions to the first and last nodes.

              IF (mpi_coords(1) .EQ. 0 .OR. mpi_coords(1) .EQ. (mpi_dims(1)-1)) THEN
                                        
                yl   = yznl(j)
                y    = yzn(j)
                yr   = yznr(j)
                dy   = yr-yl

                zl   = zznl(k)
                z    = zzn(k)
                zr   = zznr(k)
                dz   = zr-zl


                call BNDRY (fznp, nznp, bndmin, bndmax, xyzswp, p, rho, rhoe, rhopr,      &
                            vel, velt, veltt, f, b, bt, btt, x, xl, xr, dx,               &
                            y, yl, yr, dy, z, zl, zr, dz, uwi, vwi, shift, Rorb, gamma_ad)

              
              ENDIF


              ! ====================
              ! LEFT GHOST CELLS
              ! ====================

              IF (mpi_coords(1) .NE. 0) THEN

                 DO i= fznp-5, fznp-1

                    vel(i)   = physic(i,j,k)%velx
                    velt(i)  = physic(i,j,k)%vely
                    veltt(i) = physic(i,j,k)%velz

                    b(i)     = physic(i,j,k)%bxc
                    bt(i)    = physic(i,j,k)%byc
                    btt(i)   = physic(i,j,k)%bzc
                         
                    rho(i)    = physic(i,j,k)%densty
                    rhoe(i)   = physic(i,j,k)%denstye
                    rhopr(i)  = physic(i,j,k)%densty - physic(i,j,k)%denstye
                    p(i)      = physic(i,j,k)%pres
                        
                    f(i)    = physic(i,j,k)%tracer
                   
                    xl(i)   = boundPosX(i)%xznl 
                    x(i)    = boundPosX(i)%xzn  
                    xr(i)   = boundPosX(i)%xznr 
                    dx(i)   = boundPosX(i)%dx   

                 ENDDO
                
              ENDIF


              ! ====================
              ! RIGHT GHOST CELLS
              ! ====================


              IF (mpi_coords(1) .NE. (mpi_dims(1) -1)) THEN

                 DO i= nznp+1, nznp+5

                    vel(i)   = physic(i,j,k)%velx
                    velt(i)  = physic(i,j,k)%vely
                    veltt(i) = physic(i,j,k)%velz

                    b(i)     = physic(i,j,k)%bxc
                    bt(i)    = physic(i,j,k)%byc
                    btt(i)   = physic(i,j,k)%bzc
                         
                    rho(i)    = physic(i,j,k)%densty
                    rhoe(i)   = physic(i,j,k)%denstye
                    rhopr(i)  = physic(i,j,k)%densty - physic(i,j,k)%denstye
                    p(i)      = physic(i,j,k)%pres
                         
                    f(i)    = physic(i,j,k)%tracer
                    
                    xl(i)   = boundPosX(i-nznp)%xznl 
                    x(i)    = boundPosX(i-nznp)%xzn  
                    xr(i)   = boundPosX(i-nznp)%xznr 
                    dx(i)   = boundPosX(i-nznp)%dx   
                         
                 ENDDO
                   
              ENDIF


              ! Relativistic velocity and magnetic fields.
              DO i= fznp-5, nznp+5

                 w(i)      = 1.D0/DSQRT(1.D0-vel(i)*vel(i)-velt(i)*velt(i)-veltt(i)*veltt(i))
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
                 scrch1(i) =  b44/2.D0/p(i)
                 ishock(i) = 0

              ENDDO

              ! 'Detect' only works with MCMM limiter.
	      IF (limiter .EQ. 'MCMM' .OR. limiter .EQ. 'VLMM') THEN 
                  call DETECT (nznp, tid, vel, scrch1, p, ishock)
	      ENDIF



              yl   = yznl(j)
              y    = yzn(j)
              yr   = yznr(j)
              dy   = yr-yl

              zl   = zznl(k)
              z    = zzn(k)
              zr   = zznr(k)
              dz   = zr-zl


              DO i= fznp, nznp+1

                  IF (p(i) .EQ. 1.D-15) THEN
                      ishock(i) = 1
                      ishock(i+1) = 1
                  ENDIF


                  IF (BASENM(1:2) .EQ. 'S4' .AND. Y .LT. 12.D0) THEN
                      ISHOCK(I) = 1
                  ENDIF

              ENDDO


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


              DO i= fznp-1, nznp+1

                 SCRCH1(i) = velm(i)*velm(i) + veltm(i)*veltm(i) + velttm(i)*velttm(i)

                 IF (scrch1(i).GE.1.D0) THEN

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

                    SCRCH1(i) = velm(i)*velm(i) + veltm(i)*veltm(i) + velttm(i)*velttm(i)

                 END IF
                 
                 wm(i)     = 1.D0/DSQRT( 1.D0 - SCRCH1(i) )

                 SCRCH1(i) = velp(i)*velp(i) + veltp(i)*veltp(i) + velttp(i)*velttp(i)

                 IF (scrch1(i).GE.1.D0) THEN

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

                    SCRCH1(i) = velp(i)*velp(i) + veltp(i)*veltp(i) + velttp(i)*velttp(i)

                 END IF

                 wp(i)     = 1.D0/DSQRT( 1.D0 - SCRCH1(i) ) 

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

                 ! In the interphase:
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

                 br(i)    = bl(i)
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


              ! Equation of state for unknown thermodynamical quantities.
              IF (eqs .EQ. 2) THEN

                 call EOS_SYNGE (1, tid, fznp, nznp+1, me, mp, rhol, xel, teml, ul,   &
                                 pl, gaml, hl, csl, dpdrhl, dpdul)	
                 call EOS_SYNGE (1, tid, fznp, nznp+1, me, mp, rhor, xer, temr, ur,   &
                                 pr, gamr, hr, csr, dpdrhr, dpdur)

              ELSE

                 call EOS (1, fznp, nznp+1, gamma_ad, ul, rhol, pl, dpdrhl, dpdul, hl, csl)   
                 call EOS (1, fznp, nznp+1, gamma_ad, ur, rhor, pr, dpdrhr, dpdur, hr, csr) 

              ENDIF



              !	----------
              !	GETTING THE EXTREME FAST MAGNETOSONIC EIGENVALUES
              !	----------

              DO i=fznp, nznp+1
	
	 		
                 ! Coefficients for the left interface
	
                 CALF2  = B44L(I)/(RHOL(I)*HL(I) + B44L(I))
                 OMEGA2 = CSL(I)**2*(1.D0 - CALF2) + CALF2
                 AUX1   = OMEGA2/WL(I)**2
                 AUX2   = CSL(I)**2/(RHOL(I)*HL(I) + B44L(I))/WL(I)**4

                 CA0    = VELL(I)**4*(1.D0 - OMEGA2) - VELL(I)**2*AUX1 + B41L(I)**2*AUX2
                 CA1    = -4.D0*VELL(I)**3*(1.D0 - OMEGA2) + 2.D0*VELL(I)*AUX1 -        &
                          2.D0*B40L(I)*B41L(I)*AUX2
                 CA2    = 6.D0*VELL(I)**2*(1.D0 - OMEGA2) - (1.D0 - VELL(I)**2)*AUX1 +  &
                          (B40L(I)**2 - B41L(I)**2)*AUX2
                 CA3    = -4.D0*VELL(I)*(1.D0 - OMEGA2) - 2*VELL(I)*AUX1 +              &
                          2.D0*B40L(I)*B41L(I)*AUX2
                 CA4    = 1.D0 - OMEGA2*(VELL(I)**2 + VELtL(I)**2 + VELttL(I)**2) -     &
                          B40L(I)**2*AUX2


                 ! Calculate the extreme fast magnetosonic eigenvalues for the left interface.
                 call FAST (ca0,ca1,ca2,ca3,ca4,lamb1xl,lamb4xl)
		 

                 ! Coefficients for the right interface

	
                 CALF2  = B44R(I)/(RHOR(I)*HR(I) + B44R(I))
                 OMEGA2 = CSR(I)**2*(1.D0 - CALF2) + CALF2
                 AUX1   = OMEGA2/WR(I)**2
                 AUX2   = CSR(I)**2/(RHOR(I)*HR(I) + B44R(I))/WR(I)**4

                 CA0    = VELR(I)**4*(1.D0 - OMEGA2) - VELR(I)**2*AUX1 + B41R(I)**2*AUX2
                 CA1    = -4.D0*VELR(I)**3*(1.D0 - OMEGA2) + 2.D0*VELR(I)*AUX1 -        &
                           2.D0*B40R(I)*B41R(I)*AUX2
                 CA2    = 6.D0*VELR(I)**2*(1.D0 - OMEGA2) - (1.D0 - VELR(I)**2)*AUX1 +  &
                          (B40R(I)**2 - B41R(I)**2)*AUX2
                 CA3    = -4.D0*VELR(I)*(1.D0 - OMEGA2) - 2*VELR(I)*AUX1 +              &
                           2.D0*B40R(I)*B41R(I)*AUX2
                 CA4    = 1.D0 - OMEGA2*(VELR(I)**2 + VELtR(I)**2 + VELttR(I)**2) -     &
                           B40R(I)**2*AUX2


                 ! Calculate the extreme fast magnetosonic eigenvalues for the right interface.
                 call FAST (ca0,ca1,ca2,ca3,ca4,lamb1xr,lamb4xr)

                 ! Upper and lower bounds of the wavespeeds at the interfaces.  
                 lambda = MAX(ABS(lamb4xl),ABS(lamb4xr),ABS(lamb1xl),ABS(lamb1xr)) 

                 ! Maximum 1/tstep.
                 dtestmax = lambda/(xr(i)-xl(i)) 

                 IF (dtestmax .GT. dtcc_omp(tid)) THEN

                    ict(tid)   = i
                    jct(tid)   = j
                    kct(tid)   = k
                    dtcc_omp(tid) = dtestmax
                    dtt(tid) = dtestmax
                    soundt(tid) = 0.5D0*(csr(i)+csl(i+1))

                 ENDIF

                 dtcc = MAX (dtcc, dtcc_omp(tid))

              ENDDO    
           ENDDO
        ENDDO
!$OMP END DO

        ! Synchronization: Each thread waits for all the other threads at the end of the sweep.
!$OMP BARRIER


        ! =============================
        !           Y-SWEEP
        ! =============================

        IF (nsdim .GE. 2) THEN

           ! YSweep own parameters.

           xyzswp  = 2
           bndmin = bndmny
           bndmax = bndmxy
           fzn    = 1
           nzn    = ny
           fznp   = fzn
           nznp   = nzn
           nxp    = nx
           nzp    = nz
           

!$OMP DO SCHEDULE (static)
           DO k= 1, nzp
              DO i= 1, nxp

                DO j = fznp, nznp+1
                   bs(j) = CT(i,j,k)%by
                ENDDO
                 
                 ! Take the primitive variables at the cell centers.
                 DO j= fznp, nznp

                    vel(j)  = physic(i,j,k)%vely
                    velt(j) = physic(i,j,k)%velz
                    veltt(j)= physic(i,j,k)%velx

                    b(j)    = physic(i,j,k)%byc
                    bt(j)   = physic(i,j,k)%bzc
                    btt(j)  = physic(i,j,k)%bxc
                   
                    rho(j)  = physic(i,j,k)%densty
                    rhoe(j) = physic(i,j,k)%denstye
                    rhopr(j)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                    p(j)    = physic(i,j,k)%pres
                   
                    f(j) = trac(i,j,k)
                   
                    xl(j)   = yznl(j)
                    x(j)    = yzn(j)
                    xr(j)   = yznr(j)
                    dx(j)   = xr(j)-xl(j)
                    
                 ENDDO


                 ! We need to provide boundary conditions to the first and last nodes. 

                 IF (mpi_coords(2) .EQ. 0  .OR.    &
                     mpi_coords(2) .EQ. (mpi_dims(2)-1)) THEN 

                   yl   = zznl(k)
                   y    = zzn(k)
                   yr   = zznr(k)
                   dy   = yr-yl

                   zl   = xznl(i)
                   z    = xzn(i)
                   zr   = xznr(i)
                   dz   = zr-zl


                   call BNDRY (fznp, nznp, bndmin, bndmax, xyzswp, p, rho, rhoe, rhopr,      &
                               vel, velt, veltt, f, b, bt, btt, x, xl, xr, dx,               &
                               y, yl, yr, dy, z, zl, zr, dz, uwi, vwi, shift, Rorb, gamma_ad)



		      IF (mpi_coords(2) .EQ. 0) THEN
                      IF (basenm(1:2) .EQ. 'MH' .OR. basenm(1:2) .EQ. 'MC') THEN
                         DO j= fznp-5, fznp-1

                            IF (physic(i,j,k)%bxc**2+physic(i,j,k)%bzc**2 .NE. 0.D0) THEN

                            b(j)    = 0.D0
                            bt(j)   = physic(i,j,k)%bzc
                            btt(j)  = physic(i,j,k)%bxc

                            ENDIF

                         ENDDO
                      ENDIF
		      ENDIF

              
                 ENDIF


                 ! ====================                
                 ! LEFT GHOST CELLS
                 ! ====================


                 IF (mpi_coords(2) .NE. 0) THEN  

                    DO j= fznp-5, fznp-1

                       vel(j)  = physic(i,j,k)%vely
                       velt(j) = physic(i,j,k)%velz
                       veltt(j)= physic(i,j,k)%velx

                       b(j)    = physic(i,j,k)%byc
                       bt(j)   = physic(i,j,k)%bzc
                       btt(j)  = physic(i,j,k)%bxc
                   
                       rho(j)  = physic(i,j,k)%densty
                       rhoe(j) = physic(i,j,k)%denstye
                       rhopr(j)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                       p(j)    = physic(i,j,k)%pres

                       f(j) = physic(i,j,k)%tracer
                       
                       xl(j)   = boundPosY(j)%yznl
                       x(j)    = boundPosY(j)%yzn
                       xr(j)   = boundPosY(j)%yznr
                       dx(j)   = boundPosY(j)%dx

                    ENDDO

                 ENDIF


                 ! ====================
                 ! RIGHT GHOST CELLS
                 ! ====================

                 IF (mpi_coords(2) .NE. mpi_dims(2)-1) THEN

                    DO j= nznp+1, nznp+5

                       vel(j)  = physic(i,j,k)%vely
                       velt(j) = physic(i,j,k)%velz
                       veltt(j)= physic(i,j,k)%velx

                       b(j)    = physic(i,j,k)%byc
                       bt(j)   = physic(i,j,k)%bzc
                       btt(j)  = physic(i,j,k)%bxc
                   
                       rho(j)  = physic(i,j,k)%densty
                       rhoe(j) = physic(i,j,k)%denstye
                       rhopr(j)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                       p(j)    = physic(i,j,k)%pres

                       f(j) = physic(i,j,k)%tracer
                   
                       xl(j)   = boundPosY(j-nznp)%yznl
                       x(j)    = boundPosY(j-nznp)%yzn
                       xr(j)   = boundPosY(j-nznp)%yznr
                       dx(j)   = boundPosY(j-nznp)%dx

                    ENDDO

                 ENDIF  

		                              
                 ! Relativistic velocity and magnetic fields.
                 DO j= fznp-5, nznp+5

                    w(j)      = 1.D0/DSQRT(1.D0-vel(j)*vel(j)-velt(j)*velt(j)-veltt(j)*veltt(j))
                    wvel(j)   = w(j)*vel(j)
                    wvelt(j)  = w(j)*velt(j)
                    wveltt(j) = w(j)*veltt(j)
                    b40  = w(j)*(vel(j)*b(j) + velt(j)*bt(j) + veltt(j)*btt(j))  
                    b41  = b(j)/w(j) + b40*vel(j)
                    b42  = bt(j)/w(j) + b40*velt(j)
                    b43  = btt(j)/w(j) + b40*veltt(j)

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    ! Beta factor (magnetic pressure fraction) for shock detection.
                    scrch1(j) =  b44/2.D0/p(j)
                    ishock(j) = 0

                 ENDDO

                 ! 'Detect' only works with MCMM limiter.
                 IF (limiter .EQ. 'MCMM' .OR. limiter .EQ. 'VLMM') THEN
                     call DETECT (nznp, tid, vel, scrch1, p, ishock)
		 ENDIF



                 DO j= fznp, nznp+1

                     IF (p(j) .EQ. 1.D-15) THEN
                         ishock(j) = 1
                	 ishock(j+1) = 1
            	     ENDIF


                     IF (BASENM(1:2) .EQ. 'S4' .AND. X(J) .LT. 12.D0) THEN
                	ISHOCK(J) = 1
                     ENDIF

                 ENDDO


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


                 DO j= fznp-1, nznp+1

                    SCRCH1(j) = velm(j)*velm(j) + veltm(j)*veltm(j) + velttm(j)*velttm(j)

                    IF (scrch1(j).GE.1.D0) THEN

                       RHOM(j) = RHO(j)
                       RHOEM(j)= RHOE(j)
                       RHOPM(j)= RHOPR(j)
                       PM(j)   = P(j)
                       FM(j)   = F(j)

                       BTM(j)  = BT(j)
                       BTTM(j) = BTT(j)

                       VELM(j)   = VEL(j)
                       VELTM(j)  = VELT(j)
                       VELTTM(j) = VELTT(j)

                       SCRCH1(j) = velm(j)*velm(j) + veltm(j)*veltm(j) + velttm(j)*velttm(j)

                    ENDIF

                    wm(j)     = 1.D0/DSQRT( 1.D0 - SCRCH1(j) )

                    SCRCH1(j) = velp(j)*velp(j) + veltp(j)*veltp(j) + velttp(j)*velttp(j)

                    IF (scrch1(j).GE.1.D0) THEN

                       RHOP(j) = RHO(j)
                       RHOEP(j)= RHOE(j)
                       RHOPP(j)= RHOPR(j)
                       PP(j)   = P(j)
                       FP(j)   = F(j)

                       BTP(j)  = BT(j)
                       BTTP(j) = BTT(j)

                       VELP(j)   = VEL(j)
                       VELTP(j)  = VELT(j)
                       VELTTP(j) = VELTT(j)

                       SCRCH1(j) = velp(j)*velp(j) + veltp(j)*veltp(j) + velttp(j)*velttp(j)

                    ENDIF

                    wp(j)     = 1.D0/DSQRT( 1.D0 - SCRCH1(j) )

                 ENDDO  
      

                 DO j= fznp, nznp+1
                   
                    !!------------------!------------------!-------
                    !!                  !                  !
                    !!                  !                  !
                    !!U1M    U 1     U1P!U2M     U 2    U2P!U3M
                    !!         *         !         *        !
                    !!                  !                  !
                    !!                  !                  !
                    !!------------------!------------------!--------
                    !1                  2                  3

                    ! In the interphase:
                    ! UL (j) = UP (j-1)
                    ! UR (j) = UM (j)

                    ! Left states:
                    rhol(j)  = rhop(j-1)
                    rhoel(j) = rhoep(j-1)
                    rhoprl(j)= rhopp(j-1)
                    pl(j)    = pp(j-1)
                    tt1      = pl(j) / rhol(j)
                    
                    vell(j)  = velp(j-1)
                    veltl(j) = veltp(j-1)
                    velttl(j)= velttp(j-1)
                    wl(j)    = wp(j-1)
                    xel(j)   = rhoel(j)/rhol(j)
                    teml(j)  = pl(j) / (rhol(j)*( xel(j)/me + (1.D0-xel(j))/mp )) 
   
                    bl(j)    = bs(j)
                    btl(j)   = btp(j-1)
                    bttl(j)  = bttp(j-1)
                   
                    b40l(j)  = wl(j)*(vell(j)*bl(j) + veltl(j)*btl(j) + velttl(j)*bttl(j))  
                    b41l(j)  = bl(j)/wl(j) + b40l(j)*vell(j)
                    b42l(j)  = btl(j)/wl(j) + b40l(j)*veltl(j)
                    b43l(j)  = bttl(j)/wl(j) + b40l(j)*velttl(j)
           
                    ! Assuming the metric to be (-,+,+,+).
                    b44l(j) = -b40l(j)*b40l(j)+b41l(j)*b41l(j)+b42l(j)*b42l(j)+b43l(j)*b43l(j)

                    fl(j) = fp(j-1)
           

                    ! Right states:
                    rhor(j)  = rhom(j)
                    rhoer(j) = rhoem(j)
                    rhoprr(j)= rhopm(j) 
                    pr(j)    = pm(j)
                    tt3      = pr(j) / rhor(j)

                    velr(j)  = velm(j)
                    veltr(j) = veltm(j)
                    velttr(j) = velttm(j)
                    wr(j)    = wm(j)
                    xer(j)   = rhoer(j)/rhor(j)
                    temr(j)  = pr(j) / (rhor(j)*( xer(j)/me + (1.D0-xer(j))/mp ))

                    br(j)    = bl(j)
                    btr(j)   = btm(j)
                    bttr(j)  = bttm(j)
                   
                    b40r(j)  = wr(j)*(velr(j)*br(j) + veltr(j)*btr(j) + velttr(j)*bttr(j))  
                    b41r(j)  = br(j)/wr(j) + b40r(j)*velr(j)
                    b42r(j)  = btr(j)/wr(j) + b40r(j)*veltr(j)
                    b43r(j)  = bttr(j)/wr(j) + b40r(j)*velttr(j)
                   
                    ! Assuming the metric to be (-,+,+,+).
                    b44r(j) = -b40r(j)*b40r(j)+b41r(j)*b41r(j)+b42r(j)*b42r(j)+b43r(j)*b43r(j)

                    fr(j) = fm(j)
           
                 ENDDO


                 ! Equation of state for unknown thermodynamical quantities.
                 IF (eqs .EQ. 2) THEN

                    call EOS_SYNGE (1, tid, fznp, nznp+1, me, mp, rhol, xel, teml, ul,   &
                                    pl, gaml, hl, csl, dpdrhl, dpdul)	
                    call EOS_SYNGE (1, tid, fznp, nznp+1, me, mp, rhor, xer, temr, ur,   &
                                    pr, gamr, hr, csr, dpdrhr, dpdur)

                 ELSE

                    call EOS (1, fznp, nznp+1, gamma_ad, ul, rhol, pl, dpdrhl, dpdul, hl, csl)   
                    call EOS (1, fznp, nznp+1, gamma_ad, ur, rhor, pr, dpdrhr, dpdur, hr, csr) 

                 ENDIF

   
                 !	----------
                 !	GETTING THE EXTREME FAST MAGNETOSONIC EIGENVALUES
                 !	----------

                 DO j=fznp, nznp+1
                    
                    ! Coefficients for the left interface

                    CALF2  = B44L(j)/(RHOL(j)*HL(j) + B44L(j))
                    OMEGA2 = CSL(j)**2*(1.D0 - CALF2) + CALF2
                    AUX1   = OMEGA2/WL(j)**2
                    AUX2   = CSL(j)**2/(RHOL(j)*HL(j) + B44L(j))/WL(j)**4

                    CA0    = VELL(j)**4*(1.D0 - OMEGA2) - VELL(j)**2*AUX1 + B41L(j)**2*AUX2
                    CA1    = -4.D0*VELL(j)**3*(1.D0 - OMEGA2) + 2.D0*VELL(j)*AUX1 -         &
                              2.D0*B40L(j)*B41L(j)*AUX2
                    CA2    = 6.D0*VELL(j)**2*(1.D0 - OMEGA2) - (1.D0 - VELL(j)**2)*AUX1 +   &
                             (B40L(j)**2 - B41L(j)**2)*AUX2
                    CA3    = -4.D0*VELL(j)*(1.D0 - OMEGA2) - 2*VELL(j)*AUX1 +               &
                              2.D0*B40L(j)*B41L(j)*AUX2
                    CA4    = 1.D0 - OMEGA2*(VELL(j)**2 + VELtL(j)**2 + VELttL(j)**2) -      &
                             B40L(j)**2*AUX2

                    ! Calculate the extreme fast magnetosonic eigenvalues for the left interface.
                    call FAST (ca0,ca1,ca2,ca3,ca4,lamb1yl,lamb4yl)
		 

                    ! Coefficients for the right interface

	
                    CALF2  = B44R(j)/(RHOR(j)*HR(j) + B44R(j))
                    OMEGA2 = CSR(j)**2*(1.D0 - CALF2) + CALF2
                    AUX1   = OMEGA2/WR(j)**2
                    AUX2   = CSR(j)**2/(RHOR(j)*HR(j) + B44R(j))/WR(j)**4

                    CA0    = VELR(j)**4*(1.D0 - OMEGA2) - VELR(j)**2*AUX1 + B41R(j)**2*AUX2
                    CA1    = -4.D0*VELR(j)**3*(1.D0 - OMEGA2) + 2.D0*VELR(j)*AUX1 -   &
                              2.D0*B40R(j)*B41R(j)*AUX2
                    CA2    = 6.D0*VELR(j)**2*(1.D0 - OMEGA2) - (1.D0 - VELR(j)**2)*AUX1 +  &
                            (B40R(j)**2 - B41R(j)**2)*AUX2
                    CA3    = -4.D0*VELR(j)*(1.D0 - OMEGA2) - 2*VELR(j)*AUX1 +   &
                              2.D0*B40R(j)*B41R(j)*AUX2
                    CA4    = 1.D0 - OMEGA2*(VELR(j)**2 + VELtR(j)**2 + VELttR(j)**2) -  &
                             B40R(j)**2*AUX2


                    ! Calculate the extreme fast magnetosonic eigenvalues for the right interface.
                    call FAST (ca0,ca1,ca2,ca3,ca4,lamb1yr,lamb4yr)

                    ! Upper and lower bounds of the wavespeeds at the interfaces.  
                    lambda = MAX(ABS(lamb4yl),ABS(lamb4yr),ABS(lamb1yl),ABS(lamb1yr)) 

                    ! Maximum tstep.
                    dtestmax = lambda/(xr(j)-xl(j))

                    IF (dtestmax .GT. dtcc_omp(tid)) THEN

                       ict(tid)   = i
                       jct(tid)   = j
                       kct(tid)   = k
                       dtcc_omp(tid) = dtestmax
                       dtt(tid) = dtestmax
                       soundt(tid) = 0.5D0*(csr(j)+csl(j+1))

                    ENDIF

                    dtcc = MAX (dtcc, dtcc_omp(tid))
                   
                 ENDDO    
              ENDDO
           ENDDO
!$OMP END DO

	ENDIF

        ! Synchronization: Each thread waits for all the other threads at the end of the sweep.
!$OMP BARRIER

        ! =============================
        !           Z-SWEEP
        ! =============================

        IF (nsdim .EQ. 3) THEN

           ! ZSweep own parameters.
           xyzswp  = 3
           bndmin = bndmnz
           bndmax = bndmxz
           fzn    = 1
           nzn    = nz
           fznp   = fzn
           nznp   = nzn
           nxp    = nx
           nyp    = ny

           
!$OMP DO SCHEDULE (static)

           DO j= 1, nyp
              DO i= 1, nxp

                 DO k=fznp, nznp+1
                    bs(k) = CT(i,j,k)%bz
                 ENDDO
                 
                 ! Take the primitive variables at the cell centers.
                 DO k= fznp, nznp

                    vel(k)  = physic(i,j,k)%velz
                    velt(k) = physic(i,j,k)%velx
                    veltt(k)= physic(i,j,k)%vely

                    b(k)    = physic(i,j,k)%bzc
                    bt(k)   = physic(i,j,k)%bxc
                    btt(k)  = physic(i,j,k)%byc
                   
                    rho(k)  = physic(i,j,k)%densty
                    rhoe(k) = physic(i,j,k)%denstye
                    rhopr(k)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                    p(k)    = physic(i,j,k)%pres

                    f(k)    = trac(i,j,k)
                   
                    xl(k)   = zznl(k)
                    x(k)    = zzn(k)
                    xr(k)   = zznr(k)
                    dx(k)   = xr(k)-xl(k)
                    
                 ENDDO


                 ! We need to provide boundary conditions to the first and last nodes.

                 IF (mpi_coords(3) .EQ. 0 .OR.    &
                     mpi_coords(3) .EQ. (mpi_dims(3)-1)) THEN
 

                    yl   = xznl(i)
                    y    = xzn(i)
                    yr   = xznr(i)
                    dy   = yr-yl

                    zl   = yznl(j)
                    z    = yzn(j)
                    zr   = yznr(j)
                    dz   = zr-zl


                    call BNDRY (fznp, nznp, bndmin, bndmax, xyzswp, p, rho, rhoe, rhopr,      &
                                vel, velt, veltt, f, b, bt, btt, x, xl, xr, dx,               &
                                y, yl, yr, dy, z, zl, zr, dz, uwi, vwi, shift, Rorb, gamma_ad)

              
                 ENDIF


                 ! ====================
                 ! LEFT GHOST CELLS
                 ! ====================


                 IF (mpi_coords(3) .NE. 0) THEN  

                    DO k= fznp-5, fznp-1

                       vel(k)  = physic(i,j,k)%velz
                       velt(k) = physic(i,j,k)%velx
                       veltt(k)= physic(i,j,k)%vely
                       
                       b(k)    = physic(i,j,k)%bzc
                       bt(k)   = physic(i,j,k)%bxc
                       btt(k)  = physic(i,j,k)%byc
                   
                       rho(k)  = physic(i,j,k)%densty
                       rhoe(k) = physic(i,j,k)%denstye
                       rhopr(k)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                       p(k)    = physic(i,j,k)%pres

                       f(k) = physic(i,j,k)%tracer
                   
                       xl(k)   = boundPosZ(k)%zznl
                       x(k)    = boundPosZ(k)%zzn
                       xr(k)   = boundPosZ(k)%zznr
                       dx(k)   = boundPosZ(k)%dx
                       
                    ENDDO

                 ENDIF


                 ! ====================
                 ! RIGHT GHOST CELLS
                 ! ====================

                 IF (mpi_coords(3) .NE. mpi_dims(3)-1) THEN

                    DO k= nznp+1, nznp+5

                       vel(k)  = physic(i,j,k)%velz
                       velt(k) = physic(i,j,k)%velx
                       veltt(k)= physic(i,j,k)%vely

                       b(k)    = physic(i,j,k)%bzc
                       bt(k)   = physic(i,j,k)%bxc
                       btt(k)  = physic(i,j,k)%byc
                   
                       rho(k)  = physic(i,j,k)%densty
                       rhoe(k) = physic(i,j,k)%denstye
                       rhopr(k)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                       p(k)    = physic(i,j,k)%pres

                       f(k) = physic(i,j,k)%tracer
                   
                       xl(k)   = boundPosZ(k-nznp)%zznl
                       x(k)    = boundPosZ(k-nznp)%zzn
                       xr(k)   = boundPosZ(k-nznp)%zznr
                       dx(k)   = boundPosZ(k-nznp)%dx

                    ENDDO

                 ENDIF   

            
                 ! Relativistic velocity and magnetic fields.
                 DO k= fznp-5, nznp+5

                    w(k)      = 1.D0/DSQRT(1.D0-vel(k)*vel(k)-velt(k)*velt(k)-veltt(k)*veltt(k))
                    wvel(k)   = w(k)*vel(k)
                    wvelt(k)  = w(k)*velt(k)
                    wveltt(k) = w(k)*veltt(k)
                    b40  = w(k)*(vel(k)*b(k) + velt(k)*bt(k) + veltt(k)*btt(k))  
                    b41  = b(k)/w(k) + b40*vel(k)
                    b42  = bt(k)/w(k) + b40*velt(k)
                    b43  = btt(k)/w(k) + b40*veltt(k)

                    ! Assuming the Minkowski metric to be (-,+,+,+)
                    b44 = -b40*b40+b41*b41+b42*b42+b43*b43

                    ! Beta factor (magnetic pressure fraction) for shock detection.
                    scrch1(k) =  b44/2.D0/p(k)
                    ishock(k) = 0

                 ENDDO

                 ! 'Detect' only works with MCMM limiter.
                 IF (limiter .EQ. 'MCMM' .OR. limiter .EQ. 'VLMM') THEN
                     call DETECT (nznp, tid, vel, scrch1, p, ishock)
		 ENDIF


                    yl   = xznl(i)
                    y    = xzn(i)
                    yr   = xznr(i)
                    dy   = yr-yl

                    zl   = yznl(j)
                    z    = yzn(j)
                    zr   = yznr(j)
                    dz   = zr-zl


                 DO k= fznp, nznp+1

                     IF (p(k) .EQ. 1.D-15) THEN
                	ishock(k) = 1
                	ishock(k+1) = 1
            	     ENDIF


                     IF (BASENM(1:2) .EQ. 'S4' .AND. Z .LT. 12.D0) THEN
                        ISHOCK(K) = 1
                     ENDIF

                 ENDDO


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


                 DO k= fznp-1, nznp+1

                    SCRCH1(k) = velm(k)*velm(k) + veltm(k)*veltm(k) + velttm(k)*velttm(k)

                    IF (scrch1(k).GE.1.D0) THEN

                       RHOM(k) = RHO(k)
                       RHOEM(k)= RHOE(k)
                       RHOPM(k)= RHOPR(k)
                       PM(k)   = P(k)
                       FM(k)   = F(k)

                       BTM(k)  = BT(k)
                       BTTM(k) = BTT(k)

                       VELM(k)   = VEL(k)
                       VELTM(k)  = VELT(k)
                       VELTTM(k) = VELTT(k)

                       SCRCH1(k) = velm(k)*velm(k) + veltm(k)*veltm(k) + velttm(k)*velttm(k)

                    END IF

                    wm(k)     = 1.D0/DSQRT( 1.D0 - SCRCH1(k) )

                    SCRCH1(k) = velp(k)*velp(k) + veltp(k)*veltp(k) + velttp(k)*velttp(k)

                    IF (scrch1(k).GE.1.D0) THEN

                       RHOP(k) = RHO(k)
                       RHOEP(k)= RHOE(k)
                       RHOPP(k)= RHOPR(k)
                       PP(k)   = P(k)
                       FP(k)   = F(k)

                       BTP(k)  = BT(k)
                       BTTP(k) = BTT(k)

                       VELP(k)   = VEL(k)
                       VELTP(k)  = VELT(k)
                       VELTTP(k) = VELTT(k)

                       SCRCH1(k) = velp(k)*velp(k) + veltp(k)*veltp(k) + velttp(k)*velttp(k)

                    END IF

                    wp(k)     = 1.D0/DSQRT( 1.D0 - SCRCH1(k) )

                 ENDDO

                 DO k= fznp, nznp+1
                      
                    !!------------------!------------------!-------
                    !!                  !                  !
                    !!                  !                  !
                    !!U1M    U 1     U1P!U2M     U 2    U2P!U3M
                    !!        *         !         *        !
                    !!                  !                  !
                    !!                  !                  !
                    !!------------------!------------------!--------
                    !1                  2                  3

                    ! In the interphase:
                    ! UL (j) = UP (j-1)
                    ! UR (j) = UM (j)

                    ! Left states:
                    rhol(k)  = rhop(k-1)
                    rhoel(k) = rhoep(k-1)
                    rhoprl(k)= rhopp(k-1)
                    pl(k)    = pp(k-1)
                    tt1      = pl(k) / rhol(k)

                    vell(k)  = velp(k-1)
                    veltl(k) = veltp(k-1)
                    velttl(k)= velttp(k-1)
                    wl(k)    = wp(k-1)
                    xel(k)   = rhoel(k)/rhol(k)
                    teml(k)  = pl(k) / (rhol(k)*( xel(k)/me + (1.D0-xel(k))/mp )) 
   
                    bl(k)    = bs(k)
                    btl(k)   = btp(k-1)
                    bttl(k)  = bttp(k-1)
                    
                    b40l(k)  = wl(k)*(vell(k)*bl(k) + veltl(k)*btl(k) + velttl(k)*bttl(k))  
                    b41l(k)  = bl(k)/wl(k) + b40l(k)*vell(k)
                    b42l(k)  = btl(k)/wl(k) + b40l(k)*veltl(k)
                    b43l(k)  = bttl(k)/wl(k) + b40l(k)*velttl(k)
           
                    ! Assuming the metric to be (-,+,+,+).
                    b44l(k) = -b40l(k)*b40l(k)+b41l(k)*b41l(k)+b42l(k)*b42l(k)+b43l(k)*b43l(k)

                    fl(k) = fp(k-1)
           

                    ! Right states:
                    rhor(k)  = rhom(k)
                    rhoer(k) = rhoem(k)
                    rhoprr(k)= rhopm(k) 
                    pr(k)    = pm(k)
                    tt3      = pr(k) / rhor(k)

                    velr(k)  = velm(k)
                    veltr(k) = veltm(k)
                    velttr(k) = velttm(k)
                    wr(k)    = wm(k)
                    xer(k)   = rhoer(k)/rhor(k)
                    temr(k)  = pr(k) / (rhor(k)*( xer(k)/me + (1.D0-xer(k))/mp ))

                    br(k)    = bl(k)
                    btr(k)   = btm(k)
                    bttr(k)  = bttm(k)
                   
                    b40r(k)  = wr(k)*(velr(k)*br(k) + veltr(k)*btr(k) + velttr(k)*bttr(k))  
                    b41r(k)  = br(k)/wr(k) + b40r(k)*velr(k)
                    b42r(k)  = btr(k)/wr(k) + b40r(k)*veltr(k)
                    b43r(k)  = bttr(k)/wr(k) + b40r(k)*velttr(k)
                   
                    ! Assuming the metric to be (-,+,+,+).
                    b44r(k) = -b40r(k)*b40r(k)+b41r(k)*b41r(k)+b42r(k)*b42r(k)+b43r(k)*b43r(k)

                    fr(k) = fm(k)
           
                 ENDDO


                 ! Equation of state for unknown thermodynamical quantities.
                 IF (eqs .EQ. 2) THEN

                    call EOS_SYNGE (1, tid, fznp, nznp+1, me, mp, rhol, xel, teml, ul,   &
                                    pl, gaml, hl, csl, dpdrhl, dpdul)	
                    call EOS_SYNGE (1, tid, fznp, nznp+1, me, mp, rhor, xer, temr, ur,   &
                                    pr, gamr, hr, csr, dpdrhr, dpdur)

                 ELSE

                    call EOS (1, fznp, nznp+1, gamma_ad, ul, rhol, pl, dpdrhl, dpdul, hl, csl)   
                    call EOS (1, fznp, nznp+1, gamma_ad, ur, rhor, pr, dpdrhr, dpdur, hr, csr) 

                 ENDIF

   
                 !	----------
                 !	GETTING THE EXTREME FAST MAGNETOSONIC EIGENVALUES
                 !	----------

                 DO k = fznp, nznp+1
				
                    ! Coefficients for the left interface


                    CALF2  = B44L(K)/(RHOL(K)*HL(K) + B44L(K))
                    OMEGA2 = CSL(K)**2*(1.D0 - CALF2) + CALF2
                    AUX1   = OMEGA2/WL(K)**2
                    AUX2   = CSL(K)**2/(RHOL(K)*HL(K) + B44L(K))/WL(K)**4

                    CA0    = VELL(K)**4*(1.D0 - OMEGA2) - VELL(K)**2*AUX1 + B41L(K)**2*AUX2
                    CA1    = -4.D0*VELL(K)**3*(1.D0 - OMEGA2) + 2.D0*VELL(K)*AUX1 -         &
                              2.D0*B40L(K)*B41L(K)*AUX2
                    CA2    = 6.D0*VELL(K)**2*(1.D0 - OMEGA2) - (1.D0 - VELL(K)**2)*AUX1 +   &
                             (B40L(K)**2 - B41L(K)**2)*AUX2
                    CA3    = -4.D0*VELL(K)*(1.D0 - OMEGA2) - 2*VELL(K)*AUX1 +               &
                              2.D0*B40L(K)*B41L(K)*AUX2
                    CA4    = 1.D0 - OMEGA2*(VELL(K)**2 + VELtL(K)**2 + VELttL(K)**2) -      &
                             B40L(K)**2*AUX2


                    ! Calculate the extreme fast magnetosonic eigenvalues for the left interface.
                    call FAST (ca0,ca1,ca2,ca3,ca4,lamb1zl,lamb4zl)
		 

                    ! Coefficients for the right interface
	
                    CALF2  = B44R(K)/(RHOR(K)*HR(K) + B44R(K))
                    OMEGA2 = CSR(K)**2*(1.D0 - CALF2) + CALF2
                    AUX1   = OMEGA2/WR(K)**2
                    AUX2   = CSR(K)**2/(RHOR(K)*HR(K) + B44R(K))/WR(K)**4

                    CA0    = VELR(K)**4*(1.D0 - OMEGA2) - VELR(K)**2*AUX1 + B41R(K)**2*AUX2
                    CA1    = -4.D0*VELR(K)**3*(1.D0 - OMEGA2) + 2.D0*VELR(K)*AUX1 -         &
                              2.D0*B40R(K)*B41R(K)*AUX2
                    CA2    = 6.D0*VELR(K)**2*(1.D0 - OMEGA2) - (1.D0 - VELR(K)**2)*AUX1 +   &
                             (B40R(K)**2 - B41R(K)**2)*AUX2
                    CA3    = -4.D0*VELR(K)*(1.D0 - OMEGA2) - 2*VELR(K)*AUX1 +               &
                              2.D0*B40R(K)*B41R(K)*AUX2
                    CA4    = 1.D0 - OMEGA2*(VELR(K)**2 + VELtR(K)**2 + VELttR(K)**2) -      &
                             B40R(K)**2*AUX2

                    ! Calculate the extreme fast magnetosonic eigenvalues for the right interface.
                    call FAST (ca0,ca1,ca2,ca3,ca4,lamb1zr,lamb4zr)

                    ! Upper and lower bounds of the wavespeeds at the interfaces.  
                    lambda = MAX(ABS(lamb4zl),ABS(lamb4zr),ABS(lamb1zl),ABS(lamb1zr)) 

                    ! Maximum 1/tstep.
                    dtestmax = lambda / (xr(k)-xl(k))

                    IF (dtestmax .GT. dtcc_omp(tid)) THEN
                      
                       ict(tid)   = i
                       jct(tid)   = j
                       kct(tid)   = k
                       dtcc_omp(tid) = dtestmax
                       dtt(tid) = dtestmax
                       soundt(tid) = 0.5D0*(csr(k)+csl(k+1))

                    ENDIF

                    dtcc = MAX (dtcc, dtcc_omp(tid))

                 ENDDO    
              ENDDO
           ENDDO
!$OMP END DO

!$OMP BARRIER

!$OMP DO SCHEDULE (static)
           DO j= 1, nyp
              DO i= 1, nxp
                 DO k=fznp, nznp

		    ! THERMAL COOLING
		    ! ..................................................................

                    IF (basenm(1:2) .EQ. 'S4') THEN

                        xe  =  physic(i,j,k)%denstye/physic(i,j,k)%densty
                        tem =  physic(i,j,k)%pres*m_0*cl**2/ (physic(i,j,k)%densty*    &
                               (xe/me + (1.D0-xe)/mp))/kb

                        IF (tem .LT. 1.D4) THEN

                            dedt = 0.D0
                            dtco  = 1.D0

                        ELSE IF (tem .GE. 1.D4 .AND. tem .LT. 1.D5) THEN

                            dedt = (physic(i,j,k)%densty*rho_0a*                     &
                                   (xe/me + (1.D0-xe)/mp))**2/m_0**2 * 7.D-27*tem
                            dedt = dedt/rho_0a/9.D20/3.D10*R_b
                            dtco = 0.1D0*physic(i,j,k)%energy/dedt

                        ELSE IF (tem .GE. 1.D5 .AND. tem .LT. 4.D7) THEN

                            dedt = (physic(i,j,k)%densty*rho_0a*                   &
                                   (xe/me + (1.D0-xe)/mp))**2/m_0**2 *              &
                                    7.D-19*tem**(-0.6D0)
                            dedt = dedt/rho_0a/9.D20/3.D10*R_b
                            dtco = 0.1D0*physic(i,j,k)%energy/dedt

                        ELSE IF (tem .GE. 4.D7) THEN

                            dedt = (physic(i,j,k)%densty*rho_0a*                     &
                                   (xe/me + (1.D0-xe)/mp))**2/m_0**2 *               & 
                                    3.D-27*tem**(0.5D0)
                            dedt = dedt/rho_0a/9.D20/3.D10*R_b
                            dtco = 0.1D0*physic(i,j,k)%energy/dedt

                        ENDIF


                        IF (dtco .LT. dtcmin) THEN
                            dtcmin = dtco
                        ENDIF

                        IF (dtco .LT. dtc_omp(tid)) THEN
                            dtc_omp(tid)  = dtco
                        ENDIF
            
                        dtc = MIN(dtc, dtc_omp(tid))
                 
                    ENDIF

                 ENDDO    
              ENDDO
           ENDDO
!$OMP END DO

	ENDIF	
     
!$OMP END PARALLEL

!!!$OMP MASTER 
         
        DO kk = 1, np
           
           ! Search for the minimum and its position in the array.
           IF (dtt(kk) .EQ. dtcc) THEN
              
              ic = ict(kk)
              jc = jct(kk)
              kc = kct(kk)
	
              v  = DSQRT (physic (ic,jc,kc)%velx*physic(ic,jc,kc)%velx+        &
                          physic (ic,jc,kc)%vely*physic(ic,jc,kc)%vely+        & 
                          physic (ic,jc,kc)%velz*physic(ic,jc,kc)%velz)
              
              sound = soundt(kk)   

	      EXIT
            
           ENDIF

        ENDDO

!!!$OMP END MASTER

!!!$OMP END PARALLEL

	IF (basenm(1:2) .NE. 'S4') dtc = 1.D0

        ! Courant Condition.
        IF (nsdim .EQ. 1) THEN
           dt = DMIN1 (cfl/dtcc, 1.20*olddt)
        ELSE
           dt = MIN(cfl/dtcc/dsqrt(DBLE(nsdim)), 1.02*olddt, dtc)
        ENDIF

 

        ! Global minimum time step.
        IF (mpi_dims(1)*mpi_dims(2)*mpi_dims(3) .GT. 1) THEN
           call REDUCEMINDT2 (ic, jc, kc, dt, sound, v, pmin)
        ENDIF



        ! dt is written on the screen every 'itstep' timesteps (input).
        IF (MOD(nstep, itstp) .EQ. 0) THEN

           IF (mpi_dims(1)*mpi_dims(2)*mpi_dims(3) .GT. 1) THEN

              ! Since all nodes have dtmin information, we can use the first node
              ! to write the results.

              IF (cartrank .EQ. 0) THEN
                 write (6,                                                  &
                 '(i5,2x,i4,2x,e9.2,3(2x,i5),e11.2,2x,e11.2,2x,e11.2)')     &
                  nstep, pmin, dt, ic, jc, kc, sound, v, timej
              ENDIF

           ELSE
              write (6,'(i5,2x,e9.2,3(2x,i5),e11.2,2x,e11.2,2x,e11.2)')     &
                       nstep, dt, ic, jc, kc, sound, v, timej
           ENDIF

        ENDIF

        !--------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE TSTEP_3D
! =============================================================================== 
