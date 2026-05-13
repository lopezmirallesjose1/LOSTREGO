
!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************


! NAME:
!       SWEEPS (2D VERSION)

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
!       Runge-Kutta time evolution with Constrained Transport (CT).

! ACADEMIC REFERENCE:
!       "Efficient Implementation of Essentially Non-oscillatory Shock-
!        Capturing Schemes." C.W.Shu & S.Osher. 1988. JCP 77,439-471
!       "An HLLC Riemann solver for relativistic flows-II.Magnetohydrodynamics."
!        A.Mignone & G.Bodo. 2006. MNRAS. 368, 1040-1054.
!       "Relativistic MHD simulations of extragalctic jets."
!        T.Leismann et al. 2005. A&A. 436, 503-526.
!        See also: "MHD relativista numerica: Aplicaciones en Relatividad Espacial y General."
!        L.A.Ruiz, 2008.


! ===============================================================================        
SUBROUTINE SWEEPS_2D (basenm, limiter, irk, rkin, iars, iccv,                            & 
                      igeomx, igeomy, igeomz, bndmnx, bndmxx, bndmny, bndmxy,            &
                      bndmnz, bndmxz, nx, mnx5, mnx6, nsdim, ny, mny5, mny6,             &
                      eqs, nz, mnz5, mnz6, stelar, cool, gamma_ad, alpha_R, beta_R, dt,  &
                      me, mp, xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,         &
                      boundPosX, boundPosY, boundPosZ, grav, physic, CT, uwi, vwi, shift, Rorb)
! ===============================================================================  


        USE tipos
        USE workarrays
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(4) :: basenm, limiter

        integer(is) :: rkin, iars, iccv, igeomx, igeomy, igeomz,                &
                       bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                       nx, mnx5, nsdim, ny, mny5, eqs, nz, mnz5,  &
                       stelar, cool, mnx6, mny6, mnz6, irk

        real(rs) :: gamma_ad, alpha_R, beta_R, dt, me, mp                 

        ! sizes mnx1, mny1, mnz1
        real(rs) :: xznl(*), xzn(*), xznr(*),   &
                    yznl(*), yzn(*), yznr(*),   &
                    zznl(*), zzn(*), zznr(*)

        type(TboundPosX) :: boundPosX(-4:5)
        type(TboundPosY) :: boundPosY(-4:5)
        type(TboundPosZ) :: boundPosZ(-4:5)

	real(rs) :: uwi, vwi, shift, Rorb

        !------------------------------------------------------ Output Variables

        real(rs)   :: grav(-4:mny5)

        type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)

        !------------------------------------------------------- Local Variables

        integer(is) :: nthreads, curth

        integer(is) :: ii,i, j, k, fznp, nznp, nxp, nyp, nzp,     &
                       bndmin, bndmax, igeom, fzn, nzn, xyzswp

        integer :: tid

# ifdef OPENMP
        integer :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
# endif

        real(rs)   :: ar, br, cr, dr, er, aux, alpha1, alpha2, alpha3, alpha4, dtdx, dtdy, dtdz

        real(rs), pointer :: vel(:), velt(:), veltt(:),        &
                             rho(:), rhoe(:), rhopr(:), p(:),  &
                             b(:), bt(:), btt(:),              &
                             x(:), xl(:), xr(:), dx(:),        &
                             area(:), areal(:), arear(:),      &
                             dvol(:), bs(:)

        real(rs), pointer :: f(:), fict(:)

        real(rs), pointer :: rflx(:), reflx(:), mflx(:), mtflx(:),  &
                             mttflx(:), eflx(:), fflx(:), btflx(:), & 
                             bttflx(:)


        real(rs) :: y, yl, yr, dy,        &
                    z, zl, zr, dz


        real(rs), pointer :: trac(:,:,:)

        !real(rs), pointer :: fbyx(:,:,:), fbzx(:,:,:),        &
        !                     fbxy(:,:,:), fbzy(:,:,:),        &
        !                     fbxz(:,:,:), fbyz(:,:,:)


        !real(rs), pointer :: omega_x(:,:,:), omega_y(:,:,:), omega_z(:,:,:)

        character(32) :: message1

	type(ConsTransp) :: CT(-4:mnx6, -4:mny6, -4:mnz6)
     
	real(rs) :: b2, bct2                        

        real(rs) :: SCRCH1, SCRCH2, SCRCH3, SCRCH4, SCRCH5, SCRCH6, FLUXC
        real(rs) :: SCRCH11, SCRCH12, SCRCH21, SCRCH22, SCRCH41, SCRCH42, SCRCH51, SCRCH52
        
        !-------------------------------------------------------------------------- BEGIN

        nthreads = omp_th
	
        ! We define two sets of magnetic field vectors:

        !    1. The staggered interface magnetic field, B.
        !    2. The zone centered magnetic field Bc (xi,yj,zk) defined as:

        !       Bxc(i) = 0.5*(Bxi + Bxi+1)
        !       Byc(j) = 0.5*(Byj + Byj+1)
        !       Bzc(k) = 0.5*(Bzk + Bzk+1)

        ! The zone centered vector is computed at the end of the sweep
        ! by a simple interpolation from the staggered field components (CT).

        ! The zone centered vector Bc is required to set the boundary conditions
        ! and to calculate the source terms and the fluxes of the other variables.



!$OMP PARALLEL DEFAULT(shared)                                                     & 
!$OMP PRIVATE(curth,tid,i,j,k,fznp,nznp,nxp,nyp,nzp,bndmin,bndmax,igeom,fzn,nzn)   &
!$OMP PRIVATE(ar,br,cr,aux,xyzswp)                                                 &
!$OMP PRIVATE(vel,velt,veltt,rho,rhoe,rhopr,p,b,bt,btt,bs,x,xl,xr,dx)              &
!$OMP PRIVATE(y,yl,yr,dy,z,zl,zr,dz)                                               &
!$OMP PRIVATE(area,areal,arear,dvol,f,fict)                                        &
!$OMP PRIVATE(rflx,reflx,mflx,mtflx,mttflx,eflx,fflx,btflx,bttflx,trac)              


# ifdef OPENMP
        call my_omp_get_thread_num (curth)
        tid = curth +1
# else
        tid = 1
# endif


        ! -4:mn6 vectors
        rflx   => wa(tid)%wa6sweeps100
        reflx  => wa(tid)%wa6sweeps101
        mflx   => wa(tid)%wa6sweeps102
        mtflx  => wa(tid)%wa6sweeps103
        mttflx => wa(tid)%wa6sweeps104
        eflx   => wa(tid)%wa6sweeps105
        fflx   => wa(tid)%wa6sweeps106
        btflx  => wa(tid)%wa6sweeps107
        bttflx => wa(tid)%wa6sweeps108


        ! -4:mn5 vectors
        vel   => wa(tid)%wa5sweeps100
        velt  => wa(tid)%wa5sweeps101
        veltt => wa(tid)%wa5sweeps102 
        rho   => wa(tid)%wa5sweeps103
        rhoe  => wa(tid)%wa5sweeps104
        rhopr => wa(tid)%wa5sweeps105
        p     => wa(tid)%wa5sweeps106
        x     => wa(tid)%wa5sweeps107
        xl    => wa(tid)%wa5sweeps108
        xr    => wa(tid)%wa5sweeps109
        dx    => wa(tid)%wa5sweeps110
        area  => wa(tid)%wa5sweeps111
        areal => wa(tid)%wa5sweeps112
        arear => wa(tid)%wa5sweeps113
        dvol  => wa(tid)%wa5sweeps114

        f    => wa(tid)%wa5sweeps115
        fict => wa(tid)%wa5sweeps116

        b    => wa(tid)%wa5sweeps117
        bt   => wa(tid)%wa5sweeps118
        btt  => wa(tid)%wa5sweeps119
	bs   => wa(tid)%wa5sweeps120

        ! Tridimensional working arrays
          
        trac =>wa(tid)%wa3sweeps100
        
  
       ! =============================
       !           X-SWEEP
       ! =============================

       ! Runge-Kutta parameters.
       ! rkin goes from 1 to irk (Runge-Kutta order).
       ! Runge-Kutta main loop is defined in the main program.


       IF (irk .LE. 3) THEN	

          IF (rkin .EQ. 1) THEN

             ar = 1.D0
             br = 0.D0
             cr = 1.D0

          ELSE IF (rkin .EQ. 2) THEN

             ar = beta_R/alpha_R
             br = 1.D0/alpha_R
             cr = br

          ELSE

             ar = 1.D0/3.D0
             br = 2.D0/3.D0
             cr = br

          ENDIF
	   
       ENDIF

       ! XSweep own parameters.
       xyzswp  = 1
       bndmin = bndmnx
       bndmax = bndmxx
       igeom  = igeomx
       fzn    = 1
       nzn    = nx
       fznp   = fzn
       nznp   = nzn
       nyp    = ny
       nzp    = nz


       ! Take the primitive variables at the cell centers.

       DO k= 1, nzp

!$OMP DO SCHEDULE(static)
          DO j= 1, nyp

	     DO i=fznp, nznp+1
		bs(i) = CT(i,j,k)%bx
	     ENDDO

             DO i= fznp, nznp
		
                vel(i)   = physic(i,j,k)%velx
                velt(i)  = physic(i,j,k)%vely
                veltt(i) = physic(i,j,k)%velz 

                b(i)     = physic(i,j,k)%bxc2
                bt(i)    = physic(i,j,k)%byc2
                btt(i)   = physic(i,j,k)%bzc2
                      
                rho(i)  = physic(i,j,k)%densty
                rhoe(i) = physic(i,j,k)%denstye
                rhopr(i)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                p(i)    = physic(i,j,k)%pres
                   
                trac(i,j,k) = physic(i,j,k)%tracer2
                f(i) = physic(i,j,k)%tracer2


                xl(i)   = xznl(i)
                x(i)    = xzn(i)
                xr(i)   = xznr(i)
                dx(i)   = xr(i)-xl(i)
	       
                fict(i) = 0.0
                
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




		   IF (mpi_coords(1) .EQ. 0) THEN

                      DO i= fznp-5, fznp-1

                         physic(i,j,k)%velx  =  vel(i)
                         physic(i,j,k)%vely  =  velt(i)
                         physic(i,j,k)%velz  =  veltt(i)

                         physic(i,j,k)%bxc2   =  b(i)
                         physic(i,j,k)%byc2   =  bt(i)
                         physic(i,j,k)%bzc2   =  btt(i)

                         physic(i,j,k)%densty  = rho(i)
                         physic(i,j,k)%denstye = rhoe(i)
                         physic(i,j,k)%pres    = p(i)

                         physic(i,j,k)%tracer  = f(i)

		      ENDDO

		   ENDIF


		   IF (mpi_coords(1) .EQ. mpi_dims(1)-1) THEN

                      DO i= nznp+1, nznp+5

                         physic(i,j,k)%velx  =  vel(i)
                         physic(i,j,k)%vely  =  velt(i)
                         physic(i,j,k)%velz  =  veltt(i)

                         physic(i,j,k)%bxc2   =  b(i)
                         physic(i,j,k)%byc2   =  bt(i)
                         physic(i,j,k)%bzc2   =  btt(i)

                         physic(i,j,k)%densty  = rho(i)
                         physic(i,j,k)%denstye = rhoe(i)
                         physic(i,j,k)%pres    = p(i)

                         physic(i,j,k)%tracer  = f(i)

		      ENDDO

		   ENDIF


                   
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

                   fict(i) = 0.0
                         
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

                   fict(i) = 0.0
                         
                ENDDO
                   
             ENDIF

             ! Calculate the sections of the cells in the normal direction.

             call GEOM (fznp, nznp, igeom, dx, area, areal, arear, dvol)  
		
             ii=j

                yl   = yznl(j)
                y    = yzn(j)
                yr   = yznr(j)
                dy   = yr-yl

                zl   = zznl(k)
                z    = zzn(k)
                zr   = zznr(k)
                dz   = zr-zl

             ! Calculate numerical fluxes.
             call MHYDROW3 (ii,tid,xyzswp,basenm, limiter, fznp, nznp, iars,            &
                            eqs, gamma_ad, me, mp, x, xl, xr, areal,                    &
                            p, rho, rhoe, rhopr, vel, velt, veltt, f, b, bs, bt, btt,   &
                            rflx, reflx, mflx, mtflx, mttflx,                           &
                            eflx, btflx, bttflx, fflx,y)
            
             !========================
             ! Runge-Kutta substep 
             !========================

             ! Advance of conserved variables.

             DO i= fznp, nznp

                aux = dt / dx(i)
		

                physic(i,j,k)%tracer = ar * physic(i,j,k)%odenstw     &
                                          * physic(i,j,k)%otracer     &
                                     + br * physic(i,j,k)%denstw      &
                                          * physic(i,j,k)%tracer      &
                                     - cr * (fflx(i+1)-fflx(i))*aux

                physic(i,j,k)%denstw = ar * physic(i,j,k)%odenstw     &
                                     + br * physic(i,j,k)%denstw      &
                                     - cr * (rflx(i+1)-rflx(i))*aux

                physic(i,j,k)%denstew = ar * physic(i,j,k)%odenstew   &
                                      + br * physic(i,j,k)%denstew    &
                                      - cr * (reflx(i+1)-reflx(i))*aux   

                physic(i,j,k)%momenx = ar * physic(i,j,k)%omomenx     &
                                     + br * physic(i,j,k)%momenx      &
                                     - cr * (mflx(i+1)-mflx(i))*aux 

                physic(i,j,k)%momeny = ar * physic(i,j,k)%omomeny     &
                                     + br * physic(i,j,k)%momeny      &
                                     - cr * (mtflx(i+1)-mtflx(i))*aux 

                physic(i,j,k)%momenz = ar * physic(i,j,k)%omomenz     &
                                     + br * physic(i,j,k)%momenz      &
                                     - cr * (mttflx(i+1)-mttflx(i))*aux 

                physic(i,j,k)%energy = ar * physic(i,j,k)%oenergy     &
                                     + br * physic(i,j,k)%energy      &
                                     - cr * (eflx(i+1)-eflx(i))*aux 

 
                ! Magnetic field.                     
                                     
                IF (nsdim .EQ. 2) THEN

                   physic(i,j,k)%bzc    = ar * physic(i,j,k)%obzc        &
                                        + br * physic(i,j,k)%bzc         &
                                        - cr * (bttflx(i+1)-bttflx(i))*aux 
	           

      	           IF (iccv .NE. 0) THEN


                      physic(i,j,k)%bxc    = ar * physic(i,j,k)%obxc        &
                                           + br * physic(i,j,k)%bxc         
                                          

                      physic(i,j,k)%byc    = ar * physic(i,j,k)%obyc        &
                                           + br * physic(i,j,k)%byc         &
                                           - cr * (btflx(i+1)-btflx(i))*aux


	      	   ENDIF


                ELSE IF (nsdim .EQ. 1) THEN

	           physic(i,j,k)%bxc    = ar * physic(i,j,k)%obxc        &
                                        + br * physic(i,j,k)%bxc

                   physic(i,j,k)%bzc    = ar * physic(i,j,k)%obzc        &
                                        + br * physic(i,j,k)%bzc         &
                                        - cr * (bttflx(i+1)-bttflx(i))*aux

                   physic(i,j,k)%byc    = ar * physic(i,j,k)%obyc        &
                                        + br * physic(i,j,k)%byc         &
                                        - cr * (btflx(i+1)-btflx(i))*aux            


                ENDIF


           
                ! Normalizing Tracer for one-dimensional simulations...

                IF (nsdim .EQ. 1) THEN
                   physic(i,j,k)%tracer = physic(i,j,k)%tracer / physic(i,j,k)%denstw
                ENDIF


             ENDDO 

             ! CT-fluxes.
             IF (nsdim .GT. 1) THEN

                DO i = fznp, nznp + 1

                   CT(i,j,k)%fbyx = btflx(i)
                   CT(i,j,k)%fbzx = bttflx(i)
                   CT(i,j,k)%rflxx = rflx(i)
                   
                ENDDO

             ENDIF

	    
          ENDDO
!$OMP END DO 

       ENDDO

       ! Synchronization: Each thread waits for all the other threads at the end of the sweep.
!$OMP BARRIER

       ! =============================
       !           Y-SWEEP
       ! =============================

       ! Runge-Kutta parameters.
       ! rkin goes from 1 to irk (Runge-Kutta order).
       ! Runge-Kutta Main Loop is defined in the main program.

       IF (irk .LE. 3) THEN


          IF (rkin .EQ. 1) THEN

	     ar = 0.D0
             br = 1.D0
             cr = 1.D0

          ELSE IF (rkin .EQ. 2) THEN

	     ar = 0.D0
             br = 1.D0
             cr = 1.D0/alpha_R

          ELSE

	     ar = 0.D0
             br = 1.D0
             cr = 2.D0/3.D0

          ENDIF


       ENDIF

       IF (nsdim .GT. 1) THEN

          ! YSweep own parameters.

          xyzswp  = 2
          bndmin = bndmny
          bndmax = bndmxy
          igeom  = igeomy
          fzn    = 1
          nzn    = ny
          fznp   = fzn
          nznp   = nzn
          nxp    = nx
          nzp    = nz
 
          ! Take the primitive variables at the cell centers. 
 

          DO k= 1, nzp

!$OMP DO SCHEDULE(static)
             DO i=1,nxp

	        DO j = fznp, nznp+1
                   bs(j) = CT(i,j,k)%by
                ENDDO

                DO j= fznp, nznp

                   vel(j)  = physic(i,j,k)%vely
                   velt(j) = physic(i,j,k)%velz
                   veltt(j)= physic(i,j,k)%velx


                   b(j)    = physic(i,j,k)%byc2
                   bt(j)   = physic(i,j,k)%bzc2
                   btt(j)  = physic(i,j,k)%bxc2

                   
                   rho(j)  = physic(i,j,k)%densty
                   rhoe(j) = physic(i,j,k)%denstye
                   rhopr(j)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                   p(j)    = physic(i,j,k)%pres

                   f(j) = physic(i,j,k)%tracer2
                   
                   xl(j)   = yznl(j)
                   x(j)    = yzn(j)
                   xr(j)   = yznr(j)
                   dx(j)   = xr(j)-xl(j)
                   fict(j) = 0.0

                ENDDO


                ! We need to provide boundary conditions to the first and last nodes. 

                IF (mpi_coords(2) .EQ. 0 .OR. mpi_coords(2) .EQ. (mpi_dims(2)-1)) THEN
 
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

                      DO j= fznp-5, fznp-1

                         physic(i,j,k)%vely  =  vel(j)
                         physic(i,j,k)%velz  =  velt(j)
                         physic(i,j,k)%velx  =  veltt(j)

                         physic(i,j,k)%byc2   =  b(j)
                         physic(i,j,k)%bzc2   =  bt(j)
                         physic(i,j,k)%bxc2   =  btt(j)

                         physic(i,j,k)%densty  = rho(j)
                         physic(i,j,k)%denstye = rhoe(j)
                         physic(i,j,k)%pres    = p(j)

                         physic(i,j,k)%tracer  = f(j)

		      ENDDO

		   ENDIF


		   IF (mpi_coords(2) .EQ. mpi_dims(2)-1) THEN

                      DO j= nznp+1, nznp+5

                         physic(i,j,k)%vely  =  vel(j)
                         physic(i,j,k)%velz  =  velt(j)
                         physic(i,j,k)%velx  =  veltt(j)

                         physic(i,j,k)%byc2   =  b(j)
                         physic(i,j,k)%bzc2   =  bt(j)
                         physic(i,j,k)%bxc2   =  btt(j)

                         physic(i,j,k)%densty  = rho(j)
                         physic(i,j,k)%denstye = rhoe(j)
                         physic(i,j,k)%pres    = p(j)

                         physic(i,j,k)%tracer  = f(j)

		      ENDDO

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

                      fict(j) = 0.0

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

                      fict(j) = 0.0

                   ENDDO

                ENDIF   


                ! Calculate the sections of the cells in the normal direction.

                call GEOM (fznp, nznp, igeom, dx, area, areal, arear, dvol) 
                
                ! Calculate numerical fluxes.

	        ii=i

                call MHYDROW3 (ii,tid, xyzswp,basenm,limiter, fznp, nznp, iars,            &
                              eqs, gamma_ad, me, mp, x, xl, xr, areal,                     &
                              p, rho, rhoe, rhopr, vel, velt, veltt, f, b, bs, bt, btt,    &
                              rflx, reflx, mflx, mtflx, mttflx,                            &
                              eflx, btflx, bttflx, fflx,y) 

               
                !========================
                ! Runge-Kutta substep 
                !========================

                ! Advance of conserved variables.
                   

                DO j= fznp, nznp

                   aux = dt / dx(j)

                   physic(i,j,k)%tracer = ar * physic(i,j,k)%otracer     &
                                        + br * physic(i,j,k)%tracer      &
                                        - cr * (fflx(j+1)-fflx(j))*aux

                   physic(i,j,k)%denstw = ar * physic(i,j,k)%odenstw     &
                                        + br * physic(i,j,k)%denstw      &
                                        - cr * (rflx(j+1)-rflx(j))*aux

                   physic(i,j,k)%denstew = ar * physic(i,j,k)%odenstew   &
                                         + br * physic(i,j,k)%denstew    &
                                         - cr * (reflx(j+1)-reflx(j))*aux   

                   physic(i,j,k)%momenx = ar * physic(i,j,k)%omomenx     &
                                        + br * physic(i,j,k)%momenx      &
                                        - cr * (mttflx(j+1)-mttflx(j))*aux 

                   physic(i,j,k)%momeny = ar * physic(i,j,k)%omomeny     &
                                        + br * physic(i,j,k)%momeny      &
                                        - cr * (mflx(j+1)-mflx(j))*aux 

                   physic(i,j,k)%momenz = ar * physic(i,j,k)%omomenz     &
                                        + br * physic(i,j,k)%momenz      &
                                        - cr * (mtflx(j+1)-mtflx(j))*aux 

                   physic(i,j,k)%energy = ar * physic(i,j,k)%oenergy     &
                                        + br * physic(i,j,k)%energy      &
                                        - cr * (eflx(j+1)-eflx(j))*aux


                 
                   ! Magnetic field.


		   IF (iccv .NE. 0) THEN

                      physic(i,j,k)%bxc    = br * physic(i,j,k)%bxc         &
                                           - cr * (bttflx(j+1)-bttflx(j))*aux

                      physic(i,j,k)%byc    = br * physic(i,j,k)%byc

                                           
		   ENDIF 

                 
                   IF (nsdim .NE. 3) THEN

                      physic(i,j,k)%bzc    = br * physic(i,j,k)%bzc         &
                                           - cr * (btflx(j+1)-btflx(j))*aux
			
                   ENDIF


                   ! Normalizing Tracer for two-dimensional simulations...

                   IF (nsdim .EQ. 2) THEN

                      physic(i,j,k)%tracer = physic(i,j,k)%tracer/    &
                                             physic(i,j,k)%denstw


                   ENDIF

                ENDDO

                !CT- fluxes.
                
                 DO j=fznp, nznp+1

                   CT(i,j,k)%fbxy = bttflx(j)
                   CT(i,j,k)%fbzy = btflx(j)
                   CT(i,j,k)%rflxy= rflx(j)
                   
                ENDDO
             ENDDO

!$OMP END DO
          ENDDO


       ENDIF

!$OMP END PARALLEL


       ! ===================================
       !          CONSTRAINED TRANSPORT
       ! ===================================

       !!! Evolution of the magnetic field for 2D and 3D simulations.    

       !                                                      |
       !                                                      |
       !                                                      |
       !                                                      |
       !                                      |       O       |
       !                                      |               |
       !        -------------/-|--------------|-------------/-|------------------------------/|--------
       !                   /   |      ^       |           /   |                            /  |
       !                 /     |      |       | Ωz,i+1/2j+1/2 |                          /    |
       !               /--------------O-------|-------O       |        O               /      |
       !             / |       |   By,j+1/2   |     / |       ^                      /        | 
       !           /   |       |              |   /   |       |                    /          |
       !         /     |       |              | /     |       |                  /            | 
       !       /       o---> Bx,i-1/2         /       O---> Bx,i+1/2           /              | 
       !  ----|--------|--------------x-------|-------|-------|--------x------|               |
       !      |    Fxi-1/2,j  /-----i,j,k-----|       |       /-----i+1,j,k---|-------------- /-------
       !      |     ---|--> /         ^       |       |     / |               |             /
       !      |        |  /   ^       |       v       |   /   |               |           /
       !      |        |/-----|-------o-------|-------| /     |               |         /
       !      |       /       |    By,j-1/2   |       /       |               |       /
       !      |     /     Fyi,j-1/2           |     /         |               |     /              ^ y
       !      |   /                           |   /           |               |   /                |
       !      | /                             | /             |               | /                  |
       !   ---|-------------------------------|---------------|---------------|--------------------o------->      
       !                                      |               |                                  /         x
       !                                      |               |                              z /
       !                                      |
       !                                      |
       !                                      |
       !                                      |


       !!! Nomenclature:

       ! fbxy: the flux of Bx in the y-direction. 
       ! fbxz: the flux of Bx in the z-direction. 

       ! fbyx: the flux of By in the x-direction.
       ! fbyz: the flux of By in the z-direction.

       ! fbzx: the flux of Bz in the x-direction.
       ! fbzy: the flux of Bz in the y-direction.


       IF (nsdim .EQ. 2) THEN

          !!! Boundary conditions for the staggered fluxes.

          ! ------------
          ! LEFT X-BOUNDARY
          ! ------------

          ! Reflecting Boundary.

          IF (bndmnx .EQ. 1) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

             DO  k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO  j=1,ny+1
                   CT(0,j,k)%fbxy = - CT(1,j,k)%fbxy
                ENDDO
!$OMP END DO
             ENDDO
!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmnx .EQ. 2 .OR. bndmnx .EQ. 4) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO j=1,ny+1
                   CT(0,j,k)%fbxy = CT(1,j,k)%fbxy
                ENDDO
!$OMP END DO
             ENDDO


             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO j=1,ny+1
                   CT(0,j,k)%rflxy = CT(1,j,k)%rflxy
                ENDDO
!$OMP END DO
             ENDDO
             
!$OMP END PARALLEL

          ! Periodic Boundary.

          ELSE IF (bndmnx .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO j=1,ny+1
                   CT(0,j,k)%fbxy = CT(nx,j,k)%fbxy
                ENDDO
!$OMP END DO
             ENDDO

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO j=1,ny+1
                   CT(0,j,k)%rflxy = CT(nx,j,k)%rflxy
                ENDDO
!$OMP END DO
             ENDDO


!$OMP END PARALLEL


          ELSE

             message1 = 'Error in boundaries of B fluxes.'
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)

          ENDIF

          ! ------------
          ! RIGHT X-BOUNDARY
          ! ------------


          ! Reflecting Boundary.

          IF (bndmxx .EQ. 1) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

             DO  k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO  j=1,ny+1
                   CT(nx+1,j,k)%fbxy = - CT(nx,j,k)%fbxy
                ENDDO
!$OMP END DO
             ENDDO
!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmxx .EQ. 2 .OR. bndmxx .EQ. 9) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO j=1,ny+1
                   CT(nx+1,j,k)%fbxy = CT(nx,j,k)%fbxy
                ENDDO
!$OMP END DO
             ENDDO

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO j=1,ny+1
                   CT(nx+1,j,k)%rflxy = CT(nx,j,k)%rflxy
                ENDDO
!$OMP END DO
             ENDDO

             
!$OMP END PARALLEL

          ! Periodic Boundary.

          ELSE IF (bndmxx .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO j=1,ny+1
                   CT(nx+1,j,k)%fbxy = CT(1,j,k)%fbxy
                ENDDO
!$OMP END DO
             ENDDO

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO j=1,ny+1
                   CT(nx+1,j,k)%rflxy = CT(1,j,k)%rflxy
                ENDDO
!$OMP END DO
             ENDDO



!$OMP END PARALLEL

          ELSE
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)
          ENDIF


          ! ------------
          ! LEFT Y-BOUNDARY
          ! ------------

          ! Reflecting Boundary.

          IF (bndmny .EQ. 1) THEN
 
!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

             DO  k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO  i=1,nx+1
                   CT(i,0,k)%fbyx = - CT(i,1,k)%fbyx
                ENDDO
!$OMP END DO
             ENDDO
!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmny .EQ. 2 .OR. bndmny .EQ. 7) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)

                DO i=1,nx+1
                   CT(i,0,k)%fbyx = CT(i,1,k)%fbyx
                ENDDO
!$OMP END DO
             ENDDO

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)

                DO i=1,nx+1
                   CT(i,0,k)%rflxx = CT(i,1,k)%rflxx
                ENDDO
!$OMP END DO
             ENDDO
             
!$OMP END PARALLEL

         ! Periodic Boundary.

          ELSE IF (bndmny .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)

                DO i=1,nx+1
                   CT(i,0,k)%fbyx = CT(i,ny,k)%fbyx
                ENDDO
!$OMP END DO
             ENDDO


             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)

                DO i=1,nx+1
                   CT(i,0,k)%rflxx = CT(i,ny,k)%rflxx
                ENDDO
!$OMP END DO
             ENDDO


!$OMP END PARALLEL


          ELSE

             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)


          ENDIF

          ! ------------
          ! RIGHT Y-BOUNDARY
          ! ------------


          ! Reflecting Boundary.

          IF (bndmxy .EQ. 1) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

             DO  k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO  i=1,nx+1
                   CT(i,ny+1,nz)%fbyx = - CT(i,ny,nz)%fbyx
                ENDDO
!$OMP END DO
             ENDDO
!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmxy .EQ. 2) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO i=1,nx+1
                   CT(i,ny+1,k)%fbyx = CT(i,ny,k)%fbyx
                ENDDO
!$OMP END DO
             ENDDO

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO i=1,nx+1
                   CT(i,ny+1,k)%rflxx = CT(i,ny,k)%rflxx
                ENDDO
!$OMP END DO
             ENDDO
             
!$OMP END PARALLEL

          ! Periodic Boundary.

          ELSE IF (bndmxy .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO i=1,nx+1
                   CT(i,ny+1,k)%fbyx = CT(i,1,k)%fbyx
                ENDDO
!$OMP END DO
             ENDDO

             DO k=1,nz
!$OMP DO SCHEDULE (STATIC)
                DO i=1,nx+1
                   CT(i,ny+1,k)%rflxx = CT(i,1,k)%rflxx
                ENDDO
!$OMP END DO
             ENDDO

!$OMP END PARALLEL


          ELSE
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)
          ENDIF

       ENDIF


       IF (nsdim .EQ. 2) THEN

!$OMP PARALLEL DEFAULT(shared)         &
!$OMP PRIVATE(i,j,k,ar,br,cr)          &
!$OMP PRIVATE(FLUXC, SCRCH1, SCRCH2, SCRCH3, SCRCH4, SCRCH5, SCRCH6) &
!$OMP PRIVATE(SCRCH11, SCRCH12, SCRCH21, SCRCH22, SCRCH41, SCRCH42, SCRCH51, SCRCH52)




       IF (rkin .EQ. 1) THEN

          ar = 1.D0
          br = 0.D0
          cr = 1.D0

       ELSE IF (rkin .EQ. 2) THEN

          ar = beta_R/alpha_R
          br = 1.D0/alpha_R
          cr = br

       ELSE

          ar = 1.D0/3.D0
          br = 2.D0/3.D0
          cr = br

       ENDIF


       physic(0,0,1)%velx = (physic(0,1,1)%velx+physic(1,0,1)%velx)/2.D0
       physic(0,0,1)%vely = (physic(0,1,1)%vely+physic(1,0,1)%vely)/2.D0
       physic(0,0,1)%bxc2 = (physic(0,1,1)%bxc2+physic(1,0,1)%bxc2)/2.D0
       physic(0,0,1)%byc2 = (physic(0,1,1)%byc2+physic(1,0,1)%byc2)/2.D0
       
       physic(nx+1,0,1)%velx = (physic(nx,0,1)%velx+physic(nx+1,1,1)%velx)/2.D0
       physic(nx+1,0,1)%vely = (physic(nx,0,1)%vely+physic(nx+1,1,1)%vely)/2.D0
       physic(nx+1,0,1)%bxc2 = (physic(nx,0,1)%bxc2+physic(nx+1,1,1)%bxc2)/2.D0
       physic(nx+1,0,1)%byc2 = (physic(nx,0,1)%byc2+physic(nx+1,1,1)%byc2)/2.D0

       physic(0,ny+1,1)%velx = (physic(0,ny,1)%velx+physic(1,ny+1,1)%velx)/2.D0
       physic(0,ny+1,1)%vely = (physic(0,ny,1)%vely+physic(1,ny+1,1)%vely)/2.D0
       physic(0,ny+1,1)%bxc2 = (physic(0,ny,1)%bxc2+physic(1,ny+1,1)%bxc2)/2.D0
       physic(0,ny+1,1)%byc2 = (physic(0,ny,1)%byc2+physic(1,ny+1,1)%byc2)/2.D0

       physic(nx+1,ny+1,1)%velx = (physic(nx,ny+1,1)%velx+physic(nx+1,ny,1)%velx)/2.D0
       physic(nx+1,ny+1,1)%vely = (physic(nx,ny+1,1)%vely+physic(nx+1,ny,1)%vely)/2.D0
       physic(nx+1,ny+1,1)%bxc2 = (physic(nx,ny+1,1)%bxc2+physic(nx+1,ny,1)%bxc2)/2.D0
       physic(nx+1,ny+1,1)%byc2 = (physic(nx,ny+1,1)%byc2+physic(nx+1,ny,1)%byc2)/2.D0

       
          !!! Electromotive forces Ω in the corners.

          DO k=1, nz
!$OMP DO SCHEDULE(static)
             DO j= 1,  ny+1
                DO i= 1, nx+1


                   
                   CT(i,j,k)%omega_z = 0.25D0*(CT(i,j-1,k)%fbyx + CT(i,j,k)%fbyx - &
                                               CT(i-1,j,k)%fbxy-CT(i,j,k)%fbxy)



                   ! ====================
		   ! CONTACT-CT (UPWIND)
                   ! ====================


		   ! SCRCH1;
		   ! ...................................		

		   IF (CT(I,J-1,K)%RFLXX .GT. 1.d-15) THEN

                        FLUXC  = -(physic(i-1,j-1,k)%vely*physic(i-1,j-1,k)%bxc2-physic(i-1,j-1,k)%byc2*physic(i-1,j-1,k)%velx)
			SCRCH1 = 2.D0/(YZNR(J)-YZNL(J))*(-CT(I-1,J,K)%fbxy-FLUXC)

		   ELSE IF (CT(I,J-1,K)%RFLXX .LT. -1.D-15) THEN

                        FLUXC  = -(physic(i,j-1,k)%vely*physic(i,j-1,k)%bxc2-physic(i,j-1,k)%byc2*physic(i,j-1,k)%velx)                
                        SCRCH1 = 2.D0/(YZNR(J)-YZNL(J))*(-CT(I,J,K)%fbxy-FLUXC)

	           ELSE

                      SCRCH11 = 2.D0/(YZNR(J)-YZNL(J))*(-CT(I-1,J,K)%fbxy+(physic(i-1,j-1,k)%vely*physic(i-1,j-1,k)%bxc2-  &
                                                                           physic(i-1,j-1,k)%byc2*physic(i-1,j-1,k)%velx))	
                      SCRCH12 = 2.D0/(YZNR(J)-YZNL(J))*(-CT(I,J,K)%fbxy+(physic(i,j-1,k)%vely*physic(i,j-1,k)%bxc2-        &
                                                                         physic(i,j-1,k)%byc2*physic(i,j-1,k)%velx))
		      SCRCH1 = 0.5D0*(SCRCH11+SCRCH12)			

		   ENDIF


                   ! SCRCH2;
                   ! ...................................


		   IF (CT(I,J-1,K)%RFLXX .GT. 1.D-15) THEN

                        FLUXC  = -(physic(i-1,j,k)%vely*physic(i-1,j,k)%bxc2-physic(i-1,j,k)%byc2*physic(i-1,j,k)%velx)
			SCRCH2 = 2.D0/(YZNR(J)-YZNL(J))*(FLUXC+CT(I-1,J,K)%fbxy)

		   ELSE IF (CT(I,J-1,K)%RFLXX .LT. -1.D-15) THEN

                        FLUXC  = -(physic(i,j,k)%vely*physic(i,j,k)%bxc2-physic(i,j,k)%byc2*physic(i,j,k)%velx)                      
                        SCRCH2 = 2.D0/(YZNR(J)-YZNL(J))*(FLUXC+CT(I,J,K)%fbxy)

	           ELSE

                        SCRCH21 = 2.D0/(YZNR(J)-YZNL(J))*(-(physic(i-1,j,k)%vely*physic(i-1,j,k)%bxc2-                &
                                                            physic(i-1,j,k)%byc2*physic(i-1,j,k)%velx)+CT(I-1,J,K)%fbxy)	
                        SCRCH22 = 2.D0/(YZNR(J)-YZNL(J))*(-(physic(i,j,k)%vely*physic(i,j,k)%bxc2-                    &
                                                            physic(i,j,k)%byc2*physic(i,j,k)%velx)+CT(I,J,K)%fbxy)
			SCRCH2 = 0.5D0*(SCRCH21+SCRCH22)			

		   ENDIF


		   ! SCRCH4;
		   ! ...................................		

		   IF (CT(I-1,J,K)%RFLXY .GT. 1.D-15) THEN

                        FLUXC  = (physic(i-1,j-1,k)%velx*physic(i-1,j-1,k)%byc2-physic(i-1,j-1,k)%bxc2*physic(i-1,j-1,k)%vely)
			SCRCH4 = 2.D0/(XZNR(I)-XZNL(I))*(CT(I,J-1,K)%fbyx-FLUXC)

		   ELSE IF (CT(I-1,J,K)%RFLXY .LT. -1.D-15) THEN

                        FLUXC  = (physic(i-1,j,k)%velx*physic(i-1,j,k)%byc2-physic(i-1,j,k)%bxc2*physic(i-1,j,k)%vely)                
                        SCRCH4 = 2.D0/(XZNR(I)-XZNL(I))*(CT(I,J,K)%fbyx-FLUXC)

	           ELSE

                      SCRCH41 = 2.D0/(XZNR(I)-XZNL(I))*(CT(I,J-1,K)%fbyx-(physic(i-1,j-1,k)%velx*physic(i-1,j-1,k)%byc2-  &
                                                                           physic(i-1,j-1,k)%bxc2*physic(i-1,j-1,k)%vely))	
                      SCRCH42 = 2.D0/(XZNR(I)-XZNL(I))*(CT(I,J,K)%fbyx-(physic(i-1,j,k)%velx*physic(i-1,j,k)%byc2-        &
                                                                        physic(i-1,j,k)%bxc2*physic(i-1,j,k)%vely))
		      SCRCH4 = 0.5D0*(SCRCH41+SCRCH42)			

		   ENDIF


                   ! SCRCH5;
                   ! ...................................


		   IF (CT(I-1,J,K)%RFLXY .GT. 1.D-15) THEN

                        FLUXC  = (physic(i,j-1,k)%velx*physic(i,j-1,k)%byc2-physic(i,j-1,k)%bxc2*physic(i,j-1,k)%vely)
			SCRCH5 = 2.D0/(XZNR(I)-XZNL(I))*(FLUXC-CT(I,J-1,K)%fbyx)

		   ELSE IF (CT(I-1,J,K)%RFLXY .LT. -1.D-15) THEN

                        FLUXC  = (physic(i,j,k)%velx*physic(i,j,k)%byc2-physic(i,j,k)%bxc2*physic(i,j,k)%vely)                      
                        SCRCH5 = 2.D0/(XZNR(I)-XZNL(I))*(FLUXC-CT(I,J,K)%fbyx)

	           ELSE

                        SCRCH51 = 2.D0/(XZNR(I)-XZNL(I))*((physic(i,j-1,k)%velx*physic(i,j-1,k)%byc2-                &
                                                            physic(i,j-1,k)%bxc2*physic(i,j-1,k)%vely)-CT(I,J-1,K)%fbyx)	
                        SCRCH52 = 2.D0/(XZNR(I)-XZNL(I))*((physic(i,j,k)%velx*physic(i,j,k)%byc2-                    &
                                                            physic(i,j,k)%bxc2*physic(i,j,k)%vely)-CT(I,J,K)%fbyx)
			SCRCH5 = 0.5D0*(SCRCH51+SCRCH52)			

		   ENDIF



                   
		   SCRCH3 = SCRCH1 - SCRCH2		
		   SCRCH6 = SCRCH4 - SCRCH5

		   CT(i,j,k)%omega_z = CT(i,j,k)%omega_z + SCRCH3*(YZNR(J)-YZNL(J))/8.D0 + &
                                                           SCRCH6*(XZNR(I)-XZNL(I))/8.D0
                   
                   
                ENDDO
             ENDDO
!$OMP END DO
          ENDDO


          !!! CT magnetic field evolution.
     
          DO k= 1, nz
!$OMP DO SCHEDULE(static)
             DO j= 1, ny
                DO i=1, nx

                   CT(i,j,k)%bx = ar*CT(i,j,k)%obx + br*CT(i,j,k)%bx   &
                   + cr*dt*(CT(i,j+1,k)%omega_z-CT(i,j,k)%omega_z)/(yznr(j)-yznl(j))

                   CT(i,j,k)%by = ar*CT(i,j,k)%oby + br*CT(i,j,k)%by   &
                   + cr*dt*(CT(i,j,k)%omega_z-CT(i+1,j,k)%omega_z)/(xznr(i)-xznl(i))

                ENDDO
             ENDDO
!$OMP END DO
          ENDDO


          DO k= 1, nz
!$OMP DO SCHEDULE(static)
             DO j= 1, ny
                CT(nx+1,j,k)%bx = ar*CT(nx+1,j,k)%obx + br*CT(nx+1,j,k)%bx  &
                + cr*dt*(CT(nx+1,j+1,k)%omega_z-CT(nx+1,j,k)%omega_z)/(yznr(j)-yznl(j))
             ENDDO
!$OMP END DO
          ENDDO


          DO k= 1, nz
!$OMP DO SCHEDULE(static)
             DO i= 1, nx
                CT(i,ny+1,k)%by = ar*CT(i,ny+1,k)%oby + br*CT(i,ny+1,k)%by  &
                + cr*dt*(CT(i,ny+1,k)%omega_z-CT(i+1,ny+1,k)%omega_z)/(xznr(i)-xznl(i))
             ENDDO
!$OMP END DO
          ENDDO

!$OMP END PARALLEL

!        ENDIF
                


!$OMP PARALLEL DEFAULT(shared)        &
!$OMP PRIVATE(i,j,k)

        DO k=1, nz
!$OMP DO SCHEDULE(static)
           DO j=1,ny
              DO i=1,nx

                 physic(i,j,k)%bxct = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                    
                 IF (nsdim .EQ. 2) THEN

                    physic(i,j,k)%byct = (CT(i,j,k)%by + CT(i,j+1,k)%by)/2.D0 

                 ELSE IF (nsdim .EQ. 3) THEN

                    physic(i,j,k)%byc = (CT(i,j,k)%by + CT(i,j+1,k)%by)*0.5 
                    physic(i,j,k)%bzc = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)*0.5

                 ENDIF
       
              ENDDO
           ENDDO
!$OMP END DO
        ENDDO

!$OMP END PARALLEL


        IF (iccv .NE. 0) THEN
       
!$OMP PARALLEL DEFAULT (shared) &
!$OMP PRIVATE (i,j,k,b2,bct2)

           DO k=1, nz
!$OMP DO SCHEDULE (static)
              DO j=1,ny
                 DO i=1,nx

	           IF (iccv .EQ. 2) physic(i,j,k)%bxzn = physic(i,j,k)%bxc
		   IF (iccv .EQ. 2) physic(i,j,k)%byzn = physic(i,j,k)%byc

                   b2 = physic(i,j,k)%bxc**2 + physic(i,j,k)%byc**2
                   bct2 = physic(i,j,k)%bxct**2 + physic(i,j,k)%byct**2

                   physic(i,j,k)%energy = physic(i,j,k)%energy - 0.5D0*(b2-bct2)

                ENDDO
              ENDDO
!$OMP END DO
           ENDDO

!$OMP END PARALLEL

        ENDIF


!$OMP PARALLEL DEFAULT (shared)  &
!$OMP PRIVATE (i,j,k)

        DO k=1,nz
!$OMP DO SCHEDULE (static)
           DO j=1,ny
              DO i=1, nx
 
	         physic(i,j,k)%bxc = physic(i,j,k)%bxct
	
	         IF (nsdim .EQ. 2) THEN
                    physic(i,j,k)%byc = physic(i,j,k)%byct
 		 ENDIF

             ENDDO
           ENDDO
!$OMP END DO
        ENDDO
!$OMP END PARALLEL

ENDIF

        !----------------------------------------------------------------------------END

! =============================================================================== 
END SUBROUTINE SWEEPS_2D
! =============================================================================== 

 
