!*********************************
!     3D-RMHD CODE: LOSTREGO
!*********************************

      
! NAME:
!       SWEEPS (3D VERSION)

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
!        See also: "MHD relativista numerica: Aplicaciones en Relatividad Especial y General."
!        Luis Antón Ruiz, 2008.


! ===============================================================================      
SUBROUTINE SWEEPS_3D (basenm, limiter, irk, rkin, iars, iccv,                            & 
                      igeomx, igeomy, igeomz, bndmnx, bndmxx, bndmny, bndmxy,            &
                      bndmnz, bndmxz, nx, mnx5, mnx6, nsdim, ny, mny5, mny6,             &
                      eqs, nz, mnz5, mnz6, stelar, cool, gamma_ad, alpha_R, beta_R, dt,  &
                      me, mp, xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr,         &
                      boundPosX, boundPosY, boundPosZ, grav, physic, CT, uwi, vwi, shift,&
                      Rorb, kb, rho_0a, R_b)
! ===============================================================================


        USE tipos
        USE workarrays
        implicit none

# include "types.h"

        !------------------------------------------------------- Input Variables

        character(4) :: basenm, limiter

        integer(is) :: rkin, iars, iccv, igeomx, igeomy, igeomz,        &
                       bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,  &
                       nx, mnx5, nsdim, ny, mny5, eqs, nz, mnz5,        &
                       stelar, cool, mnx6, mny6, mnz6, irk

        real(rs) :: gamma_ad, alpha_R, beta_R, dt, me, mp, uwi, vwi, shift, Rorb                 

        ! sizes mnx1, mny1, mnz1
        real(rs) :: xznl(*), xzn(*), xznr(*),   &
                    yznl(*), yzn(*), yznr(*),   &
                    zznl(*), zzn(*), zznr(*)

        type(TboundPosX) :: boundPosX(-4:5)
        type(TboundPosY) :: boundPosY(-4:5)
        type(TboundPosZ) :: boundPosZ(-4:5)

	! COOLING
	real(rs) :: kb, rho_0a, R_b

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

        real(rs)   :: ar, br, cr, aux, alpha1, alpha2, alpha3, alpha4, dtdx, dtdy, dtdz

        real(rs), pointer :: vel(:), velt(:), veltt(:),        &
                             rho(:), rhoe(:), rhopr(:), p(:),  &
                             b(:), bt(:), btt(:),              &
                             x(:), xl(:), xr(:), dx(:),        &
                             area(:), areal(:), arear(:),      &
                             dvol(:), bs(:)

	real(rs) :: y, yl, yr, dy,        &
		    z, zl, zr, dz

        real(rs), pointer :: f(:), fict(:)

        real(rs), pointer :: rflx(:), reflx(:), mflx(:), mtflx(:),  &
                             mttflx(:), eflx(:), fflx(:), btflx(:), & 
                             bttflx(:)

        real(rs), pointer :: trac(:,:,:)

        character(32) :: message1

	type(ConsTransp) :: CT(-4:mnx6, -4:mny6, -4:mnz6)
     
	real(rs) :: b2, bct2, div(nx-1,ny-1), tt2, tt3, theta1, theta2                        

	real(rs) :: RBM, BMJ, WB, BPHI1, VELB, RHOB, PB, DIVE

        real(rs) :: SCRCH1, SCRCH2, SCRCH3, SCRCH4, SCRCH5, SCRCH6, FLUXC
        real(rs) :: SCRCH11, SCRCH12, SCRCH21, SCRCH22, SCRCH41, SCRCH42, SCRCH51, SCRCH52

        real(rs) :: VELX, VELY, VELZ, R, dbx_dx, dby_dy, dbz_dz, deltax1, deltax2, ss1, ss2

	! COOLING
	real(rs) :: xe_s, meff, tem, dedt

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

        ! The zone centered field vector Bc is required to set the boundary conditions
        ! and to calculate the numerical fluxes.


!$OMP PARALLEL DEFAULT(shared)                                                     & 
!$OMP PRIVATE(curth,tid,i,j,k,fznp,nznp,nxp,nyp,nzp,bndmin,bndmax,igeom,fzn,nzn)   &
!$OMP PRIVATE(ar,br,cr,aux,xyzswp,trac)                                            &
!$OMP PRIVATE(vel,velt,veltt,rho,rhoe,rhopr,p,b,bt,btt,bs,x,xl,xr,dx)              &
!$OMP PRIVATE(y,yl,yr,dy,z,zl,zr,dz)                                               &
!$OMP PRIVATE(area,areal,arear,dvol,f,fict)                                        &
!$OMP PRIVATE(rflx,reflx,mflx,mtflx,mttflx,eflx,fflx,btflx,bttflx)                 &
!$OMP PRIVATE(xe_s,meff,tem,dedt)


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

!$OMP DO SCHEDULE(static)
       DO k= 1, nzp
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

                
                trac(i,j,k) = physic(i,j,k)%tracer
                f(i) = physic(i,j,k)%tracer

                xl(i)   = xznl(i)
                x(i)    = xzn(i)
                xr(i)   = xznr(i)
                dx(i)   = xr(i)-xl(i)

                fict(i) = 0.D0
                
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

                         physic(i,j,k)%tracer2  = f(i)

		      ENDDO

		   ENDIF

		   IF (mpi_coords(1) .EQ. (mpi_dims(1)-1)) THEN

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

                         physic(i,j,k)%tracer2  = f(i)

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
                         
                   physic(i,j,k)%bxc2 = physic(i,j,k)%bxc
                   physic(i,j,k)%byc2 = physic(i,j,k)%byc
                   physic(i,j,k)%bzc2 = physic(i,j,k)%bzc

                   rho(i)    = physic(i,j,k)%densty
                   rhoe(i)   = physic(i,j,k)%denstye
                   rhopr(i)  = physic(i,j,k)%densty - physic(i,j,k)%denstye
                   p(i)      = physic(i,j,k)%pres
               
                   f(i)    = physic(i,j,k)%tracer
                   
                   xl(i)   = boundPosX(i)%xznl 
                   x(i)    = boundPosX(i)%xzn  
                   xr(i)   = boundPosX(i)%xznr 
                   dx(i)   = boundPosX(i)%dx   

                   fict(i) = 0.D0
                         
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
                         
                   physic(i,j,k)%bxc2 = physic(i,j,k)%bxc
                   physic(i,j,k)%byc2 = physic(i,j,k)%byc
                   physic(i,j,k)%bzc2 = physic(i,j,k)%bzc


                   rho(i)    = physic(i,j,k)%densty
                   rhoe(i)   = physic(i,j,k)%denstye
                   rhopr(i)  = physic(i,j,k)%densty - physic(i,j,k)%denstye
                   p(i)      = physic(i,j,k)%pres
                         
                   f(i)    = physic(i,j,k)%tracer
                   
                   xl(i)   = boundPosX(i-nznp)%xznl 
                   x(i)    = boundPosX(i-nznp)%xzn  
                   xr(i)   = boundPosX(i-nznp)%xznr 
                   dx(i)   = boundPosX(i-nznp)%dx   

                   fict(i) = 0.D0
                         
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
             call MHYDROW3 (ii,tid,xyzswp,basenm, limiter, fznp, nznp, iars,           &
                            eqs, gamma_ad, me, mp, x, xl, xr, areal,                   &
                            p, rho, rhoe, rhopr, vel, velt, veltt, f, b, bs, bt, btt,  &
                            rflx, reflx, mflx, mtflx, mttflx,                          &
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


                ELSE IF (nsdim .EQ. 1 .OR. (nsdim .EQ. 3 .AND. iccv .NE. 0)) THEN

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
		   CT(i,j,k)%rflxx= rflx(i)

                ENDDO

             ENDIF

          ENDDO
       ENDDO
!$OMP END DO



       ! Synchronization: Each thread waits for all the other threads at the end of the sweep.
!$OMP BARRIER

       ! =============================
       !           Y-SWEEP
       ! =============================

       ! Runge-Kutta parameters.
       ! rkin goes from 1 to irk (Runge-Kutta order).
       ! Runge-Kutta Main Loop is defined in the main program.


       IF (rkin .EQ. 1) THEN

          br = 1.D0
          cr = 1.D0

       ELSE IF (rkin .EQ. 2) THEN

          br = 1.D0
          cr = 1.D0/alpha_R

       ELSE

          br = 1.D0
          cr = 2.D0/3.D0

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
 
!$OMP DO SCHEDULE(static)
          DO k= 1, nzp
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
                   fict(j) = 0.D0

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

                      IF (basenm(1:2) .EQ. 'MH' .OR. basenm(1:2) .EQ. 'MC') THEN
                         DO j= fznp-5, fznp-1

			    IF (physic(i,j,k)%bxc**2+physic(i,j,k)%bzc**2 .NE. 0.D0) THEN

                            b(j)    = 0.D0
                            bt(j)   = physic(i,j,k)%bzc
                            btt(j)  = physic(i,j,k)%bxc

			    ENDIF

                         ENDDO
                      ENDIF
                      
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

                         physic(i,j,k)%tracer2  = f(j)

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

                         physic(i,j,k)%tracer2  = f(j)

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
                   
                      physic(i,j,k)%bxc2 = physic(i,j,k)%bxc
                      physic(i,j,k)%byc2 = physic(i,j,k)%byc
                      physic(i,j,k)%bzc2 = physic(i,j,k)%bzc

                      rho(j)  = physic(i,j,k)%densty
                      rhoe(j) = physic(i,j,k)%denstye
                      rhopr(j)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                      p(j)    = physic(i,j,k)%pres

                      f(j) = physic(i,j,k)%tracer
                   
                      xl(j)   = boundPosY(j)%yznl
                      x(j)    = boundPosY(j)%yzn
                      xr(j)   = boundPosY(j)%yznr
                      dx(j)   = boundPosY(j)%dx

                      fict(j) = 0.D0

                   ENDDO

                ENDIF


                ! ====================
                ! RIGHT GHOST CELLS
                ! ====================

                IF (mpi_coords(2) .NE. (mpi_dims(2)-1)) THEN

                   DO j= nznp+1, nznp+5

                      vel(j)  = physic(i,j,k)%vely
                      velt(j) = physic(i,j,k)%velz
                      veltt(j)= physic(i,j,k)%velx

                      b(j)    = physic(i,j,k)%byc
                      bt(j)   = physic(i,j,k)%bzc
                      btt(j)  = physic(i,j,k)%bxc
                   
                      physic(i,j,k)%bxc2 = physic(i,j,k)%bxc
                      physic(i,j,k)%byc2 = physic(i,j,k)%byc
                      physic(i,j,k)%bzc2 = physic(i,j,k)%bzc


                      rho(j)  = physic(i,j,k)%densty
                      rhoe(j) = physic(i,j,k)%denstye
                      rhopr(j)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                      p(j)    = physic(i,j,k)%pres

                      f(j) = physic(i,j,k)%tracer
                   
                      xl(j)   = boundPosY(j-nznp)%yznl
                      x(j)    = boundPosY(j-nznp)%yzn
                      xr(j)   = boundPosY(j-nznp)%yznr
                      dx(j)   = boundPosY(j-nznp)%dx

                      fict(j) = 0.D0

                   ENDDO

                ENDIF   


                ! Calculate the sections of the cells in the normal direction.

                call GEOM (fznp, nznp, igeom, dx, area, areal, arear, dvol) 
                
                ! Calculate numerical fluxes.

	        ii=i

                yl   = zznl(k)
                y    = zzn(k)
                yr   = zznr(k)
                dy   = yr-yl

                zl   = xznl(i)
                z    = xzn(i)
                zr   = xznr(i)
                dz   = zr-zl


                call MHYDROW3 (ii,tid, xyzswp,basenm,limiter, fznp, nznp, iars,          &
                              eqs, gamma_ad, me, mp, x, xl, xr, areal,                   &
                              p, rho, rhoe, rhopr, vel, velt, veltt, f, b, bs, bt, btt,  &
                              rflx, reflx, mflx, mtflx, mttflx,                          &
                              eflx, btflx, bttflx, fflx,y) 


               
                !========================
                ! Runge-Kutta substep 
                !========================

                ! Advance of conserved variables.
                   

                DO j= fznp, nznp

                   aux = dt / dx(j)

                   physic(i,j,k)%tracer = br * physic(i,j,k)%tracer      &
                                        - cr * (fflx(j+1)-fflx(j))*aux

                   physic(i,j,k)%denstw = br * physic(i,j,k)%denstw      &
                                        - cr * (rflx(j+1)-rflx(j))*aux

                   physic(i,j,k)%denstew = br * physic(i,j,k)%denstew    &
                                         - cr * (reflx(j+1)-reflx(j))*aux   

                   physic(i,j,k)%momenx = br * physic(i,j,k)%momenx      &
                                        - cr * (mttflx(j+1)-mttflx(j))*aux 

                   physic(i,j,k)%momeny = br * physic(i,j,k)%momeny      &
                                        - cr * (mflx(j+1)-mflx(j))*aux 

                   physic(i,j,k)%momenz = br * physic(i,j,k)%momenz      &
                                        - cr * (mtflx(j+1)-mtflx(j))*aux 

                   physic(i,j,k)%energy = br * physic(i,j,k)%energy      &
                                        - cr * (eflx(j+1)-eflx(j))*aux



                   ! Magnetic field.


		   IF (iccv .NE. 0) THEN

                      physic(i,j,k)%bxc    = br * physic(i,j,k)%bxc         &
                                           - cr * (bttflx(j+1)-bttflx(j))*aux

                      physic(i,j,k)%byc    = br * physic(i,j,k)%byc         
                                           

		   ENDIF 

                 
                   IF (nsdim .NE. 3 .OR. (nsdim .EQ. 3 .AND. iccv .NE. 0)) THEN
                     
                      physic(i,j,k)%bzc    = br * physic(i,j,k)%bzc       &
                                           - cr * (btflx(j+1)-btflx(j))*aux
			
                   ENDIF


                   ! Normalizing Tracer for two-dimensional simulations...

                   IF (nsdim .EQ. 2) THEN

                      physic(i,j,k)%tracer = physic(i,j,k)%tracer/    &
                                             physic(i,j,k)%denstw

                   ENDIF

                ENDDO
                
                ! CT-fluxes.

                 DO j=fznp, nznp+1

                   CT(i,j,k)%fbxy = bttflx(j)
                   CT(i,j,k)%fbzy = btflx(j)
                   CT(i,j,k)%rflxy= rflx(j)


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

      
       ! Runge-Kutta parameters.
       ! rkin goes from 1 to irk (Runge-Kutta order).
       ! Runge-Kutta Main Loop is defined in the main program.


       IF (rkin .EQ. 1) THEN

          br = 1.D0
          cr = 1.D0

       ELSE IF (rkin .EQ. 2) THEN

          br = 1.D0
          cr = 1.D0/alpha_R

       ELSE

          br = 1.D0
          cr = 2.D0/3.D0

       ENDIF


       IF (nsdim .GT. 2) THEN

          ! ZSweep own parameters.

          xyzswp  = 3
          bndmin = bndmnz
          bndmax = bndmxz
          igeom  = igeomz
          fzn    = 1
          nzn    = nz
          fznp   = fzn
          nznp   = nzn
          nxp    = nx
          nyp    = ny
 

          ! Take the primitive variables at the cell centers. 
 
!$OMP DO SCHEDULE(static)

          DO j= 1, nyp
             DO i=1,nxp

                DO k=fznp, nznp+1
                   bs(k) = CT(i,j,k)%bz
                ENDDO

                DO k= fznp, nznp

                   vel(k)  = physic(i,j,k)%velz
                   velt(k) = physic(i,j,k)%velx
                   veltt(k)= physic(i,j,k)%vely

                   b(k)    = physic(i,j,k)%bzc2
                   bt(k)   = physic(i,j,k)%bxc2
                   btt(k)  = physic(i,j,k)%byc2
                   
                   rho(k)  = physic(i,j,k)%densty
                   rhoe(k) = physic(i,j,k)%denstye
                   rhopr(k)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                   p(k)    = physic(i,j,k)%pres

                   f(k) = physic(i,j,k)%tracer2
                   
                   xl(k)   = zznl(k)
                   x(k)    = zzn(k)
                   xr(k)   = zznr(k)
                   dx(k)   = xr(k)-xl(k)

                   fict(k) = 0.D0

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



		   IF (mpi_coords(3) .EQ. 0) THEN

                      DO k= fznp-5, fznp-1

                         physic(i,j,k)%velz  =  vel(k)
                         physic(i,j,k)%velx  =  velt(k)
                         physic(i,j,k)%vely  =  veltt(k)

                         physic(i,j,k)%bzc2   =  b(k)
                         physic(i,j,k)%bxc2   =  bt(k)
                         physic(i,j,k)%byc2   =  btt(k)

                         physic(i,j,k)%densty  = rho(k)
                         physic(i,j,k)%denstye = rhoe(k)
                         physic(i,j,k)%pres    = p(k)

                         physic(i,j,k)%tracer2  = f(k)

		      ENDDO

		   ENDIF


		   IF (mpi_coords(3) .EQ. mpi_dims(3)-1) THEN

                      DO k= nznp+1, nznp+5

                         physic(i,j,k)%velz  =  vel(k)
                         physic(i,j,k)%velx  =  velt(k)
                         physic(i,j,k)%vely  =  veltt(k)

                         physic(i,j,k)%bzc2   =  b(k)
                         physic(i,j,k)%bxc2   =  bt(k)
                         physic(i,j,k)%byc2   =  btt(k)

                         physic(i,j,k)%densty  = rho(k)
                         physic(i,j,k)%denstye = rhoe(k)
                         physic(i,j,k)%pres    = p(k)

                         physic(i,j,k)%tracer2  = f(k)

		      ENDDO

		   ENDIF
                                
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
                   
                      physic(i,j,k)%bxc2 = physic(i,j,k)%bxc
                      physic(i,j,k)%byc2 = physic(i,j,k)%byc
                      physic(i,j,k)%bzc2 = physic(i,j,k)%bzc


                      rho(k)  = physic(i,j,k)%densty
                      rhoe(k) = physic(i,j,k)%denstye
                      rhopr(k)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                      p(k)    = physic(i,j,k)%pres

                      f(k) = physic(i,j,k)%tracer
                   
                      xl(k)   = boundPosZ(k)%zznl
                      x(k)    = boundPosZ(k)%zzn
                      xr(k)   = boundPosZ(k)%zznr
                      dx(k)   = boundPosZ(k)%dx

                      fict(k) = 0.D0

                   ENDDO

                ENDIF


                ! ====================
                ! RIGHT GHOST CELLS
                ! ====================


                IF (mpi_coords(3) .NE. (mpi_dims(3)-1)) THEN

                   DO k= nznp+1, nznp+5

                      vel(k)  = physic(i,j,k)%velz
                      velt(k) = physic(i,j,k)%velx
                      veltt(k)= physic(i,j,k)%vely

                      b(k)    = physic(i,j,k)%bzc
                      bt(k)   = physic(i,j,k)%bxc
                      btt(k)  = physic(i,j,k)%byc
         
                      physic(i,j,k)%bxc2 = physic(i,j,k)%bxc
                      physic(i,j,k)%byc2 = physic(i,j,k)%byc
                      physic(i,j,k)%bzc2 = physic(i,j,k)%bzc

          
                      rho(k)  = physic(i,j,k)%densty
                      rhoe(k) = physic(i,j,k)%denstye
                      rhopr(k)= physic(i,j,k)%densty - physic(i,j,k)%denstye
                      p(k)    = physic(i,j,k)%pres

                      f(k) = physic(i,j,k)%tracer
                   
                      xl(k)   = boundPosZ(k-nznp)%zznl
                      x(k)    = boundPosZ(k-nznp)%zzn
                      xr(k)   = boundPosZ(k-nznp)%zznr
                      dx(k)   = boundPosZ(k-nznp)%dx

                      fict(k) = 0.D0

                   ENDDO

                ENDIF   


                ! Calculate the sections of the cells in the normal direction.

                call GEOM (fznp, nznp, igeom, dx, area, areal, arear, dvol) 
              
                ! Calculate numerical fluxes.

		ii=k

                 yl   = xznl(i)
                 y    = xzn(i)
                 yr   = xznr(i)
                 dy   = yr-yl

                 zl   = yznl(j)
                 z    = yzn(j)
                 zr   = yznr(j)
                 dz   = zr-zl


                call MHYDROW3 (ii,tid,xyzswp,basenm, limiter, fznp, nznp, iars,             &
                               eqs, gamma_ad, me, mp, x, xl, xr, areal,                     &
                               p, rho, rhoe, rhopr, vel, velt, veltt, f, b, bs, bt, btt,    &
                               rflx, reflx, mflx, mtflx, mttflx,                            &
                               eflx, btflx, bttflx, fflx, z)


                !========================
                ! Runge-Kutta substep 
                !========================

                ! Advance of conserved variables.
                   

                DO k= fznp, nznp

                   aux = dt / dx(k)

                   physic(i,j,k)%tracer = br * physic(i,j,k)%tracer      &
                                        - cr * (fflx(k+1)-fflx(k))*aux

                   physic(i,j,k)%denstw = br * physic(i,j,k)%denstw      &
                                        - cr * (rflx(k+1)-rflx(k))*aux

                   physic(i,j,k)%denstew = br * physic(i,j,k)%denstew    &
                                         - cr * (reflx(k+1)-reflx(k))*aux   

                   physic(i,j,k)%momenx = br * physic(i,j,k)%momenx      &
                                        - cr * (mtflx(k+1)-mtflx(k))*aux 

                   physic(i,j,k)%momeny = br * physic(i,j,k)%momeny      &
                                        - cr * (mttflx(k+1)-mttflx(k))*aux 

                   physic(i,j,k)%momenz = br * physic(i,j,k)%momenz      &
                                        - cr * (mflx(k+1)-mflx(k))*aux 

                   physic(i,j,k)%energy = br * physic(i,j,k)%energy      &
                                        - cr * (eflx(k+1)-eflx(k))*aux 



		   ! THERMAL COOLING
		   ! ---------------------------------------------------------------------

		   IF (BASENM(1:2) .EQ. 'S4') THEN

                      xe_s  =  physic(i,j,k)%denstye/physic(i,j,k)%densty
                      meff= 1.D0/9.1095D-28 * (xe_s + (1.D0-xe_s)/1836.11D0)
                      tem = physic(i,j,k)%pres/physic(i,j,k)%densty*9.D20/meff/kb
                  
                      IF (tem .LT. 1.D4) THEN
                         dedt = 0.D0
                      ELSE IF (tem .GE. 1.D4 .AND. tem .LT. 1.D5) THEN
                         dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 7.D-27*tem 
                      ELSE IF (tem .GE. 1.D5 .AND. tem .LT. 4.D7) THEN
                         dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 7.D-19*tem**(-0.6D0)
                      ELSE IF (tem .GE. 4.D7) THEN
                         dedt = (physic(i,j,k)%densty*rho_0a*meff)**2 * 3.D-27*tem**(0.5D0)
                      ENDIF
                  
                      dedt = dedt/rho_0a/9.D20/3.D10*R_b
                      physic(i,j,k)%energy = physic(i,j,k)%energy - cr*dedt*dt
                     
                      IF (physic(i,j,k)%energy .LE. 0.D0) THEN

                    	   WRITE(6,*) 'SWEEPS_3D: MINIMUM ENERGY REACHED AT POINT:'
                    	   print*, i,j,k
                    	   WRITE(6,*) '             RESETING COOLING TERMS'

                           physic(i,j,k)%energy = physic(i,j,k)%energy + cr*dedt*dt

		      ENDIF
			

                   ENDIF

                   ! ---------------------------------------------------------------------

                   
                   ! Magnetic field.


                   IF (iccv .NE. 0) THEN

                      physic(i,j,k)%bxc    = br * physic(i,j,k)%bxc         &
                                           - cr * (btflx(k+1)-btflx(k))*aux

                      physic(i,j,k)%byc    = br * physic(i,j,k)%byc         &
                                           - cr * (bttflx(k+1)-bttflx(k))*aux

                      physic(i,j,k)%bzc    = br * physic(i,j,k)%bzc         
                                           

                   ENDIF
         

                   ! Normalizing Tracer for three-dimensional simulations.
             
                   physic(i,j,k)%tracer = physic(i,j,k)%tracer/    &
                                          physic(i,j,k)%denstw

                ENDDO
                
                ! CT-fluxes.

                DO k= fznp, nznp +1

                   CT(i,j,k)%fbxz = btflx(k)
                   CT(i,j,k)%fbyz = bttflx(k)
                   CT(i,j,k)%rflxz= rflx(k)


                ENDDO
             ENDDO
          ENDDO

!$OMP END DO

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

       ! rflxx : density flux in the x-direction.
       ! rflxy : density flux in the y-direction.
       ! rflxz : density flux in the z-direction.

       
          !!! Boundary conditions for the staggered fluxes.

          ! ------------
          ! LEFT X-BOUNDARY
          ! ------------

          ! Reflecting Boundary.

	  IF (mpi_coords(1) .EQ. 0) THEN

          IF (bndmnx .EQ. 1) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz
                DO  j=1,ny+1
                   CT(0,j,k)%fbxy = - CT(1,j,k)%fbxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz+1
                DO  j=1,ny
                   CT(0,j,k)%fbxz = - CT(1,j,k)%fbxz
                ENDDO
             ENDDO
!OMP END DO

!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmnx .EQ. 2) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO j=1,ny+1
                   CT(0,j,k)%fbxy = CT(1,j,k)%fbxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO j=1,ny
                   CT(0,j,k)%fbxz = CT(1,j,k)%fbxz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)             
             DO k=1,nz
                DO j=1,ny+1
                   CT(0,j,k)%rflxy = CT(1,j,k)%rflxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)             
             DO k=1,nz+1
                DO j=1,ny
                   CT(0,j,k)%rflxz = CT(1,j,k)%rflxz
                ENDDO
             ENDDO
!$OMP END DO
             
!$OMP END PARALLEL


          ! Periodic Boundary.

          ELSE IF (bndmnx .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO j=1,ny+1
                   CT(0,j,k)%fbxy = CT(nx,j,k)%fbxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO j=1,ny
                   CT(0,j,k)%fbxz = CT(nx,j,k)%fbxz
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL



          ELSE

             message1 = 'Error in boundaries of B fluxes.'
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)

          ENDIF

	  ENDIF

          ! ------------
          ! RIGHT X-BOUNDARY
          ! ------------

	  IF (mpi_coords(1) .EQ. (mpi_dims(1)-1)) THEN

          ! Reflecting Boundary.

          IF (bndmxx .EQ. 1) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz
                DO  j=1,ny+1
                   CT(nx+1,j,k)%fbxy = - CT(nx,j,k)%fbxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz+1
                DO  j=1,ny
                   CT(nx+1,j,k)%fbxz = - CT(nx,j,k)%fbxz
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmxx .EQ. 2) THEN


!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO j=1,ny+1
                   CT(nx+1,j,k)%fbxy = CT(nx,j,k)%fbxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO j=1,ny
                   CT(nx+1,j,k)%fbxz = CT(nx,j,k)%fbxz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO j=1,ny+1
                   CT(nx+1,j,k)%rflxy = CT(nx,j,k)%rflxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO j=1,ny
                   CT(nx+1,j,k)%rflxz = CT(nx,j,k)%rflxz
                ENDDO
             ENDDO
!$OMP END DO


             
!$OMP END PARALLEL



          ! Periodic Boundary.

          ELSE IF (bndmxx .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K)

!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz
                DO  j=1,ny+1
                   CT(nx+1,j,k)%fbxy =  CT(1,j,k)%fbxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz+1
                DO  j=1,ny
                   CT(nx+1,j,k)%fbxz =  CT(1,j,k)%fbxz
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL


          ELSE IF (bndmxx .EQ. 9) THEN


!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (J,K,THETA1,THETA2,VELX,VELY,VELZ,SCRCH1)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO j=1,ny+1
                   CT(nx+1,j,k)%fbxy = CT(nx,j,k)%fbxy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO j=1,ny
                   CT(nx+1,j,k)%fbxz = CT(nx,j,k)%fbxz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO j=1,ny+1

                       
                    ! INCLINATION
                    theta1 = DATAN(zzn(k)/       &
                             DSQRT((yznl(j)+shift)**2+(Rorb-(xzn(nx)+xznr(nx)-xznl(nx)))**2))

	            ! AZIMUTHAL ANGLE.
                    theta2 = DATAN((yznl(j)+shift)/(Rorb-(xzn(nx)+xznr(nx)-xznl(nx))))
		
                    VELX = -vwi*DCOS(theta1)*DCOS(theta2)
                    VELY =  vwi*DCOS(theta1)*DSIN(theta2)
                    VELZ =  vwi*DSIN(theta1)

			!VELX = 0.D0
			!VELY = 0.D0
			!VELZ = 0.D0

                    SCRCH1 = 1.D0/DSQRT(1.D0-VELX**2-VELY**2-VELZ**2)
                  
                    CT(nx+1,j,k)%rflxy = physic(nx,j,k)%densty*SCRCH1*VELY


                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO j=1,ny

                    ! INCLINATION
                    theta1 = DATAN(zznl(k)/       &
                             DSQRT((yzn(j)+shift)**2+(Rorb-(xzn(nx)+xznr(nx)-xznl(nx)))**2))

	            ! AZIMUTHAL ANGLE.
                    theta2 = DATAN((yzn(j)+shift)/(Rorb-(xzn(nx)+xznr(nx)-xznl(nx))))
		
                    VELX = -vwi*DCOS(theta1)*DCOS(theta2)
                    VELY =  vwi*DCOS(theta1)*DSIN(theta2)
                    VELZ =  vwi*DSIN(theta1)

		        !VELX = 0.D0
                        !VELY = 0.D0
                        !VELZ = 0.D0


                    SCRCH1 = 1.D0/DSQRT(1.D0-VELX**2-VELY**2-VELZ**2)
                  
                    CT(nx+1,j,k)%rflxz = physic(nx,j,k)%densty*SCRCH1*VELZ
                   


                ENDDO
             ENDDO
!$OMP END DO

             
!$OMP END PARALLEL


          
          ELSE

             message1 = 'Error in boundaries of B fluxes.'
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)

          ENDIF

	  ENDIF

          ! ------------
          ! LEFT Y-BOUNDARY
          ! ------------

	  IF (mpi_coords(2) .EQ. 0) THEN

          ! Reflecting Boundary.

          IF (bndmny .EQ. 1) THEN
 
!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz
                DO  i=1,nx+1
                   CT(i,0,k)%fbyx = - CT(i,1,k)%fbyx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz+1
                DO  i=1,nx
                   CT(i,0,k)%fbyz = - CT(i,1,k)%fbyz
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmny .EQ. 2) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,0,k)%fbyx = CT(i,1,k)%fbyx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,0,k)%fbyz = CT(i,1,k)%fbyz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER             
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,0,k)%rflxx = CT(i,1,k)%rflxx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,0,k)%rflxz = CT(i,1,k)%rflxz
                ENDDO
             ENDDO
!$OMP END DO
             
!$OMP END PARALLEL

          ! Periodic Boundary.

          ELSE IF (bndmny .EQ. 3) THEN


!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,0,k)%fbyx = CT(i,ny,k)%fbyx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,0,k)%fbyz = CT(i,ny,k)%fbyz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,0,k)%rflxx = CT(i,ny,k)%rflxx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,0,k)%rflxz = CT(i,ny,k)%rflxz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP END PARALLEL


          ! Periodic Boundary.

          ELSE IF (bndmny .EQ. 5) THEN


!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,0,k)%fbyx = CT(i,10,k)%fbyx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,0,k)%fbyz = CT(i,10,k)%fbyz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,0,k)%rflxx = CT(i,10,k)%rflxx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,0,k)%rflxz = CT(i,10,k)%rflxz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP END PARALLEL



          ! +++++++++++++++++++++++++++++++++++++++++++++++++++
          ELSE IF (bndmny .EQ. 7 .OR. bndmny .EQ. 8) THEN


!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K,R)           &
!$OMP PRIVATE (THETA1,THETA2,SCRCH1,VELX,VELY,VELZ)
             
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1

                    R = DSQRT(XZNL(I)**2+ZZN(K)**2)

		    IF (R .LE. 1.5D0) THEN

                        CT(i,0,k)%fbyx =-CT(i,0,k)%bx*0.55D0/DCOSH(R**4)

		    ELSE IF  (R .GT. 1.5D0 .AND. R .LE. 2.D0) THEN

		        IF (physic(i,1,k)%pres .GT. 1.D-8) THEN

                            CT(i,0,k)%fbyx =-CT(i,0,k)%bx*(-physic(i,1,k)%vely)

			ELSE


                            ! INCLINATION
                            theta1 = DATAN(zzn(k)/       &
                                     DSQRT((yzn(1)+shift-(yznr(1)-yznl(1)))**2+(Rorb-xznl(i))**2))

                            ! AZIMUTHAL ANGLE.
                       	    theta2 = DATAN((yzn(1)+shift-(yznr(1)-yznl(1)))/(Rorb-xznl(i)))


                            VELX = -vwi*DCOS(theta1)*DCOS(theta2)
                            VELY =  vwi*DCOS(theta1)*DSIN(theta2)
                            VELZ =  vwi*DSIN(theta1)

                            CT(i,0,k)%fbyx =-CT(i,0,k)%bx*VELY

			ENDIF


		    ELSE


                        IF (physic(i,1,k)%pres .GT. 1.D-8) THEN

                            CT(i,0,k)%fbyx =-CT(i,1,k)%bx*(-physic(i,1,k)%vely)

                        ELSE


                            ! INCLINATION
                            theta1 = DATAN(zzn(k)/       &
                                     DSQRT((yzn(1)+shift-(yznr(1)-yznl(1)))**2+(Rorb-xznl(i))**2))

                            ! AZIMUTHAL ANGLE.
                            theta2 = DATAN((yzn(1)+shift-(yznr(1)-yznl(1)))/(Rorb-xznl(i)))


                            VELX = -vwi*DCOS(theta1)*DCOS(theta2)
                            VELY =  vwi*DCOS(theta1)*DSIN(theta2)
                            VELZ =  vwi*DSIN(theta1)

                            CT(i,0,k)%fbyx =-CT(i,1,k)%bx*VELY                 

                        ENDIF



		    ENDIF


                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx

                    R = DSQRT(XZN(I)**2+ZZNL(K)**2)


                    IF (R .LE. 1.5D0) THEN

                        CT(i,0,k)%fbyz =-CT(i,0,k)%bz*0.55D0/DCOSH(R**4)

                    ELSE IF  (R .GT. 1.5D0 .AND. R .LE. 2.D0) THEN

                        IF (physic(i,1,k)%pres .GT. 1.D-8) THEN

                            CT(i,0,k)%fbyz =-CT(i,0,k)%bz*(-physic(i,1,k)%vely)

                        ELSE


                            ! INCLINATION
                            theta1 = DATAN(zznl(k)/	 &
                                     DSQRT((yzn(1)+shift-(yznr(1)-yznl(1)))**2+(Rorb-xzn(i))**2))

                            ! AZIMUTHAL ANGLE.
                            theta2 = DATAN((yzn(1)+shift-(yznr(1)-yznl(1)))/(Rorb-xzn(i)))


                            VELX = -vwi*DCOS(theta1)*DCOS(theta2)
                            VELY =  vwi*DCOS(theta1)*DSIN(theta2)
                            VELZ =  vwi*DSIN(theta1)

                            CT(i,0,k)%fbyz =-CT(i,0,k)%bz*VELY

                        ENDIF


                    ELSE


                        IF (physic(i,1,k)%pres .GT. 1.D-8) THEN

                            CT(i,0,k)%fbyz =-CT(i,1,k)%bz*(-physic(i,1,k)%vely)

                        ELSE


                            ! INCLINATION
                            theta1 = DATAN(zznl(k)/	 &
                                     DSQRT((yzn(1)+shift-(yznr(1)-yznl(1)))**2+(Rorb-xzn(i))**2))

                            ! AZIMUTHAL ANGLE.
                            theta2 = DATAN((yzn(1)+shift-(yznr(1)-yznl(1)))/(Rorb-xzn(i)))


                            VELX = -vwi*DCOS(theta1)*DCOS(theta2)
                            VELY =  vwi*DCOS(theta1)*DSIN(theta2)
                            VELZ =  vwi*DSIN(theta1)

                            CT(i,0,k)%fbyz =-CT(i,1,k)%bz*VELY

                        ENDIF



                    ENDIF

                          
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1

                    R = DSQRT(XZNL(I)**2+ZZN(K)**2)

                    IF (R .LE. 1.D0) THEN

                       CT(i,0,k)%rflxx = 0.D0
                       
                    ELSE

		       IF (R .LE. 1.5D0) THEN


                           !IF (0.55D0/DCOSH(R**4) .GT. vwi) THEN

                              VELY   =  0.55D0/DCOSH(R**4)
                              VELZ   = 0.D0
                              VELX   = 0.D0

                           !ELSE
 
                           !   VELY   =  physic(i,1,k)%vely
                           !   VELZ   =  physic(i,1,k)%velz
                           !   VELX   =  physic(i,1,k)%velx

                           !ENDIF

			   SCRCH1 = 1.D0/DSQRT(1.D0-VELX**2-VELY**2-VELZ**2)

                           CT(i,0,k)%rflxx = 1.D0*SCRCH1*VELX
                           !CT(i,0,k)%rflxx = ((8.8D0-1.D0)/DCOSH(R**4)+1.D0)*SCRCH1*VELX


		       ELSE IF (R .GT. 1.5D0 .AND. physic(i,1,k)%pres .GT. 1.D-8) THEN


                              VELY   =  -physic(i,1,k)%vely
                              VELZ   =  physic(i,1,k)%velz
                              VELX   =  physic(i,1,k)%velx

                           SCRCH1 = 1.D0/DSQRT(1.D0-VELX**2-VELY**2-VELZ**2)

                           CT(i,0,k)%rflxx = physic(i,1,k)%densty*SCRCH1*VELX


		       ELSE
                       
                       ! INCLINATION
                       theta1 = DATAN(zzn(k)/       &
                                DSQRT((yzn(1)+shift-(yznr(1)-yznl(1)))**2+(Rorb-xznl(i))**2))

	               ! AZIMUTHAL ANGLE.
                       theta2 = DATAN((yzn(1)+shift-(yznr(1)-yznl(1)))/(Rorb-xznl(i)))
	

                        !VELX = 0.D0
                        !VELY = 0.D0
                        !VELZ = 0.D0

	
                       VELX = -vwi*DCOS(theta1)*DCOS(theta2)
                       VELY =  vwi*DCOS(theta1)*DSIN(theta2)
                       VELZ =  vwi*DSIN(theta1)

                       SCRCH1 = 1.D0/DSQRT(1.D0-VELX**2-VELY**2-VELZ**2)
                  
                       CT(i,0,k)%rflxx = physic(i,1,k)%densty*SCRCH1*VELX



		       ENDIF


                    ENDIF
                    

                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx

                    R = DSQRT(XZN(I)**2+ZZNL(K)**2)

                    IF (R .LE. 1.D0) THEN

                       CT(i,0,k)%rflxz = 0.D0
                       
                    ELSE
                       
                       IF (R .LE. 1.5D0) THEN


                           !IF (0.55D0/DCOSH(R**4) .GT. vwi) THEN

                              VELY   =  0.55D0/DCOSH(R**4)
                              VELZ   = 0.D0
                              VELX   = 0.D0

                           !ELSE

                           !   VELY   =  physic(i,1,k)%vely
                           !   VELZ   =  physic(i,1,k)%velz
                           !   VELX   =  physic(i,1,k)%velx

                           !ENDIF

                           SCRCH1 = 1.D0/DSQRT(1.D0-VELX**2-VELY**2-VELZ**2)

                           CT(i,0,k)%rflxz = 1.D0*SCRCH1*VELZ


                       ELSE IF (R .GT. 1.5D0 .AND. physic(i,1,k)%pres .GT. 1.D-8) THEN


                              VELY   =  -physic(i,1,k)%vely
                              VELZ   =  physic(i,1,k)%velz
                              VELX   =  physic(i,1,k)%velx

                           SCRCH1 = 1.D0/DSQRT(1.D0-VELX**2-VELY**2-VELZ**2)

                           CT(i,0,k)%rflxz = physic(i,1,k)%densty*SCRCH1*VELZ
                           !CT(i,0,k)%rflxz = ((8.8D0-1.D0)/DCOSH(R**4)+1.D0)*SCRCH1*VELZ


		       ELSE

                       ! INCLINATION
                       theta1 = DATAN(zznl(k)/       &
                                DSQRT((yzn(1)+shift-(yznr(1)-yznl(1)))**2+(Rorb-xzn(i))**2))

	               ! AZIMUTHAL ANGLE.
                       theta2 = DATAN((yzn(1)+shift-(yznr(1)-yznl(1)))/(Rorb-xzn(i)))
	
	
                        !VELX = 0.D0
                        !VELY = 0.D0
                        !VELZ = 0.D0

                       VELX = -vwi*DCOS(theta1)*DCOS(theta2)
                       VELY =  vwi*DCOS(theta1)*DSIN(theta2)
                       VELZ =  vwi*DSIN(theta1)

                       SCRCH1 = 1.D0/DSQRT(1.D0-VELX**2-VELY**2-VELZ**2)
                  
                       CT(i,0,k)%rflxz = physic(i,1,k)%densty*SCRCH1*VELZ

		
		       ENDIF


                    ENDIF                          


                ENDDO
             ENDDO
!$OMP END DO


             
!$OMP END PARALLEL


         

          ELSE
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++


             message1 = 'Error in boundaries of B fluxes.'
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)

          ENDIF

	  ENDIF

          ! ------------
          ! RIGHT Y-BOUNDARY
          ! ------------

	  IF (mpi_coords(2) .EQ. mpi_dims(2)-1) THEN


          ! Reflecting Boundary.

          IF (bndmxy .EQ. 1) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz
                DO  i=1,nx+1
                   CT(i,ny+1,nz)%fbyx = - CT(i,ny,nz)%fbyx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  k=1,nz+1
                DO  i=1,nx
                   CT(i,ny+1,nz)%fbyz = - CT(i,ny,nz)%fbyz
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmxy .EQ. 2) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,ny+1,k)%fbyx = CT(i,ny,k)%fbyx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,ny+1,k)%fbyz = CT(i,ny,k)%fbyz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,ny+1,k)%rflxx = CT(i,ny,k)%rflxx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,ny+1,k)%rflxz = CT(i,ny,k)%rflxz
                ENDDO
             ENDDO
!$OMP END DO
             

!$OMP END PARALLEL

          ELSE IF (bndmxy .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,K)

!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,ny+1,k)%fbyx = CT(i,1,k)%fbyx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,ny+1,k)%fbyz = CT(i,1,k)%fbyz
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz
                DO i=1,nx+1
                   CT(i,ny+1,k)%rflxx = CT(i,1,k)%rflxx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO k=1,nz+1
                DO i=1,nx
                   CT(i,ny+1,k)%rflxz = CT(i,1,k)%rflxz
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL


          ELSE

             message1 = 'Error in boundaries of B fluxes.'
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)

          ENDIF

	  ENDIF

          ! ------------
          ! LEFT Z-BOUNDARY
          ! ------------

	  IF (mpi_coords(3) .EQ. 0) THEN

          ! Reflecting Boundary.

          IF (bndmnz .EQ. 1) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,J)

!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny
                DO  i=1,nx+1
                   CT(i,j,0)%fbzx = - CT(i,j,1)%fbzx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny+1
                DO  i=1,nx
                   CT(i,j,0)%fbzy = - CT(i,j,1)%fbzy
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmnz .EQ. 2) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,J)

!$OMP DO SCHEDULE(STATIC)
             DO  j=1,ny
                DO  i=1,nx+1
                   CT(i,j,0)%fbzx =  CT(i,j,1)%fbzx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)
             DO  j=1,ny+1
                DO  i=1,nx
                   CT(i,j,0)%fbzy =  CT(i,j,1)%fbzy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)
             DO  j=1,ny
                DO  i=1,nx+1
                   CT(i,j,0)%rflxx =  CT(i,j,1)%rflxx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)
             DO  j=1,ny+1
                DO  i=1,nx
                   CT(i,j,0)%rflxy =  CT(i,j,1)%rflxy
                ENDDO
             ENDDO
!$OMP END DO
             
!$OMP END PARALLEL


          ! Periodic Boundary.

          ELSE IF (bndmnz .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,J)

!$OMP DO SCHEDULE(STATIC)
             DO  j=1,ny
                DO  i=1,nx+1
                   CT(i,j,0)%fbzx =  CT(i,j,nz)%fbzx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)
             DO  j=1,ny+1
                DO  i=1,nx
                   CT(i,j,0)%fbzy =  CT(i,j,nz)%fbzy
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL


          ELSE

             message1 = 'Error in boundaries of B fluxes.'
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)

          ENDIF

	  ENDIF

          ! ------------
          ! RIGHT Z-BOUNDARY
          ! ------------

	  IF (mpi_coords(3) .EQ. mpi_dims(3)-1) THEN


          ! Reflecting Boundary.

          IF (bndmxz .EQ. 1) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,J)

!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny
                DO  i=1,nx+1
                   CT(i,j,nz+1)%fbzx = - CT(i,j,nz)%fbzx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny+1
                DO  i=1,nx
                   CT(i,j,nz+1)%fbzy = - CT(i,j,nz)%fbzy
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL

          ! Zero-gradient Boundary.

          ELSE IF (bndmxz .EQ. 2) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,J)

!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny
                DO  i=1,nx+1
                   CT(i,j,nz+1)%fbzx =  CT(i,j,nz)%fbzx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny+1
                DO  i=1,nx
                   CT(i,j,nz+1)%fbzy =  CT(i,j,nz)%fbzy
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny
                DO  i=1,nx+1
                   CT(i,j,nz+1)%rflxx =  CT(i,j,nz)%rflxx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny+1
                DO  i=1,nx
                   CT(i,j,nz+1)%rflxy =  CT(i,j,nz)%rflxy
                ENDDO
             ENDDO
!$OMP END DO


!$OMP END PARALLEL


          ! Periodic Boundary.

          ELSE IF (bndmxz .EQ. 3) THEN

!$OMP PARALLEL DEFAULT (SHARED) &
!$OMP PRIVATE (I,J)

!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny
                DO  i=1,nx+1
                   CT(i,j,nz+1)%fbzx =  CT(i,j,1)%fbzx
                ENDDO
             ENDDO
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE (STATIC)
             DO  j=1,ny+1
                DO  i=1,nx
                   CT(i,j,nz+1)%fbzy =  CT(i,j,1)%fbzy
                ENDDO
             ENDDO
!$OMP END DO

!$OMP END PARALLEL


          ELSE

             message1 = 'Error in boundaries of B fluxes.'
             call ERROR ('SWEEPS', message1, __FILE__, __LINE__)

          ENDIF

	  ENDIF 


	  IF (mpi_dims(1)*mpi_dims(2)*mpi_dims(3) .NE. 1) THEN

              call PAR_BARRIER()

              call INTERCAMBIABOUNDCT (CT, nx, ny, nz, mnx6, mny6, mnz6)

	      call PAR_BARRIER ()

              call INTERCAMBIABOUNDCTT (CT, nx, ny, nz, mnx6, mny6, mnz6)

	  ENDIF


!        =======================================================================


          ! CORNER 0-0-K
          ! ----------------------------

          IF (mpi_coords(1) .EQ. 0 .AND. mpi_coords(2) .EQ. 0) THEN

             DO k= 1, nz

                physic(0,0,k)%velx = (physic(1,0,k)%velx+physic(0,1,k)%velx)/2.D0
                physic(0,0,k)%vely = (physic(1,0,k)%vely+physic(0,1,k)%vely)/2.D0
                physic(0,0,k)%velz = (physic(1,0,k)%velz+physic(0,1,k)%velz)/2.D0

                physic(0,0,k)%bxc2 = (physic(1,0,k)%bxc2+physic(0,1,k)%bxc2)/2.D0
                physic(0,0,k)%byc2 = (physic(1,0,k)%byc2+physic(0,1,k)%byc2)/2.D0
                physic(0,0,k)%bzc2 = (physic(1,0,k)%bzc2+physic(0,1,k)%bzc2)/2.D0

             ENDDO

          ENDIF

          IF (mpi_coords(1) .EQ. mpi_dims(1)-1 .AND. mpi_coords(2) .EQ. 0) THEN

             DO k= 1, nz

                physic(nx+1,0,k)%velx = (physic(nx,0,k)%velx+physic(nx+1,1,k)%velx)/2.D0
                physic(nx+1,0,k)%vely = (physic(nx,0,k)%vely+physic(nx+1,1,k)%vely)/2.D0
                physic(nx+1,0,k)%velz = (physic(nx,0,k)%velz+physic(nx+1,1,k)%velz)/2.D0

                physic(nx+1,0,k)%bxc2 = (physic(nx,0,k)%bxc2+physic(nx+1,1,k)%bxc2)/2.D0
                physic(nx+1,0,k)%byc2 = (physic(nx,0,k)%byc2+physic(nx+1,1,k)%byc2)/2.D0
                physic(nx+1,0,k)%bzc2 = (physic(nx,0,k)%bzc2+physic(nx+1,1,k)%bzc2)/2.D0

             ENDDO

          ENDIF



          IF (mpi_coords(1) .EQ. 0 .AND. mpi_coords(2) .EQ. mpi_dims(2)-1) THEN

             DO k= 1, nz

                physic(0,ny+1,k)%velx = (physic(0,ny,k)%velx+physic(1,ny+1,k)%velx)/2.D0
                physic(0,ny+1,k)%vely = (physic(0,ny,k)%vely+physic(1,ny+1,k)%vely)/2.D0
                physic(0,ny+1,k)%velz = (physic(0,ny,k)%velz+physic(1,ny+1,k)%velz)/2.D0

                physic(0,ny+1,k)%bxc2 = (physic(0,ny,k)%bxc2+physic(1,ny+1,k)%bxc2)/2.D0
                physic(0,ny+1,k)%byc2 = (physic(0,ny,k)%byc2+physic(1,ny+1,k)%byc2)/2.D0
                physic(0,ny+1,k)%bzc2 = (physic(0,ny,k)%bzc2+physic(1,ny+1,k)%bzc2)/2.D0

             ENDDO

          ENDIF



          IF (mpi_coords(1) .EQ. mpi_dims(1)-1 .AND. mpi_coords(2) .EQ. mpi_dims(2)-1) THEN

             DO k= 1, nz

                physic(nx+1,ny+1,k)%velx = (physic(nx+1,ny,k)%velx+physic(nx,ny+1,k)%velx)/2.D0
                physic(nx+1,ny+1,k)%vely = (physic(nx+1,ny,k)%vely+physic(nx,ny+1,k)%vely)/2.D0
                physic(nx+1,ny+1,k)%velz = (physic(nx+1,ny,k)%velz+physic(nx,ny+1,k)%velz)/2.D0

                physic(nx+1,ny+1,k)%bxc2 = (physic(nx+1,ny,k)%bxc2+physic(nx,ny+1,k)%bxc2)/2.D0
                physic(nx+1,ny+1,k)%byc2 = (physic(nx+1,ny,k)%byc2+physic(nx,ny+1,k)%byc2)/2.D0
                physic(nx+1,ny+1,k)%bzc2 = (physic(nx+1,ny,k)%bzc2+physic(nx,ny+1,k)%bzc2)/2.D0

             ENDDO

          ENDIF


          !IF (mpi_coords(1) .EQ. 0) THEN

	  !   IF (bndmnx .EQ. 2) THEN

          !   DO k= 1, nz

          !      physic(0,0,k)%velx = physic(1,0,k)%velx
          !      physic(0,0,k)%vely = physic(1,0,k)%vely
          !      physic(0,0,k)%velz = physic(1,0,k)%velz

          !      physic(0,0,k)%bxc2 = physic(1,0,k)%bxc2
          !      physic(0,0,k)%byc2 = physic(1,0,k)%byc2
          !      physic(0,0,k)%bzc2 = physic(1,0,k)%bzc2

          !   ENDDO


          !   DO k= 1, nz

          !      physic(0,ny+1,k)%velx =  physic(1,ny+1,k)%velx
          !      physic(0,ny+1,k)%vely =  physic(1,ny+1,k)%vely
          !      physic(0,ny+1,k)%velz =  physic(1,ny+1,k)%velz

          !      physic(0,ny+1,k)%bxc2 = physic(1,ny+1,k)%bxc2
          !      physic(0,ny+1,k)%byc2 = physic(1,ny+1,k)%byc2
          !      physic(0,ny+1,k)%bzc2 = physic(1,ny+1,k)%bzc2

          !   ENDDO

	  !   ENDIF

          !ENDIF


          !IF (mpi_coords(1) .EQ. mpi_dims(1)-1) THEN


	  !   IF (bndmxx .EQ. 2 .OR. bndmxx .EQ. 9) THEN

          !   DO k= 1, nz

          !      physic(nx+1,0,k)%velx = physic(nx,0,k)%velx
          !      physic(nx+1,0,k)%vely = physic(nx,0,k)%vely
          !      physic(nx+1,0,k)%velz = physic(nx,0,k)%velz

          !      physic(nx+1,0,k)%bxc2 = physic(nx,0,k)%bxc2
          !      physic(nx+1,0,k)%byc2 = physic(nx,0,k)%byc2
          !      physic(nx+1,0,k)%bzc2 = physic(nx,0,k)%bzc2

          !   ENDDO


          !   DO k= 1, nz

          !      physic(nx+1,ny+1,k)%velx =  physic(nx,ny+1,k)%velx     
          !      physic(nx+1,ny+1,k)%vely =  physic(nx,ny+1,k)%vely     
          !      physic(nx+1,ny+1,k)%velz =  physic(nx,ny+1,k)%velz     

          !      physic(nx+1,ny+1,k)%bxc2 = physic(nx,ny+1,k)%bxc2    
          !      physic(nx+1,ny+1,k)%byc2 = physic(nx,ny+1,k)%byc2     
          !      physic(nx+1,ny+1,k)%bzc2 = physic(nx,ny+1,k)%bzc2     

          !   ENDDO

	  !   ENDIF

          !ENDIF


	  
          call PAR_BARRIER()

          call INTERCAMBIABOUNDFIELDK (physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy)

          call PAR_BARRIER()


          ! CORNER I-0-0
          ! ----------------------------


          IF (mpi_coords(2) .EQ. 0 .AND. mpi_coords(3) .EQ. 0) THEN

             DO i= 1, nx

                physic(i,0,0)%velx = (physic(i,1,0)%velx+physic(i,0,1)%velx)/2.D0
                physic(i,0,0)%vely = (physic(i,1,0)%vely+physic(i,0,1)%vely)/2.D0
                physic(i,0,0)%velz = (physic(i,1,0)%velz+physic(i,0,1)%velz)/2.D0

                physic(i,0,0)%bxc2 = (physic(i,1,0)%bxc2+physic(i,0,1)%bxc2)/2.D0
                physic(i,0,0)%byc2 = (physic(i,1,0)%byc2+physic(i,0,1)%byc2)/2.D0
                physic(i,0,0)%bzc2 = (physic(i,1,0)%bzc2+physic(i,0,1)%bzc2)/2.D0

             ENDDO

          ENDIF

          IF (mpi_coords(2) .EQ. mpi_dims(2)-1 .AND. mpi_coords(3) .EQ. 0) THEN

             DO i= 1, nx

                physic(i,ny+1,0)%velx = (physic(i,ny,0)%velx+physic(i,ny+1,1)%velx)/2.D0
                physic(i,ny+1,0)%vely = (physic(i,ny,0)%vely+physic(i,ny+1,1)%vely)/2.D0
                physic(i,ny+1,0)%velz = (physic(i,ny,0)%velz+physic(i,ny+1,1)%velz)/2.D0

                physic(i,ny+1,0)%bxc2 = (physic(i,ny,0)%bxc2+physic(i,ny+1,1)%bxc2)/2.D0
                physic(i,ny+1,0)%byc2 = (physic(i,ny,0)%byc2+physic(i,ny+1,1)%byc2)/2.D0
                physic(i,ny+1,0)%bzc2 = (physic(i,ny,0)%bzc2+physic(i,ny+1,1)%bzc2)/2.D0

             ENDDO

          ENDIF


          IF (mpi_coords(2) .EQ. 0 .AND. mpi_coords(3) .EQ. mpi_dims(3)-1) THEN

             DO i= 1, nx

                physic(i,0,nz+1)%velx = (physic(i,0,nz)%velx+physic(i,1,nz+1)%velx)/2.D0
                physic(i,0,nz+1)%vely = (physic(i,0,nz)%vely+physic(i,1,nz+1)%vely)/2.D0
                physic(i,0,nz+1)%velz = (physic(i,0,nz)%velz+physic(i,1,nz+1)%velz)/2.D0

                physic(i,0,nz+1)%bxc2 = (physic(i,0,nz)%bxc2+physic(i,1,nz+1)%bxc2)/2.D0
                physic(i,0,nz+1)%byc2 = (physic(i,0,nz)%byc2+physic(i,1,nz+1)%byc2)/2.D0
                physic(i,0,nz+1)%bzc2 = (physic(i,0,nz)%bzc2+physic(i,1,nz+1)%bzc2)/2.D0

             ENDDO

          ENDIF


          IF (mpi_coords(2) .EQ. mpi_dims(2)-1 .AND. mpi_coords(3) .EQ. mpi_dims(3)-1) THEN

             DO i= 1, nx

                physic(i,ny+1,nz+1)%velx = (physic(i,ny+1,nz)%velx+physic(i,ny,nz+1)%velx)/2.D0
                physic(i,ny+1,nz+1)%vely = (physic(i,ny+1,nz)%vely+physic(i,ny,nz+1)%vely)/2.D0
                physic(i,ny+1,nz+1)%velz = (physic(i,ny+1,nz)%velz+physic(i,ny,nz+1)%velz)/2.D0

                physic(i,ny+1,nz+1)%bxc2 = (physic(i,ny+1,nz)%bxc2+physic(i,ny,nz+1)%bxc2)/2.D0
                physic(i,ny+1,nz+1)%byc2 = (physic(i,ny+1,nz)%byc2+physic(i,ny,nz+1)%byc2)/2.D0
                physic(i,ny+1,nz+1)%bzc2 = (physic(i,ny+1,nz)%bzc2+physic(i,ny,nz+1)%bzc2)/2.D0

             ENDDO

          ENDIF



          !IF (mpi_coords(2) .EQ. 0) THEN


	  !   IF (bndmny .EQ. 2) THEN


          !   DO i= 1, nx


          !     physic(i,0,0)%velx = physic(i,1,0)%velx
          !     physic(i,0,0)%vely = physic(i,1,0)%vely
          !     physic(i,0,0)%velz = physic(i,1,0)%velz


          !     physic(i,0,0)%bxc2 = physic(i,1,0)%bxc2
          !     physic(i,0,0)%byc2 = physic(i,1,0)%byc2
          !     physic(i,0,0)%bzc2 = physic(i,1,0)%bzc2


          !   ENDDO


          !   DO i= 1, nx

          !      physic(i,0,nz+1)%velx =  physic(i,1,nz+1)%velx
          !      physic(i,0,nz+1)%vely =  physic(i,1,nz+1)%vely
          !      physic(i,0,nz+1)%velz =  physic(i,1,nz+1)%velz

          !      physic(i,0,nz+1)%bxc2 = physic(i,1,nz+1)%bxc2
          !      physic(i,0,nz+1)%byc2 = physic(i,1,nz+1)%byc2
          !      physic(i,0,nz+1)%bzc2 = physic(i,1,nz+1)%bzc2

          !   ENDDO


	   !  ELSE IF (bndmny .EQ. 7) THEN


           !  DO i= 1, nx

           !    R = DSQRT(XZN(I)**2+(ZZN(1)-(ZZNR(1)-ZZNL(1)))**2)

           !    IF (R .LE. 1.D0) THEN

           !        physic(i,0,0)%velx = 0.D0
           !        physic(i,0,0)%vely = 0.55D0
           !        physic(i,0,0)%velz = 0.D0

           !    ELSE

           !        physic(i,0,0)%velx = physic(i,1,0)%velx
           !        physic(i,0,0)%vely = physic(i,1,0)%vely
           !        physic(i,0,0)%velz = physic(i,1,0)%velz

	   !    ENDIF


           !     physic(i,0,0)%bxc2 = physic(i,0,0)%bxc
           !     physic(i,0,0)%byc2 = 0.D0
           !     physic(i,0,0)%bzc2 = physic(i,0,0)%bzc

           !  ENDDO



           !  DO i= 1, nx

           !    R = DSQRT(XZN(I)**2+(ZZN(NZ)+(ZZNR(NZ)-ZZNL(NZ)))**2)

           !    IF (R .LE. 1.D0) THEN

           !        physic(i,0,NZ+1)%velx = 0.D0
           !        physic(i,0,NZ+1)%vely = 0.55D0
           !        physic(i,0,NZ+1)%velz = 0.D0

           !    ELSE

           !        physic(i,0,NZ+1)%velx = physic(i,1,NZ+1)%velx
           !        physic(i,0,NZ+1)%vely = physic(i,1,NZ+1)%vely
           !        physic(i,0,NZ+1)%velz = physic(i,1,NZ+1)%velz

           !    ENDIF


           !     physic(i,0,nz+1)%bxc2 = physic(i,0,NZ+1)%bxc 
           !     physic(i,0,nz+1)%byc2 = 0.D0              
           !     physic(i,0,nz+1)%bzc2 = physic(i,0,NZ+1)%bzc 

           !  ENDDO


	   !  ENDIF

          !ENDIF


          !IF (mpi_coords(2) .EQ. mpi_dims(2)-1) THEN

	  !   IF (bndmxy .EQ. 2) THEN

          !   DO i= 1, nx

                !physic(i,ny+1,0)%velx = physic(i,ny,0)%velx
                !physic(i,ny+1,0)%vely = physic(i,ny,0)%vely
                !physic(i,ny+1,0)%velz = physic(i,ny,0)%velz

                !physic(i,ny+1,0)%bxc2 = physic(i,ny,0)%bxc2
                !physic(i,ny+1,0)%byc2 = physic(i,ny,0)%byc2
                !physic(i,ny+1,0)%bzc2 = physic(i,ny,0)%bzc2

             !ENDDO


           !  DO i=1, nx

           !     physic(i,ny+1,nz+1)%velx =  physic(i,ny,nz+1)%velx
           !     physic(i,ny+1,nz+1)%vely =  physic(i,ny,nz+1)%vely
           !     physic(i,ny+1,nz+1)%velz =  physic(i,ny,nz+1)%velz

           !     physic(i,ny+1,nz+1)%bxc2 = physic(i,ny,nz+1)%bxc2
           !     physic(i,ny+1,nz+1)%byc2 = physic(i,ny,nz+1)%byc2
           !     physic(i,ny+1,nz+1)%bzc2 = physic(i,ny,nz+1)%bzc2

           !  ENDDO

	   !  ENDIF

          !ENDIF


	  call PAR_BARRIER()

          call INTERCAMBIABOUNDFIELDI (physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy)

          call PAR_BARRIER()


          ! CORNER 0-J-0
          ! ----------------------------


          IF (mpi_coords(3) .EQ. 0 .AND. mpi_coords(1) .EQ. 0) THEN

             DO j= 1, ny

                physic(0,j,0)%velx = (physic(0,j,1)%velx+physic(1,j,0)%velx)/2.D0
                physic(0,j,0)%vely = (physic(0,j,1)%vely+physic(1,j,0)%vely)/2.D0
                physic(0,j,0)%velz = (physic(0,j,1)%velz+physic(1,j,0)%velz)/2.D0

                physic(0,j,0)%bxc2 = (physic(0,j,1)%bxc2+physic(1,j,0)%bxc2)/2.D0
                physic(0,j,0)%byc2 = (physic(0,j,1)%byc2+physic(1,j,0)%byc2)/2.D0
                physic(0,j,0)%bzc2 = (physic(0,j,1)%bzc2+physic(1,j,0)%bzc2)/2.D0

             ENDDO

          ENDIF



          IF (mpi_coords(3) .EQ. mpi_dims(3)-1 .AND. mpi_coords(1) .EQ. 0) THEN

             DO j= 1, ny

                physic(0,j,nz+1)%velx = (physic(0,j,nz)%velx+physic(1,j,nz+1)%velx)/2.D0
                physic(0,j,nz+1)%vely = (physic(0,j,nz)%vely+physic(1,j,nz+1)%vely)/2.D0
                physic(0,j,nz+1)%velz = (physic(0,j,nz)%velz+physic(1,j,nz+1)%velz)/2.D0

                physic(0,j,nz+1)%bxc2 = (physic(0,j,nz)%bxc2+physic(1,j,nz+1)%bxc2)/2.D0
                physic(0,j,nz+1)%byc2 = (physic(0,j,nz)%byc2+physic(1,j,nz+1)%byc2)/2.D0
                physic(0,j,nz+1)%bzc2 = (physic(0,j,nz)%bzc2+physic(1,j,nz+1)%bzc2)/2.D0

             ENDDO

          ENDIF


          IF (mpi_coords(3) .EQ. 0 .AND. mpi_coords(1) .EQ. mpi_dims(1)-1) THEN

             DO j= 1, ny

                physic(nx+1,j,0)%velx = (physic(nx,j,0)%velx+physic(nx+1,j,1)%velx)/2.D0
                physic(nx+1,j,0)%vely = (physic(nx,j,0)%vely+physic(nx+1,j,1)%vely)/2.D0
                physic(nx+1,j,0)%velz = (physic(nx,j,0)%velz+physic(nx+1,j,1)%velz)/2.D0

                physic(nx+1,j,0)%bxc2 = (physic(nx,j,0)%bxc2+physic(nx+1,j,1)%bxc2)/2.D0
                physic(nx+1,j,0)%byc2 = (physic(nx,j,0)%byc2+physic(nx+1,j,1)%byc2)/2.D0
                physic(nx+1,j,0)%bzc2 = (physic(nx,j,0)%bzc2+physic(nx+1,j,1)%bzc2)/2.D0

             ENDDO

          ENDIF


          IF (mpi_coords(3) .EQ. mpi_dims(3)-1 .AND. mpi_coords(1) .EQ. mpi_dims(1)-1) THEN

             DO j= 1, ny

                physic(nx+1,j,nz+1)%velx = (physic(nx,j,nz+1)%velx+physic(nx+1,j,nz)%velx)/2.D0
                physic(nx+1,j,nz+1)%vely = (physic(nx,j,nz+1)%vely+physic(nx+1,j,nz)%vely)/2.D0
                physic(nx+1,j,nz+1)%velz = (physic(nx,j,nz+1)%velz+physic(nx+1,j,nz)%velz)/2.D0

                physic(nx+1,j,nz+1)%bxc2 = (physic(nx,j,nz+1)%bxc2+physic(nx+1,j,nz)%bxc2)/2.D0
                physic(nx+1,j,nz+1)%byc2 = (physic(nx,j,nz+1)%byc2+physic(nx+1,j,nz)%byc2)/2.D0
                physic(nx+1,j,nz+1)%bzc2 = (physic(nx,j,nz+1)%bzc2+physic(nx+1,j,nz)%bzc2)/2.D0

             ENDDO

          ENDIF



          !IF (mpi_coords(3) .EQ. 0) THEN

	  !   IF (bndmnz .EQ. 2) THEN

          !   DO j= 1, ny

          !      physic(0,j,0)%velx = physic(0,j,1)%velx
          !      physic(0,j,0)%vely = physic(0,j,1)%vely
          !      physic(0,j,0)%velz = physic(0,j,1)%velz

          !      physic(0,j,0)%bxc2 = physic(0,j,1)%bxc2
          !      physic(0,j,0)%byc2 = physic(0,j,1)%byc2
          !      physic(0,j,0)%bzc2 = physic(0,j,1)%bzc2

          !   ENDDO


         !    DO j= 1, ny

         !       physic(nx+1,j,0)%velx =  physic(nx+1,j,1)%velx
         !       physic(nx+1,j,0)%vely =  physic(nx+1,j,1)%vely
         !       physic(nx+1,j,0)%velz =  physic(nx+1,j,1)%velz

         !       physic(nx+1,j,0)%bxc2 = physic(nx+1,j,1)%bxc2
         !       physic(nx+1,j,0)%byc2 = physic(nx+1,j,1)%byc2
         !       physic(nx+1,j,0)%bzc2 = physic(nx+1,j,1)%bzc2

         !    ENDDO

	 !    ENDIF

         ! ENDIF


         ! IF (mpi_coords(3) .EQ. mpi_dims(3)-1) THEN

	 !    IF (bndmxz .EQ. 2) THEN

         !    DO j= 1, ny

         !       physic(0,j,nz+1)%velx = physic(0,j,nz)%velx
         !       physic(0,j,nz+1)%vely = physic(0,j,nz)%vely
         !       physic(0,j,nz+1)%velz = physic(0,j,nz)%velz

         !       physic(0,j,nz+1)%bxc2 = physic(0,j,nz)%bxc2
         !       physic(0,j,nz+1)%byc2 = physic(0,j,nz)%byc2
         !       physic(0,j,nz+1)%bzc2 = physic(0,j,nz)%bzc2

         !    ENDDO


         !    DO j=1, ny

         !       physic(nx+1,j,nz+1)%velx =  physic(nx+1,j,nz)%velx
         !       physic(nx+1,j,nz+1)%vely =  physic(nx+1,j,nz)%vely
         !       physic(nx+1,j,nz+1)%velz =  physic(nx+1,j,nz)%velz

         !       physic(nx+1,j,nz+1)%bxc2 = physic(nx+1,j,nz)%bxc2
         !       physic(nx+1,j,nz+1)%byc2 = physic(nx+1,j,nz)%byc2
         !       physic(nx+1,j,nz+1)%bzc2 = physic(nx+1,j,nz)%bzc2

         !    ENDDO

	 !    ENDIF

         ! ENDIF


	  call PAR_BARRIER()

          call INTERCAMBIABOUNDFIELDJ (physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy)

          call PAR_BARRIER()



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



          
          ! OMEGA - Z .............................

                             
          !!! Electromotive forces Ω in the corners.

!$OMP DO SCHEDULE(static)
          DO k=1, nz
             DO j= 1,  ny+1
                DO i= 1, nx+1

                   ! ====================
		   ! ARITHMETIC AVERAGING.
                   ! ====================
                   
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
          ENDDO
!$OMP END DO

!$OMP BARRIER

!$OMP DO SCHEDULE(static)
          DO k=1, nz+1
             DO j= 1, ny
                DO i= 1, nx+1


                   ! ====================
		   ! ARITHMETIC AVERAGING.
                   ! ====================
                  
                   CT(i,j,k)%omega_y = 0.25D0*(CT(i-1,j,k)%fbxz + CT(i,j,k)%fbxz - &
                                               CT(i,j,k-1)%fbzx-CT(i,j,k)%fbzx)



             
                   ! ====================
		   ! CONTACT-CT (UPWIND)
                   ! ====================


		   ! SCRCH1;
		   ! ...................................		

		   IF (CT(I-1,J,K)%RFLXZ .GT. 1.d-15) THEN

                        FLUXC  = -(physic(i-1,j,k-1)%velx*physic(i-1,j,k-1)%bzc2-physic(i-1,j,k-1)%bxc2*physic(i-1,j,k-1)%velz)
			SCRCH1 = 2.D0/(XZNR(I)-XZNL(I))*(-CT(I,J,K-1)%fbzx-FLUXC)

		   ELSE IF (CT(I-1,J,K)%RFLXZ .LT. -1.D-15) THEN

                        FLUXC  = -(physic(i-1,j,k)%velx*physic(i-1,j,k)%bzc2-physic(i-1,j,k)%bxc2*physic(i-1,j,k)%velz)                
                        SCRCH1 = 2.D0/(XZNR(I)-XZNL(I))*(-CT(I,J,K)%fbzx-FLUXC)

	           ELSE

                      SCRCH11 = 2.D0/(XZNR(I)-XZNL(I))*(-CT(I,J,K-1)%fbzx+(physic(i-1,j,k-1)%velx*physic(i-1,j,k-1)%bzc2-  &
                                                                           physic(i-1,j,k-1)%bxc2*physic(i-1,j,k-1)%velz))	
                      SCRCH12 = 2.D0/(XZNR(I)-XZNL(I))*(-CT(I,J,K)%fbzx+(physic(i-1,j,k)%velx*physic(i-1,j,k)%bzc2-        &
                                                                         physic(i-1,j,k)%bxc2*physic(i-1,j,k)%velz))
		      SCRCH1 = 0.5D0*(SCRCH11+SCRCH12)			

		   ENDIF

                   
                   ! SCRCH2;
                   ! ...................................


		   IF (CT(I-1,J,K)%RFLXZ .GT. 1.D-15) THEN

                        FLUXC  = -(physic(i,j,k-1)%velx*physic(i,j,k-1)%bzc2-physic(i,j,k-1)%bxc2*physic(i,j,k-1)%velz)
			SCRCH2 = 2.D0/(XZNR(I)-XZNL(I))*(FLUXC+CT(I,J,K-1)%fbzx)

		   ELSE IF (CT(I-1,J,K)%RFLXZ .LT. -1.D-15) THEN

                        FLUXC  = -(physic(i,j,k)%velx*physic(i,j,k)%bzc2-physic(i,j,k)%bxc2*physic(i,j,k)%velz)                      
                        SCRCH2 = 2.D0/(XZNR(I)-XZNL(I))*(FLUXC+CT(I,J,K)%fbzx)

	           ELSE

                        SCRCH21 = 2.D0/(XZNR(I)-XZNL(I))*(-(physic(i,j,k-1)%velx*physic(i,j,k-1)%bzc2-                &
                                                            physic(i,j,k-1)%bxc2*physic(i,j,k-1)%velz)+CT(I,J,K-1)%fbzx)	
                        SCRCH22 = 2.D0/(XZNR(I)-XZNL(I))*(-(physic(i,j,k)%velx*physic(i,j,k)%bzc2-                    &
                                                            physic(i,j,k)%bxc2*physic(i,j,k)%velz)+CT(I,J,K)%fbzx)
			SCRCH2 = 0.5D0*(SCRCH21+SCRCH22)			

		   ENDIF


                   
		   ! SCRCH4;
		   ! ...................................		

		   IF (CT(I,J,K-1)%RFLXX .GT. 1.D-15) THEN

                        FLUXC  = (physic(i-1,j,k-1)%velz*physic(i-1,j,k-1)%bxc2-physic(i-1,j,k-1)%bzc2*physic(i-1,j,k-1)%velx)
			SCRCH4 = 2.D0/(ZZNR(K)-ZZNL(K))*(CT(I-1,J,K)%fbxz-FLUXC)

		   ELSE IF (CT(I,J,K-1)%RFLXX .LT. -1.D-15) THEN

                        FLUXC  = (physic(i,j,k-1)%velz*physic(i,j,k-1)%bxc2-physic(i,j,k-1)%bzc2*physic(i,j,k-1)%velx)                
                        SCRCH4 = 2.D0/(ZZNR(K)-ZZNL(K))*(CT(I,J,K)%fbxz-FLUXC)

	           ELSE

                      SCRCH41 = 2.D0/(ZZNR(K)-ZZNL(K))*(CT(I-1,J,K)%fbxz-(physic(i-1,j,k-1)%velz*physic(i-1,j,k-1)%bxc2-  &
                                                                           physic(i-1,j,k-1)%bzc2*physic(i-1,j,k-1)%velx))	
                      SCRCH42 = 2.D0/(ZZNR(K)-ZZNL(K))*(CT(I,J,K)%fbxz-(physic(i,j,k-1)%velz*physic(i,j,k-1)%bxc2-        &
                                                                        physic(i,j,k-1)%bzc2*physic(i,j,k-1)%velx))
		      SCRCH4 = 0.5D0*(SCRCH41+SCRCH42)			

		   ENDIF

                   
                   ! SCRCH5;
                   ! ...................................


		   IF (CT(I,J,K-1)%RFLXX .GT. 1.D-15) THEN

                        FLUXC  = (physic(i-1,j,k)%velz*physic(i-1,j,k)%bxc2-physic(i-1,j,k)%bzc2*physic(i-1,j,k)%velx)
			SCRCH5 = 2.D0/(ZZNR(K)-ZZNL(K))*(FLUXC-CT(I-1,J,K)%fbxz)

		   ELSE IF (CT(I,J,K-1)%RFLXX .LT. -1.D-15) THEN

                        FLUXC  = (physic(i,j,k)%velz*physic(i,j,k)%bxc2-physic(i,j,k)%bzc2*physic(i,j,k)%velx)                      
                        SCRCH5 = 2.D0/(ZZNR(K)-ZZNL(K))*(FLUXC-CT(I,J,K)%fbxz)

	           ELSE

                        SCRCH51 = 2.D0/(ZZNR(K)-ZZNL(K))*((physic(i-1,j,k)%velz*physic(i-1,j,k)%bxc2-                &
                                                            physic(i-1,j,k)%bzc2*physic(i-1,j,k)%velx)-CT(I-1,J,K)%fbxz)	
                        SCRCH52 = 2.D0/(ZZNR(K)-ZZNL(K))*((physic(i,j,k)%velz*physic(i,j,k)%bxc2-                    &
                                                            physic(i,j,k)%bzc2*physic(i,j,k)%velx)-CT(I,J,K)%fbxz)
			SCRCH5 = 0.5D0*(SCRCH51+SCRCH52)			

		   ENDIF


                   SCRCH3 = SCRCH1 - SCRCH2		
		   SCRCH6 = SCRCH4 - SCRCH5




		   CT(i,j,k)%omega_y = CT(i,j,k)%omega_y + SCRCH3*(XZNR(I)-XZNL(I))/8.D0 + &
                                                           SCRCH6*(ZZNR(K)-ZZNL(K))/8.D0


                   
                ENDDO
             ENDDO
          ENDDO
!$OMP END DO

!$OMP BARRIER

!$OMP DO SCHEDULE(static)
          DO k=1, nz+1
             DO j= 1, ny+1
                DO i= 1, nx

                   ! ====================
		   ! ARITHMETIC AVERAGING.
                   ! ====================
                   
                   CT(i,j,k)%omega_x = 0.25D0*(CT(i,j,k-1)%fbzy + CT(i,j,k)%fbzy - &
                                               CT(i,j-1,k)%fbyz-CT(i,j,k)%fbyz)


                             
                   ! ====================
		   ! CONTACT-CT (UPWIND)
                   ! ====================


		   ! SCRCH1;
		   ! ...................................		

		   IF (CT(I,J,K-1)%RFLXY .GT. 1.d-15) THEN

                        FLUXC  = -(physic(i,j-1,k-1)%velz*physic(i,j-1,k-1)%byc2-physic(i,j-1,k-1)%bzc2*physic(i,j-1,k-1)%vely)
			SCRCH1 = 2.D0/(ZZNR(k)-ZZNL(K))*(-CT(I,J-1,K)%fbyz-FLUXC)

		   ELSE IF (CT(I,J,K-1)%RFLXY .LT. -1.D-15) THEN

                        FLUXC  = -(physic(i,j,k-1)%velz*physic(i,j,k-1)%byc2-physic(i,j,k-1)%bzc2*physic(i,j,k-1)%vely)                
                        SCRCH1 = 2.D0/(ZZNR(K)-ZZNL(K))*(-CT(I,J,K)%fbyz-FLUXC)

	           ELSE

                      SCRCH11 = 2.D0/(ZZNR(K)-ZZNL(K))*(-CT(I,J-1,K)%fbyz+(physic(i,j-1,k-1)%velz*physic(i,j-1,k-1)%byc2-  &
                                                                           physic(i,j-1,k-1)%bzc2*physic(i,j-1,k-1)%vely))	
                      SCRCH12 = 2.D0/(ZZNR(K)-ZZNL(K))*(-CT(I,J,K)%fbyz+(physic(i,j,k-1)%velz*physic(i,j,k-1)%byc2-        &
                                                                         physic(i,j,k-1)%bzc2*physic(i,j,k-1)%vely))
		      SCRCH1 = 0.5D0*(SCRCH11+SCRCH12)			

		   ENDIF

                   
                   ! SCRCH2;
                   ! ...................................


		   IF (CT(I,J,K-1)%RFLXY .GT. 1.D-15) THEN

                        FLUXC  = -(physic(i,j-1,k)%velz*physic(i,j-1,k)%byc2-physic(i,j-1,k)%bzc2*physic(i,j-1,k)%vely)
			SCRCH2 = 2.D0/(ZZNR(K)-ZZNL(K))*(FLUXC+CT(I,J-1,K)%fbyz)

		   ELSE IF (CT(I,J,K-1)%RFLXY .LT. -1.D-15) THEN

                        FLUXC  = -(physic(i,j,k)%velz*physic(i,j,k)%byc2-physic(i,j,k)%bzc2*physic(i,j,k)%vely)                      
                        SCRCH2 = 2.D0/(ZZNR(K)-ZZNL(K))*(FLUXC+CT(I,J,K)%fbyz)

	           ELSE

                        SCRCH21 = 2.D0/(ZZNR(K)-ZZNL(K))*(-(physic(i,j-1,k)%velz*physic(i,j-1,k)%byc2-                &
                                                            physic(i,j-1,k)%bzc2*physic(i,j-1,k)%vely)+CT(I,J-1,K)%fbyz)	
                        SCRCH22 = 2.D0/(ZZNR(K)-ZZNL(K))*(-(physic(i,j,k)%velz*physic(i,j,k)%byc2-                    &
                                                            physic(i,j,k)%bzc2*physic(i,j,k)%vely)+CT(I,J,K)%fbyz)
			SCRCH2 = 0.5D0*(SCRCH21+SCRCH22)			

		   ENDIF


                   
		   ! SCRCH4;
		   ! ...................................		

		   IF (CT(I,J-1,K)%RFLXZ .GT. 1.D-15) THEN

                        FLUXC  = (physic(i,j-1,k-1)%vely*physic(i,j-1,k-1)%bzc2-physic(i,j-1,k-1)%byc2*physic(i,j-1,k-1)%velz)
			SCRCH4 = 2.D0/(YZNR(J)-YZNL(J))*(CT(I,J,K-1)%fbzy-FLUXC)

		   ELSE IF (CT(I,J-1,K)%RFLXZ .LT. -1.D-15) THEN

                        FLUXC  = (physic(i,j-1,k)%vely*physic(i,j-1,k)%bzc2-physic(i,j-1,k)%byc2*physic(i,j-1,k)%velz)                
                        SCRCH4 = 2.D0/(YZNR(J)-YZNL(J))*(CT(I,J,K)%fbzy-FLUXC)

	           ELSE

                      SCRCH41 = 2.D0/(YZNR(J)-YZNL(J))*(CT(I,J,K-1)%fbzy-(physic(i,j-1,k-1)%vely*physic(i,j-1,k-1)%bzc2-  &
                                                                           physic(i,j-1,k-1)%byc2*physic(i,j-1,k-1)%velz))	
                      SCRCH42 = 2.D0/(YZNR(J)-YZNL(J))*(CT(I,J,K)%fbzy-(physic(i,j-1,k)%vely*physic(i,j-1,k)%bzc2-        &
                                                                        physic(i,j-1,k)%byc2*physic(i,j-1,k)%velz))
		      SCRCH4 = 0.5D0*(SCRCH41+SCRCH42)			

		   ENDIF

                   
                   ! SCRCH5;
                   ! ...................................


		   IF (CT(I,J-1,K)%RFLXZ .GT. 1.D-15) THEN

                        FLUXC  = (physic(i,j,k-1)%vely*physic(i,j,k-1)%bzc2-physic(i,j,k-1)%byc2*physic(i,j,k-1)%velz)
			SCRCH5 = 2.D0/(YZNR(J)-YZNL(J))*(FLUXC-CT(I,J,K-1)%fbzy)

		   ELSE IF (CT(I,J-1,K)%RFLXZ .LT. -1.D-15) THEN

                        FLUXC  = (physic(i,j,k)%vely*physic(i,j,k)%bzc2-physic(i,j,k)%byc2*physic(i,j,k)%velz)                      
                        SCRCH5 = 2.D0/(YZNR(J)-YZNL(J))*(FLUXC-CT(I,J,K)%fbzy)

	           ELSE

                        SCRCH51 = 2.D0/(YZNR(J)-YZNL(J))*((physic(i,j,k-1)%vely*physic(i,j,k-1)%bzc2-                &
                                                            physic(i,j,k-1)%byc2*physic(i,j,k-1)%velz)-CT(I,J,K-1)%fbzy)	
                        SCRCH52 = 2.D0/(YZNR(J)-YZNL(J))*((physic(i,j,k)%vely*physic(i,j,k)%bzc2-                    &
                                                            physic(i,j,k)%byc2*physic(i,j,k)%velz)-CT(I,J,K)%fbzy)
			SCRCH5 = 0.5D0*(SCRCH51+SCRCH52)			

		   ENDIF


                   SCRCH3 = SCRCH1 - SCRCH2		
		   SCRCH6 = SCRCH4 - SCRCH5




		   CT(i,j,k)%omega_x = CT(i,j,k)%omega_x + SCRCH3*(ZZNR(K)-ZZNL(K))/8.D0 + &
                                                           SCRCH6*(YZNR(J)-YZNL(J))/8.D0

                                      
                ENDDO
             ENDDO
          ENDDO
!$OMP END DO

!$OMP BARRIER
 
          !!! CT magnetic field evolution.
     
!$OMP DO SCHEDULE(static)
          DO k= 1, nz
             DO j= 1, ny
                DO i=1, nx

                   CT(i,j,k)%bx = ar*CT(i,j,k)%obx + br*CT(i,j,k)%bx                   &
                   + cr*dt*(CT(i,j+1,k)%omega_z-CT(i,j,k)%omega_z)/(yznr(j)-yznl(j))   &
                   + cr*dt*(CT(i,j,k)%omega_y-CT(i,j,k+1)%omega_y)/(zznr(k)-zznl(k))


                   CT(i,j,k)%by = ar*CT(i,j,k)%oby + br*CT(i,j,k)%by                   &
                   + cr*dt*(CT(i,j,k)%omega_z-CT(i+1,j,k)%omega_z)/(xznr(i)-xznl(i))   &
                   + cr*dt*(CT(i,j,k+1)%omega_x-CT(i,j,k)%omega_x)/(zznr(k)-zznl(k))


                   CT(i,j,k)%bz = ar*CT(i,j,k)%obz + br*CT(i,j,k)%bz +                 &
                   cr*dt*(CT(i+1,j,k)%omega_y-CT(i,j,k)%omega_y)/(xznr(i)-xznl(i))   + &
                   cr*dt*(CT(i,j,k)%omega_x-CT(i,j+1,k)%omega_x)/(yznr(j)-yznl(j))



                ENDDO
             ENDDO
          ENDDO
!$OMP END DO


!$OMP BARRIER


!$OMP DO SCHEDULE(static)
          DO k= 1, nz
             DO j= 1, ny

                CT(nx+1,j,k)%bx = ar*CT(nx+1,j,k)%obx + br*CT(nx+1,j,k)%bx                &
                + cr*dt*(CT(nx+1,j+1,k)%omega_z-CT(nx+1,j,k)%omega_z)/(yznr(j)-yznl(j))   &
                + cr*dt*(CT(nx+1,j,k)%omega_y-CT(nx+1,j,k+1)%omega_y)/(zznr(k)-zznl(k))

             ENDDO
          ENDDO
!$OMP END DO

!$OMP BARRIER

!$OMP DO SCHEDULE(static)
          DO k= 1, nz
             DO i= 1, nx

                CT(i,ny+1,k)%by = ar*CT(i,ny+1,k)%oby + br*CT(i,ny+1,k)%by                &
                + cr*dt*(CT(i,ny+1,k)%omega_z-CT(i+1,ny+1,k)%omega_z)/(xznr(i)-xznl(i))   &
                + cr*dt*(CT(i,ny+1,k+1)%omega_x-CT(i,ny+1,k)%omega_x)/(zznr(k)-zznl(k))


             ENDDO
          ENDDO
!$OMP END DO

!$OMP BARRIER

!$OMP DO SCHEDULE(static)
          DO j= 1, ny
             DO i= 1, nx

                CT(i,j,nz+1)%bz = ar*CT(i,j,nz+1)%obz + br*CT(i,j,nz+1)%bz +              &
                cr*dt*(CT(i+1,j,nz+1)%omega_y-CT(i,j,nz+1)%omega_y)/(xznr(i)-xznl(i))+    &
                cr*dt*(CT(i,j,nz+1)%omega_x-CT(i,j+1,nz+1)%omega_x)/(yznr(j)-yznl(j))

             ENDDO
          ENDDO
!$OMP END DO

!$OMP END PARALLEL

                
!$OMP PARALLEL DEFAULT(shared)        &
!$OMP PRIVATE(i,j,k)

!$OMP DO SCHEDULE(static)
        DO k=1, nz
           DO j=1,ny
              DO i=1,nx



                 physic(i,j,k)%bxct = (CT(i,j,k)%bx + CT(i+1,j,k)%bx)/2.D0
                    
                 IF (nsdim .EQ. 2) THEN

                    physic(i,j,k)%byct = (CT(i,j,k)%by + CT(i,j+1,k)%by)/2.D0

                 ELSE IF (nsdim .EQ. 3) THEN

                    physic(i,j,k)%byct = (CT(i,j,k)%by + CT(i,j+1,k)%by)/2.D0 
                    physic(i,j,k)%bzct = (CT(i,j,k)%bz + CT(i,j,k+1)%bz)/2.D0

                 ENDIF
       

              ENDDO
           ENDDO
        ENDDO
!$OMP END DO

!$OMP END PARALLEL


        IF (iccv .NE. 0) THEN
       
!$OMP PARALLEL DEFAULT (shared) &
!$OMP PRIVATE (i,j,k,b2,bct2)

!$OMP DO SCHEDULE (static)
           DO k=1, nz
              DO j=1,ny
                 DO i=1,nx

	           IF (iccv .EQ. 2) physic(i,j,k)%bxzn = physic(i,j,k)%bxc
		   IF (iccv .EQ. 2) physic(i,j,k)%byzn = physic(i,j,k)%byc
                   IF (iccv .EQ. 2) physic(i,j,k)%bzzn = physic(i,j,k)%bzc

                   b2 = physic(i,j,k)%bxc**2 + physic(i,j,k)%byc**2 + physic(i,j,k)%bzc**2
                   bct2 = physic(i,j,k)%bxct**2 + physic(i,j,k)%byct**2 + physic(i,j,k)%bzct**2

                   physic(i,j,k)%energy = physic(i,j,k)%energy - 0.5D0*(b2-bct2)

                ENDDO
              ENDDO
           ENDDO
!$OMP END DO

!$OMP END PARALLEL

        ENDIF


!$OMP PARALLEL DEFAULT (shared)  &
!$OMP PRIVATE (i,j,k)

!$OMP DO SCHEDULE (static)
        DO k=1,nz
           DO j=1,ny
              DO i=1, nx
 
	         physic(i,j,k)%bxc = physic(i,j,k)%bxct
	
	         IF (nsdim .EQ. 2) THEN

                    physic(i,j,k)%byc = physic(i,j,k)%byct

 		 ELSE IF (nsdim .EQ. 3) THEN

                    physic(i,j,k)%byc = physic(i,j,k)%byct
                    physic(i,j,k)%bzc = physic(i,j,k)%bzct

		 ENDIF

             ENDDO
           ENDDO
        ENDDO
!$OMP END DO

!$OMP END PARALLEL


        !DO k=1,nz-1
        !   DO j=1,ny-1
        !      DO i=1, nx-1


	!	DIVE = (CT(i+1,j,k)%bx-CT(i,j,k)%bx)/(XZNR(I)-XZNL(I))+&
        !              (CT(i,j+1,k)%by-CT(i,j,k)%by)/(YZNR(J)-YZNL(J))+&			
        !              (CT(i,j,k+1)%bz-CT(i,j,k)%bz)/(ZZNR(K)-ZZNL(K))
 !		IF (DIVE .GT. 1.D-6) print*, 'divergence', I, j, K
 !            ENDDO
 !          ENDDO
 !       ENDDO





        !----------------------------------------------------------------------------END

! =============================================================================== 
END SUBROUTINE SWEEPS_3D
! =============================================================================== 


