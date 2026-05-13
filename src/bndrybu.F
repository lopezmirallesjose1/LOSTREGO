!*********************************
!     3D-RMHD CODE: LOSTREGO
!*********************************
      
! NAME:
!       BNDRY

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
!       Provide boundary conditions
!       in the direction of the sweep.


! ===============================================================================
SUBROUTINE BNDRY (fznp, nznp, bndmin, bndmax, xyzswp, p, rho, rhoe, rhopr,      &
                  vel, velt, veltt, f, b, bt, btt, x, xl, xr, dx,               &
                  y, yl, yr, dy, z, zl, zr, dz, uwi, vwi, shift, Rorb, gamma_ad)
! ===============================================================================                                                        

                  
        USE parallel
        USE tipos
        
        implicit none

# include "types.h"


        !     Boundary conditions.
        !     bndm.. = 1 ===> reflecting boundary                 (min+max)
        !     bndm.. = 2 ===> flow out boundary (zero gradient)   (min+max)
        !     bndm.. = 3 ===> periodic boundary                   (min+max)
	!     bndm.. = 7 ===> homogenous wind
	!     bndm.. = 8 ===> pertubed jet + fixed                (min)


        !------------------------------------------------------- Input Variables

        integer(is):: fznp, nznp

        integer(is):: bndmin, bndmax, xyzswp

	real(rs)   :: y, yl, yr, dy, z, zl, zr, dz, uwi, vwi, shift, Rorb, gamma_ad

        !------------------------------------------------------ Output Variables

        ! size -4:mn5
        real(rs)   :: p(-4:*), rho(-4:*), rhoe(-4:*), rhopr(-4:*), vel(-4:*),     &
                      velt(-4:*), veltt(-4:*), f(-4:*), b(-4:*), bt(-4:*), btt(-4:*)

        real(rs)   :: x(-4:*), xl(-4:*), xr(-4:*), dx(-4:*)

        !------------------------------------------------------- Local Variables

        integer(is):: i

	real(rs) :: theta1, theta2

	real(rs) :: RHOB, VELB, PB, TT1, GAMM1, RBM, BMJ, WB, BPHI1, scrch1 
 
        !--------------------------------------------------------------------------- BEGIN

        
        IF ((xyzswp .EQ. 1 .AND. mpi_coords(1) .EQ. 0) .OR. &
            (xyzswp .EQ. 2 .AND. mpi_coords(2) .EQ. 0) .OR. &
            (xyzswp .EQ. 3 .AND. mpi_coords(3) .EQ. 0)) THEN

            
            !    =================
            !    LEFT BOUNDARY
            !    =================


            !! REFLECTING BOUNDARY (BND .. 1)
	     ! -------------------------------------------------

            IF (bndmin .EQ. 1) THEN

               DO i= -4, 0

                  rho(fznp-1+i)   = rho(fznp-i)
                  rhoe(fznp-1+i)  = rhoe(fznp-i)
                  rhopr(fznp-1+i) = rho(fznp-1+i) - rhoe(fznp-1+i)
                  vel(fznp-1+i)   = -vel(fznp-i)
                  velt(fznp-1+i)  = velt(fznp-i)
                  veltt(fznp-1+i) = veltt(fznp-i)
                  p(fznp-1+i)     = p(fznp-i)
                  f(fznp-1+i)     = f(fznp-i)
                  b(fznp-1+i)     = b(fznp-i)
                  bt(fznp-1+i)    = bt(fznp-i)
                  btt(fznp-1+i)   = btt(fznp-i)
                  dx(fznp-1+i)    = dx(fznp-i)

               ENDDO

               DO i=fznp-1, fznp-5, -1

                  xl(i) = xl(i+1) - dx(i)
                  xr(i) = xr(i+1) - dx(i+1)
                  x(i)  = 0.5D0*(xl(i) + xr(i))

               ENDDO


            !! FLOW-OUT BOUNDARY (BND .. 2)
             ! -------------------------------------------------


            ELSE IF (bndmin .EQ. 2) THEN

               DO i= fznp-5, fznp-1

                  rho(i)   =  rho(fznp)
                  rhoe(i)  =  rhoe(fznp)
                  rhopr(i) =  rho(i) - rhoe(i)
                  vel(i)   =  vel(fznp)
                  velt(i)  =  velt(fznp)
                  veltt(i) =  veltt(fznp)
                  p(i)     =  p(fznp)
                  f(i)     =  f(fznp)
                  b(i)     =  b(fznp)
                  bt(i)    =  bt(fznp)
                  btt(i)   =  btt(fznp) 
                  dx(i)    =  dx(fznp)

               ENDDO

               DO i=fznp-1, fznp-5, -1

                  xl(i) = xl(i+1) - dx(i)
                  xr(i) = xr(i+1) - dx(i+1)
                  x(i)  = 0.5D0*(xl(i) + xr(i))

               ENDDO


            !! PERIODIC BOUNDARY (BND .. 3)
             ! -------------------------------------------------


            ELSE IF (bndmin .EQ. 3) THEN

               DO i= -4, 0

                  rho   (fznp-1+i) = rho  (nznp+i)
                  rhoe  (fznp-1+i) = rhoe  (nznp+i)
                  rhopr (fznp-1+i) = rho  (fznp-1+i) - rhoe (fznp-1+i)
                  vel   (fznp-1+i) = vel  (nznp+i)
                  velt  (fznp-1+i) = velt (nznp+i)
                  veltt (fznp-1+i) = veltt (nznp+i)
                  p     (fznp-1+i) = p    (nznp+i)
                  f     (fznp-1+i) = f    (nznp+i)
                  b     (fznp-1+i) = b    (nznp+i)
                  bt    (fznp-1+i) = bt   (nznp+i)
                  btt   (fznp-1+i) = btt  (nznp+i)
                  dx    (fznp-1+i) = dx   (nznp+i)

               ENDDO

               DO i= fznp-1, fznp-5, -1

                  xl(i) = xl(i+1) - dx(i)
                  xr(i) = xr(i+1) - dx(i+1)
                  x(i)  = 0.5D0*(xl(i) + xr(i))

               ENDDO



            !! INFLOW SLAB JET 2D (BND .. 4)
             ! -------------------------------------------------


            ELSE IF (bndmin .EQ. 4) THEN


               RHOB = 0.1D0
               VELB = 0.998D0
               PB   = 0.01D0

                DO i = -4, 0
                   dx(i) = dx (fznp)
                ENDDO


                DO i=fznp-1, fznp-5, -1

                   xl(i) = xl(i+1) - dx(i)
                   xr(i) = xr(i+1) - dx(i+1)
                   x(i)  = 0.5D0*(xl(i) + xr(i))

                ENDDO


               DO i= -4, 0

                  tt1 = ABS(Y)

                  IF (tt1 .LE. 1.D0) THEN

                     rho(i)  = rhob
                     rhoe(i) = rhob
                     rhopr(i)= rho(i)-rhoe(i)
                     vel(i)  = velb
                     velt(i) = 0.D0
                     veltt(i)= 0.D0
                     p(i)    = pb
                     f(i)    = 1.D0

                     b(i)     = 1.D0
                     bt(i)    = 0.D0
                     btt(i)   = 0.D0


                  ELSE


                     rho(fznp-1+i)   = rho(fznp-i)
                     rhoe(fznp-1+i)  = rhoe(fznp-i)
                     rhopr(fznp-1+i) = rho(fznp-1+i) - rhoe(fznp-1+i)
                     vel(fznp-1+i)   = -vel(fznp-i)
                     velt(fznp-1+i)  = velt(fznp-i)
                     veltt(fznp-1+i) = veltt(fznp-i)
                     p(fznp-1+i)     = p(fznp-i)
                     f(fznp-1+i)     = f(fznp-i)
                     b(fznp-1+i)     = b(fznp-i)
                     bt(fznp-1+i)    = bt(fznp-i)
                     btt(fznp-1+i)   = btt(fznp-i)

	          ENDIF

               ENDDO



            !! INFLOW TEST + REFLECTING (BND .. 7)
             ! -------------------------------------------------

             ELSE IF (bndmin .EQ. 7) THEN


                DO i = -4, 0
                   dx(i) = dx (fznp)
                ENDDO


                DO i=fznp-1, fznp-5, -1

                   xl(i) = xl(i+1) - dx(i)
                   xr(i) = xr(i+1) - dx(i+1)
                   x(i)  = 0.5D0*(xl(i) + xr(i))

                ENDDO

                DO i= -4, 0


                        rho(i)  = rho(fznp)
                        rhoe(i) = rho(fznp)
                        rhopr(i)=rho(i)-rhoe(i)

                        ! INCLINATION
                        theta1 = DATAN(Y/	&
                                 DSQRT((x(i)+shift)**2+(Rorb-Z)**2))


                        ! AZIMUTHAL ANGLE.
                        theta2 = DATAN((x(i)+shift)/(Rorb-Z))

                        vel(i)  =  vwi*DCOS(theta1)*DSIN(theta2)
                        velt(i) =  vwi*DSIN(theta1)
                        veltt(i)= -vwi*DCOS(theta1)*DCOS(theta2)

                        !vel(i) = 0.D0
                        !velt(i) = 0.D0
                        !veltt(i) = 0.D0


                        p(i)    = (gamma_ad-1.D0)*uwi*1.D0

                        f(i)    = 0.D0


		 ENDDO



            !! INFLOW PERTURBED JET + FIXED (BND .. 8)
             ! -------------------------------------------------


            ELSE IF (bndmin .EQ. 8) THEN

               RHOB = 8.8D0
               VELB = 0.55D0
                PB  = 0.00263D0
               RBM  = 0.37D0
                BMJ = 0.106812D0
                WB  = 1.D0/DSQRT(1.D0-VELB**2)

                DO i = -4, 0
                   dx(i) = dx (fznp)
                ENDDO


                DO i=fznp-1, fznp-5, -1

                   xl(i) = xl(i+1) - dx(i)
                   xr(i) = xr(i+1) - dx(i+1)
                   x(i)  = 0.5D0*(xl(i) + xr(i))

                ENDDO


               DO i= -4, 0

                  b(i)     =  b(fznp)
                  bt(i)    =  bt(fznp)
                  btt(i)   =  btt(fznp)

                  tt1 = DSQRT(Z*Z + Y*Y)

                  IF (tt1 .LE. 1.D0) THEN

		     BPHI1 = (2.D0*BMJ*(1.D0/RBM))/(1.D0+(1.D0/RBM)**2)

             	     rho(i)  = (rhob-1.D0)/DCOSH(TT1**4)+1.D0
             	     rhoe(i) = (rhob-1.D0)/DCOSH(TT1**4)+1.D0
		     rhopr(i)= rho(i)-rhoe(i)
            	     vel(i)  = velb/DCOSH(tt1**4)
            	     velt(i) = 0.D0
            	     veltt(i)= 0.D0
            	     p(i)    = 2.D0*(BMJ/(WB*(1.D0+(TT1/RBM)**2)))**2+PB-BPHI1**2/(2.D0*WB**2)*(1.D0+RBM**2)
            	     f(i)    = 1.D0

		
		  ELSE

		     IF (tt1 .LE. 1.5D0) THEN

                        rho(i)   =  1.D0
                  	rhoe(i)  =  1.D0
                  	rhopr(i) =  rho(i) - rhoe(i)

		  	!IF (velb/DCOSH(tt1**4) .GT. vwi) THEN
		  			
                  	   vel(i)   = velb/DCOSH(tt1**4)
		  	   velt(i)  = 0.D0
		  	   veltt(i) = 0.D0

	                !ELSE

		  	!   vel(i)   =  vel(fznp)
                   	!   velt(i)  =  velt(fznp)
                  	!   veltt(i) =  veltt(fznp)

		        !ENDIF

                  	p(i)     =  p(fznp)
                  	f(i)     =  0.D0


		     ELSE IF (tt1 .GT. 1.5D0 .AND. p(fznp) .GT. 1.D-8) THEN


                        rho(fznp-1+i)   = rho(fznp-i)
                        rhoe(fznp-1+i)  = rhoe(fznp-i)
                        rhopr(fznp-1+i) = rho(fznp-1+i) - rhoe(fznp-1+i)
                        vel(fznp-1+i)   = -vel(fznp-i)
                        velt(fznp-1+i)  = velt(fznp-i)
                        veltt(fznp-1+i) = veltt(fznp-i)
                        p(fznp-1+i)     = p(fznp-i)
                        f(fznp-1+i)     = 0.D0

		     ELSE

             	        rho(i)  = rho(fznp)
		        rhoe(i) = rho(fznp)
		        rhopr(i)=rho(i)-rhoe(i)
			
		        ! INCLINATION
                        theta1 = DATAN(Y/       &
                                 DSQRT((x(i)+shift)**2+(Rorb-Z)**2))


                        ! AZIMUTHAL ANGLE.
                        theta2 = DATAN((x(i)+shift)/(Rorb-Z))

               	        vel(i)  =  vwi*DCOS(theta1)*DSIN(theta2)   
            	        velt(i) =  vwi*DSIN(theta1)
            	        veltt(i)= -vwi*DCOS(theta1)*DCOS(theta2)
	
			!vel(i) = 0.D0
			!velt(i) = 0.D0
			!veltt(i) = 0.D0


	                p(i)    = (gamma_ad-1.D0)*uwi*1.D0

               	        f(i)    = 0.D0

		     ENDIF

      	    	  ENDIF 

    	       ENDDO  


            ELSE 

               IF (xyzswp .EQ. 1) THEN
                  call ERROR ('BNDRY','Unexpected x-min boundary', __FILE__, __LINE__)
               ELSE IF (xyzswp .EQ. 2) THEN
                  call ERROR ('BNDRY','Unexpected y-min boundary', __FILE__, __LINE__)
               ELSE
                  call ERROR ('BNDRY','Unexpected z-min boundary', __FILE__, __LINE__)
               ENDIF

            ENDIF

        ENDIF


        IF ((xyzswp .EQ. 1 .AND. mpi_coords(1) .EQ. mpi_dims(1)-1) .OR.    &
            (xyzswp .EQ. 2 .AND. mpi_coords(2) .EQ. mpi_dims(2)-1) .OR.    &
            (xyzswp .EQ. 3 .AND. mpi_coords(3) .EQ. mpi_dims(3)-1)) THEN


            !    =================
            !    RIGHT BOUNDARY
            !    =================

            !! REFLECTING BOUNDARY (BND .. 1)
             ! -------------------------------------------------


            IF (bndmax .EQ. 1) THEN

               DO i= 1, 5

                  rho(nznp+i)   = rho(nznp+1-i)
                  rhoe(nznp+i)  = rhoe(nznp+1-i)
                  rhopr(nznp+i) = rho(nznp+i) - rhoe(nznp+i)
                  vel(nznp+i)   = -vel(nznp+1-i)
                  velt(nznp+i)  = velt(nznp+1-i)
                  veltt(nznp+i) = veltt(nznp+1-i)
                  p(nznp+i)     = p(nznp+1-i)
                  f(nznp+i)     = f(nznp+1-i)
                  b(nznp+i)     = -b(nznp+1-i)
                  bt(nznp+i)    = bt(nznp+1-i)
                  btt(nznp+i)   = btt(nznp+1-i)
                  dx(nznp+i)    = dx(nznp+1-i)

               ENDDO

               DO i=nznp+1, nznp+5

                  xl(i) = xl(i-1) + dx(i-1)
                  xr(i) = xr(i-1) + dx(i)
                  x(i)  = 0.5D0*(xl(i) + xr(i))

               ENDDO


            !! FLOW-OUT (ZERO GRADIENT) BOUNDARY (BND .. 2)
             ! -------------------------------------------------


            ELSE IF (bndmax .EQ. 2) THEN

               DO i= nznp+1, nznp+5

                  rho(i)   =  rho(nznp)
                  rhoe(i)  =  rhoe(nznp)
                  rhopr(i) =  rho(i) - rhoe(i)
                  vel(i)   =  vel(nznp)
                  velt(i)  =  velt(nznp)
                  veltt(i) =  veltt(nznp)
                  p(i)     =  p(nznp)
                  f(i)     =  f(nznp)
                  b(i)     =  b(nznp)
                  bt(i)    =  bt(nznp)
                  btt(i)   =  btt(nznp) 
                  dx(i)    =  dx(nznp)

               ENDDO

               DO i=nznp+1, nznp+5

                  xl(i) = xl(i-1) + dx(i-1)
                  xr(i) = xr(i-1) + dx(i)
                  x(i)  = 0.5D0*(xl(i) + xr(i))

               ENDDO


            !! PERIODIC BOUNDARY (BND .. 3)
             ! -------------------------------------------------


            ELSE IF (bndmax .EQ. 3) THEN

               DO i= 1, 5

                  rho   (nznp+i) = rho  (i)
                  rhoe  (nznp+i) = rhoe (i)
                  rhopr (nznp+i) = rho  (nznp+i) - rhoe (nznp+i)
                  vel   (nznp+i) = vel  (i)
                  velt  (nznp+i) = velt (i)
                  veltt (nznp+i) = veltt (i)
                  p     (nznp+i) = p    (i)
                  f     (nznp+i) = f    (i)
                  b     (nznp+i) = b    (i)
                  bt    (nznp+i) = bt   (i)
                  btt   (nznp+i) = btt  (i)
                  dx    (nznp+i) = dx   (i)

               ENDDO

               DO i= nznp+1, nznp+5

                  xl(i) = xl(i-1) + dx(i-1)
                  xr(i) = xr(i-1) + dx(i)
                  x(i)  = 0.5D0*(xl(i) + xr(i))

               ENDDO



            !! WIND INFLOW BOUNDARY (BND .. 9)
             ! -------------------------------------------------	    

       	     ELSE IF (bndmax .EQ. 9) THEN

	        DO i = nznp+1, nznp+5

                   xl(i) = xl(i-1) + dx(nznp)
       		   xr(i) = xr(i-1) + dx(nznp)
       		   x(i)  = 0.5D0*(xl(i) + xr(i))

                   ! INCLINATION
                   theta1 = DATAN(z/       &
                            DSQRT((y+shift)**2+(Rorb-x(i))**2))

                   ! AZIMUTHAL ANGLE.
                   theta2 = DATAN((y+shift)/(Rorb-x(i)))

                   vel(i)   = -vwi*DCOS(theta1)*DCOS(theta2)
                   velt(i)  =  vwi*DCOS(theta1)*DSIN(theta2)
                   veltt(i) =  vwi*DSIN(theta1)


       		   rho(i)  = rho(nznp)
		   rhoe(i) = rhoe(nznp)
		   rhopr(i) =  rho(i) - rhoe(i)


       		   p(i)    = (gamma_ad - 1.D0)*uwi*1.D0

      		   f(i)    =  0.D0
         	   dx(i)   = dx(nznp)

		   b (i)   = 0.D0
	           bt (i)  = 0.D0
		   btt(i)  = 0.D0

               ENDDO

            ELSE 

               IF (xyzswp .EQ. 1) THEN
                  call ERROR ('BNDRY','Unexpected x-max boundary', __FILE__, __LINE__ )
               ELSE IF (xyzswp .EQ. 2) THEN
                  call ERROR ('BNDR+Y','Unexpected y-max boundary', __FILE__, __LINE__ )
               ELSE
                  call ERROR ('BNDRY','Unexpected z-max boundary', __FILE__, __LINE__ )
               ENDIF

            ENDIF

	
        ENDIF
 
	RETURN

        !----------------------------------------------------------------------------- END

! ===============================================================================        
END SUBROUTINE BNDRY
! ===============================================================================
