!*********************************
!     3D-RMHD CODE: LOSTREGO
!*********************************

! NAME:
!       PLM (slope-limiter method)

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
!       Linear interpolation for cell reconstruction.

! ACADEMIC REFERENCE:
!       "An HLLC Riemann solver for relativistic flows-II, Magnetohydrodynamics."
!        A. Mignone & G. Bodo. 2006. Mon. Not. Astron. Soc. 368, 1040-1054
!        See also: "Relativistic Hydrodynamics". L.Rezzolla & O. Zanotti (Chap. 9.3)


	
! ===============================================================================
SUBROUTINE PLM (LIMITER, ISHOCK, N, X, XL, XR, UM, U, UP)
! ===============================================================================
        
	USE workarrays
	implicit none
#include "types.h"

        !-------------------------------------------------------Input Variables

        character(4) :: LIMITER

        !size -4:mn5
	integer(is) :: N, ISHOCK(-4:*)

	real(rs)   :: U(-4:*), X(-4:*), XL(-4:*), XR(-4:*)

        !------------------------------------------------------ Output Variables

        !size -4:mn5
	real(rs)   :: UM(-4:*), UP(-4:*)     
	  
        !------------------------------------------------------ Local Variables

	integer(is) ::  I   

	real(rs) :: DELTAX1, DELTAX2, SS1, SS2, SS3,  &
	            sign, signSS1, signSS2, SLOPE

	real(rs) :: ALPHA, EPSILN 
		    

        !------------------------------------------------------------------------- BEGIN



        DO I= 0, N+1


	   ! Delta variables
	   DELTAX1  = X(I) - X(I-1)
	   DELTAX2  = X(I+1) - X(I)

	   ! Slope variables
	   SS1  = (U(I) - U(I-1))/DELTAX1
	   SS2  = (U(I+1) - U(I))/DELTAX2
	   SS3  = (U(I+1) - U(I-1))/(DELTAX1 + DELTAX2)

           ! For strong shocks (ISHOCK=1), the MC/VANL prescriptions can produce numerical oscillations.
           ! A weak form of flattening is introduced replacing the limiter by the MINMOD slope.	   
           ! The limiter MCMM/VLMM behaves like MINMOD if a shock is detected, and like MC/VANL otherwise.
	   
	   IF (LIMITER .EQ. 'MMOD' .OR. (LIMITER .EQ. 'MCMM' .AND. ISHOCK(I) .EQ. 1)     &
                                   .OR. (LIMITER .EQ. 'VLMM' .AND. ISHOCK(I) .EQ. 1)) THEN


	      IF ((SS1*SS2) .GT. 0.D0) THEN
		 ! We are imposing the sign of SS1 to be the same as SS2.
		 SLOPE = DMIN1(ABS(SS1),ABS(SS2))*(SS1/ABS(SS1))
	      ELSE
		 SLOPE = 0.D0
	      ENDIF


	   ELSE IF (LIMITER .EQ. 'MC00' .OR. (LIMITER .EQ. 'MCMM' .AND. ISHOCK(I) .EQ. 0)) THEN

	      IF (SS1 .EQ. 0.D0) THEN
		 signSS1 = 0.D0
	      ELSE
		 signSS1 = SS1/ABS(SS1)
	      END IF

	      IF (SS2 .EQ. 0.D0) THEN
		 signSS2 = 0.D0
	      ELSE
		 signSS2 = SS2/ABS(SS2)
	      ENDIF

	      sign = 0.5D0*(signSS1 + signSS2)
	      SLOPE = sign*DMIN1(2.D0*ABS(SS1),2.D0*ABS(SS2),ABS(SS3))


           ELSE IF (LIMITER.EQ.'VANL'.OR.(LIMITER.EQ.'VLMM'.AND.ISHOCK(I).EQ.0)) THEN

              IF (SS3.NE.0.D0) THEN
                 SLOPE = MAX(0.D0,SS1*SS2)/SS3
              ELSE
                 SLOPE = 0.D0
              ENDIF



           ELSE IF (LIMITER .EQ. 'GODU') THEN	

	       SLOPE = 0.D0

	   ENDIF


           ! Linear interpolation to obtain edge variables.
	   UM(I) = U(I) + SLOPE*(XL(I) - X(I)) 
	   UP(I) = U(I) + SLOPE*(XR(I) - X(I))    

	ENDDO
	
	RETURN

        !------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE PLM
! ===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


! NAME:
!       MP 

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
!       MP Reconstruction (MP3, MP5, MP7, MP9) in planar coordinates
!       and homogeneous grid.

! ACADEMIC REFERENCE:
!       Suresh & Huynh, JCP, vol. 136, 83 (1997)



! ===============================================================================
SUBROUTINE MPREC (tid, LIMITER, ISHOCK, N, X, XL, XR, UM, U, UP)
! ===============================================================================


        USE workarrays
        implicit none
#include "types.h"


        !-------------------------------------------------------Input Variables

        character(4) :: LIMITER

        !size -4:mn5
        integer(is) :: N, ISHOCK(-4:*), tid

        real(rs)   :: U(-4:*), X(-4:*), XL(-4:*), XR(-4:*)

        !------------------------------------------------------ Output Variables

        !size -4:mn5
        real(rs)   :: UM(-4:*), UP(-4:*)

        !------------------------------------------------------ Local Variables

	integer(is) :: I

	real(rs) :: SS, SSP, SLOPEP, SLOPEM, ALPHA, EPSILN,  &
	 	    UUL, UAV, UMD, ULC, UMN, UMX, SLOPE

	real(rs), pointer :: UMPP(:), UMPM(:), D(:), DM4P(:), DM4M(:) 

        !------------------------------------------------------------------------- BEGIN


        ! -4:mn6 vectors

        UMPP     => wa(tid)%wa6MP100
        UMPM     => wa(tid)%wa6MP101
        D        => wa(tid)%wa6MP102
        DM4P     => wa(tid)%wa6MP103
        DM4M     => wa(tid)%wa6MP104



        ALPHA  = 4.D0
        EPSILN = 1.D-10


        !     ORIGINAL INTERFACE VALUES ACCORDING TO S&H97, EQ. (2.1)
        !     -------------------------------------------------------
      
      	IF (LIMITER.EQ.'MP3_') THEN
            DO I = 0, N+1
                UP(I) = (- 1.D0*U(I-1) + 5.D0*U(I) + 2.D0*U(I+1))/6.D0
                UM(I) = (- 1.D0*U(I+1) + 5.D0*U(I) + 2.D0*U(I-1))/6.D0
            ENDDO
        ELSE IF (LIMITER.EQ.'MP5_') THEN
            DO I = 0, N+1
                UP(I) = (2.D0*U(I-2) - 13.D0*U(I-1) + 47.D0*U(I) + 27.D0*U(I+1) - 3.D0*U(I+2))/60.D0
                UM(I) = (2.D0*U(I+2) - 13.D0*U(I+1) + 47.D0*U(I) + 27.D0*U(I-1) - 3.D0*U(I-2))/60.D0
            ENDDO
        ELSE IF (LIMITER.EQ.'MP7_') THEN
            DO I = 0, N+1
                UP(I) = (-3.D0*U(I-3) + 25.D0*U(I-2) - 101.D0*U(I-1) +                 &
                        319.D0*U(I) + 214.D0*U(I+1) - 38.D0*U(I+2) + 4.D0*U(I+3))/420.D0
                UM(I) = (-3.D0*U(I+3) + 25.D0*U(I+2) - 101.D0*U(I+1) +                 & 
                        319.D0*U(I) + 214.D0*U(I-1) - 38.D0*U(I-2) + 4.D0*U(I-3))/420.D0
            ENDDO
        ELSE IF (LIMITER.EQ.'MP9_') THEN
            DO I = 0, N+1
                UP(I) = (4.D0*U(I-4) - 41.D0*U(I-3) + 199.D0*U(I-2) -   &
                        641.D0*U(I-1) + 1879.D0*U(I) + 1375.D0*U(I+1) - 305.D0*U(I+2) + 55.D0*U(I+3) -5.D0*U(I+4))/2520.D0
                UM(I) = (4.D0*U(I+4) - 41.D0*U(I+3) + 199.D0*U(I+2) -   &
                        641.D0*U(I+1) + 1879.D0*U(I) + 1375.D0*U(I-1) - 305.D0*U(I-2) + 55.D0*U(I-3) -5.D0*U(I-4))/2520.D0
            ENDDO
        ENDIF


        !     MP VALUES ACCORDING TO S&H97, EQ. (2.12)
        !     ----------------------------------------

      	DO I = 0, N+1
	   SS  = U(I)   - U(I-1)
	   SSP = U(I+1) - U(I)

           IF ((SS*SSP) .GT. 0.D0) THEN
	      SLOPEP = DMIN1(ABS(ALPHA*SS),ABS(SSP))*SS/ABS(SS)
	      SLOPEM = DMIN1(ABS(SS),ABS(ALPHA*SSP))*SS/ABS(SS)
	   ELSE
	      SLOPEP = 0.D0
              SLOPEM = 0.D0
	   ENDIF

           UMPP(I) = U(I) + SLOPEP
           UMPM(I) = U(I) - SLOPEM
        ENDDO


        !     MODIFICATION OF THE ORIGINAL INTERFACE VALUES ACCORDING TO 
        !     THE ALGORITHM DESCRIBED AT THE END OF SECT. 2.3, S&H97
        !     ----------------------------------------------------------


        DO I = -1, N+2
            D(I) = U(I-1) + U(I+1) - 2.D0*U(I)
        ENDDO

        DO I = -1, N+1
            IF (D(I)*(4.D0*D(I) - D(I+1))  .GT. 0.D0 .AND.      &
                D(I)*(4.D0*D(I+1) - D(I))  .GT. 0.D0 .AND.      &
                D(I)*D(I+1)                .GT. 0.D0) THEN

                DM4P(I) = DMIN1(ABS(4.D0*D(I) - D(I+1)),ABS(4.D0*D(I+1) - D(I)),  &
                                ABS(D(I)),ABS(D(I+1)))*D(I)/ABS(D(I))
	 
            ELSE
                DM4P(I) = 0.D0
	    ENDIF

            DM4M(I+1) = DM4P(I)
        ENDDO

        DO I = 0, N+1
           IF ((UP(I) - U(I))*(UP(I) - UMPP(I)) .GT. EPSILN) THEN
              UUL    = U(I) + ALPHA*(U(I) - U(I-1))
              UAV    = (U(I) + U(I+1))/2.D0
              UMD    = UAV - DM4P(I)/2.D0
              ULC    = U(I) + (U(I) - U(I-1))/2.D0 + 4.D0*DM4P(I-1)/3.D0
              UMN    = DMAX1(DMIN1(U(I),U(I+1),UMD), DMIN1(U(I),UUL,ULC))
              UMX    = DMIN1(DMAX1(U(I),U(I+1),UMD), DMAX1(U(I),UUL,ULC))

              IF ((UMN - UP(I))*(UMX - UP(I)) .GT. 0.D0) THEN
	         SLOPEP = DMIN1(ABS(UMN - UP(I)),ABS(UMX - UP(I)))*(UMN - UP(I))/ABS(UMN - UP(I))
	      ELSE
	         SLOPEP = 0.D0
	      ENDIF

              UP(I) = UP(I) + SLOPEP
           ENDIF
        ENDDO

        DO I = 0, N+1
           IF ((UM(I) - U(I))*(UM(I) - UMPM(I)) .GT. EPSILN) THEN
              UUL    = U(I) + ALPHA*(U(I) - U(I+1))
              UAV    = (U(I) + U(I-1))/2.D0
              UMD    = UAV - DM4M(I)/2.D0
              ULC    = U(I) + (U(I) - U(I+1))/2.D0 + 4.D0*DM4M(I+1)/3.D0
              UMN    = DMAX1(DMIN1(U(I),U(I-1),UMD), DMIN1(U(I),UUL,ULC))
              UMX    = DMIN1(DMAX1(U(I),U(I-1),UMD), DMAX1(U(I),UUL,ULC))

              IF ((UMN - UM(I))*(UMX - UM(I)) .GT. 0.D0) THEN
	         SLOPEM = DMIN1(ABS(UMN - UM(I)),ABS(UMX - UM(I)))*(UMN - UM(I))/ABS(UMN - UM(I))
	      ELSE
	         SLOPEM = 0.D0
	      ENDIF

              UM(I) = UM(I) + SLOPEM
           ENDIF
        ENDDO


        DO I = 0, N+1
           IF (ISHOCK(I).EQ.1) THEN
	       SS  = U(I)   - U(I-1)
	       SSP = U(I+1) - U(I)

	       IF ((SS*SSP) .GT. 0.D0) THEN
	           SLOPE = DMIN1(ABS(SS),ABS(SSP))*SS/ABS(SS)
	       ELSE
	           SLOPE = 0.D0
	       ENDIF
               UM(I) = U(I) - SLOPE/2.D0
               UP(I) = U(I) + SLOPE/2.D0
           ENDIF
        ENDDO

        RETURN

! ===============================================================================
END SUBROUTINE MPREC
! ===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


! NAME:
!       DETECT

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
!       Shock detection to switch the limiter.

! ACADEMIC REFERENCE:
!       "An HLLC Riemann solver for relativistic flows-II, Magnetohydrodynamics 
!        (Appendix A)." A. Mignone & G. Bodo. 2006. Mon. Not. Astron. Soc. 368, 1040-1054.
!        See also: "Relativistic Hydrodynamics". L.Rezzolla & O. Zanotti (Chap. 9.3)


! ===============================================================================
SUBROUTINE DETECT (N, tid, VEL, BETA, P, ISHOCK)
! ===============================================================================

	USE workarrays
	implicit none
#include "types.h"

        !------------------------------------------------------- Input Variables

        ! size -4:mn5
        integer(is) :: N, tid

        real(rs)   :: VEL(-4:*), P(-4:*), BETA(-4:*)	 

        !------------------------------------------------------ Output Variables

        ! size -4:mn5
        integer(is)   :: ISHOCK(-4:*)   	  

        !------------------------------------------------------- Local Variables
	
        integer(is) ::  I

        real(rs)  :: EPSI, EPSIB

        real(rs), pointer :: xi(:) 

        !------------------------------------------------------------------------- BEGIN
 
	xi => wa(tid)%wa5plm100

	! The threshold (EPSI) is established attending to particular needs.
	! An adittional constrained (epsiB) for the magnetic pressure 
	! is provided according to Mignone & Bodo (2006). 
	EPSI = 1.D0
	EPSIB = 1.D0

	! First condition: pressure
	DO I = 0, N+2

	   ! Initialization
	   ISHOCK(I) = 0     
	   xi(I) = 0.D0

	   IF (BETA(I+1) .NE. 0.D0 .AND. BETA(I-1) .NE. 0.D0) THEN

	   IF (DABS(P(I+1) - P(I-1))/MIN(P(I+1),P(I-1)) .LE. EPSI .AND.       &
	   DABS(P(I+1)*BETA(I+1) - P(I-1)*BETA(I-1))/MIN(P(I+1)*BETA(I+1),    & 
           P(I-1)*BETA(I-1)) .LE. EPSIB) xi(I) = 1.D0

	   ELSE

           IF (DABS(P(I+1) - P(I-1))/MIN(P(I+1),P(I-1)) .LE. EPSI) xi(I) = 1.D0

	   ENDIF

	ENDDO

	! Second condition: velocity gradient.
	DO I = 0, N+1
	   IF (VEL(I-1) .GT. VEL(I+1) .AND. xi(I) .EQ. 0.D0) ISHOCK(I) = 1
	ENDDO

	RETURN
        
        !------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE DETECT
! ===============================================================================
