!*********************************
!     3D-RMHD CODE: LOSTREGO
!*********************************


! NAME:
!       GRID 

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
!       Compute the grid coordinates.

! ACADEMIC REFERENCE:
!       See references in SUBROUTINE INIT.


! ===============================================================================       
SUBROUTINE GRID (basenm, nx, ny, nz, nsdim, gridlx, gridly, gridlz,                   &
                 xCPU, yCPU, zCPU, bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,    &
                 delx, dely, delz, xznl, xzn, xznr, yznl, yzn, yznr,                  &
                 zznl, zzn, zznr)
! =============================================================================== 


        USE tipos
        USE parallel
        USE memoria
        implicit none
# include "types.h"
      
        !------------------------------------------------------- Input Variables

        character(4) :: basenm

        integer(is)  :: nx, ny ,nz, nsdim

        real(rs)       :: gridlx, gridly, gridlz
        
        integer(is)    :: xCPU, yCPU, zCPU
      
        !------------------------------------------------------ Output Variables
             
        integer(is):: bndmnx, bndmxx, bndmny, bndmxy, bndmxz, bndmnz

        real(rs)   :: delx, dely, delz
        
        real(rs)   :: xznl(*), xzn(*), xznr(*), &
                      yznl(*), yzn(*), yznr(*), &
                      zznl(*), zzn(*), zznr(*)

        !------------------------------------------------------- Local Variables

        character(40) :: message1
        
        integer(is) :: I, ii, N

	real(rs) :: R

        !-------------------------------------------------------------------------- BEGIN
        
        !     =====================
        !     Type of coordinate
        !     =====================

        !     igeom = 0 ===> planar (Cartesian) geometry
 
        !     =====================
        !     Boundary conditions
        !     =====================

        !     bndm.. = 1 ===> reflecting boundary (min+max)
        !     bndm.. = 2 ===> flow out boundary   (min+max)
        !     bndm.. = 3 ===> periodic boundary   (min+max)

        !     ----------------------
        !     1D Cartesian RMHD Test
        !     ----------------------

        !     BALSARA'S TEST :: 'BA'
        !     (BALSARA, 2001, ApJ, vol. 132, 83)

        !     GENERIC ALFVEN TEST :: 'GR'
        !     (GIACOMAZZO & REZZOLLA, 2006, JFM, vol. 562, 223)

        !     KOMISSAROV'S SHOCK TUBE :: 'KO'
        !     (KOMISSAROV, 1999, MNRAS, vol. 303, 343)

        IF (BASENM(1:2) .EQ. 'BA' .OR. BASENM(1:2) .EQ. 'GR' .OR.    &
            BASENM(1:2) .EQ. 'KO' .OR. BASENM .EQ. 'CONT' .OR.       &
            BASENM .EQ. 'DZBU') THEN

           BNDMNX = 2
           BNDMXX = 2

        !     LARGE AMPLITUDE CP ALFVEN WAVE TEST :: 'DZAW'
        !     (DEL ZANNA ET AL. 2007, A&A, vol. 473, 11)

        ELSE IF (BASENM .EQ. 'DZAW') THEN

           BNDMNX = 3
           BNDMXX = 3
           
        !     ----------------------
        !     2D Cartesian RMHD Test
        !     ----------------------

        !     KOMISSAROV'S CYLINDRICAL MAGNETIZED BLAST WAVE :: 'CMBW'
        !     (KOMISSAROV, 1999, MNRAS, vol. 303, 343)

        !     RELATIVISTIC ROTOR :: 'RROT'
        !     (DEL ZANNA ET AL. 2003, A&A, vol. 400, 397)

        ELSE IF (BASENM .EQ. 'CMBW' .OR. BASENM .EQ. 'RROT' .OR.    &
                 BASENM .EQ. 'SYBW') THEN

           BNDMNX = 2
           BNDMXX = 2

           BNDMNY = 2
           BNDMXY = 2

        ELSE IF (BASENM .EQ. 'IMPL') THEN

           BNDMNX = 2
           BNDMXX = 2

           BNDMNY = 2
           BNDMXY = 2           


        ELSE IF (BASENM .EQ. 'SLKO') THEN

           BNDMNX = 4
           BNDMXX = 2

           BNDMNY = 2
           BNDMXY = 2


        ELSE IF (BASENM .EQ. 'KH2D' .OR. BASENM .EQ. 'ORTA') THEN

           BNDMNX = 3
           BNDMXX = 3

           BNDMNY = 3
           BNDMXY = 3


        !     ----------------------
        !     3D Cartesian RMHD Test
        !     ----------------------

        !     KOMISSAROV'S SPHERICAL MAGNETIZED BLAST WAVE :: 'SMBW'
        !     (KOMISSAROV, 1999, MNRAS, vol. 303, 343)

        ELSE IF (BASENM .EQ. 'SMBW' .OR. BASENM .EQ. '3ROT' .OR. BASENM .EQ. 'SC3D') THEN

           BNDMNX = 2
           BNDMXX = 2

           BNDMNY = 2
           BNDMXY = 2

	   BNDMNZ = 2
	   BNDMXZ = 2


        ELSE IF (BASENM .EQ. 'KH3D') THEN

           BNDMNX = 2
           BNDMXX = 2

           BNDMNY = 2
           BNDMXY = 2

           BNDMNZ = 2
           BNDMXZ = 2


        ELSE IF (BASENM .EQ. 'LOOP') THEN

           BNDMNX = 3
           BNDMXX = 3

           BNDMNY = 3
           BNDMXY = 3

           BNDMNZ = 3
           BNDMXZ = 3



        ELSE IF (BASENM .EQ. 'MI10') THEN

           BNDMNX = 2
           BNDMXX = 2

           BNDMNY = 8
           BNDMXY = 2

           BNDMNZ = 2
           BNDMXZ = 2


        ELSE IF (BASENM(1:2) .EQ. 'S4') THEN

           BNDMNX = 2
           BNDMXX = 2

           BNDMNY = 2
           BNDMXY = 2

           BNDMNZ = 2
           BNDMXZ = 2


	!     ====================================
        !     MICROQUASAR JET WITH HOMOGENOUS WIND 
	!     ====================================

        ELSE IF ((BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'H') .OR.  &
                 (BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'C')) THEN


           BNDMNX = 2
           BNDMXX = 9

           BNDMNY = 8
           BNDMXY = 2

           BNDMNZ = 2
           BNDMXZ = 2

         
        ELSE

           call ERROR ('GRIDBU','Check boundary conditions.' , __FILE__, __LINE__)

        ENDIF


        ! Check periodic boundary condition.

        IF ((BNDMNX.EQ.3.AND.BNDMXX.NE.3) .OR.   &
            (BNDMNY.EQ.3.AND.BNDMXY.NE.3) .OR.   &
            (BNDMNZ.EQ.3.AND.BNDMXZ.NE.3)) THEN

           call ERROR ('GRIDBU','Incorrect periodic boundary.' , __FILE__, __LINE__)

        ENDIF



        !     =================
        !       SET UP X-GRID
        !     =================

        IF (xCPU .NE. 1) THEN

 
           ! gridlx :: length of the x-grid (global).
           !     nx :: number of cells in the x-direction (local after DECOMPOSE_MPI_DOMAIN).
    
           delx    = gridlx/(REAL(nx)*xCPU)


	   IF (BASENM .EQ. 'SC3D') THEN

              DO ii= 1, nx+1
                 xznl(ii) = REAL(iniTramo(1,mpi_coords(1))+ii-2)*delx 
              ENDDO


           ELSE IF (BASENM .EQ. 'S433') THEN

                 DO ii= 1, nx+1

                    IF (iniTramo(1,mpi_coords(1))+ii-2 .GT. 10 .AND. iniTramo(1,mpi_coords(1))+ii-2 .LE. 410) THEN
                        N = 0
                    ELSE
                        N = iniTramo(1,mpi_coords(1))+ii-2 - 210
                    ENDIF

                    R = 1.2821069420985605

                    xznl(ii) = REAL(iniTramo(1,mpi_coords(1))+ii-2-10)*0.1D0 -200*0.1D0

                    IF (N .LT. 0) THEN
                        xznl(ii) = -200*0.1D0 - R*(1.D0-R**abs(200+N))/(1.D0-R)*0.1D0
                    ELSE IF (N .GT. 0) THEN
                        xznl(ii) = 200*0.1D0 + R*(1.D0-R**(N-200))/(1.D0-R)*0.1D0
                    ENDIF

                 ENDDO


           ELSE IF (BASENM .EQ. 'S430' .OR. BASENM .EQ. 'S431' .OR. BASENM .EQ. 'S432') THEN

                 DO ii= 1, nx+1

                    IF (iniTramo(1,mpi_coords(1))+ii-2 .GT. 150 .AND. iniTramo(1,mpi_coords(1))+ii-2 .LE. 1150) THEN
                        N = 0
                    ELSE
                        N = iniTramo(1,mpi_coords(1))+ii-2 - 650
                    ENDIF

                    R = 1.019553802767D0

                    xznl(ii) = REAL(iniTramo(1,mpi_coords(1))+ii-2-150)*0.1D0 -500*0.1D0

                    IF (N .LT. 0) THEN
                        xznl(ii) = -500*0.1D0 - R*(1.D0-R**abs(500+N))/(1.D0-R)*0.1D0
                    ELSE IF (N .GT. 0) THEN
                        xznl(ii) = 500*0.1D0 + R*(1.D0-R**(N-500))/(1.D0-R)*0.1D0
                    ENDIF

                 ENDDO


	   ELSE

              DO ii= 1, nx+1
                 xznl(ii) = REAL(iniTramo(1,mpi_coords(1))+ii-2)*delx - gridlx/2.D0
              ENDDO

	   ENDIF


           DO ii= 1, nx
              xznr(ii) = xznl(ii+1)
           ENDDO

	   xznr(nx+1) = xznr(nx) + delx

           DO ii= 1, nx
              xzn(ii) = 0.5D0*(xznr(ii)+xznl(ii))
           ENDDO


        ELSE

           IF (BASENM(1:2).EQ.'BA'.OR.BASENM.EQ.'RROT'.OR.BASENM .EQ. '3ROT' .OR. &
               BASENM .EQ. 'CONT' .OR. BASENM .EQ. 'KH3D' .OR. BASENM .EQ. 'KH2D') THEN
              XZNL(1) = -0.5D0
           ELSE IF (BASENM.EQ.'GR2_'.OR.BASENM(1:2).EQ.'KO') THEN
              XZNL(1) = -2.D0
           ELSE IF (BASENM.EQ.'CMBW' .OR. BASENM .EQ. 'SMBW') THEN
              XZNL(1) = -6.D0
           ELSE IF (BASENM .EQ. 'IMPL' .OR. BASENM .EQ. 'DZBU' .OR. BASENM .EQ. 'SC3D') THEN
              XZNL(1) = 0.D0
	   ELSE IF (BASENM .EQ. 'SLKO' .OR. BASENM .EQ. 'ORTA') THEN
	      XZNL(1) = 0.D0
	   ELSE IF ((BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'H') .OR.  &
                    (BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'C') .OR. BASENM .EQ. 'S433') THEN
	      XZNL(1) = -10.D0
	   ELSE IF (BASENM .EQ. 'MI10') THEN
	      XZNL(1) = -8.D0
	   ELSE IF (BASENM .EQ. 'LOOP') THEN
              XZNL(1) = -5.D0
	   ELSE IF (BASENM .EQ. 'SYBW') THEN
              XZNL(1) = -20.D0
           ELSE
              XZNL(1) = 0.D0
           END IF

           DELX = GRIDLX/NX

           DO  I=2,NX+1
              XZNL(I) = XZNL(I-1) + DELX
           ENDDO

           DO  I=1,NX+1
              XZNR(I) = XZNL(I) + DELX
           ENDDO

           DO  I=1,NX
              XZN(I)  = 0.5D0*(XZNL(I) + XZNR(I))
           ENDDO

           YZNL(1) = 0.D0
           YZNR(1) = 0.D0
           YZN(1)  = 0.D0

           ZZNL(1) = 0.D0
           ZZNR(1) = 0.D0
           ZZN(1)  = 0.D0

        ENDIF

     
        !     =================
        !       SET UP Y-GRID
        !     =================

        IF (nsdim .GE. 2) THEN


           IF (yCPU .NE. 1) THEN

              ! gridly :: length of the y-grid (global).
              !     ny :: number of cells in the y-direction (local after DECOMPOSE_MPI_DOMAIN).

              dely = gridly/(REAL(ny)*yCPU)

	      IF ((BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'H')  .OR.   &
                  (BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'C')  .OR.   &
                   BASENM .EQ. 'LOOP') THEN
             
                 DO ii= 1, ny+1
                    yznl(ii) = REAL(iniTramo(2,mpi_coords(2))+ii-2)*dely
                 ENDDO


              ELSE IF (BASENM .EQ. 'S433') THEN

                 DO ii= 1, ny+1

                    IF (iniTramo(2,mpi_coords(2))+ii-2 .LE. 400) THEN
                        N = 0
                    ELSE
                        N = iniTramo(2,mpi_coords(2))+ii-2 - 400
                    ENDIF

                    R = 1.137074296153224

                    yznl(ii) = REAL(iniTramo(2,mpi_coords(2))+ii-2)*0.1D0

                    IF (N .NE. 0) THEN
                        yznl(ii) = 400*0.1D0 + R*(1.D0-R**N)/(1.D0-R)*0.1D0
                    ENDIF

                 ENDDO


              ELSE IF (BASENM .EQ. 'S430' .OR. BASENM .EQ. 'S431' .OR. BASENM .EQ. 'S432') THEN

                 DO ii= 1, ny+1

                    IF (iniTramo(2,mpi_coords(2))+ii-2 .LE. 1000) THEN
                        N = 0
                    ELSE
                        N = iniTramo(2,mpi_coords(2))+ii-2 - 1000
                    ENDIF

                    R = 1.0174627181D0

                    yznl(ii) = REAL(iniTramo(2,mpi_coords(2))+ii-2)*0.1D0

                    IF (N .NE. 0) THEN
                        yznl(ii) = 1000*0.1D0 + R*(1.D0-R**N)/(1.D0-R)*0.1D0
                    ENDIF

                 ENDDO


	      ELSE

                 DO ii= 1, ny+1
                    yznl(ii) = REAL(iniTramo(2,mpi_coords(2))+ii-2)*dely - gridly/2.D0
                 ENDDO

	      ENDIF  
             
              DO ii= 1, ny
                 yznr(ii) = yznl(ii+1)
              ENDDO

              yznr(ny+1) = yznr(ny) + dely

              DO ii= 1, ny
                 yzn(ii) = 0.5D0*(yznr(ii)+yznl(ii))
              ENDDO

           ELSE

              IF (BASENM.EQ.'CMBW' .OR. BASENM .EQ. 'SMBW') THEN
                 YZNL(1) = -6.D0
              ELSE IF (BASENM.EQ.'RROT' .OR. BASENM .EQ. '3ROT' .OR. BASENM .EQ. 'SC3D') THEN
                 YZNL(1) = -0.5D0
              ELSE IF (BASENM .EQ. 'IMPL') THEN
                 YZNL(1) = -0.5D0
              ELSE IF (BASENM .EQ. 'LOOP' .OR. BASENM .EQ. 'ORTA') THEN
                 YZNL(1) = 0.D0
	      ELSE IF (BASENM .EQ. 'SLKO') THEN
		 YZNL(1) = -14.D0
	      ELSE IF (BASENM .EQ. 'KH3D' .OR. BASENM .EQ. 'KH2D') THEN
		 YZNL(1) = -1.D0
	      ELSE IF ((BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'H') .OR.  &
                       (BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'C') .OR.  &
                        BASENM .EQ. 'MI10') THEN
		 YZNL(1) = 0.D0
              ELSE
                 YZNL(1) = 0.D0
              END IF

              DELY = GRIDLY/NY

              DO I=2,NY+1
                 YZNL(I) = YZNL(I-1) + DELY
              ENDDO

              DO  I=1,NY+1
                 YZNR(I) = YZNL(I) + DELY
              ENDDO

              DO I=1,NY
                 YZN(I)  = 0.5D0*(YZNL(I) + YZNR(I))
              ENDDO

              IF (nsdim .NE. 3) THEN

                 ZZNL(1) = 0.D0
                 ZZNR(1) = 0.D0
                 ZZN(1)  = 0.D0

              ENDIF

           ENDIF

	ENDIF
  
        IF (nsdim .EQ. 3) THEN

           IF (zCPU .NE. 1) THEN

              ! gridlz :: length of the z-grid (global).
              !     nz :: number of cells in the z-direction (local after DECOMPOSE_MPI_DOMAIN).

              delz = gridlz/(REAL(nz)*zCPU)


              IF (BASENM .EQ. 'S433') THEN

                 DO ii= 1, nz+1

                    IF (iniTramo(3,mpi_coords(3))+ii-2 .GT. 10 .AND. iniTramo(3,mpi_coords(3))+ii-2 .LE. 410) THEN
                        N = 0
                    ELSE
                        N = iniTramo(3,mpi_coords(3))+ii-2 - 210
                    ENDIF

                    R = 1.2821069420985605

                    zznl(ii) = REAL(iniTramo(3,mpi_coords(3))+ii-2-10)*0.1D0 -200*0.1D0

                    IF (N .LT. 0) THEN
                        zznl(ii) = -200*0.1D0 - R*(1.D0-R**abs(200+N))/(1.D0-R)*0.1D0
                    ELSE IF (N .GT. 0) THEN
                        zznl(ii) = 200*0.1D0 + R*(1.D0-R**(N-200))/(1.D0-R)*0.1D0
                    ENDIF

                 ENDDO


              ELSE IF (BASENM .EQ. 'S430' .OR. BASENM .EQ. 'S431' .OR. BASENM .EQ. 'S432') THEN

                 DO ii= 1, nz+1

                    IF (iniTramo(3,mpi_coords(3))+ii-2 .GT. 150 .AND. iniTramo(3,mpi_coords(3))+ii-2 .LE. 1150) THEN
                        N = 0
                    ELSE
                        N = iniTramo(3,mpi_coords(3))+ii-2 - 650
                    ENDIF

                    R = 1.019553802767D0

                    zznl(ii) = REAL(iniTramo(3,mpi_coords(3))+ii-2-150)*0.1D0 -500*0.1D0

                    IF (N .LT. 0) THEN
                        zznl(ii) = -500*0.1D0 - R*(1.D0-R**abs(500+N))/(1.D0-R)*0.1D0
                    ELSE IF (N .GT. 0) THEN
                        zznl(ii) = 500*0.1D0 + R*(1.D0-R**(N-500))/(1.D0-R)*0.1D0
                    ENDIF

                 ENDDO

	      ELSE	             

                  DO ii= 1, nz+1
                     zznl(ii) = REAL(iniTramo(3,mpi_coords(3))+ii-2)*delz - gridlz/2.D0
                  ENDDO  

	      ENDIF

             
              DO ii= 1, nz
                 zznr(ii) = zznl(ii+1)
              ENDDO

              zznr(nz+1) = zznr(nz) + delz

              DO ii= 1, nz
                 zzn(ii) = 0.5D0*(zznr(ii)+zznl(ii))
              ENDDO

           ELSE

              IF (BASENM .EQ. 'SMBW') THEN
                 ZZNL(1) = -6.D0
	      ELSE IF( BASENM .EQ. '3ROT' .OR. BASENM .EQ. 'SC3D' .OR. BASENM .EQ. 'KH3D') THEN	
	         ZZNL(1) = -0.5D0
              ELSE IF (BASENM .EQ. 'MI10') THEN
                 ZZNL(1) = -8.D0
              ELSE IF (BASENM .EQ. 'LOOP') THEN
                 ZZNL(1) = -5.D0
	      ELSE IF ((BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'H') .OR.   &
                       (BASENM(1:1) .EQ. 'M' .AND. BASENM(2:2) .EQ. 'C') .OR. BASENM .EQ. 'S433') THEN
		 ZZNL(1) = -10.D0
	      ENDIF

              DELZ = GRIDLZ/NZ

              DO I=2,NZ+1
                 ZZNL(I) = ZZNL(I-1) + DELZ
              ENDDO

              DO  I=1,NZ+1
                 ZZNR(I) = ZZNL(I) + DELZ
              ENDDO

              DO I=1,NZ
                 ZZN(I)  = 0.5D0*(ZZNL(I) + ZZNR(I))
              ENDDO

           ENDIF

        ENDIF

	RETURN

        !-------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE GRID
! =============================================================================== 
