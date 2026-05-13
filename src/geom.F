
!*********************************
!     3D-RMHD CODE: LOSTREGO
!*********************************

! NAME:
!       GEOM

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
!       Calculate the sections of the cells and the cell volumes which are 
!       perpendicular to the direction of the sweep, in Cartesian
!       coordinates. Other geometries will be implemented in the future.

      
! ACADEMIC REFERENCE:


! ===============================================================================      
SUBROUTINE GEOM (fznp, nznp, igeom, dx, area, areal, arear, dvol)
! =============================================================================== 

        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        integer(is) :: fznp, nznp, igeom

        ! size -4:mn5
        real(rs)    :: dx(-4:*)

        !------------------------------------------------------ Output Variables

        ! size -4:mn5
        real(rs)   :: area(-4:*), areal(-4:*), arear(-4:*), dvol(-4:*)
 
        !------------------------------------------------------- Local Variables

        integer(is):: i, ig

        character(46) :: message1

        !-------------------------------------------------------------------------- BEGIN

        IF (igeom .NE. 0) THEN

           message1 = 'This coordinate system is not implemented yet.'
           call ERROR ('GEOM', message1, __FILE__, __LINE__)

        ENDIF

        ig = igeom + 1

        ! CARTESIAN COORDINATES

        IF (ig .EQ. 1) THEN

           DO i= fznp, nznp+1

              areal(i) = 1.D0
              arear(i) = 1.D0
              area(i)  = 1.D0
              dvol(i)  = dx(i)

           ENDDO

        ENDIF

        !-------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE GEOM
! =============================================================================== 
