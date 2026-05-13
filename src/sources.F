!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************


! NAME:
!       SOURCES

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
!       Compute source terms (pseudogravity, stelar mass, cooling...) inside each sweep.

! ACADEMIC REFERENCE:


! ===============================================================================       
SUBROUTINE SOURCES (i, j, k, xyzswp, stelar, cool, basenm, grav, physic, cr, dt, &
                    mnx5, mny5, mnz5)
! =============================================================================== 

	USE tipos
	USE workarrays
	USE parallel

        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        integer (is) :: i, j, k, xyzswp, stelar, cool
	integer (is) :: mnx5, mny5, mnz5

        character(4) :: basenm

        real(rs)     :: grav (-4:mny5), cr, dt

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)

        !------------------------------------------------------- Local Variables

        character(44) :: message1
        character(36) :: message2
        character(27) :: message3

        !-------------------------------------------------------------------------- BEGIN


        IF (xyzswp .EQ. 1) THEN

           ! Pseudogravity is only applied with Jet simulations.          
           IF (basenm (1:1) .EQ. 'J') THEN

              message1 = 'PSGRAV: Jet simulations are not implemented.'
              call ERROR ('SOURCES', message1, __FILE__, __LINE__)

              ! For future implementations....
              physic(i,j,k)%momenx = physic(i,j,k)%momenx         &
                                   + cr*dt*physic(i,j,k)%gravx

              physic(i,j,k)%energy = physic(i,j,k)%energy         &
                                   + cr*dt*physic(i,j,k)%gravx    &
                                   * physic(i,j,k)%velx

           ENDIF


        ELSE IF (xyzswp .EQ. 2) THEN

           ! Pseudogravity is only applied with Jet simulations. 

           IF (basenm (1:1) .EQ. 'J' .AND. basenm(2:2) .EQ. '3') THEN

              call ERROR ('SOURCES', message1, __FILE__, __LINE__)

              physic(i,j,k)%momeny = physic(i,j,k)%momeny         &
                                   + cr*dt*physic(i,j,k)%gravy

              physic(i,j,k)%energy = physic(i,j,k)%energy         &
                                   + cr*dt*physic(i,j,k)%gravy    &
                                   * physic(i,j,k)%vely


           ELSE IF (basenm (1:1) .EQ. 'J' .AND. basenm(2:2) .NE. '3') THEN  

              call ERROR ('SOURCES', message1, __FILE__, __LINE__)

              physic(i,j,k)%momeny = physic(i,j,k)%momeny         &
                                   + cr*dt*grav(j)

              physic(i,j,k)%energy = physic(i,j,k)%energy         &
                                   + cr*dt*grav(j)                &
                                   * physic(i,j,k)%vely 

           ENDIF

           message2 = 'Stelar mass load is not implemented.'                      
           IF (stelar .EQ. 1) call ERROR ('SOURCES', message2, __FILE__, __LINE__)

           message3 = 'Cooling is not implemented.'
           IF (cool .EQ. 1)   call ERROR ('SOURCES', message3, __FILE__, __LINE__)


        ELSE

           ! Pseudogravity is only applied with Jet simulations. 

           IF (basenm (1:1) .EQ. 'J' .AND. basenm(2:2) .EQ. '3') THEN

              physic(i,j,k)%momenz = physic(i,j,k)%momenz         &
                                   + cr*dt*physic(i,j,k)%gravz

              physic(i,j,k)%energy = physic(i,j,k)%energy         &
                                   + cr*dt*physic(i,j,k)%gravz    &
                                   * physic(i,j,k)%velz

           ENDIF

           IF (stelar .EQ. 1) call ERROR ('SOURCES', message2, __FILE__, __LINE__)
           IF (cool .EQ. 1)   call ERROR ('SOURCES', message3, __FILE__, __LINE__)

        ENDIF


        !-------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE SOURCES
! =============================================================================== 
