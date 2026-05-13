
!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************

! NAME:
!       ERROR 

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
!       Debug: Print an error message.

! ACADEMIC REFERENCE:



! ===============================================================================      
SUBROUTINE ERROR (routine, message, fich, line)
! ===============================================================================

        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(*) :: routine
        character(*) :: message
        character(*) :: fich
        integer(is)  :: line

        !------------------------------------------------------- Output Variables

        ! ...

        !------------------------------------------------------- Local Variables

        ! ...

        !------------------------------------------------------------------------- BEGIN

        write(*,*) '*************************************************************'
        write(*,*) 'ERROR'
        write(*,*) ' '
        write(*,*) routine, ': ', message
        write(*,*) ' '
        write(*,*) '           File: ', fich
        write(*,*) '           Line: ', line
        write(*,*) ' '
        write(*,*) '*************************************************************'

        STOP 'ERROR'

        !--------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE ERROR
! ===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! NAME:
!       WARNING

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
!       Show an advertisement when a variable is updated.

! ACADEMIC REFERENCE:



! ===============================================================================      
SUBROUTINE WARNING (routine, message, fich, line)
! ===============================================================================

        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(*) :: routine
        character(*) :: message
        character(*) :: fich
        integer(is)  :: line

        !------------------------------------------------------- Output Variables

        ! ...

        !------------------------------------------------------- Local Variables

        ! ...

        !------------------------------------------------------------------------- BEGIN

        write(*,*) '*************************************************************'
        write(*,*) 'WARNING'
        write(*,*) ' '
        write(*,*) routine, ': ', message
        write(*,*) ' '
        write(*,*) '           File: ', fich
        write(*,*) '           Line: ', line
        write(*,*) ' '
        write(*,*) '*************************************************************'

        ! The program does not stop after a WARNING message.

        !--------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE WARNING
! ===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


! NAME:
!       int2ch4 

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
!       Convert a number into an array of characteres.

! ACADEMIC REFERENCE:



! ===============================================================================
FUNCTION int2ch4 (number)
! ===============================================================================      

        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        integer(is)  :: number

        !------------------------------------------------------- Output Variables

        character(4) :: int2ch4

        !------------------------------------------------------- Local Variables

        integer(is)  :: tt, cero, jj
        character    :: cc (4)

        !------------------------------------------------------------------------- BEGIN

        IF (number .GT. 9999 .OR. number .LT. 0) THEN
           call ERROR ('int2ch4', 'The number is out of range', __FILE__, __LINE__)
        ELSE

           cero = ICHAR ('0')  ! 48 is the ascii code of zero
           tt   = number

           DO jj = 4, 1, -1

              cc(jj) = CHAR (cero + MOD(tt,10))
              tt = tt/10

           ENDDO 

           int2ch4 = cc(1)//cc(2)//cc(3)//cc(4)

        ENDIF

        return

        !--------------------------------------------------------------------------- END

! ===============================================================================
END FUNCTION int2ch4
! ===============================================================================


