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
!       Construct new filnames for output and restart files.

! ACADEMIC REFERENCE:


! =============================================================================== 
SUBROUTINE FILNAM (chd, basenm, datfil, rstfil, outfil, tdlfil, suffixd, suffix)
! =============================================================================== 

        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(4):: basenm
        character(1):: chd

        !------------------------------------------------------ Output Variables

        character(8):: datfil
        character(7):: rstfil, outfil, tdlfil
        character(3):: suffixd
        character(2):: suffix

        !------------------------------------------------------- Local Variables

        character(1):: sf1, sf2, sf3

        integer(is) :: isf1, isf2, isf3

        !------------------------------------------------------------------------- BEGIN

        IF (chd .NE. 'd') THEN

           ! Suffix will take values aa, ab, ac...az, ba, bb, bc....bz...
           IF (suffix(2:2) .EQ. 'z' .OR. suffix(2:2) .EQ. 'Z') THEN

              sf1  = suffix(1:1)
              isf1 = ICHAR(sf1)
              sf2  = suffix(2:2)
              isf2 = ICHAR(sf2)

              isf1 = isf1 + 1
              isf2 = isf2 - 25  !Restart z --> a
              suffix(1:1) = CHAR(isf1)
              suffix(2:2) = CHAR(isf2)
              
           ELSE

              sf2  = suffix(2:2)
              isf2 = ICHAR(sf2)
              isf2 = isf2 + 1
              suffix(2:2) = CHAR(isf2)

           ENDIF
        
              ! Construct restart/output file names.
              rstfil = basenm // 'R' // suffix ! Restart file
              outfil = basenm // 'O' // suffix ! Output file
           
        ELSE

           ! Suffix will take values aaa, aab, aac...aaz, aba, abb, abc....
           IF (suffixd(3:3) .EQ. 'z' .OR. suffixd(3:3) .EQ. 'Z') THEN

              sf1  = suffixd(1:1)
              isf1 = ICHAR(sf1)
              sf2  = suffixd(2:2)
              isf2 = ICHAR(sf2)
              sf3  = suffixd(3:3)
              isf3 = ICHAR(sf3)

              IF (suffixd(2:2) .EQ. 'z' .OR. suffixd(2:2) .EQ. 'Z') THEN
                    
                 isf1 = isf1 + 1
                 isf2 = isf2 - 25
                 isf3 = isf3 - 25
                 suffixd(1:1) = CHAR(isf1)
                 suffixd(2:2) = CHAR(isf2)
                 suffixd(3:3) = CHAR(isf3)

              ELSE

                 isf2 = isf2 + 1
                 isf3 = isf3 - 25
                 suffixd(2:2) = CHAR(isf2)
                 suffixd(3:3) = CHAR(isf3)

              ENDIF

           ELSE

              sf3  = suffixd(3:3)
              isf3 = ICHAR(sf3)
              isf3 = isf3 + 1
              suffixd(3:3) = CHAR(isf3)
              
           ENDIF

           tdlfil = basenm // suffixd
           datfil = basenm // 'D' // suffixd
      
        ENDIF
      
        !--------------------------------------------------------------------------- END

! =============================================================================== 
END SUBROUTINE FILNAM
! =============================================================================== 
