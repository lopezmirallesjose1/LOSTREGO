!*********************************
!     3D-RMHD CODE: LOSTREGO
!*********************************

! NAME:
!       PLTOUT

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
!       Generate Output (Fortran formatted) files for 1D tests.

! ACADEMIC REFERENCE:



! ===============================================================================
SUBROUTINE PLTOUT (path, basenm, nsdim, nx, ny, nz,        &
                   mnx5, mny5, mnz5, xzn, physic,          &
                   rstfil, outfil, tdlfil, datfil,         &
                   nout1, tout1, suffix, suffixd,          &
	           limiter, nend, irstrt, nrstrt,          &
	           nout, itstp, iars, irk, niter,          &
		   igeomx, igeomy, igeomz, eqs,            &
		   tmax, trstrt, tout, cfl, dtini,         &
                   small, smlrho, smallp, smallu,          &
		   gridlx, gridly, gridlz)
! ===============================================================================

                   
        USE math_constants
        USE tipos
        USE workarrays
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        character(256):: path
        character(4)  :: basenm

        integer(is):: nx, ny, nz, mnx5, mny5, mnz5, nsdim

        real(rs)   :: xzn(*)
 
        type(Tphysic):: physic (-4:mnx5,-4:mny5,-4:mnz5)

        character(4):: limiter

        integer(is) :: nend, irstrt, nrstrt, nout, itstp, iars,        &
                       irk,  niter,                                    &                           
                       igeomx, igeomy, igeomz, eqs

        real(rs)    :: tmax, trstrt, tout, cfl, dtini,                 &
                       small, smlrho, smallp, smallu,                  &
                       gridlx, gridly, gridlz                          

        !------------------------------------------------------ Output Variables

        character(7):: rstfil, outfil, tdlfil
        character(8):: datfil
        character(2):: suffix
        character(3):: suffixd

        integer(is) :: nout1
        real(rs) :: tout1

        !------------------------------------------------------- Local Variables

        character(256) :: fich

        character(4), external :: int2ch4

        integer(is)   :: i, j, k, ierr

        character(33) :: message1


        !-------------------------------------------------------------------------- BEGIN

        !  ========================
        !     1-DIMENSIONAL TEST
        !  ========================

        IF (basenm(1:2) .EQ. 'BA' .OR. basenm(1:2) .EQ. 'GR' .OR.   &
            basenm(1:2) .EQ. 'KO' .OR. basenm .EQ. 'DZAW' .OR. BASENM .EQ. 'CONT' .OR. &
	    basenm .EQ. 'DZBU') THEN

            fich = TRIM(path)//outfil//'_'//int2ch4(cartrank)
            open (12, file= fich, form= 'formatted', iostat= ierr)

            IF (ierr .NE. 0) THEN
               call ERROR ('PLTOUT', 'Error opening file', __FILE__, __LINE__)
            ENDIF

	    
	    !write(12,*) suffix
	    !write(12,*) suffixd
	    !write(12,*) basenm
            !write(12,*) limiter
	    write(12,*) nend
            write(12,*) irstrt
	    write(12,*) nrstrt
	    write(12,*) nout
	    write(12,*) itstp
            write(12,*) iars
            write(12,*) irk
	    write(12,*) niter
            write(12,*) nsdim
            write(12,*) nx
            write(12,*) ny
            write(12,*) nz
            write(12,*) igeomx
            write(12,*) igeomy
            write(12,*) igeomz
            write(12,*) eqs
	    write(12,*) tmax
            write(12,*) trstrt
	    write(12,*) tout
            write(12,*) cfl
            write(12,*) dtini
            write(12,*) small 
            write(12,*) smlrho
            write(12,*) smallp
            write(12,*) smallu
            write(12,*) gridlx
            write(12,*) gridly
            write(12,*) gridlz

      
            ! Write physical variables.
            
            DO i = 1, nx

               write (12,'(11(f16.8,1x))') xzn(i), physic(i,1,1)%pres,                            &
                                                   physic(i,1,1)%densty, physic(i,1,1)%denstye,   &
                                                   physic(i,1,1)%velx, physic(i,1,1)%vely,        &
                                                   physic(i,1,1)%velz, physic(i,1,1)%bxc,         &
                                                   physic(i,1,1)%byc, physic(i,1,1)%bzc,          & 
						   physic(i,1,1)%tracer
                                                 

            ENDDO
	       
                                       
            close (12)

        ELSE

           message1 = 'PLTOUT only works with 1D tests.'
           call ERROR ('PLTOUT', message1, __FILE__, __LINE__)

        ENDIF

        nout1=0
        tout1=0.0

        call FILNAM ('a', basenm, datfil, rstfil, outfil, tdlfil, suffixd, suffix)


        !---------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE PLTOUT
! ===============================================================================
