

!===============================================================================
! NAME              :: INTERCAMBIABOUNDCT
! GLOBAL VARIABLES  :: mpi_coords, mpi_dims, xright, xleft, yright, yleft, zright, zleft.
! DESCRIPTION       :: Send and Receive directives for the interchange of boundary staggered
!                      field between adjacent nodes.
!===============================================================================

     
SUBROUTINE INTERCAMBIABOUNDCT (CT, nx, ny, nz, mnx6, mny6, mnz6)

        !   ================================
        !     EXCHANGE RIGHT SIDE (Node i)
        !   ================================ 

        !                                       ....
        !                                       |  |
        !      Node i-1     |         Node i    || v    Node i+1     |
        ! ------------------|-----------------|--|--|----------------|----------
        !                   |                 |--|--|                |
        !                   |                 |--|--|                |  
        !                   |                 |--|--|                |   BOUNDARY
        !                   |                 |--|--|                |  CONDITIONS
        !                   |                 |--|--|                |
        !                   |                 |--|--|                |
        ! ------------------x-----------------|--x--|----------------x-----------
        !                                       ^  |                 
        !                                       |  |
        !                                       *--*
        !                                      
        ! PHILOSOPHY:
        ! i-Node SENDS the left boundary physical variables to i+1-Node.
        ! i+1-Node SENDS the right boundary physical variables to i-Node.


        !   ================================
        !     EXCHANGE LEFT SIDE (Node i)
        !   ================================ 

        !                  ....                    
        !                  |  |                  
        !      Node i-1    || v      Node i    |     |   Node i+1     |
        ! ---------------|--|--|---------------|--|--|----------------|----------
        !                |--|--|               |--|--|                |
        !                |--|--|               |--|--|                |  
        !                |--|--|               |--|--|                |  
        !                |--|--|               |--|--|                | 
        !                |--|--|               |--|--|                |
        !                |--|--|               |--|--|                |
        ! ---------------|--x--|---------------|--x--|----------------x-----------
        !                 ^  |                 
        !                 |  |
        !                 *--*
        !                                      
        ! PHILOSOPHY:
        ! i-Node RECEIVES the left boundary physical variables from i-1-Node.
        ! i-1-Node RECEIVES the right boundary physical variables from i-Node.

        USE tipos
        USE parallel
        USE workarrays
        USE memoria

	implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

# include "types.h"


        !------------------------------------------------------- Input Variables

        integer(is) :: nx, ny, nz, mnx6, mny6, mnz6

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
 
        !------------------------------------------------------- Local Variables

# ifdef PARALELO

            integer(is) :: lon, ierr, req_cnt 
            integer(is) :: status (MPI_STATUS_SIZE)
            integer(is) :: arr_request (12)
      
            !type(TboundPhysic), pointer :: boundsX (:,:,:), boundsY (:,:,:), boundsZ (:,:,:)    
            integer ::  arr_status (MPI_STATUS_SIZE, 12)

            
            !--------------------------------------------------------------------------- BEGIN
      

            ! Send / Receive counter.
            req_cnt = 0
 
# define bptagxleft 52345
# define bptagxright 52346
# define bptagyleft 52347
# define bptagyright 52348
# define bptagzleft 52349 
# define bptagzright 52350

            ! ---------    
            !    MEMORY
            ! --------- 

            ! Periodic boundary
            !allocate (boundsX(5,ny,nz), stat= ierr)
            !IF (ierr .NE. 0) THEN
            !   call ERROR ('INTERCAMBIABOUNDPHYSIC',                              &
            !               'Allocation of bounds', __FILE__, __LINE__) 
            !ENDIF

            !call MEMPUSH (5*ny*nz, rs, 'boundsX')
 
            !allocate (boundsY(nx,5,nz), stat= ierr)
            !IF (ierr .NE. 0) THEN
            !   call ERROR ('intercambiaBoundPhysic',                              &
            !               'Allocation of bounds', __FILE__, __LINE__) 
            !ENDIF

            !call MEMPUSH (nx*5*nz, rs, 'boundsY')
 
            !allocate (boundsZ(nx,ny,5), stat= ierr) 
            !IF (ierr .NE. 0) THEN
            !   call ERROR ('intercambiaBoundPhysic',                              &
            !               'Allocation of bounds', __FILE__, __LINE__) 
            !ENDIF

            !call MEMPUSH (nx*ny*5, rs, 'boundsZ')

            ! Non-periodic boundary (defined in modulos.F).
            call GETBOUNDSMEMORYSPACECT (nx, ny, nz)



            ! --------------------------
            ! NON-PERIODIC BOUNDARY
            ! --------------------------

       
               !   ================================
               !     EXCHANGE RIGHT SIDE (Node i)
               !   ================================ 
           
               !-------
               ! X-dir
               !-------
                    
               ! Last node doesn't change information from the right side.
               IF (mpi_coords(1) .NE. (mpi_dims(1)-1)) THEN
              
                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = 2*(ny+1)*nz*2
         
                  call FROMCTTOBOUNDX (nx-1, nx, ny, nz, mnx6, mny6, mnz6, CTboundsXs1xy, CT)


                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTboundsXs1xy, lon, MPI_DOUBLE_PRECISION, xright, bptagxright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)
        
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &                  
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTboundsXr1xy, lon, MPI_DOUBLE_PRECISION, xright,                 &
                                  bptagxleft, CART_WORLD,  arr_request(req_cnt), ierr)
              
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF        
              
               ENDIF
        
               !-------
               ! Y-dir
               !-------
                    
               ! Last node doesn't change information from the right side.
               IF (mpi_coords(2) .NE. (mpi_dims(2)-1)) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = nx*2*(nz+1)*2

                  call FROMCTTOBOUNDY (ny-1, ny, nx, nz, mnx6, mny6, mnz6, CTboundsYs1yz, CT)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTboundsYs1yz, lon, MPI_DOUBLE_PRECISION, yright,                 &
                                  bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTboundsYr1yz, lon, MPI_DOUBLE_PRECISION, yright,                 &
                                  bptagyleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      & 
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF
          
               ENDIF

               !-------
               ! Z-dir
               !-------
                    
               ! Last node doesn't change information from the right side.
               IF (mpi_coords(3) .NE. (mpi_dims(3)-1)) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = (nx+1)*ny*2*2

                  call FROMCTTOBOUNDZ (nz-1, nz, nx, ny, mnx6, mny6, mnz6, CTboundsZs1zx, CT)
         
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTboundsZs1zx, lon, MPI_DOUBLE_PRECISION, zright,                 &
                                  bptagzright, CART_WORLD, arr_request(req_cnt), ierr )

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ( 'INTERCAMBIABOUNDPHYSIC',                                     &
                                  'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTboundsZr1zx, lon, MPI_DOUBLE_PRECISION, zright,                 &
                                  bptagzleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF                 

               ENDIF


               !   ================================
               !     EXCHANGE LEFT SIDE (Node i)
               !   ================================ 

               !-------
               ! X-dir
               !-------
 
               ! First node doesn't change information from the left side.
               IF (mpi_coords(1) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = 2*(ny+1)*nz*2

                  call FROMCTTOBOUNDX (1, 2, ny, nz, mnx6, mny6, mnz6, CTboundsXs2xy, CT)
         
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTboundsXs2xy, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                                  bptagxleft, CART_WORLD, arr_request(req_cnt), ierr)
 
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTboundsXr2xy, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                                  bptagxright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF         
          
               ENDIF


               !-------
               ! Y-dir
               !-------

               ! First node doesn't change information from the left side.
               IF (mpi_coords(2) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = nx*2*(nz+1)*2

                  call FROMCTTOBOUNDY (1, 2, nx, nz, mnx6, mny6, mnz6, CTboundsYs2yz, CT)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTboundsYs2yz, lon, MPI_DOUBLE_PRECISION, yleft,                  &
                                  bptagyleft, CART_WORLD, arr_request(req_cnt),ierr )
                  
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
                  
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTboundsYr2yz, lon, MPI_DOUBLE_PRECISION, yleft,                  &
                                  bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF      

               ENDIF


               !-------
               ! Z-dir
               !-------

               ! First node doesn't change information from the left side.
               IF (mpi_coords(3) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = (nx+1)*ny*2*2

                  call FROMCTTOBOUNDZ (1, 2, nx, ny, mnx6, mny6, mnz6, CTboundsZs2zx, CT)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTboundsZs2zx, lon, MPI_DOUBLE_PRECISION, zleft,                  &
                                  bptagzleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTboundsZr2zx, lon, MPI_DOUBLE_PRECISION, zleft,                  &
                                  bptagzright, CART_WORLD, arr_request(req_cnt), ierr)

          
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv',  __FILE__, __LINE__)
                  ENDIF
                  
               ENDIF

               ! Wait for all the communications to finish.
               call MPI_WAITALL (req_cnt, arr_request, arr_status, ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',                                         &
                              'Error in MPI_WAITALL', __FILE__, __LINE__)
               ENDIF


               IF (mpi_coords(1) .NE. (mpi_dims(1)-1)) THEN
                  call FROMBOUNDTOCTX (nx+1, nx+2, ny, nz, mnx6, mny6, mnz6, CTboundsXr1xy, CT) 
               ENDIF

               IF (mpi_coords(2) .NE. (mpi_dims(2)-1)) THEN
                  call FROMBOUNDTOCTY (ny+1, ny+2, nx, nz, mnx6, mny6, mnz6, CTboundsYr1yz, CT) 
               ENDIF

               IF (mpi_coords(3) .NE. (mpi_dims(3)-1)) THEN
                  call FROMBOUNDTOCTZ (nz+1, nz+2, nx, ny, mnx6, mny6, mnz6, CTboundsZr1zx, CT) 
               ENDIF

               IF (mpi_coords(1) .GT. 0) THEN
                  call FROMBOUNDTOCTX (-1, 0, ny, nz, mnx6, mny6, mnz6, CTboundsXr2xy, CT)
               ENDIF

               IF (mpi_coords(2) .GT. 0) THEN
                  call FROMBOUNDTOCTY (-1, 0, nx, nz, mnx6, mny6, mnz6, CTboundsYr2yz, CT)
               ENDIF

               IF (mpi_coords(3) .GT. 0) THEN
                  call FROMBOUNDTOCTZ (-1, 0, nx, ny, mnx6, mny6, mnz6, CTboundsZr2zx, CT)
               ENDIF



            ! ---------    
            !    MEMORY
            ! --------- 

            ! Periodic boundary.

            !deallocate (boundsX)
            !call MEMPOP ('boundsX')
            !deallocate (boundsY)
            !call MEMPOP ('boundsY')
            !deallocate (boundsZ)
            !call MEMPOP ('boundsZ')

            ! Non-periodic boundary.

            call FREEBOUNDSMEMORYSPACECT ()

# endif

        !--------------------------------------------------------------------------- END

END SUBROUTINE INTERCAMBIABOUNDCT


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOCTX (ini, fin, ny, nz, mnx6, mny6, mnz6, boundsCT, CT)

        USE tipos
        USE parallel
        USE workarrays
        USE memoria
        implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

#include <types.h>

        !------------------------------------------------------- Input Variables

        integer(is) :: ini, fin, ny, nz, mnx6, mny6, mnz6
        type(TboundCTxy) :: boundsCT (2,ny+1,nz)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN


        DO kk= 1, nz
           DO jj= 1, ny+1
              DO ii= ini, fin
                 
                 CT(ii,jj,kk)%fbxy     = boundsCT(ii-ini+1,jj,kk)%fbxy
                 CT(ii,jj,kk)%rflxy    = boundsCT(ii-ini+1,jj,kk)%rflxy
                 
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOCTX

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Send y-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOCTY (ini, fin, nx, nz, mnx6, mny6, mnz6, boundsCT, CT)

        USE tipos
        USE parallel
        USE workarrays
        USE memoria
        implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

#include <types.h>

        !------------------------------------------------------- Input Variables

        integer(is) :: ini, fin, nx, nz, mnx6, mny6, mnz6
        type(TboundCTyz) :: boundsCT (nx,2,nz+1)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz+1
           DO jj= ini, fin
              DO ii= 1, nx

                 CT(ii,jj,kk)%fbyz     = boundsCT(ii,jj-ini+1,kk)%fbyz
                 CT(ii,jj,kk)%rflxz    = boundsCT(ii,jj-ini+1,kk)%rflxz

              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOCTY


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYZ
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save z-boundary information into physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOCTZ (ini, fin, nx, ny, mnx6, mny6, mnz6, boundsCT, CT)

        USE tipos
        USE parallel
        USE workarrays
        USE memoria
        implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

#include <types.h>

        !------------------------------------------------------- Input Variables

        integer(is) :: ini, fin, nx, ny, mnx6, mny6, mnz6
        type(TboundCTzx) :: boundsCT (nx+1,ny,2)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= ini, fin
           DO jj= 1, ny
              DO ii= 1, nx+1
                 
                 CT(ii,jj,kk)%fbzx     = boundsCT(ii,jj,kk-ini+1)%fbzx
                 CT(ii,jj,kk)%rflxx    = boundsCT(ii,jj,kk-ini+1)%rflxx
                 
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOCTZ

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMCTTOBOUNDX (ini, fin, ny, nz, mnx6, mny6, mnz6, boundsCT, CT)

        USE tipos
        USE parallel
        USE workarrays
        USE memoria
        implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

#include <types.h>

        !------------------------------------------------------- Input Variables

        integer(is) :: ini, fin, ny, nz, mnx6, mny6, mnz6
        type(TboundCTxy) :: boundsCT (2,ny+1,nz)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN


        DO kk= 1, nz
           DO jj= 1, ny+1
              DO ii= ini, fin

                 boundsCT(ii-ini+1,jj,kk)%fbxy     = CT(ii,jj,kk)%fbxy
                 boundsCT(ii-ini+1,jj,kk)%rflxy     = CT(ii,jj,kk)%rflxy
                 
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMCTTOBOUNDX


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMCTTOBOUNDY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save y-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMCTTOBOUNDY (ini, fin, nx, nz, mnx6, mny6, mnz6, boundsCT, CT)

        USE tipos
        USE parallel
        USE workarrays
        USE memoria
        implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

#include <types.h>

        !------------------------------------------------------- Input Variables

        integer(is) :: ini, fin, nx, nz, mnx6, mny6, mnz6
        type(TboundCTyz) :: boundsCT (nx,2,nz+1)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN


        DO kk= 1, nz+1
           DO jj= ini, fin
              DO ii= 1, nx

                 boundsCT(ii,jj-ini+1,kk)%fbyz     = CT(ii,jj,kk)%fbyz
                 boundsCT(ii,jj-ini+1,kk)%rflxz     = CT(ii,jj,kk)%rflxz

                 
              ENDDO
           ENDDO
        ENDDO




        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMCTTOBOUNDY


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDZ
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save z-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMCTTOBOUNDZ (ini, fin, nx, ny, mnx6, mny6, mnz6, boundsCT, CT)

        USE tipos
        USE parallel
        USE workarrays
        USE memoria
        implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

#include <types.h>

        !------------------------------------------------------- Input Variables

        integer(is) :: ini, fin, nx, ny, mnx6, mny6, mnz6
        type(TboundCTzx) :: boundsCT (nx+1,ny,2)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= ini, fin
           DO jj= 1, ny
              DO ii= 1, nx+1
                 
                 boundsCT(ii,jj,kk-ini+1)%fbzx     = CT(ii,jj,kk)%fbzx
                 boundsCT(ii,jj,kk-ini+1)%rflxx     = CT(ii,jj,kk)%rflxx
                 
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMCTTOBOUNDZ


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

