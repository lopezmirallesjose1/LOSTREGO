

!===============================================================================
! NAME              :: INTERCAMBIABOUNDCT
! GLOBAL VARIABLES  :: mpi_coords, mpi_dims, xright, xleft, yright, yleft, zright, zleft.
! DESCRIPTION       :: Send and Receive directives for the interchange of boundary staggered
!                      field between adjacent nodes.
!===============================================================================

     
SUBROUTINE INTERCAMBIABOUNDCTT (CT, nx, ny, nz, mnx6, mny6, mnz6)

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
 
# define bptagxleft 42345
# define bptagxright 42346
# define bptagyleft 42347
# define bptagyright 42348
# define bptagzleft 42349 
# define bptagzright 42350

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
            call GETBOUNDSMEMORYSPACECTT (nx, ny, nz)
 


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
                  lon = 2*ny*(nz+1)*2
         
                  call FROMCTTOBOUNDXT (nx-1, nx, ny, nz, mnx6, mny6, mnz6, CTTboundsXs1xz, CT)
                   
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTTboundsXs1xz, lon, MPI_DOUBLE_PRECISION, xright, bptagxright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)
        
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &                  
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTTboundsXr1xz, lon, MPI_DOUBLE_PRECISION, xright,                 &
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
                  lon = (nx+1)*2*nz*2

                  call FROMCTTOBOUNDYT (ny-1, ny, nx, nz, mnx6, mny6, mnz6, CTTboundsYs1yx, CT)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTTboundsYs1yx, lon, MPI_DOUBLE_PRECISION, yright,                 &
                                  bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTTboundsYr1yx, lon, MPI_DOUBLE_PRECISION, yright,                 &
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
                  lon = nx*(ny+1)*2*2

                  call FROMCTTOBOUNDZT (nz-1, nz, nx, ny, mnx6, mny6, mnz6, CTTboundsZs1zy, CT)
         
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTTboundsZs1zy, lon, MPI_DOUBLE_PRECISION, zright,                 &
                                  bptagzright, CART_WORLD, arr_request(req_cnt), ierr )

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ( 'INTERCAMBIABOUNDPHYSIC',                                     &
                                  'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTTboundsZr1zy, lon, MPI_DOUBLE_PRECISION, zright,                 &
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
                  lon = 2*ny*(nz+1)*2

                  call FROMCTTOBOUNDXT (1, 2, ny, nz, mnx6, mny6, mnz6, CTTboundsXs2xz, CT)
         
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTTboundsXs2xz, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                                  bptagxleft, CART_WORLD, arr_request(req_cnt), ierr)
 
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTTboundsXr2xz, lon, MPI_DOUBLE_PRECISION, xleft,                  &
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
                  lon = (nx+1)*2*nz*2

                  call FROMCTTOBOUNDYT (1, 2, nx, nz, mnx6, mny6, mnz6, CTTboundsYs2yx, CT)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTTboundsYs2yx, lon, MPI_DOUBLE_PRECISION, yleft,                  &
                                  bptagyleft, CART_WORLD, arr_request(req_cnt),ierr )
                  
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
                  
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTTboundsYr2yx, lon, MPI_DOUBLE_PRECISION, yleft,                  &
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
                  lon = nx*(ny+1)*2*2

                  call FROMCTTOBOUNDZT (1, 2, nx, ny, mnx6, mny6, mnz6, CTTboundsZs2zy, CT)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (CTTboundsZs2zy, lon, MPI_DOUBLE_PRECISION, zleft,                  &
                                  bptagzleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (CTTboundsZr2zy, lon, MPI_DOUBLE_PRECISION, zleft,                  &
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
                  call FROMBOUNDTOCTTX (nx+1, nx+2, ny, nz, mnx6, mny6, mnz6, CTTboundsXr1xz, CT) 
               ENDIF

               IF (mpi_coords(2) .NE. (mpi_dims(2)-1)) THEN
                  call FROMBOUNDTOCTTY (ny+1, ny+2, nx, nz, mnx6, mny6, mnz6, CTTboundsYr1yx, CT) 
               ENDIF

               IF (mpi_coords(3) .NE. (mpi_dims(3)-1)) THEN
                  call FROMBOUNDTOCTTZ (nz+1, nz+2, nx, ny, mnx6, mny6, mnz6, CTTboundsZr1zy, CT) 
               ENDIF

               IF (mpi_coords(1) .GT. 0) THEN
                  call FROMBOUNDTOCTTX (-1, 0, ny, nz, mnx6, mny6, mnz6, CTTboundsXr2xz, CT)
               ENDIF

               IF (mpi_coords(2) .GT. 0) THEN
                  call FROMBOUNDTOCTTY (-1, 0, nx, nz, mnx6, mny6, mnz6, CTTboundsYr2yx, CT)
               ENDIF

               IF (mpi_coords(3) .GT. 0) THEN
                  call FROMBOUNDTOCTTZ (-1, 0, nx, ny, mnx6, mny6, mnz6, CTTboundsZr2zy, CT)
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

            call FREEBOUNDSMEMORYSPACECTT ()

# endif

        !--------------------------------------------------------------------------- END

END SUBROUTINE INTERCAMBIABOUNDCTT


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOCTTX (ini, fin, ny, nz, mnx6, mny6, mnz6, boundsCT, CT)

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
        type(TboundCTxz) :: boundsCT (2,ny,nz+1)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz+1
           DO jj= 1, ny
              DO ii= ini, fin

                 CT(ii,jj,kk)%fbxz     = boundsCT(ii-ini+1,jj,kk)%fbxz
                 CT(ii,jj,kk)%rflxz     = boundsCT(ii-ini+1,jj,kk)%rflxz

              ENDDO
           ENDDO
        ENDDO


        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOCTTX

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Send y-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOCTTY (ini, fin, nx, nz, mnx6, mny6, mnz6, boundsCT, CT)

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
        type(TboundCTyx) :: boundsCT (nx+1,2,nz)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
           DO jj= ini, fin
              DO ii= 1, nx+1
                 
                 CT(ii,jj,kk)%fbyx     = boundsCT(ii,jj-ini+1,kk)%fbyx
                 CT(ii,jj,kk)%rflxx     = boundsCT(ii,jj-ini+1,kk)%rflxx
                 

              ENDDO
           ENDDO
        ENDDO


        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOCTTY


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYZ
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save z-boundary information into physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOCTTZ (ini, fin, nx, ny, mnx6, mny6, mnz6, boundsCT, CT)

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
        type(TboundCTzy) :: boundsCT (nx,ny+1,2)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN


        DO kk= ini, fin
           DO jj= 1, ny+1
              DO ii= 1, nx

                 CT(ii,jj,kk)%fbzy     = boundsCT(ii,jj,kk-ini+1)%fbzy
                 CT(ii,jj,kk)%rflxy     = boundsCT(ii,jj,kk-ini+1)%rflxy

                 
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOCTTZ

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMCTTOBOUNDXT (ini, fin, ny, nz, mnx6, mny6, mnz6, boundsCT, CT)

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
        type(TboundCTxz) :: boundsCT (2,ny,nz+1)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz+1
           DO jj= 1, ny
              DO ii= ini, fin

                 boundsCT(ii-ini+1,jj,kk)%fbxz     = CT(ii,jj,kk)%fbxz
                 boundsCT(ii-ini+1,jj,kk)%rflxz     = CT(ii,jj,kk)%rflxz

                 
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMCTTOBOUNDXT


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMCTTOBOUNDY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save y-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMCTTOBOUNDYT (ini, fin, nx, nz, mnx6, mny6, mnz6, boundsCT, CT)

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
        type(TboundCTyx) :: boundsCT (nx+1,2,nz)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
           DO jj= ini, fin
              DO ii= 1, nx+1

                 boundsCT(ii,jj-ini+1,kk)%fbyx     = CT(ii,jj,kk)%fbyx
                 boundsCT(ii,jj-ini+1,kk)%rflxx     = CT(ii,jj,kk)%rflxx

                 
              ENDDO
           ENDDO
        ENDDO


        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMCTTOBOUNDYT


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDZ
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save z-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMCTTOBOUNDZT (ini, fin, nx, ny, mnx6, mny6, mnz6, boundsCT, CT)

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
        type(TboundCTzy) :: boundsCT (nx,ny+1,2)

        !------------------------------------------------------- Output Variables

        type(ConsTransp) :: CT (-4:mnx6, -4:mny6, -4:mnz6)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= ini, fin
           DO jj= 1, ny+1
              DO ii= 1, nx

                 boundsCT(ii,jj,kk-ini+1)%fbzy     = CT(ii,jj,kk)%fbzy
                 boundsCT(ii,jj,kk-ini+1)%rflxy     = CT(ii,jj,kk)%rflxy
                 
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMCTTOBOUNDZT


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

