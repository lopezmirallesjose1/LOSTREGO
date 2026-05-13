
!==============================================
!==============================================

!       HYBRID OMP + MPI PARALLELIZATION

!==============================================
!==============================================

      
!===============================================================================
! NAME              :: INTERCAMBIABOUNDFIELD
! GLOBAL VARIABLES  :: mpi_coords, mpi_dims, xright, xleft, yright, yleft, zright, zleft.
! DESCRIPTION       :: Send and Receive directives for the interchange of boundary physical
!                      variables between adjacent nodes.
!===============================================================================

     
SUBROUTINE INTERCAMBIABOUNDFIELDI (physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy)

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

        integer(is) :: nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
 
        !------------------------------------------------------- Local Variables

# ifdef PARALELO

            integer(is) :: lon, ierr, req_cnt 
            integer(is) :: status (MPI_STATUS_SIZE)
            integer(is) :: arr_request (12)
      
            !type(TboundPhysic), pointer :: boundsX (:,:,:), boundsY (:,:,:), boundsZ (:,:,:)    
            integer ::  arr_status (MPI_STATUS_SIZE, 12)


            ! TboundPhysic arrays used in intercambiaBoundPhysic
            type(TboundField) :: boundsYs1i(2,nx), boundsZs1i(2,nx), boundsZs1iR(nx)
            type(TboundField) :: boundsYr2i(2,nx), boundsZr2i(2,nx), boundsZr2iR(nx)
            type(TboundField) :: boundsYr1i(2,nx), boundsZs2i(2,nx), boundsZs2iR(nx)
            type(TboundField) :: boundsYs2i(2,nx), boundsZr1i(2,nx), boundsZr1iR(nx)

            
            !--------------------------------------------------------------------------- BEGIN
      

            ! Send / Receive counter.
            req_cnt = 0
 
# define bptagxleft 72345
# define bptagxright 72346
# define bptagyleft 72347
# define bptagyright 72348
# define bptagzleft 72349 
# define bptagzright 72350

            
            ! ---------    
            !    MEMORY
            ! --------- 

            ! Non-periodic boundary (defined in modulos.F).
            !call GETBOUNDSMEMORYSPACE (nx, ny, nz)
 

        
               !   ================================
               !     EXCHANGE RIGHT SIDE (Node i)
               !   ================================ 
           
               !-------
               ! Y-dir
               !-------
                    
               ! Last node doesn't change information from the right side.
               IF (mpi_coords(2) .NE. (mpi_dims(2)-1)) THEN
              
                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = nx*2*1*6
         
                  call IFROMFIELDTOBOUNDY (ny, ny, nz, nx, mnx5, mny5, mnz5, boundsYs1i, physic)
                   
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsYs1i, lon, MPI_DOUBLE_PRECISION, yright, bptagyright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)
        
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &                  
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsYr1i, lon, MPI_DOUBLE_PRECISION, yright,                 &
                                  bptagyleft, CART_WORLD,  arr_request(req_cnt), ierr)
              
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
                  lon = nx*2*1*6

                  call IFROMFIELDTOBOUNDZ (nz, nz, ny, nx, mnx5, mny5, mnz5, boundsZs1i, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsZs1i, lon, MPI_DOUBLE_PRECISION, zright, bptagzright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsZr1i, lon, MPI_DOUBLE_PRECISION, zright,                 &
                                  bptagzleft, CART_WORLD,  arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF



               !call PAR_BARRIER()

               !-------
               ! Z-dir (L)
               !-------
                    
               ! Last node doesn't change information from the right side.
               !IF (mpi_coords(3) .NE. (mpi_dims(3)-1) .AND. mpi_coords(2) .EQ. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
               !   lon = nx*1*1*6

               !   call IFROMFIELDTOBOUNDZ (nz, nz, 0, nx, mnx5, mny5, mnz5, boundsZs1iL, physic)
          
               !   req_cnt = req_cnt + 1
               !   call MPI_ISEND (boundsZs1iL, lon, MPI_DOUBLE_PRECISION, zright,               &
               !                   bptagzright, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Send', __FILE__, __LINE__)
               !   ENDIF
          
               !   req_cnt = req_cnt + 1
               !   call MPI_IRECV (boundsZr1iL, lon, MPI_DOUBLE_PRECISION, zright,                 &
               !                   bptagzleft, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      & 
               !                  'Error in MPI_Recv', __FILE__, __LINE__)
               !   ENDIF
          
               !ENDIF

               !call PAR_BARRIER()

               !-------
               ! Z-dir (R)
               !-------

               ! Last node doesn't change information from the right side.
               !IF (mpi_coords(3) .NE. (mpi_dims(3)-1) .AND. mpi_coords(2) .EQ. mpi_dims(2)-1) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
               !   lon = nx*1*1*6

               !   call IFROMFIELDTOBOUNDZ (nz, nz, ny+1, nx, mnx5, mny5, mnz5, boundsZs1iR, physic)

               !   req_cnt = req_cnt + 1
               !   call MPI_ISEND (boundsZs1iR, lon, MPI_DOUBLE_PRECISION, zright,               &
               !                   bptagzright, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Send', __FILE__, __LINE__)
               !   ENDIF

               !   req_cnt = req_cnt + 1
               !   call MPI_IRECV (boundsZr1iR, lon, MPI_DOUBLE_PRECISION, zright,                 &
               !                   bptagzleft, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Recv', __FILE__, __LINE__)
               !   ENDIF

               !ENDIF


	       !call PAR_BARRIER()


               !   ================================
               !     EXCHANGE LEFT SIDE (Node i)
               !   ================================


               IF (mpi_coords(2) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = nx*2*1*6

                  call IFROMFIELDTOBOUNDY (1, 1, nz, nx, mnx5, mny5, mnz5, boundsYs2i, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsYs2i, lon, MPI_DOUBLE_PRECISION, yleft,                  &
                                  bptagyleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsYr2i, lon, MPI_DOUBLE_PRECISION, yleft,                  &
                                  bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF


               IF (mpi_coords(3) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = nx*2*1*6

                  call IFROMFIELDTOBOUNDZ (1, 1, ny, nx, mnx5, mny5, mnz5, boundsZs2i, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsZs2i, lon, MPI_DOUBLE_PRECISION, zleft,                  &
                                  bptagzleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsZr2i, lon, MPI_DOUBLE_PRECISION, zleft,                  &
                                  bptagzright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF


               !call PAR_BARRIER()

               !-------
               ! Z-dir (L)
               !-------


               !IF (mpi_coords(3) .GT. 0 .AND. mpi_coords(2) .EQ. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
               !   lon = nx*1*1*6

               !   call IFROMFIELDTOBOUNDZ (1, 1, 0, nx, mnx5, mny5, mnz5, boundsZs2iL, physic)

               !   req_cnt = req_cnt + 1
               !   call MPI_ISEND (boundsZs2iL, lon, MPI_DOUBLE_PRECISION, zleft,                  &
               !                   bptagzleft, CART_WORLD, arr_request(req_cnt),ierr )

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Send', __FILE__, __LINE__)
               !   ENDIF

               !   req_cnt = req_cnt + 1
               !   call MPI_IRECV (boundsZr2iL, lon, MPI_DOUBLE_PRECISION, zleft,                  &
               !                   bptagzright, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Recv', __FILE__, __LINE__)
               !   ENDIF

               !ENDIF
               !call PAR_BARRIER()

               
               !-------
               ! Z-dir (R)
               !-------


               !IF (mpi_coords(3) .GT. 0 .AND. mpi_coords(2) .EQ. mpi_dims(2)-1) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
               !   lon = nx*1*1*6

               !   call IFROMFIELDTOBOUNDZ (1, 1, ny+1, nx, mnx5, mny5, mnz5, boundsZs2iR, physic)

               !   req_cnt = req_cnt + 1
               !   call MPI_ISEND (boundsZs2iR, lon, MPI_DOUBLE_PRECISION, zleft,                  &
               !                   bptagzleft, CART_WORLD, arr_request(req_cnt),ierr )

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Send', __FILE__, __LINE__)
               !   ENDIF

               !   req_cnt = req_cnt + 1
               !   call MPI_IRECV (boundsZr2iR, lon, MPI_DOUBLE_PRECISION, zleft,                  &
               !                   bptagzright, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Recv', __FILE__, __LINE__)
               !   ENDIF

               !ENDIF             




               ! Wait for all the communications to finish.
               call MPI_WAITALL (req_cnt, arr_request, arr_status, ierr)

	       !call PAR_BARRIER()

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',                                         &
                              'Error in MPI_WAITALL', __FILE__, __LINE__)
               ENDIF




               IF (mpi_coords(2) .NE. (mpi_dims(2)-1)) THEN
                  call IFROMBOUNDTOFIELDY (ny+1, ny+1, nz, nx, mnx5, mny5, mnz5, boundsYr1i, physic)
               ENDIF


               IF (mpi_coords(3) .NE. (mpi_dims(3)-1)) THEN
                  call IFROMBOUNDTOFIELDZ (nz+1, nz+1, ny, nx, mnx5, mny5, mnz5, boundsZr1i, physic)
               ENDIF


               !call PAR_BARRIER()


               IF (mpi_coords(2) .GT. 0) THEN
                  call IFROMBOUNDTOFIELDY (0, 0, nz, nx, mnx5, mny5, mnz5, boundsYr2i, physic)
               ENDIF

               IF (mpi_coords(3) .GT. 0) THEN
                  call IFROMBOUNDTOFIELDZ (0, 0, ny, nx, mnx5, mny5, mnz5, boundsZr2i, physic)
               ENDIF


               !call PAR_BARRIER()


               ! LEFT Z

               !IF (mpi_coords(3) .NE. mpi_dims(3)-1 .AND. mpi_coords(2) .EQ. 0) THEN
               !   call IFROMBOUNDTOFIELDZ (nz+1, nz+1, 0, nx, mnx5, mny5, mnz5, boundsZr1iL, physic)
               !ENDIF

               !call PAR_BARRIER()


               !IF (mpi_coords(3) .GT. 0 .AND. mpi_coords(2) .EQ. 0) THEN
               !   call IFROMBOUNDTOFIELDZ (0, 0, 0, nx, mnx5, mny5, mnz5, boundsZr2iL, physic)
               !ENDIF

               !call PAR_BARRIER()


               ! RIGHT Z

               !IF (mpi_coords(3) .NE. mpi_dims(3)-1 .AND. mpi_coords(2) .EQ. mpi_dims(2)-1) THEN
               !   call IFROMBOUNDTOFIELDZ (nz+1, nz+1, ny+1, nx, mnx5, mny5, mnz5, boundsZr1iR, physic)
               !ENDIF

               !call PAR_BARRIER()


               !IF (mpi_coords(3) .GT. 0 .AND. mpi_coords(2) .EQ. mpi_dims(2)-1) THEN
               !   call IFROMBOUNDTOFIELDZ (0, 0, ny+1, nx, mnx5, mny5, mnz5, boundsZr2iR, physic)
               !ENDIF


            ! ---------    
            !    MEMORY
            ! --------- 


            ! Non-periodic boundary.

            !call FREEBOUNDSMEMORYSPACE ()

# endif

        !--------------------------------------------------------------------------- END

END SUBROUTINE INTERCAMBIABOUNDFIELDI



     


      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE IFROMBOUNDTOFIELDY (ini, fin, nz, nx, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, nz, nx, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,nx)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO ii= 1, nx
                 
              physic(ii,ini,0)%velx     = bounds(1,ii)%velx
              physic(ii,ini,0)%vely     = bounds(1,ii)%vely
              physic(ii,ini,0)%velz     = bounds(1,ii)%velz
              physic(ii,ini,0)%bxc2     = bounds(1,ii)%bxc
              physic(ii,ini,0)%byc2     = bounds(1,ii)%byc
              physic(ii,ini,0)%bzc2     = bounds(1,ii)%bzc
              

              physic(ii,ini,nz+1)%velx     = bounds(2,ii)%velx
              physic(ii,ini,nz+1)%vely     = bounds(2,ii)%vely
              physic(ii,ini,nz+1)%velz     = bounds(2,ii)%velz
              physic(ii,ini,nz+1)%bxc2     = bounds(2,ii)%bxc
              physic(ii,ini,nz+1)%byc2     = bounds(2,ii)%byc
              physic(ii,ini,nz+1)%bzc2     = bounds(2,ii)%bzc



        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE IFROMBOUNDTOFIELDY

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Send y-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE IFROMBOUNDTOFIELDZ (ini, fin, ny, nx, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, ny, nx, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,nx)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO ii= 1, nx
                 
              physic(ii,0,ini)%velx     = bounds(1,ii)%velx
              physic(ii,0,ini)%vely     = bounds(1,ii)%vely
              physic(ii,0,ini)%velz     = bounds(1,ii)%velz
              physic(ii,0,ini)%bxc2     = bounds(1,ii)%bxc
              physic(ii,0,ini)%byc2     = bounds(1,ii)%byc
              physic(ii,0,ini)%bzc2     = bounds(1,ii)%bzc


              physic(ii,ny+1,ini)%velx     = bounds(2,ii)%velx
              physic(ii,ny+1,ini)%vely     = bounds(2,ii)%vely
              physic(ii,ny+1,ini)%velz     = bounds(2,ii)%velz
              physic(ii,ny+1,ini)%bxc2     = bounds(2,ii)%bxc
              physic(ii,ny+1,ini)%byc2     = bounds(2,ii)%byc
              physic(ii,ny+1,ini)%bzc2     = bounds(2,ii)%bzc
     
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE IFROMBOUNDTOFIELDZ


      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE IFROMFIELDTOBOUNDY (ini, fin, nz, nx, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, nz, nx, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,nx)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO ii= 1, nx
           
              bounds(1,ii)%velx    = physic(ii,ini,0)%velx    
              bounds(1,ii)%vely    = physic(ii,ini,0)%vely    
              bounds(1,ii)%velz    = physic(ii,ini,0)%velz
              bounds(1,ii)%bxc     = physic(ii,ini,0)%bxc2
              bounds(1,ii)%byc     = physic(ii,ini,0)%byc2
              bounds(1,ii)%bzc     = physic(ii,ini,0)%bzc2

              bounds(2,ii)%velx    = physic(ii,ini,nz+1)%velx
              bounds(2,ii)%vely    = physic(ii,ini,nz+1)%vely
              bounds(2,ii)%velz    = physic(ii,ini,nz+1)%velz
              bounds(2,ii)%bxc     = physic(ii,ini,nz+1)%bxc2
              bounds(2,ii)%byc     = physic(ii,ini,nz+1)%byc2
              bounds(2,ii)%bzc     = physic(ii,ini,nz+1)%bzc2


        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE IFROMFIELDTOBOUNDY


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save y-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE IFROMFIELDTOBOUNDZ (ini, fin, ny, nx, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, ny, nx, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,nx)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO ii= 1, nx
                 
              bounds(1,ii)%velx    = physic(ii,0,ini)%velx
              bounds(1,ii)%vely    = physic(ii,0,ini)%vely
              bounds(1,ii)%velz    = physic(ii,0,ini)%velz
              bounds(1,ii)%bxc     = physic(ii,0,ini)%bxc2
              bounds(1,ii)%byc     = physic(ii,0,ini)%byc2
              bounds(1,ii)%bzc     = physic(ii,0,ini)%bzc2

              bounds(2,ii)%velx    = physic(ii,ny+1,ini)%velx
              bounds(2,ii)%vely    = physic(ii,ny+1,ini)%vely
              bounds(2,ii)%velz    = physic(ii,ny+1,ini)%velz
              bounds(2,ii)%bxc     = physic(ii,ny+1,ini)%bxc2
              bounds(2,ii)%byc     = physic(ii,ny+1,ini)%byc2
              bounds(2,ii)%bzc     = physic(ii,ny+1,ini)%bzc2
      
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE IFROMFIELDTOBOUNDZ


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

