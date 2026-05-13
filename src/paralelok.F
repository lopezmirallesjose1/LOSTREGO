
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

     
SUBROUTINE INTERCAMBIABOUNDFIELDK (physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy)

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
            type(TboundField) :: boundsXs1k(2,nz), boundsYs1k(2,nz), boundsYs1kR(nz)
            type(TboundField) :: boundsXr2k(2,nz), boundsYr2k(2,nz), boundsYr2kR(nz)
            type(TboundField) :: boundsXs2k(2,nz), boundsYs2k(2,nz), boundsYs2kR(nz)
            type(TboundField) :: boundsXr1k(2,nz), boundsYr1k(2,nz), boundsYr1kR(nz)            

            !--------------------------------------------------------------------------- BEGIN
      

            ! Send / Receive counter.
            req_cnt = 0
 
# define bptagxleft 62345
# define bptagxright 62346
# define bptagyleft 62347
# define bptagyright 62348
# define bptagzleft 62349 
# define bptagzright 62350

            
            ! ---------    
            !    MEMORY
            ! --------- 

            ! Non-periodic boundary (defined in modulos.F).
            !call GETBOUNDSMEMORYSPACE (nx, ny, nz)
 

        
               !   ================================
               !     EXCHANGE RIGHT SIDE (Node i)
               !   ================================ 
           
               !-------
               ! X-dir
               !-------
                    
               ! Last node doesn't change information from the right side.
               IF (mpi_coords(1) .NE. (mpi_dims(1)-1)) THEN
              
                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = 2*1*nz*6
         
                  call KFROMFIELDTOBOUNDX (nx, nx, ny, nz, mnx5, mny5, mnz5, boundsXs1k, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsXs1k, lon, MPI_DOUBLE_PRECISION, xright, bptagxright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)
        
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &                  
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsXr1k, lon, MPI_DOUBLE_PRECISION, xright,                 &
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
                  lon = 2*1*nz*6

                  call KFROMFIELDTOBOUNDY (ny, ny, nx, nz, mnx5, mny5, mnz5, boundsYs1k, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsYs1k, lon, MPI_DOUBLE_PRECISION, yright, bptagyright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF


                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsYr1k, lon, MPI_DOUBLE_PRECISION, yright,                 &
                                  bptagyleft, CART_WORLD,  arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF




	       !call PAR_BARRIER()
               !-------
               ! Y-dir (L)
               !-------
                    
               ! Last node doesn't change information from the right side.
               !IF (mpi_coords(2) .NE. (mpi_dims(2)-1) .AND. mpi_coords(1) .EQ. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                !  lon = 1*1*nz*6

                !  call KFROMFIELDTOBOUNDY (ny, ny, 0, nz, mnx5, mny5, mnz5, boundsYs1kL, physic)
          
                !  req_cnt = req_cnt + 1
                !  call MPI_ISEND (boundsYs1kL, lon, MPI_DOUBLE_PRECISION, yright,                 &
                !                  bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

                !  IF (ierr .NE. MPI_SUCCESS) THEN
                !     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                !                 'Error in MPI_Send', __FILE__, __LINE__)
                !  ENDIF
          
                !  req_cnt = req_cnt + 1
                !  call MPI_IRECV (boundsYr1kL, lon, MPI_DOUBLE_PRECISION, yright,                 &
                !                  bptagyleft, CART_WORLD, arr_request(req_cnt), ierr)

                !  IF (ierr .NE. MPI_SUCCESS) THEN
                !     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      & 
                !                 'Error in MPI_Recv', __FILE__, __LINE__)
                !  ENDIF
          
                !ENDIF

               !call PAR_BARRIER()

               !-------
               ! Y-dir (R)
               !-------

               ! Last node doesn't change information from the right side.
               !IF (mpi_coords(2) .NE. (mpi_dims(2)-1) .AND. mpi_coords(1) .EQ. mpi_dims(1)-1) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
               !   lon = 1*1*nz*6

               !   call KFROMFIELDTOBOUNDY (ny, ny, nx+1, nz, mnx5, mny5, mnz5, boundsYs1kR, physic)

               !   req_cnt = req_cnt + 1
               !   call MPI_ISEND (boundsYs1kR, lon, MPI_DOUBLE_PRECISION, yright,                 &
               !                   bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Send', __FILE__, __LINE__)
               !   ENDIF

               !   req_cnt = req_cnt + 1
               !   call MPI_IRECV (boundsYr1kR, lon, MPI_DOUBLE_PRECISION, yright,                 &
               !                   bptagyleft, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Recv', __FILE__, __LINE__)
               !   ENDIF

               ! ENDIF


	!	call PAR_BARRIER()


               !   ================================
               !     EXCHANGE LEFT SIDE (Node i)
               !   ================================


               IF (mpi_coords(1) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = 2*1*nz*6

                  call KFROMFIELDTOBOUNDX (1, 1, ny, nz, mnx5, mny5, mnz5, boundsXs2k, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsXs2k, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                                  bptagxleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsXr2k, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                                  bptagxright, CART_WORLD, arr_request(req_cnt), ierr)


                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF



               IF (mpi_coords(2) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = 2*1*nz*6

                  call KFROMFIELDTOBOUNDY (1, 1, nx, nz, mnx5, mny5, mnz5, boundsYs2k, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsYs2k, lon, MPI_DOUBLE_PRECISION, yleft,                  &
                                  bptagyleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsYr2k, lon, MPI_DOUBLE_PRECISION, yleft,                  &
                                  bptagyright, CART_WORLD, arr_request(req_cnt), ierr)


                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF




               !call PAR_BARRIER()

               !-------
               ! Y-dir (L)
               !-------

             !  IF (mpi_coords(2) .GT. 0 .AND. mpi_coords(1) .EQ. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
             !     lon = 1*1*nz*6

             !     call KFROMFIELDTOBOUNDY (1, 1, 0, nz, mnx5, mny5, mnz5, boundsYs2kL, physic)

             !     req_cnt = req_cnt + 1
             !     call MPI_ISEND (boundsYs2kL, lon, MPI_DOUBLE_PRECISION, yleft,                  &
             !                     bptagyleft, CART_WORLD, arr_request(req_cnt),ierr )

             !     IF (ierr .NE. MPI_SUCCESS) THEN
             !        call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
             !                    'Error in MPI_Send', __FILE__, __LINE__)
             !     ENDIF

             !     req_cnt = req_cnt + 1
             !     call MPI_IRECV (boundsYr2kL, lon, MPI_DOUBLE_PRECISION, yleft,                  &
             !                     bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

             !     IF (ierr .NE. MPI_SUCCESS) THEN
             !        call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
             !                    'Error in MPI_Recv', __FILE__, __LINE__)
             !     ENDIF

             !  ENDIF

             !  call PAR_BARRIER()


               !-------
               ! Y-dir (R)
               !-------


              ! IF (mpi_coords(2) .GT. 0 .AND. mpi_coords(1) .EQ. mpi_dims(1)-1) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
              !    lon = 1*1*nz*6

              !    call KFROMFIELDTOBOUNDY (1, 1, nx+1, nz, mnx5, mny5, mnz5, boundsYs2kR, physic)

              !    req_cnt = req_cnt + 1
              !    call MPI_ISEND (boundsYs2kR, lon, MPI_DOUBLE_PRECISION, yleft,                  &
              !                    bptagyleft, CART_WORLD, arr_request(req_cnt),ierr )

              !    IF (ierr .NE. MPI_SUCCESS) THEN
              !       call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
              !                   'Error in MPI_Send', __FILE__, __LINE__)
              !    ENDIF

              !    req_cnt = req_cnt + 1
              !    call MPI_IRECV (boundsYr2kR, lon, MPI_DOUBLE_PRECISION, yleft,                  &
              !                    bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

              !    IF (ierr .NE. MPI_SUCCESS) THEN
              !       call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
              !                   'Error in MPI_Recv', __FILE__, __LINE__)
              !    ENDIF

              ! ENDIF



               ! Wait for all the communications to finish.
               call MPI_WAITALL (req_cnt, arr_request, arr_status, ierr)


	       !call PAR_BARRIER()

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',                                         &
                              'Error in MPI_WAITALL', __FILE__, __LINE__)
               ENDIF



               IF (mpi_coords(1) .NE. (mpi_dims(1)-1)) THEN
                  call KFROMBOUNDTOFIELDX (nx+1, nx+1, ny, nz, mnx5, mny5, mnz5, boundsXr1k, physic)
               ENDIF

               IF (mpi_coords(2) .NE. (mpi_dims(2)-1)) THEN
                  call KFROMBOUNDTOFIELDY (ny+1, ny+1, nx, nz, mnx5, mny5, mnz5, boundsYr1k, physic)
               ENDIF


               !call PAR_BARRIER()


               IF (mpi_coords(1) .GT. 0) THEN
                  call KFROMBOUNDTOFIELDX (0, 0, ny, nz, mnx5, mny5, mnz5, boundsXr2k, physic)
               ENDIF

               IF (mpi_coords(2) .GT. 0) THEN
                  call KFROMBOUNDTOFIELDY (0, 0, nx, nz, mnx5, mny5, mnz5, boundsYr2k, physic)
               ENDIF

               !call PAR_BARRIER()


	       ! LEFT Y

               !IF (mpi_coords(2) .NE. mpi_dims(2)-1 .AND. mpi_coords(1) .EQ. 0) THEN
               !   call KFROMBOUNDTOFIELDY (ny+1, ny+1, 0, nz, mnx5, mny5, mnz5, boundsYr1kL, physic)
               !ENDIF

               !call PAR_BARRIER()


               !IF (mpi_coords(2) .GT. 0 .AND. mpi_coords(1) .EQ. 0) THEN
               !   call KFROMBOUNDTOFIELDY (0, 0, 0, nz, mnx5, mny5, mnz5, boundsYr2kL, physic)
               !ENDIF

               !call PAR_BARRIER()


               ! RIGHT Y

               !IF (mpi_coords(2) .NE. mpi_dims(2)-1 .AND. mpi_coords(1) .EQ. mpi_dims(1)-1) THEN
               !   call KFROMBOUNDTOFIELDY (ny+1, ny+1, nx+1, nz, mnx5, mny5, mnz5, boundsYr1kR, physic)
               !ENDIF

               !call PAR_BARRIER()


               !IF (mpi_coords(2) .GT. 0 .AND. mpi_coords(1) .EQ. mpi_dims(1)-1) THEN
               !   call KFROMBOUNDTOFIELDY (0, 0, nx+1, nz, mnx5, mny5, mnz5, boundsYr2kR, physic)
               !ENDIF


            ! ---------    
            !    MEMORY
            ! --------- 


            ! Non-periodic boundary.

            !call FREEBOUNDSMEMORYSPACE ()

# endif

        !--------------------------------------------------------------------------- END

END SUBROUTINE INTERCAMBIABOUNDFIELDK



     


      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE KFROMBOUNDTOFIELDX (ini, fin, ny, nz, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, ny, nz, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,nz)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
                 
              physic(ini,0,kk)%velx    = bounds(1,kk)%velx
              physic(ini,0,kk)%vely    = bounds(1,kk)%vely
              physic(ini,0,kk)%velz    = bounds(1,kk)%velz
              physic(ini,0,kk)%bxc2     = bounds(1,kk)%bxc
              physic(ini,0,kk)%byc2     = bounds(1,kk)%byc
              physic(ini,0,kk)%bzc2     = bounds(1,kk)%bzc
              
              physic(ini,ny+1,kk)%velx    = bounds(2,kk)%velx
              physic(ini,ny+1,kk)%vely    = bounds(2,kk)%vely
              physic(ini,ny+1,kk)%velz    = bounds(2,kk)%velz
              physic(ini,ny+1,kk)%bxc2     = bounds(2,kk)%bxc
              physic(ini,ny+1,kk)%byc2     = bounds(2,kk)%byc
              physic(ini,ny+1,kk)%bzc2     = bounds(2,kk)%bzc

        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE KFROMBOUNDTOFIELDX

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Send y-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE KFROMBOUNDTOFIELDY (ini, fin, nx, nz, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, nx, nz, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,nz)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
                 
              physic(0,ini,kk)%velx    = bounds(1,kk)%velx
              physic(0,ini,kk)%vely    = bounds(1,kk)%vely
              physic(0,ini,kk)%velz    = bounds(1,kk)%velz
              physic(0,ini,kk)%bxc2     = bounds(1,kk)%bxc
              physic(0,ini,kk)%byc2     = bounds(1,kk)%byc
              physic(0,ini,kk)%bzc2     = bounds(1,kk)%bzc

              physic(nx+1,ini,kk)%velx    = bounds(2,kk)%velx
              physic(nx+1,ini,kk)%vely    = bounds(2,kk)%vely
              physic(nx+1,ini,kk)%velz    = bounds(2,kk)%velz
              physic(nx+1,ini,kk)%bxc2     = bounds(2,kk)%bxc
              physic(nx+1,ini,kk)%byc2     = bounds(2,kk)%byc
              physic(nx+1,ini,kk)%bzc2     = bounds(2,kk)%bzc
     
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE KFROMBOUNDTOFIELDY


      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE KFROMFIELDTOBOUNDX (ini, fin, ny, nz, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, ny, nz, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,nz)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
           
              bounds(1,kk)%velx    = physic(ini,0,kk)%velx    
              bounds(1,kk)%vely    = physic(ini,0,kk)%vely    
              bounds(1,kk)%velz    = physic(ini,0,kk)%velz
              bounds(1,kk)%bxc     = physic(ini,0,kk)%bxc2
              bounds(1,kk)%byc     = physic(ini,0,kk)%byc2
              bounds(1,kk)%bzc     = physic(ini,0,kk)%bzc2


              bounds(2,kk)%velx    = physic(ini,ny+1,kk)%velx
              bounds(2,kk)%vely    = physic(ini,ny+1,kk)%vely
              bounds(2,kk)%velz    = physic(ini,ny+1,kk)%velz
              bounds(2,kk)%bxc     = physic(ini,ny+1,kk)%bxc2
              bounds(2,kk)%byc     = physic(ini,ny+1,kk)%byc2
              bounds(2,kk)%bzc     = physic(ini,ny+1,kk)%bzc2



        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE KFROMFIELDTOBOUNDX


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save y-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE KFROMFIELDTOBOUNDY (ini, fin, nx, nz, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, nx, nz, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,nz)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
                 
              bounds(1,kk)%velx    = physic(0,ini,kk)%velx
              bounds(1,kk)%vely    = physic(0,ini,kk)%vely
              bounds(1,kk)%velz    = physic(0,ini,kk)%velz
              bounds(1,kk)%bxc     = physic(0,ini,kk)%bxc2
              bounds(1,kk)%byc     = physic(0,ini,kk)%byc2
              bounds(1,kk)%bzc     = physic(0,ini,kk)%bzc2


              bounds(2,kk)%velx    = physic(nx+1,ini,kk)%velx
              bounds(2,kk)%vely    = physic(nx+1,ini,kk)%vely
              bounds(2,kk)%velz    = physic(nx+1,ini,kk)%velz
              bounds(2,kk)%bxc     = physic(nx+1,ini,kk)%bxc2
              bounds(2,kk)%byc     = physic(nx+1,ini,kk)%byc2
              bounds(2,kk)%bzc     = physic(nx+1,ini,kk)%bzc2
      
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE KFROMFIELDTOBOUNDY


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

