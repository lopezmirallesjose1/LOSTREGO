
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

     
SUBROUTINE INTERCAMBIABOUNDFIELDJ (physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy)

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


            type(TboundField) :: boundsZs1j(2,ny), boundsXs1j(2,ny), boundsXs1jR(ny)
            type(TboundField) :: boundsZr2j(2,ny), boundsXr2j(2,ny), boundsXr2jR(ny)
            type(TboundField) :: boundsZr1j(2,ny), boundsXs2j(2,ny), boundsXs2jR(ny)
            type(TboundField) :: boundsZs2j(2,ny), boundsXr1j(2,ny), boundsXr1jR(ny)



            
            !--------------------------------------------------------------------------- BEGIN
      

            ! Send / Receive counter.
            req_cnt = 0
 
# define bptagxleft 82345
# define bptagxright 82346
# define bptagyleft 82347
# define bptagyright 82348
# define bptagzleft 82349 
# define bptagzright 82350

            
            ! ---------    
            !    MEMORY
            ! --------- 

            ! Non-periodic boundary (defined in modulos.F).
            !call GETBOUNDSMEMORYSPACE (nx, ny, nz)
 

        
               !   ================================
               !     EXCHANGE RIGHT SIDE (Node i)
               !   ================================ 
           
               !-------
               ! Z-dir
               !-------
                    
               ! Last node doesn't change information from the right side.
               IF (mpi_coords(3) .NE. (mpi_dims(3)-1)) THEN
              
                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = 2*ny*1*6
         
                  call JFROMFIELDTOBOUNDZ (nz, nz, nx, ny, mnx5, mny5, mnz5, boundsZs1j, physic)
                   
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsZs1j, lon, MPI_DOUBLE_PRECISION, zright, bptagzright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)
        
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &                  
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsZr1j, lon, MPI_DOUBLE_PRECISION, zright,                 &
                                  bptagzleft, CART_WORLD,  arr_request(req_cnt), ierr)
              
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF        
              
               ENDIF



               !-------
               ! X-dir
               !-------

               ! Last node doesn't change information from the right side.
               IF (mpi_coords(1) .NE. (mpi_dims(1)-1)) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = ny*2*1*6

                  call JFROMFIELDTOBOUNDX (nx, nx, nz, ny, mnx5, mny5, mnz5, boundsXs1j, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsXs1j, lon, MPI_DOUBLE_PRECISION, xright, bptagxright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsXr1j, lon, MPI_DOUBLE_PRECISION, xright,                 &
                                  bptagxleft, CART_WORLD,  arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF










               !call PAR_BARRIER()


               !-------
               ! X-dir (L)
               !-------
                    
               ! Last node doesn't change information from the right side.
               !IF (mpi_coords(1) .NE. (mpi_dims(1)-1) .AND. mpi_coords(3) .EQ. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
               !   lon = 1*ny*1*6

               !   call JFROMFIELDTOBOUNDX (nx, nx, 0, ny, mnx5, mny5, mnz5, boundsXs1jL, physic)
          
               !   req_cnt = req_cnt + 1
               !   call MPI_ISEND (boundsXs1jL, lon, MPI_DOUBLE_PRECISION, xright,                 &
               !                   bptagxright, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Send', __FILE__, __LINE__)
               !   ENDIF
          
               !   req_cnt = req_cnt + 1
               !   call MPI_IRECV (boundsXr1jL, lon, MPI_DOUBLE_PRECISION, xright,                 &
               !                   bptagxleft, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      & 
               !                  'Error in MPI_Recv', __FILE__, __LINE__)
               !   ENDIF
          
               !ENDIF


               !call PAR_BARRIER()


               !-------
               ! X-dir (R)
               !-------

               ! Last node doesn't change information from the right side.
               !IF (mpi_coords(1) .NE. (mpi_dims(1)-1) .AND. mpi_coords(3) .EQ. mpi_dims(3)-1) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
               !   lon = 1*ny*1*6

               !   call JFROMFIELDTOBOUNDX (nx, nx, nz+1, ny, mnx5, mny5, mnz5, boundsXs1jR, physic)

               !   req_cnt = req_cnt + 1
               !   call MPI_ISEND (boundsXs1jR, lon, MPI_DOUBLE_PRECISION, xright,                 &
               !                   bptagxright, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Send', __FILE__, __LINE__)
               !   ENDIF

               !   req_cnt = req_cnt + 1
               !   call MPI_IRECV (boundsXr1jR, lon, MPI_DOUBLE_PRECISION, xright,                 &
               !                   bptagxleft, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Recv', __FILE__, __LINE__)
               !   ENDIF

               !ENDIF


		!call PAR_BARRIER()


               !   ================================
               !     EXCHANGE LEFT SIDE (Node i)
               !   ================================


               IF (mpi_coords(3) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = 2*ny*1*6

                  call JFROMFIELDTOBOUNDZ (1, 1, nx, ny, mnx5, mny5, mnz5, boundsZs2j, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsZs2j, lon, MPI_DOUBLE_PRECISION, zleft,                  &
                                  bptagzleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsZr2j, lon, MPI_DOUBLE_PRECISION, zleft,                  &
                                  bptagzright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF


               !call PAR_BARRIER()



               IF (mpi_coords(1) .GT. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = ny*2*1*6

                  call JFROMFIELDTOBOUNDX (1, 1, nz, ny, mnx5, mny5, mnz5, boundsXs2j, physic)

                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsXs2j, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                                  bptagxleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF

                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsXr2j, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                                  bptagxright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Recv', __FILE__, __LINE__)
                  ENDIF

               ENDIF




               !-------
               ! X-dir (L)
               !-------

               !IF (mpi_coords(1) .GT. 0 .AND. mpi_coords(3) .EQ. 0) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
               !   lon = 1*ny*1*6

               !   call JFROMFIELDTOBOUNDX (1, 1, 0, ny, mnx5, mny5, mnz5, boundsXs2jL, physic)

               !   req_cnt = req_cnt + 1
               !   call MPI_ISEND (boundsXs2jL, lon, MPI_DOUBLE_PRECISION, xleft,                  &
               !                   bptagxleft, CART_WORLD, arr_request(req_cnt),ierr )

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Send', __FILE__, __LINE__)
               !   ENDIF

               !   req_cnt = req_cnt + 1
               !   call MPI_IRECV (boundsXr2jL, lon, MPI_DOUBLE_PRECISION, xleft,                  &
               !                   bptagxright, CART_WORLD, arr_request(req_cnt), ierr)

               !   IF (ierr .NE. MPI_SUCCESS) THEN
               !      call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
               !                  'Error in MPI_Recv', __FILE__, __LINE__)
               !   ENDIF

               !ENDIF


               !call PAR_BARRIER()

               
               !-------
               ! X-dir (R)
               !-------

               !IF (mpi_coords(1) .GT. 0 .AND. mpi_coords(3) .EQ. mpi_dims(3)-1) THEN

                  ! Number of boundary cells (tridimensional) * number of variables.
                !  lon = 1*ny*1*6

                !  call JFROMFIELDTOBOUNDX (1, 1, nz+1, ny, mnx5, mny5, mnz5, boundsXs2jR, physic)

                !  req_cnt = req_cnt + 1
                !  call MPI_ISEND (boundsXs2jR, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                !                  bptagxleft, CART_WORLD, arr_request(req_cnt),ierr )

                !  IF (ierr .NE. MPI_SUCCESS) THEN
                !     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                !                 'Error in MPI_Send', __FILE__, __LINE__)
                !  ENDIF

                !  req_cnt = req_cnt + 1
                !  call MPI_IRECV (boundsXr2jR, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                !                  bptagxright, CART_WORLD, arr_request(req_cnt), ierr)

                !  IF (ierr .NE. MPI_SUCCESS) THEN
                !     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                !                 'Error in MPI_Recv', __FILE__, __LINE__)
                !  ENDIF

               !ENDIF             

               !call PAR_BARRIER()




               ! Wait for all the communications to finish.
               call MPI_WAITALL (req_cnt, arr_request, arr_status, ierr)

	 	!call PAR_BARRIER()

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',                                         &
                              'Error in MPI_WAITALL', __FILE__, __LINE__)
               ENDIF


               IF (mpi_coords(3) .NE. (mpi_dims(3)-1)) THEN
                  call JFROMBOUNDTOFIELDZ (nz+1, nz+1, nx, ny, mnx5, mny5, mnz5, boundsZr1j, physic)
               ENDIF

               IF (mpi_coords(1) .NE. (mpi_dims(1)-1)) THEN
                  call JFROMBOUNDTOFIELDX (nx+1, nx+1, nz, ny, mnx5, mny5, mnz5, boundsXr1j, physic)
               ENDIF



               !call PAR_BARRIER()


               IF (mpi_coords(3) .GT. 0) THEN
                  call JFROMBOUNDTOFIELDZ (0, 0, nx, ny, mnx5, mny5, mnz5, boundsZr2j, physic)
               ENDIF

               IF (mpi_coords(1) .GT. 0) THEN
                  call JFROMBOUNDTOFIELDX (0, 0, nz, ny, mnx5, mny5, mnz5, boundsXr2j, physic)
               ENDIF


               !call PAR_BARRIER()


               ! LEFT X

               !IF (mpi_coords(1) .NE. mpi_dims(1)-1 .AND. mpi_coords(3) .EQ. 0) THEN
               !   call JFROMBOUNDTOFIELDX (nx+1, nx+1, 0, ny, mnx5, mny5, mnz5, boundsXr1jL, physic)
               !ENDIF

               !call PAR_BARRIER()


               !IF (mpi_coords(1) .GT. 0 .AND. mpi_coords(3) .EQ. 0) THEN
               !   call JFROMBOUNDTOFIELDX (0, 0, 0, ny, mnx5, mny5, mnz5, boundsXr2jL, physic)
               !ENDIF

               !call PAR_BARRIER()


               ! RIGHT X

               !IF (mpi_coords(1) .NE. mpi_dims(1)-1 .AND. mpi_coords(3) .EQ. mpi_dims(3)-1) THEN
               !   call JFROMBOUNDTOFIELDX (nx+1, nx+1, nz+1, ny, mnx5, mny5, mnz5, boundsXr1jR, physic)
               !ENDIF

               !call PAR_BARRIER()


               !IF (mpi_coords(1) .GT. 0 .AND. mpi_coords(3) .EQ. mpi_dims(3)-1) THEN
               !   call JFROMBOUNDTOFIELDX (0, 0, nz+1, ny, mnx5, mny5, mnz5, boundsXr2jR, physic)
               !ENDIF



            ! ---------    
            !    MEMORY
            ! --------- 


            ! Non-periodic boundary.

            !call FREEBOUNDSMEMORYSPACE ()

# endif

        !--------------------------------------------------------------------------- END

END SUBROUTINE INTERCAMBIABOUNDFIELDJ



     


      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE JFROMBOUNDTOFIELDZ (ini, fin, nx, ny, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, nx, ny, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,ny)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO jj= 1, ny
                 
              physic(0,jj,ini)%velx     = bounds(1,jj)%velx
              physic(0,jj,ini)%vely     = bounds(1,jj)%vely
              physic(0,jj,ini)%velz     = bounds(1,jj)%velz
              physic(0,jj,ini)%bxc2     = bounds(1,jj)%bxc
              physic(0,jj,ini)%byc2     = bounds(1,jj)%byc
              physic(0,jj,ini)%bzc2     = bounds(1,jj)%bzc
              
              physic(nx+1,jj,ini)%velx     = bounds(2,jj)%velx
              physic(nx+1,jj,ini)%vely     = bounds(2,jj)%vely
              physic(nx+1,jj,ini)%velz     = bounds(2,jj)%velz
              physic(nx+1,jj,ini)%bxc2     = bounds(2,jj)%bxc
              physic(nx+1,jj,ini)%byc2     = bounds(2,jj)%byc
              physic(nx+1,jj,ini)%bzc2     = bounds(2,jj)%bzc


        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE JFROMBOUNDTOFIELDZ

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Send y-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE JFROMBOUNDTOFIELDX (ini, fin, nz, ny, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, nz, ny, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,ny)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO jj= 1, ny
                 
              physic(ini,jj,0)%velx     = bounds(1,jj)%velx
              physic(ini,jj,0)%vely     = bounds(1,jj)%vely
              physic(ini,jj,0)%velz     = bounds(1,jj)%velz
              physic(ini,jj,0)%bxc2     = bounds(1,jj)%bxc
              physic(ini,jj,0)%byc2     = bounds(1,jj)%byc
              physic(ini,jj,0)%bzc2     = bounds(1,jj)%bzc


              physic(ini,jj,nz+1)%velx     = bounds(2,jj)%velx
              physic(ini,jj,nz+1)%vely     = bounds(2,jj)%vely
              physic(ini,jj,nz+1)%velz     = bounds(2,jj)%velz
              physic(ini,jj,nz+1)%bxc2     = bounds(2,jj)%bxc
              physic(ini,jj,nz+1)%byc2     = bounds(2,jj)%byc
              physic(ini,jj,nz+1)%bzc2     = bounds(2,jj)%bzc
     
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE JFROMBOUNDTOFIELDX


      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE JFROMFIELDTOBOUNDZ (ini, fin, nx, ny, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, nx, ny, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,ny)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO jj= 1, ny
           
              bounds(1,jj)%velx    = physic(0,jj,ini)%velx    
              bounds(1,jj)%vely    = physic(0,jj,ini)%vely    
              bounds(1,jj)%velz    = physic(0,jj,ini)%velz
              bounds(1,jj)%bxc     = physic(0,jj,ini)%bxc2
              bounds(1,jj)%byc     = physic(0,jj,ini)%byc2
              bounds(1,jj)%bzc     = physic(0,jj,ini)%bzc2

              bounds(2,jj)%velx    = physic(nx+1,jj,ini)%velx
              bounds(2,jj)%vely    = physic(nx+1,jj,ini)%vely
              bounds(2,jj)%velz    = physic(nx+1,jj,ini)%velz
              bounds(2,jj)%bxc     = physic(nx+1,jj,ini)%bxc2
              bounds(2,jj)%byc     = physic(nx+1,jj,ini)%byc2
              bounds(2,jj)%bzc     = physic(nx+1,jj,ini)%bzc2


        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE JFROMFIELDTOBOUNDZ


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save y-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE JFROMFIELDTOBOUNDX (ini, fin, nz, ny, mnx5, mny5, mnz5, bounds, physic)

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

        integer(is) :: ini, fin, nz, ny, mnx5, mny5, mnz5
        type(TboundField) :: bounds (2,ny)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO jj= 1, ny
                 

              bounds(1,jj)%velx    = physic(ini,jj,0)%velx
              bounds(1,jj)%vely    = physic(ini,jj,0)%vely
              bounds(1,jj)%velz    = physic(ini,jj,0)%velz
              bounds(1,jj)%bxc     = physic(ini,jj,0)%bxc2
              bounds(1,jj)%byc     = physic(ini,jj,0)%byc2
              bounds(1,jj)%bzc     = physic(ini,jj,0)%bzc2

              bounds(2,jj)%velx    = physic(ini,jj,nz+1)%velx
              bounds(2,jj)%vely    = physic(ini,jj,nz+1)%vely
              bounds(2,jj)%velz    = physic(ini,jj,nz+1)%velz
              bounds(2,jj)%bxc     = physic(ini,jj,nz+1)%bxc2
              bounds(2,jj)%byc     = physic(ini,jj,nz+1)%byc2
              bounds(2,jj)%bzc     = physic(ini,jj,nz+1)%bzc2


        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE JFROMFIELDTOBOUNDX


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

