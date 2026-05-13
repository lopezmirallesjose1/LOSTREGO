
!==============================================
!==============================================

!       HYBRID OMP + MPI PARALLELIZATION

!==============================================
!==============================================

!===============================================================================
! NAME              :: OMP_INIT
! GLOBAL VARIABLES  :: omp_th
! DESCRIPTION       :: Set the number of threads in omp_th from OMP_NUM_TRHEADS.
!===============================================================================

!===============================================================================
SUBROUTINE OMP_INIT
!===============================================================================

        USE workarrays
        USE parallel
        USE memoria
        
        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        ! ...

        !------------------------------------------------------- Output Variables

        ! ...
   
        !------------------------------------------------------- Local Variables

        integer :: ierr 
 
        !--------------------------------------------------------------------------- BEGIN

# ifdef OPENMP

            INTERFACE

               integer FUNCTION OMP_GET_NUM_THREADS ()
                       END FUNCTION OMP_GET_NUM_THREADS

               integer FUNCTION OMP_GET_THREAD_NUM ()
                       END FUNCTION OMP_GET_THREAD_NUM

            END INTERFACE

            ! Get number of OpenMP threads.
          
!$OMP PARALLEL
            omp_th = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL   
   
# else
            omp_th = 1
# endif

        allocate (ocurx(0:omp_th-1), ocury(0:omp_th-1), ocurz(0:omp_th-1), stat=ierr)
           
        IF (ierr .NE. 0) THEN
           call ERROR ('OMP_INIT', 'Allocation of OCURX, OCURY, OCURZ', __FILE__, __LINE__)
        ENDIF

        call MEMPUSH (omp_th, is, 'ocurx')
        call MEMPUSH (omp_th, is, 'ocury')
        call MEMPUSH (omp_th, is, 'ocurz')

        return
      
        !--------------------------------------------------------------------------- END

!===============================================================================      
END SUBROUTINE OMP_INIT
!===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

           
!===============================================================================
! NAME              :: DECOMPOSE_MPI_DOMAIN
! GLOBAL VARIABLES  :: omp_th
! DESCRIPTION       :: Adapt the dimensions of the grid according to the CPU topology.
!===============================================================================

!=============================================================================== 
SUBROUTINE DECOMPOSE_MPI_DOMAIN (xCPUs, yCPUs, zCPUs, nx, ny, nz)
!=============================================================================== 

        USE parallel
        USE memoria
        USE workarrays

        implicit none
# include "types.h"

        !------------------------------------------------------- Input Variables

        integer :: xCPUs, yCPUs, zCPUs
        integer :: nx, ny, nz

        !------------------------------------------------------- Output Variables

        ! ...
   
        !------------------------------------------------------- Local Variables

        character(40) :: message1

        !--------------------------------------------------------------------------- BEGIN

      
        IF ((MOD (nx,xCPUs) .EQ. 0) .AND. (MOD (ny,yCPUs) .EQ. 0) .AND.  & 
            (MOD (nz,zCPUs) .EQ. 0) ) THEN

           nx = nx / xCPUs
           ny = ny / yCPUs
           nz = nz / zCPUs
        
        ELSE
 
           message1 = 'nx,ny or nz is not divisible by x,y,zCPU'
           call ERROR ('DECOMPOSE_MPI_DOMAIN', message1, __FILE__, __LINE__)

        ENDIF

        return

        !--------------------------------------------------------------------------- END

!=============================================================================== 
END SUBROUTINE DECOMPOSE_MPI_DOMAIN
!=============================================================================== 

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


!===============================================================================
! NAME              :: MY_OMP_GET_THREAD_NUM
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Obtain the thread number.
!===============================================================================

!===============================================================================
SUBROUTINE MY_OMP_GET_THREAD_NUM (curth)
!===============================================================================

        USE workarrays
        USE parallel
        USE memoria
        
        implicit none      
# include "types.h"

        !------------------------------------------------------- Input Variables

        ! ...

        !------------------------------------------------------- Output Variables

        integer(is) :: curth
   
        !------------------------------------------------------- Local Variables

        ! ...

        !--------------------------------------------------------------------------- BEGIN

# ifdef OPENMP

            INTERFACE

               integer FUNCTION OMP_GET_NUM_THREADS()
                       END FUNCTION OMP_GET_NUM_THREADS

               integer FUNCTION OMP_GET_THREAD_NUM()
                       END FUNCTION OMP_GET_THREAD_NUM

            END INTERFACE

            ! Get number of OpenMP thread
            curth = OMP_GET_THREAD_NUM()
  
# else
            curth = 0  
# endif
      
         return

         !--------------------------------------------------------------------------- END

!===============================================================================
END SUBROUTINE MY_OMP_GET_THREAD_NUM
!===============================================================================

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: MY_OMP_GET_NUM_THREADS
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Obtain the total number of threads.
!===============================================================================

!===============================================================================
SUBROUTINE MY_OMP_GET_NUM_THREADS (nth)
!===============================================================================

        USE workarrays
        USE parallel
        USE memoria
        
        implicit none      
# include "types.h"

        !------------------------------------------------------- Input Variables

        ! ...

        !------------------------------------------------------- Output Variables

        integer(is) :: nth
   
        !------------------------------------------------------- Local Variables

        ! ...

        !--------------------------------------------------------------------------- BEGIN

# ifdef OPENMP

            INTERFACE

               integer FUNCTION OMP_GET_NUM_THREADS()
                       END FUNCTION OMP_GET_NUM_THREADS

               integer FUNCTION OMP_GET_THREAD_NUM()
                       END FUNCTION OMP_GET_THREAD_NUM

            END INTERFACE

            ! Get total number of OMP threads.
            nth = OMP_GET_NUM_THREADS()
  
# else
            nth = 0  
# endif
      
        return

        !--------------------------------------------------------------------------- END

!===============================================================================
END SUBROUTINE MY_OMP_GET_NUM_THREADS
!===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: MPI_INIT3D
! GLOBAL VARIABLES  :: periods, xleft, xright, yleft, yright, zleft, zright,   
!                      iniTramo, finTramo, sizeTramo, mpi_dims
! DESCRIPTION       :: Construct and initialize the MPI communicator (with topology).
!===============================================================================

!===============================================================================
SUBROUTINE MPI_INIT3D (xCPUs, yCPUs, zCPUs, maxDims, max_dims)
!===============================================================================

        USE parallel
        USE memoria
        USE workarrays

	implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif
        
        
# include "types.h"


        !------------------------------------------------------- Input Variables

        integer :: xCPUs, yCPUs, zCPUs, maxDims

        !------------------------------------------------------- Output Variables

        integer(is) :: max_dims

        !------------------------------------------------------- Local Variables

        integer (is)   :: ierr
        character(41)  :: message1

        !--------------------------------------------------------------------------- BEGIN
 
# ifdef PARALELO
      
            ! Initialize MPI.
            call MPI_INIT (ierr)
            IF (ierr .NE. MPI_SUCCESS) THEN
               call ERROR ('MPI_INIT3D', 'Error in MPI_init', __FILE__, __LINE__)
            ENDIF
      
            ! Wait until all processes in the communicator (MPI_COMM_WORLD) 
            ! have reached this routine.
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)

            ! Specify whether the grid is periodic (true) or not (false) in each dimension.
            periods (1:maxDims) = .FALSE.

            ! Number of processors per dimension.
            mpi_dims(1) = xCPUs
            mpi_dims(2) = yCPUs
            mpi_dims(3) = zCPUs
     
            ! Make a new communicator attaching topology information.
            ! Output :: CART_WORLD.
            call MPI_CART_CREATE (MPI_COMM_WORLD, maxDims, mpi_dims(1:maxDims),     &
                                  periods(1:maxDims), .TRUE., CART_WORLD, ierr)
     
            IF (ierr .NE. MPI_SUCCESS) THEN
               call ERROR ('MPI_INIT3D', 'Error in MPI_CART_CREATE', __FILE__, __LINE__)   
            ENDIF

            ! ===========
            ! Get Cartesian Rank and Coordinates
            ! ===========      

            ! Output :: cartrank.
            call MPI_COMM_RANK (CART_WORLD, cartrank, ierr)
            IF (ierr .NE. MPI_SUCCESS) THEN
               call ERROR ('MPI_INIT3D', 'Error in MPI_COMM_RANK', __FILE__, __LINE__)
            ENDIF

            ! Output :: mpi_coords.
            call MPI_CART_COORDS (CART_WORLD, cartrank, maxDims, mpi_coords, ierr)
            IF (ierr .NE. MPI_SUCCESS) THEN
               call ERROR( 'MPI_INIT3D', 'Error in MPI_CART_COORDS', __FILE__, __LINE__)
            ENDIF
      
            ! One-dimensional MPI parallelization.
            IF (maxDims .EQ. 1) THEN
               mpi_coords(2:3) = 0
            ENDIF

            ! Two-dimensional MPI parallelization.
            IF (maxDims .EQ. 2) THEN
               mpi_coords(3) = 0
            ENDIF

            ! Get position of the processor's neighborhood.
            IF (xCPUs .GT. 1) THEN
         
               call MPI_CART_SHIFT (CART_WORLD, 0, 1, xleft, xright, ierr)
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('MPI_INIT3D', 'Error in MPI_CART_SHIFT', __FILE__, __LINE__)  
               ENDIF

            ELSE 
               
               xleft=0
               xright=0

            ENDIF

            IF (yCPUs .GT. 1) THEN

               call MPI_CART_SHIFT (CART_WORLD, 1, 1, yleft, yright, ierr)
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('MPI_INIT3D', 'Error in MPI_CART_SHIFT', __FILE__, __LINE__) 
               ENDIF

            ELSE 
 
               yleft=0
               yright=0

            ENDIF
      
            IF (zCPUs .GT. 1) THEN
               
               call MPI_CART_SHIFT (CART_WORLD, 2, 1, zleft, zright, ierr)
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('MPI_INIT3D', 'Error in MPI_CART_SHIFT', __FILE__, __LINE__)   
               ENDIF

            ELSE

               zleft=0
               zright=0

            ENDIF
      
# else

            xleft=0
            xright=0
            yleft=0
            yright=0
            zleft=0
            zright=0
            mpi_coords(1:3) = 0
            mpi_dims(1:3) = 1

# endif

        ! maximum of mpi_dims().
        max_dims = MAX (mpi_dims(1), mpi_dims(2), mpi_dims(3))
  
        allocate (iniTramo(1:3,0:max_dims-1), finTramo(1:3,0:max_dims-1), &
                  sizeTramo(1:3,0:max_dims-1), stat= ierr)

        IF (ierr .NE. 0) THEN

           message1 = 'ALLOCATION: iniTramo, finTramo, sizeTramo'
           call ERROR ('mpi_init3D', message1, __FILE__, __LINE__)

        ENDIF
  
        call MEMPUSH (max_dims*3, is, 'iniTramo')
        call MEMPUSH (max_dims*3, is, 'finTramo')
        call MEMPUSH (max_dims*3, is, 'sizeTramo')

        !--------------------------------------------------------------------------- END

!===============================================================================
END SUBROUTINE MPI_INIT3D
!===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: PAR_END
! GLOBAL VARIABLES  :: iniTramo, finTramo, sizeTramo
! DESCRIPTION       :: Finalize MPI parallelization.
!===============================================================================

!===============================================================================
SUBROUTINE PAR_END ()
!===============================================================================

        USE parallel
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

        ! ...

        !------------------------------------------------------- Output Variables

        ! ...

        !------------------------------------------------------- Local Variables

# ifdef PARALELO
            integer (is) :: ierr
# endif

        !--------------------------------------------------------------------------- BEGIN

        deallocate (iniTramo, finTramo, sizeTramo)
        call MEMPOP ('iniTramo')
        call MEMPOP ('finTramo')
        call MEMPOP ('sizeTramo')

# ifdef PARALELO

            call MPI_FINALIZE (ierr)
            IF (ierr .NE. MPI_SUCCESS) THEN
               call ERROR ('PAR_END', 'Error in MPI_finalize', __FILE__, __LINE__)
            endif

# endif

        !--------------------------------------------------------------------------- END

!===============================================================================
END SUBROUTINE PAR_END
!===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: PAR_BARRIER
! GLOBAL VARIABLES  :: CART_WORLD
! DESCRIPTION       :: Wait until all processes in the communicator have reached this routine.
!===============================================================================

!===============================================================================
SUBROUTINE PAR_BARRIER ( )
!===============================================================================

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

        ! ...

        !------------------------------------------------------- Output Variables

        ! ...

        !------------------------------------------------------- Local Variables

# ifdef PARALELO
            integer(is) :: ierr
# endif            

        !--------------------------------------------------------------------------- BEGIN

# ifdef PARALELO
           
            call MPI_BARRIER (CART_WORLD, ierr)
            IF (ierr .NE. 0) THEN
               call ERROR ('PAR_BARRIER', 'Error in MPI_Barrier', __FILE__, __LINE__)
            ENDIF
            
# endif

        !--------------------------------------------------------------------------- END

!===============================================================================
END SUBROUTINE PAR_BARRIER
!===============================================================================

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: DISTRIBUTE_WORK
! GLOBAL VARIABLES  :: mpi_dims, iniTramo, finTramo, sizeTramo, sizeTotal
! DESCRIPTION       :: Establish the limits for each node.
!===============================================================================

!===============================================================================
SUBROUTINE DISTRIBUTE_WORK (max_dims, longsize, longsizeX, longsizeY, longsizeZ,     &       
                            mn1, mn5, mn6, mnx1, mnx5, mnx6,                         &
                            mny1, mny5, mny6, mnz1, mnz5, mnz6)
!===============================================================================

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

        integer(is) ::  longsizeX, longsizeY, longsizeZ, max_dims

        !------------------------------------------------------- Output Variables

        integer(is) :: longsize(3)

        integer(is) :: mn1, mn5, mn6, mnx1, mnx5, mnx6,           &
                                      mny1, mny5, mny6,           &
                                      mnz1, mnz5, mnz6


        !------------------------------------------------------- Local Variables

        integer(is) :: bsize, resto, i, j, pp

        integer(is) :: bsizearray(0:max_dims-1)

        character(43) :: message1

        !--------------------------------------------------------------------------- BEGIN

        
        ! Initialize longsize.
        longsize(1)=longsizeX
        longsize(2)=longsizeY
        longsize(3)=longsizeZ


        ! Work is distributed for each dimension.
        DO i = 1, 3

           IF (mpi_dims(i) .NE. 1) THEN
		
              bsize = longsize(i) / mpi_dims(i) ! Integer / Integer
              resto = longsize(i) - bsize*mpi_dims(i)

              IF (bsize .LT. 5) THEN

                 message1 = 'Size must be .GE. than boundary constraint.'
                 call ERROR ('DISTRIBUTE_WORK', message1, __FILE__, __LINE__)

              ENDIF

              ! RESTO defines the number of nodes where we can increase one cell.
              ! This is not currently working as long as DECOMPOSE_MPI_DOMAIN prevents RESTO
              ! to be not equal zero.

              IF (resto .GT. 0) THEN

                 DO j = 0, mpi_dims(i) - 1

                    IF (j .LT. resto) THEN
                       bsizearray (j) = bsize + 1
                    ELSE
                       bsizearray (j) = bsize
                    ENDIF

                 ENDDO

              ELSE

                 DO j = 0, mpi_dims(i) - 1
                    bsizearray (j) = bsize
                 ENDDO

              ENDIF


              ! First subinterval's coordinates for each direction.
              iniTramo(i,0)  = 1
              finTramo(i,0)  = bsizearray(0)
              sizeTramo(i,0) = bsizearray(0)

              DO pp= 1, mpi_dims(i)-1

                 iniTramo(i,pp)  = finTramo(i,pp-1) + 1
                 finTramo(i,pp)  = finTramo(i,pp-1) + bsizearray(pp)
                 sizeTramo(i,pp) = bsizearray(pp)

              ENDDO

              ! Update longsize in case RESTO is not equal zero.
              sizeTotal(i) = longsize(i)
              longsize(i)  = sizeTramo(i,mpi_coords(i))

           ENDIF
                
           ! Establish local limits for the grid.

           IF (i .EQ. 1) THEN

              mnx1 = longsize(i) + 1
              mnx5 = longsize(i) + 5
              mnx6 = longsize(i) + 6

           ELSE IF (i .EQ. 2) THEN

              mny1 = longsize(i) + 1
              mny5 = longsize(i) + 5
              mny6 = longsize(i) + 6

           ELSE IF (i .EQ. 3) THEN

              mnz1 = longsize(i) + 1
              mnz5 = longsize(i) + 5
              mnz6 = longsize(i) + 6

           ENDIF

        ENDDO

        ! Maximum values.
        mn1 = MAX (mnx1, mny1, mnz1)
        mn5 = MAX (mnx5, mny5, mnz5)
        mn6 = MAX (mnx6, mny6, mnz6)
	   		
        !--------------------------------------------------------------------------- END

!===============================================================================
END SUBROUTINE DISTRIBUTE_WORK
!===============================================================================


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: INTERCAMBIABOUNDPOS
! GLOBAL VARIABLES  :: mpi_coords, mpi_dims, xright, xleft, yright, yleft, zright, zleft.
! DESCRIPTION       :: 'Send' and 'Receive' directives for boundary positions' interchange
!                      between adjacent nodes.
!===============================================================================

SUBROUTINE INTERCAMBIABOUNDPOS (nx, xzn, xznl, xznr, ny, yzn, yznl, yznr,      &
                                nz, zzn, zznl, zznr, boundposX, boundposY,     &
                                boundposZ)

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

        integer(is) :: ny, nx, nz

        real(rs) :: xzn(*), xznl(*), xznr(*)
        real(rs) :: yzn(*), yznl(*), yznr(*)
        real(rs) :: zzn(*), zznl(*), zznr(*)

        !------------------------------------------------------- Output Variables

        type(TboundPosX) :: boundposX(-4:5), boundsXaux(-4:5)
        type(TboundPosY) :: boundposY(-4:5), boundsYaux(-4:5)
        type(TboundPosZ) :: boundposZ(-4:5), boundsZaux(-4:5)   


        !------------------------------------------------------- Local Variables

        integer(is) :: ii
        real(rs)    :: dx

        integer :: req_cnt

# ifdef PARALELO

            integer(is) ::  lon, ierr
            integer(is) ::  status(MPI_STATUS_SIZE)
            integer(is) ::  arr_request(12)     
            integer     ::  arr_status(MPI_STATUS_SIZE, 12)

# endif

        !--------------------------------------------------------------------------- BEGIN

# ifdef PARALELO
 
# define tagxleft 12345
# define tagxright 12346
# define tagyleft 12347
# define tagyright 12348
# define tagzleft 12349
# define tagzright 12350
    
            ! Number of boundary cells * number of variables.
            lon = 5*4

            ! Send and Receive counter.
            req_cnt = 0

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
            ! i-Node SENDS the left boundary positions to i+1-Node.
            ! i+1-Node SENDS the right boundary positions to i-Node.

            !-------
            ! X-dir
            !-------
 
            ! REMARK: The last node doesn't change information from the right side.
            IF (mpi_coords(1) .NE. mpi_dims(1) - 1) THEN

               DO ii= -4, 0    
    
                  boundsXaux(ii)%xzn  = xzn (nx+ii)
                  boundsXaux(ii)%xznl = xznl(nx+ii)
                  boundsXaux(ii)%xznr = xznr(nx+ii)
                  boundsXaux(ii)%dx   = xznl(nx+ii+1) - xznl(nx+ii)      
 
               ENDDO 

               req_cnt = req_cnt + 1
               call MPI_ISEND (boundsXaux(-4), lon, MPI_DOUBLE_PRECISION, xright,    &
                               tagxright, CART_WORLD,  arr_request(req_cnt), ierr)
       
    
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Send',            &
                              __FILE__, __LINE__)
               ENDIF

           
               req_cnt = req_cnt + 1
               call MPI_IRECV (boundposX(1), lon, MPI_DOUBLE_PRECISION, xright,      &
                               tagxleft, CART_WORLD, arr_request(req_cnt), ierr)
           

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Recv',            &
                              __FILE__, __LINE__)
               ENDIF
         
            ENDIF

            !-------
            ! Y-dir
            !-------
 
            ! REMARK: The last node doesn't change information from the right side.
            IF (mpi_coords(2) .NE. mpi_dims(2) - 1) THEN

               DO ii= -4, 0        

                  boundsYaux(ii)%yzn  = yzn(ny+ii)
                  boundsYaux(ii)%yznl = yznl(ny+ii)
                  boundsYaux(ii)%yznr = yznr(ny+ii)
                  boundsYaux(ii)%dx   = yznl(ny+ii+1) - yznl(ny+ii) 
      
               ENDDO 

               req_cnt = req_cnt + 1
               call MPI_ISEND (boundsYaux(-4), lon, MPI_DOUBLE_PRECISION, yright,    &
                               tagyright, CART_WORLD, arr_request(req_cnt), ierr)
           
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Send',            &
                              __FILE__, __LINE__)
               ENDIF
           
               req_cnt = req_cnt + 1
               call MPI_IRECV (boundposY(1), lon, MPI_DOUBLE_PRECISION, yright,      &
                               tagyleft, CART_WORLD, arr_request(req_cnt), ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Recv',            &
                               __FILE__, __LINE__)
               ENDIF
         
            ENDIF

            !-------
            ! Z-dir
            !-------

            ! REMARK: Last node doesn't change information from the right side.
            IF (mpi_coords(3) .NE. mpi_dims(3) - 1) THEN

               DO ii= -4, 0  
      
                  boundsZaux(ii)%zzn  = zzn(nz+ii)
                  boundsZaux(ii)%zznl = zznl(nz+ii)
                  boundsZaux(ii)%zznr = zznr(nz+ii)
                  boundsZaux(ii)%dx   = zznl(nz+ii+1) - zznl(nz+ii)   
    
               ENDDO

               req_cnt = req_cnt + 1
               call MPI_ISEND (boundsZaux(-4), lon, MPI_DOUBLE_PRECISION, zright,    &
                               tagzright, CART_WORLD, arr_request(req_cnt), ierr)
           
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Send',            &
                               __FILE__, __LINE__)
               ENDIF

               req_cnt = req_cnt + 1
               call MPI_IRECV (boundposZ(1), lon, MPI_DOUBLE_PRECISION, zright,      &
                               tagzleft, CART_WORLD, arr_request(req_cnt), ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Recv',            &
                               __FILE__, __LINE__)
               ENDIF       

            ENDIF


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
            ! i-Node RECEIVES the left boundary positions from i-1-Node.
            ! i-1-Node RECEIVES the right boundary positions from i-Node.


            !-------
            ! X-dir
            !-------
 
            ! REMARK: The first node doesn't change information from the left side.
            IF (mpi_coords(1) .NE. 0) THEN

               DO ii= 1, 5
                  
                  boundsXaux(ii)%xzn  = xzn(ii)
                  boundsXaux(ii)%xznl = xznl(ii)
                  boundsXaux(ii)%xznr = xznr(ii)
                  boundsXaux(ii)%dx   = xznl(ii+1) - xznl(ii)

               ENDDO

               req_cnt = req_cnt + 1
               call MPI_ISEND (boundsXaux(1), lon, MPI_DOUBLE_PRECISION, xleft,         &
                               tagxleft, CART_WORLD, arr_request(req_cnt), ierr)
           
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Send',               &
                               __FILE__, __LINE__)
               ENDIF
           
               req_cnt = req_cnt + 1 
               call MPI_IRECV (boundposX(-4), lon, MPI_DOUBLE_PRECISION, xleft,         &
                               tagxright, CART_WORLD, arr_request(req_cnt), ierr)
           
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Recv',               &
                               __FILE__, __LINE__)
               ENDIF 
           
            ENDIF
    

            !-------
            ! Y-dir
            !-------
 
            ! First node doesn't change information from the left side.
            IF (mpi_coords(2) .NE. 0) THEN

               DO ii= 1, 5

                  boundsYaux(ii)%yzn  = yzn(ii)
                  boundsYaux(ii)%yznl = yznl(ii)
                  boundsYaux(ii)%yznr = yznr(ii)
                  boundsYaux(ii)%dx   = yznl(ii+1) - yznl(ii)

               ENDDO

               req_cnt = req_cnt + 1
               call MPI_ISEND (boundsYaux(1), lon, MPI_DOUBLE_PRECISION, yleft,         &
                               tagyleft, CART_WORLD, arr_request(req_cnt), ierr)
           
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Send',               &
                               __FILE__, __LINE__)
               ENDIF
           
               req_cnt = req_cnt + 1
               call MPI_IRECV (boundposY(-4), lon, MPI_DOUBLE_PRECISION, yleft,         &
                               tagyright, CART_WORLD, arr_request(req_cnt), ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Recv',               &
                               __FILE__, __LINE__  )
               ENDIF      

            ENDIF
    
            !-------
            ! Z-dir
            !-------
 
            ! First node doesn't change information from the left side.
            IF (mpi_coords(3) .NE. 0) THEN

               DO ii= 1, 5

                  boundsZaux(ii)%zzn  = zzn(ii)
                  boundsZaux(ii)%zznl = zznl(ii)
                  boundsZaux(ii)%zznr = zznr(ii)
                  boundsZaux(ii)%dx   = zznl(ii+1) - zznl(ii)

               ENDDO

               req_cnt = req_cnt + 1
               call MPI_ISEND (boundsZaux(1), lon, MPI_DOUBLE_PRECISION, zleft,         &
                               tagzleft, CART_WORLD, arr_request(req_cnt), ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Send',               &
                              __FILE__, __LINE__  )
               ENDIF

               req_cnt = req_cnt + 1
               call MPI_IRECV (boundposZ(-4), lon, MPI_DOUBLE_PRECISION, zleft,         &
                               tagzright, CART_WORLD, arr_request(req_cnt), ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_Recv',               &
                               __FILE__, __LINE__)
               ENDIF
 
            ENDIF
 
            ! Wait for all communication to finish.
            call MPI_WAITALL(req_cnt, arr_request, arr_status, ierr)

            IF (ierr .NE. MPI_SUCCESS) THEN
               call ERROR ('INTERCAMBIABOUNDPOS', 'Error in MPI_WAITALL',               &
                           __FILE__, __LINE__)
            ENDIF
      
# endif

        !   ==========================
        !    FILL BOUNDARY CONDITIONS
        !   ==========================
           
        !   ============
        !    LEFT
        !   ============

        !-------
        ! X-dir
        !-------

        IF (mpi_coords(1) .EQ. 0) THEN

           dx = xznl(2) - xznl(1)
           boundposX(0)%xzn  = xzn(1)  - dx
           boundposX(0)%xznl = xznl(1) - dx
           boundposX(0)%xznr = xznr(1) - dx
           boundposX(0)%dx   = dx

           DO ii= -1, -4, -1

              boundposX(ii)%xzn  = boundposX(ii+1)%xzn  - dx
              boundposX(ii)%xznl = boundposX(ii+1)%xznl - dx
              boundposX(ii)%xznr = boundposX(ii+1)%xznr - dx
              boundposX(ii)%dx   = dx

           ENDDO

        ENDIF

        !-------
        ! Y-dir
        !-------

        IF (mpi_coords(2) .EQ. 0) THEN

           dx = yznl(2) - yznl(1)
           boundposY(0)%yzn  = yzn(1)  - dx
           boundposY(0)%yznl = yznl(1) - dx
           boundposY(0)%yznr = yznr(1) - dx
           boundposY(0)%dx   = dx

           DO ii= -1, -4, -1

              boundposY(ii)%yzn  = boundposY(ii+1)%yzn  - dx
              boundposY(ii)%yznl = boundposY(ii+1)%yznl - dx
              boundposY(ii)%yznr = boundposY(ii+1)%yznr - dx
              boundposY(ii)%dx   = dx

           ENDDO

        ENDIF

        !-------
        ! Z-dir
        !-------

        IF (mpi_coords(3) .EQ. 0) THEN
           
           dx = zznl(2) - zznl(1)
           boundposZ(0)%zzn  = zzn(1)  - dx
           boundposZ(0)%zznl = zznl(1) - dx
           boundposZ(0)%zznr = zznr(1) - dx
           boundposZ(0)%dx   = dx

           DO ii= -1, -4, -1
              
              boundposZ(ii)%zzn  = boundposZ(ii+1)%zzn  - dx
              boundposZ(ii)%zznl = boundposZ(ii+1)%zznl - dx
              boundposZ(ii)%zznr = boundposZ(ii+1)%zznr - dx
              boundposZ(ii)%dx   = dx

           ENDDO

        ENDIF


        !   ============
        !    RIGHT
        !   ============

        !-------
        ! X-dir
        !-------

        IF (mpi_coords(1) .EQ. mpi_dims(1) - 1) THEN

           dx = xznl(nx+1) - xznl(nx)
           boundposX(1)%xzn  = xzn(nx)  + dx
           boundposX(1)%xznl = xznl(nx) + dx
           boundposX(1)%xznr = xznr(nx) + dx
           boundposX(1)%dx   = dx

           DO ii= 2, 5

              boundposX(ii)%xzn  = boundposX(ii-1)%xzn + dx
              boundposX(ii)%xznl = boundposX(ii-1)%xznl + dx
              boundposX(ii)%xznr = boundposX(ii-1)%xznr + dx
              boundposX(ii)%dx   = dx

           ENDDO

        ENDIF

        !-------
        ! Y-dir
        !-------

        IF (mpi_coords(2) .EQ. mpi_dims(2) - 1) THEN

           dx = yznl(ny+1) - yznl(ny)
           boundposY(1)%yzn  = yzn(ny)  + dx
           boundposY(1)%yznl = yznl(ny) + dx
           boundposY(1)%yznr = yznr(ny) + dx
           boundposY(1)%dx   = dx

           DO ii= 2, 5

              boundposY(ii)%yzn  = boundposY(ii-1)%yzn + dx
              boundposY(ii)%yznl = boundposY(ii-1)%yznl + dx
              boundposY(ii)%yznr = boundposY(ii-1)%yznr + dx
              boundposY(ii)%dx   = dx

           ENDDO

        ENDIF

        !-------
        ! Z-dir
        !-------

        IF (mpi_coords(3) .EQ. mpi_dims(3) - 1) THEN

           dx = zznl(nz+1) - zznl(nz)
           boundposZ(1)%zzn  = zzn(nz)  + dx
           boundposZ(1)%zznl = zznl(nz) + dx
           boundposZ(1)%zznr = zznr(nz) + dx
           boundposZ(1)%dx   = dx

           DO ii= 2, 5

              boundposZ(ii)%zzn  = boundposZ(ii-1)%zzn + dx
              boundposZ(ii)%zznl = boundposZ(ii-1)%zznl + dx
              boundposZ(ii)%zznr = boundposZ(ii-1)%zznr + dx
              boundposZ(ii)%dx   = dx

           ENDDO

        ENDIF

        !--------------------------------------------------------------------------- END

END SUBROUTINE INTERCAMBIABOUNDPOS


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


!===============================================================================
! NAME              :: INTERCAMBIABOUNDPHYSIC
! GLOBAL VARIABLES  :: mpi_coords, mpi_dims, xright, xleft, yright, yleft, zright, zleft.
! DESCRIPTION       :: Send and Receive directives for the interchange of boundary physical
!                      variables between adjacent nodes.
!===============================================================================

     
SUBROUTINE INTERCAMBIABOUNDPHYSIC (physic, nx, ny, nz, mnx5, mny5, mnz5, bndmny, bndmxy)

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
      
            type(TboundPhysic), pointer :: boundsX (:,:,:), boundsY (:,:,:), boundsZ (:,:,:)    
            integer ::  arr_status (MPI_STATUS_SIZE, 12)

            
            !--------------------------------------------------------------------------- BEGIN
      

            ! Send / Receive counter.
            req_cnt = 0
 
# define bptagxleft 32345
# define bptagxright 32346
# define bptagyleft 32347
# define bptagyright 32348
# define bptagzleft 32349 
# define bptagzright 32350

            ! ---------    
            !    MEMORY
            ! --------- 

            ! Periodic boundary
            allocate (boundsX(5,ny,nz), stat= ierr)
            IF (ierr .NE. 0) THEN
               call ERROR ('INTERCAMBIABOUNDPHYSIC',                              &
                           'Allocation of bounds', __FILE__, __LINE__) 
            ENDIF

            call MEMPUSH (5*ny*nz, rs, 'boundsX')
 
            allocate (boundsY(nx,5,nz), stat= ierr)
            IF (ierr .NE. 0) THEN
               call ERROR ('intercambiaBoundPhysic',                              &
                           'Allocation of bounds', __FILE__, __LINE__) 
            ENDIF

            call MEMPUSH (nx*5*nz, rs, 'boundsY')
 
            allocate (boundsZ(nx,ny,5), stat= ierr) 
            IF (ierr .NE. 0) THEN
               call ERROR ('intercambiaBoundPhysic',                              &
                           'Allocation of bounds', __FILE__, __LINE__) 
            ENDIF

            call MEMPUSH (nx*ny*5, rs, 'boundsZ')

            ! Non-periodic boundary (defined in modulos.F).
            call GETBOUNDSMEMORYSPACE (nx, ny, nz)
 

            ! --------------------------
            ! PERIODIC BOUNDARY
            ! --------------------------

            ! For periodic boundaries, the last node has to communicate with the 
            ! first one.

            ! The jet propagates in the y-direction.
            
            IF (bndmny .EQ. 3 .AND. bndmxy .EQ. 3) THEN
               
               !   ================================
               !     EXCHANGE RIGHT SIDE (Node i)
               !   ================================ 

               !-------
               ! X-dir
               !-------
          
               ! Number of boundary cells (tridimensional) * number of variables.
               lon = 5*ny*nz*11

               call FROMPHYTOBOUNDX (nx-4, nx, ny, nz, mnx5, mny5, mnz5, boundsX, physic)
	         
               call MPI_SENDRECV_REPLACE (boundsX, lon, MPI_DOUBLE_PRECISION, xright, bptagxright,     &
                                          xleft, bptagxright, CART_WORLD, status, ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',                                                & 
                              'Error in MPI_SendRecv_Replace', __FILE__, __LINE__)
               ENDIF

               call FROMBOUNDTOPHYX (-4, 0, ny, nz, mnx5, mny5, mnz5, boundsX, physic)


               !-------
               ! Y-dir
               !-------

               ! Number of boundary cells (tridimensional) * number of variables.
               lon = nx*5*nz*11

               call FROMPHYTOBOUNDY (ny-4, ny, nx, nz, mnx5, mny5, mnz5, boundsY, physic)
		
               call MPI_SENDRECV_REPLACE (boundsY, lon, MPI_DOUBLE_PRECISION, yright, bptagyright,     &
                                          yleft, bptagyright, CART_WORLD, status, ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',                                                &
                              'Error in MPI_SendRecv_Replace', __FILE__, __LINE__)
               ENDIF

               call FROMBOUNDTOPHYY (-4, 0, nx, nz, mnx5, mny5, mnz5, boundsY, physic)


               !-------
               ! Z-dir
               !-------

               ! Number of boundary cells (tridimensional) * number of variables.
               lon = nx*ny*5*11

               call FROMPHYTOBOUNDZ (nz-4, nz, nx, ny, mnx5, mny5, mnz5, boundsZ, physic)
		
               call MPI_SENDRECV_REPLACE (boundsZ, lon, MPI_DOUBLE_PRECISION, zright, bptagzright,     &
                                          zleft, bptagzright, CART_WORLD, status, ierr )

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR( 'INTERCAMBIABOUNDPHYSIC',     & 
                              'Error in MPI_SendRecv_Replace', __FILE__, __LINE__)
               ENDIF

               call FROMBOUNDTOPHYZ (-4, 0, nx, ny, mnx5, mny5, mnz5, boundsZ, physic)


               !   ================================
               !     EXCHANGE LEFT SIDE (Node i)
               !   ================================ 


               !-------
               ! X-dir
               !-------

               ! Number of boundary cells (tridimensional) * number of variables.
               lon = 5*ny*nz*11

               call FROMPHYTOBOUNDX (1, 5, ny, nz, mnx5, mny5, mnz5, boundsX, physic)

               call MPI_SENDRECV_REPLACE (boundsX, lon, MPI_DOUBLE_PRECISION, xright, bptagxleft,      &
                                          xleft, bptagxleft, CART_WORLD, status, ierr)

               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',                                                &
                              'Error in MPI_SendRecv_Replace', __FILE__, __LINE__)
               ENDIF

               call FROMBOUNDTOPHYX (nx+1, nx+5, ny, nz, mnx5, mny5, mnz5, boundsX, physic)


               !-------
               ! Y-dir
               !-------

               ! Number of boundary cells (tridimensional) * number of variables.
               lon = nx*5*nz*11

               call FROMPHYTOBOUNDY (1, 5, nx, nz, mnx5, mny5, mnz5, boundsY, physic)

               call MPI_SENDRECV_REPLACE (boundsY, lon, MPI_DOUBLE_PRECISION, yright, bptagyleft,      &
                                          yleft, bptagyleft, CART_WORLD, status, ierr )
          
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',       &
                              'Error in MPI_SendRecv_Replace', __FILE__, __LINE__)
               ENDIF

               call FROMBOUNDTOPHYY (ny+1, ny+5, nx, nz, mnx5, mny5, mnz5, boundsY, physic)


               !-------
               ! Z-dir
               !-------

               ! Number of boundary cells (tridimensional) * number of variables. 
               lon = nx*ny*5*11

               call FROMPHYTOBOUNDZ (1, 5, nx, ny, mnx5, mny5, mnz5, boundsZ, physic)
 
               call MPI_SENDRECV_REPLACE (boundsZ, lon, MPI_DOUBLE_PRECISION, zright, bptagzleft,      &
                                          zleft, bptagzleft, CART_WORLD, status, ierr )
         
               IF (ierr .NE. MPI_SUCCESS) THEN
                  call ERROR ('INTERCAMBIABOUNDPHYSIC',          &
                              'Error in MPI_SendRecv_Replace', __FILE__, __LINE__)
               ENDIF

               call FROMBOUNDTOPHYZ (nz+1, nz+5, nx, ny, mnx5, mny5, mnz5, boundsZ, physic)


            ! --------------------------
            ! NON-PERIODIC BOUNDARY
            ! --------------------------

            ELSE
    
       
               !   ================================
               !     EXCHANGE RIGHT SIDE (Node i)
               !   ================================ 
           
               !-------
               ! X-dir
               !-------
                    
               ! Last node doesn't change information from the right side.
               IF (mpi_coords(1) .NE. (mpi_dims(1)-1)) THEN
              
                  ! Number of boundary cells (tridimensional) * number of variables.
                  lon = 5*ny*nz*11
         
                  call FROMPHYTOBOUNDX (nx-4, nx, ny, nz, mnx5, mny5, mnz5, boundsXs1, physic)
                   
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsXs1, lon, MPI_DOUBLE_PRECISION, xright, bptagxright,    &
                                  CART_WORLD, arr_request(req_cnt), ierr)
        
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &                  
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsXr1, lon, MPI_DOUBLE_PRECISION, xright,                 &
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
                  lon = nx*5*nz*11

                  call FROMPHYTOBOUNDY (ny-4, ny, nx, nz, mnx5, mny5, mnz5, boundsYs1, physic)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsYs1, lon, MPI_DOUBLE_PRECISION, yright,                 &
                                  bptagyright, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
          
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsYr1, lon, MPI_DOUBLE_PRECISION, yright,                 &
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
                  lon = nx*ny*5*11

                  call FROMPHYTOBOUNDZ (nz-4, nz, nx, ny, mnx5, mny5, mnz5, boundsZs1, physic)
         
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsZs1, lon, MPI_DOUBLE_PRECISION, zright,                 &
                                  bptagzright, CART_WORLD, arr_request(req_cnt), ierr )

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ( 'INTERCAMBIABOUNDPHYSIC',                                     &
                                  'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsZr1, lon, MPI_DOUBLE_PRECISION, zright,                 &
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
                  lon = 5*ny*nz*11

                  call FROMPHYTOBOUNDX (1, 5, ny, nz, mnx5, mny5, mnz5, boundsXs2, physic)
         
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsXs2, lon, MPI_DOUBLE_PRECISION, xleft,                  &
                                  bptagxleft, CART_WORLD, arr_request(req_cnt), ierr)
 
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsXr2, lon, MPI_DOUBLE_PRECISION, xleft,                  &
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
                  lon = nx*5*nz*11

                  call FROMPHYTOBOUNDY (1, 5, nx, nz, mnx5, mny5, mnz5, boundsYs2, physic)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsYs2, lon, MPI_DOUBLE_PRECISION, yleft,                  &
                                  bptagyleft, CART_WORLD, arr_request(req_cnt),ierr )
                  
                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
                  
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsYr2, lon, MPI_DOUBLE_PRECISION, yleft,                  &
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
                  lon = nx*ny*5*11

                  call FROMPHYTOBOUNDZ (1, 5, nx, ny, mnx5, mny5, mnz5, boundsZs2, physic)
          
                  req_cnt = req_cnt + 1
                  call MPI_ISEND (boundsZs2, lon, MPI_DOUBLE_PRECISION, zleft,                  &
                                  bptagzleft, CART_WORLD, arr_request(req_cnt), ierr)

                  IF (ierr .NE. MPI_SUCCESS) THEN
                     call ERROR ('INTERCAMBIABOUNDPHYSIC',                                      &
                                 'Error in MPI_Send', __FILE__, __LINE__)
                  ENDIF
         
                  req_cnt = req_cnt + 1
                  call MPI_IRECV (boundsZr2, lon, MPI_DOUBLE_PRECISION, zleft,                  &
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
                  call FROMBOUNDTOPHYX (nx+1, nx+5, ny, nz, mnx5, mny5, mnz5, boundsXr1, physic) 
               ENDIF

               IF (mpi_coords(2) .NE. (mpi_dims(2)-1)) THEN
                  call FROMBOUNDTOPHYY (ny+1, ny+5, nx, nz, mnx5, mny5, mnz5, boundsYr1, physic) 
               ENDIF

               IF (mpi_coords(3) .NE. (mpi_dims(3)-1)) THEN
                  call FROMBOUNDTOPHYZ (nz+1, nz+5, nx, ny, mnx5, mny5, mnz5, boundsZr1, physic) 
               ENDIF

               IF (mpi_coords(1) .GT. 0) THEN
                  call FROMBOUNDTOPHYX (-4, 0, ny, nz, mnx5, mny5, mnz5, boundsXr2, physic)
               ENDIF

               IF (mpi_coords(2) .GT. 0) THEN
                  call FROMBOUNDTOPHYY (-4, 0, nx, nz, mnx5, mny5, mnz5, boundsYr2, physic)
               ENDIF

               IF (mpi_coords(3) .GT. 0) THEN
                  call FROMBOUNDTOPHYZ (-4, 0, nx, ny, mnx5, mny5, mnz5, boundsZr2, physic)
               ENDIF

            ENDIF


            ! ---------    
            !    MEMORY
            ! --------- 

            ! Periodic boundary.

            deallocate (boundsX)
            call MEMPOP ('boundsX')
            deallocate (boundsY)
            call MEMPOP ('boundsY')
            deallocate (boundsZ)
            call MEMPOP ('boundsZ')

            ! Non-periodic boundary.

            call FREEBOUNDSMEMORYSPACE ()

# endif

        !--------------------------------------------------------------------------- END

END SUBROUTINE INTERCAMBIABOUNDPHYSIC


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


!===============================================================================
! NAME              :: INTERCAMBIAWRITEPHYSIC
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: 
!===============================================================================
      
SUBROUTINE INTERCAMBIAWRITEPHYSIC (nx, ny, nz, mnx1, mny1, mnz1,    &
                                   physic, plane0, max_dims)


        USE tipos
        USE parallel
        USE workarrays

	implicit none

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif

       
# include "types.h"

        !------------------------------------------ Input Variables

        integer(is) :: nx, ny, nz, mnx1, mny1, mnz1, max_dims
        
        type(Tphysic) :: physic (0:mnx1, 0:mny1, 0:mnz1) 

        !------------------------------------------ Output Variables

        type(TplanePhysic) :: plane0 (nx,ny)

        !------------------------------------------ Local Variables

# ifdef PARALELO

            integer(is) :: ierr, ii, jj, kk
                 
            integer(is) :: memz_proc, memz(max_dims), offsetz(max_dims)

            type(TPlanePhysic), pointer :: plane (:,:)

            character(23) :: message1, message2, message3

            
            !--------------------------------------------------------------------------- BEGIN

            ! Total size (nz == 1).
            memz_proc = nx*ny*nz*11

            allocate (plane(nx,ny), stat= ierr)
            IF (ierr .NE. 0) THEN

               message1='Allocation of 2D plane.'
               call ERROR ('INTERCAMBIAWRITEPHYSIC', message1, __FILE__, __LINE__)
               
            ENDIF

            ! Sending the Z-plane data (3D --> 2D).
            DO kk= 1, nz
               DO jj= 1, ny
                  DO ii= 1, nx
                     
                     plane(ii,jj)%velx    = physic(ii,jj,kk)%velx
                     plane(ii,jj)%vely    = physic(ii,jj,kk)%vely
                     plane(ii,jj)%velz    = physic(ii,jj,kk)%velz
                     plane(ii,jj)%densty  = physic(ii,jj,kk)%densty
                     plane(ii,jj)%denstye = physic(ii,jj,kk)%denstye
                     plane(ii,jj)%pres    = physic(ii,jj,kk)%pres
                     plane(ii,jj)%eps     = physic(ii,jj,kk)%eps
                     plane(ii,jj)%bxc     = physic(ii,jj,kk)%bxc
                     plane(ii,jj)%byc     = physic(ii,jj,kk)%byc
                     plane(ii,jj)%bzc     = physic(ii,jj,kk)%bzc
                     plane(ii,jj)%tracer  = physic(ii,jj,kk)%tracer

                  ENDDO
               ENDDO
            ENDDO

            ! Gather the information from different MPI-processes and save it in destination (rank=0).
            call MPI_GATHER (memz_proc, 1, MPI_INTEGER, memz, 1,  &
                             MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      
            IF (ierr .NE. MPI_SUCCESS) THEN

               message2='Error gathering -memz-.'                  
               call ERROR ('INTERCAMBIAWRITEPHYSIC', message2, __FILE__, __LINE__)
            
            ENDIF

            ! Specify the displacement relative to recvbuf at which to place the incoming
            ! data from process i.

            IF (cartrank .EQ. 0) THEN

               offsetz(1) = 0
               DO ii= 2, max_dims-1
                  offsetz(ii) = offsetz(ii-1) + memz(ii-1)
               ENDDO
               
            ENDIF
            
            call MPI_GATHERV (plane, memz_proc, MPI_DOUBLE_PRECISION, plane0, memz,   &
                              offsetz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            
            IF (ierr .NE. MPI_SUCCESS) THEN

               message3='Error gathering -plane-'
               call ERROR ('INTERCAMBISWRITEPHYSIC', message3, __FILE__, __LINE__)

            ENDIF

            deallocate (plane)

# endif

        !--------------------------------------------------------------------------- END

END SUBROUTINE INTERCAMBIAWRITEPHYSIC

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOPHYX (ini, fin, ny, nz, mnx5, mny5, mnz5, bounds, physic)

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
        type(TboundPhysic) :: bounds (5,ny,nz)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
           DO jj= 1, ny
              DO ii= ini, fin
                 
                 physic(ii,jj,kk)%velx    = bounds(ii-ini+1,jj,kk)%velx
                 physic(ii,jj,kk)%vely    = bounds(ii-ini+1,jj,kk)%vely
                 physic(ii,jj,kk)%velz    = bounds(ii-ini+1,jj,kk)%velz
                 physic(ii,jj,kk)%bxc     = bounds(ii-ini+1,jj,kk)%bxc
                 physic(ii,jj,kk)%byc     = bounds(ii-ini+1,jj,kk)%byc
                 physic(ii,jj,kk)%bzc     = bounds(ii-ini+1,jj,kk)%bzc
                 physic(ii,jj,kk)%densty  = bounds(ii-ini+1,jj,kk)%densty
                 physic(ii,jj,kk)%denstye = bounds(ii-ini+1,jj,kk)%denstye
                 physic(ii,jj,kk)%eps     = bounds(ii-ini+1,jj,kk)%eps
                 physic(ii,jj,kk)%pres    = bounds(ii-ini+1,jj,kk)%pres
                 physic(ii,jj,kk)%tracer  = bounds(ii-ini+1,jj,kk)%tracer
      
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOPHYX

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Send y-boundary information to physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOPHYY (ini, fin, nx, nz, mnx5, mny5, mnz5, bounds, physic)

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
        type(TboundPhysic) :: bounds (nx,5,nz)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
           DO jj= ini, fin
              DO ii= 1, nx
                 
                 physic(ii,jj,kk)%velx    = bounds(ii,jj-ini+1,kk)%velx
                 physic(ii,jj,kk)%vely    = bounds(ii,jj-ini+1,kk)%vely
                 physic(ii,jj,kk)%velz    = bounds(ii,jj-ini+1,kk)%velz
                 physic(ii,jj,kk)%bxc     = bounds(ii,jj-ini+1,kk)%bxc
                 physic(ii,jj,kk)%byc     = bounds(ii,jj-ini+1,kk)%byc
                 physic(ii,jj,kk)%bzc     = bounds(ii,jj-ini+1,kk)%bzc
                 physic(ii,jj,kk)%densty  = bounds(ii,jj-ini+1,kk)%densty
                 physic(ii,jj,kk)%denstye = bounds(ii,jj-ini+1,kk)%denstye
                 physic(ii,jj,kk)%eps     = bounds(ii,jj-ini+1,kk)%eps
                 physic(ii,jj,kk)%pres    = bounds(ii,jj-ini+1,kk)%pres
                 physic(ii,jj,kk)%tracer  = bounds(ii,jj-ini+1,kk)%tracer
      
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOPHYY


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
!===============================================================================
! NAME              :: FROMBOUNDTOPHYZ
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save z-boundary information into physical variables.
!===============================================================================
      
SUBROUTINE FROMBOUNDTOPHYZ (ini, fin, nx, ny, mnx5, mny5, mnz5, bounds, physic)

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
        type(TboundPhysic) :: bounds (nx,ny,5)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= ini, fin
           DO jj= 1, ny
              DO ii= 1, nx
                 
                 physic(ii,jj,kk)%velx    = bounds(ii,jj,kk-ini+1)%velx
                 physic(ii,jj,kk)%vely    = bounds(ii,jj,kk-ini+1)%vely
                 physic(ii,jj,kk)%velz    = bounds(ii,jj,kk-ini+1)%velz
                 physic(ii,jj,kk)%bxc     = bounds(ii,jj,kk-ini+1)%bxc
                 physic(ii,jj,kk)%byc     = bounds(ii,jj,kk-ini+1)%byc
                 physic(ii,jj,kk)%bzc     = bounds(ii,jj,kk-ini+1)%bzc
                 physic(ii,jj,kk)%densty  = bounds(ii,jj,kk-ini+1)%densty
                 physic(ii,jj,kk)%denstye = bounds(ii,jj,kk-ini+1)%denstye
                 physic(ii,jj,kk)%eps     = bounds(ii,jj,kk-ini+1)%eps
                 physic(ii,jj,kk)%pres    = bounds(ii,jj,kk-ini+1)%pres
                 physic(ii,jj,kk)%tracer  = bounds(ii,jj,kk-ini+1)%tracer
      
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMBOUNDTOPHYZ

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDX
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save x-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMPHYTOBOUNDX (ini, fin, ny, nz, mnx5, mny5, mnz5, bounds, physic)

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
        type(TboundPhysic) :: bounds (5,ny,nz)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
           DO jj= 1, ny
              DO ii= ini, fin
                 
                 bounds(ii-ini+1,jj,kk)%velx    = physic(ii,jj,kk)%velx    
                 bounds(ii-ini+1,jj,kk)%vely    = physic(ii,jj,kk)%vely    
                 bounds(ii-ini+1,jj,kk)%velz    = physic(ii,jj,kk)%velz
                 bounds(ii-ini+1,jj,kk)%bxc     = physic(ii,jj,kk)%bxc
                 bounds(ii-ini+1,jj,kk)%byc     = physic(ii,jj,kk)%byc
                 bounds(ii-ini+1,jj,kk)%bzc     = physic(ii,jj,kk)%bzc
                 bounds(ii-ini+1,jj,kk)%densty  = physic(ii,jj,kk)%densty  
                 bounds(ii-ini+1,jj,kk)%denstye = physic(ii,jj,kk)%denstye  
                 bounds(ii-ini+1,jj,kk)%eps     = physic(ii,jj,kk)%eps     
                 bounds(ii-ini+1,jj,kk)%pres    = physic(ii,jj,kk)%pres     
                 bounds(ii-ini+1,jj,kk)%tracer  = physic(ii,jj,kk)%tracer 
      
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMPHYTOBOUNDX


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDY
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save y-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMPHYTOBOUNDY (ini, fin, nx, nz, mnx5, mny5, mnz5, bounds, physic)

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
        type(TboundPhysic) :: bounds (nx,5,nz)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= 1, nz
           DO jj= ini, fin
              DO ii= 1, nx
                 
                 bounds(ii,jj-ini+1,kk)%velx    = physic(ii,jj,kk)%velx    
                 bounds(ii,jj-ini+1,kk)%vely    = physic(ii,jj,kk)%vely    
                 bounds(ii,jj-ini+1,kk)%velz    = physic(ii,jj,kk)%velz
                 bounds(ii,jj-ini+1,kk)%bxc     = physic(ii,jj,kk)%bxc
                 bounds(ii,jj-ini+1,kk)%byc     = physic(ii,jj,kk)%byc
                 bounds(ii,jj-ini+1,kk)%bzc     = physic(ii,jj,kk)%bzc
                 bounds(ii,jj-ini+1,kk)%densty  = physic(ii,jj,kk)%densty  
                 bounds(ii,jj-ini+1,kk)%denstye = physic(ii,jj,kk)%denstye  
                 bounds(ii,jj-ini+1,kk)%eps     = physic(ii,jj,kk)%eps     
                 bounds(ii,jj-ini+1,kk)%pres    = physic(ii,jj,kk)%pres     
                 bounds(ii,jj-ini+1,kk)%tracer  = physic(ii,jj,kk)%tracer 
      
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMPHYTOBOUNDY


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: FROMPHYTOBOUNDZ
! GLOBAL VARIABLES  :: 
! DESCRIPTION       :: Save z-boundary information from physical variables.
!===============================================================================
      
SUBROUTINE FROMPHYTOBOUNDZ (ini, fin, nx, ny, mnx5, mny5, mnz5, bounds, physic)

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
        type(TboundPhysic) :: bounds (nx,ny,5)

        !------------------------------------------------------- Output Variables

        type(Tphysic) :: physic (-4:mnx5, -4:mny5, -4:mnz5)
   
        !------------------------------------------------------- Local Variables

        integer(is) :: ii, jj, kk
        
        !--------------------------------------------------------------------------- BEGIN

        DO kk= ini, fin
           DO jj= 1, ny
              DO ii= 1, nx
                 
                 bounds(ii,jj,kk-ini+1)%velx    = physic(ii,jj,kk)%velx    
                 bounds(ii,jj,kk-ini+1)%vely    = physic(ii,jj,kk)%vely    
                 bounds(ii,jj,kk-ini+1)%velz    = physic(ii,jj,kk)%velz
                 bounds(ii,jj,kk-ini+1)%bxc     = physic(ii,jj,kk)%bxc
                 bounds(ii,jj,kk-ini+1)%byc     = physic(ii,jj,kk)%byc
                 bounds(ii,jj,kk-ini+1)%bzc     = physic(ii,jj,kk)%bzc
                 bounds(ii,jj,kk-ini+1)%densty  = physic(ii,jj,kk)%densty  
                 bounds(ii,jj,kk-ini+1)%denstye = physic(ii,jj,kk)%denstye  
                 bounds(ii,jj,kk-ini+1)%eps     = physic(ii,jj,kk)%eps     
                 bounds(ii,jj,kk-ini+1)%pres    = physic(ii,jj,kk)%pres     
                 bounds(ii,jj,kk-ini+1)%tracer  = physic(ii,jj,kk)%tracer 
      
              ENDDO
           ENDDO
        ENDDO

        !--------------------------------------------------------------------------- END     

END SUBROUTINE FROMPHYTOBOUNDZ


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!===============================================================================
! NAME              :: REDUCEMINDT
! GLOBAL VARIABLES  :: mpi_dims
! DESCRIPTION       :: Gather dt between MPI-nodes.
!===============================================================================


SUBROUTINE REDUCEMINDT2 (ic, jc, kc, dt, sound, v, reducemindt)

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

        !------------------------------------------ Output Variables

        integer(is), intent(inout) :: ic, jc, kc

 	real(rs), intent(inout) :: dt, sound, v

	integer(is) :: reducemindt

        !------------------------------------------ Local Variables

#ifdef DEBUG_COMM
        character(256) :: fich, cadout
  	character(4),external :: int2ch4
#endif

 	character(256) :: cadout

#ifdef PARALELO
        integer(is) :: ind, maxind, ierr

  	type(TstepInfoReal)          :: stepInfoRealSource
        type(TstepInfoInt)           :: stepInfoIntSource
        type(TstepInfoReal), pointer :: stepInfoRealDest(:)
        type(TstepInfoInt),  pointer :: stepInfoIntDest(:)
#endif

        !--------------------------------------------------------------------------- BEGIN

# ifdef PARALELO

            ! TstepInfo destination vectors.
            stepInfoRealDest => wass100
            stepInfoIntDest  => wass101
            maxind = (mpi_dims(1)) * (mpi_dims(2)) * (mpi_dims(3))

            !! Information we want to gather.
            stepInfoRealSource%dt    = dt
            stepInfoRealSource%sound = sound
            stepInfoRealSource%v     = v

            stepInfoIntSource%ic = ic
            stepInfoIntSource%jc = jc
            stepInfoIntSource%kc = kc


	    !call PAR_BARRIER()

            ! Gather the information from different MPI-processes and save it in destination.
            call MPI_GATHER (stepInfoRealSource, 3, MPI_DOUBLE_PRECISION,        &
                             stepInfoRealDest,   3, MPI_DOUBLE_PRECISION,        &
                             0, CART_WORLD, ierr)


            IF (ierr .NE. MPI_SUCCESS) THEN
               call ERROR ('REDUCEMINDT', 'Error in MPI_Gather', __FILE__, __LINE__)
            ENDIF
      
  
	    !call PAR_BARRIER()


            call MPI_GATHER (stepInfoIntSource, 3, MPI_INTEGER,                  &
                             stepInfoIntDest, 3, MPI_INTEGER,                    &
                             0, CART_WORLD, ierr)

            IF (ierr .NE. MPI_SUCCESS) THEN
               call ERROR ('REDUCEMINDT', 'Error in MPI_Gather', __FILE__, __LINE__)
            ENDIF

	    call PAR_BARRIER()


            ! dt has been gathered to the first CPU.
            IF (cartrank .EQ. 0) THEN

	       
               REDUCEMINDT = 0
               
               DO ind= 1, maxind-1
         
                  IF ( stepInfoRealDest(ind)%dt .LT. dt ) THEN

                     dt = stepInfoRealDest(ind)%dt
                     REDUCEMINDT = ind            
 
                  ENDIF

               ENDDO
   
               IF (REDUCEMINDT .EQ. 0) THEN

                  ic = stepInfoIntSource%ic
                  jc = stepInfoIntSource%jc
                  kc = stepInfoIntSource%kc
                  sound = stepInfoRealSource%sound
                  v     = stepInfoRealSource%v
                  dt    = stepInfoRealSource%dt

               ELSE

                  ic = stepInfoIntDest(reduceMinDT)%ic
                  jc = stepInfoIntDest(reduceMinDT)%jc
                  kc = stepInfoIntDest(reduceMinDT)%kc
                  sound = stepInfoRealDest(reduceMinDT)%sound
                  v     = stepInfoRealDest(reduceMinDT)%v

               ENDIF
               
            ENDIF
      
            !call MPI_BARRIER (CART_WORLD, ierr)
            !IF (ierr .NE. 0) THEN
            !   call ERROR ('PAR_BARRIER', 'Error in MPI_BARRIER', __FILE__, __LINE__)
            !ENDIF
             
	    call PAR_BARRIER()


            ! Send the information to all MPI-nodes.
            call MPI_BCAST (dt, 1, MPI_DOUBLE_PRECISION, 0, CART_WORLD, ierr)
            call MPI_BCAST (reducemindt, 1, MPI_INT, 0, CART_WORLD, ierr)


# endif


        !--------------------------------------------------------------------------- END

END SUBROUTINE REDUCEMINDT2
