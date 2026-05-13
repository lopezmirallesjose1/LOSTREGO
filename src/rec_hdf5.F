

!*********************************
!     3D-RMHD CODE: L�STREGO
!*********************************


#ifdef HDF


! =============================================================================== 
SUBROUTINE RRHDF5_INIT(H_openr, H_opend)
! =============================================================================== 

        USE hdf5
        USE math_constants
        USE physic_constants
        USE parallel
        USE tipos
        USE workarrays

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif
      
        IMPLICIT NONE


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif
      
# include "types.h"

#ifdef PARALELO
        integer  ierror
#endif
      
        integer :: error, H_openr, H_opend

        call h5open_f(error)

        H_openr = 0
        H_opend = 0

        RETURN

! ===============================================================================
END SUBROUTINE RRHDF5_INIT 
! ===============================================================================


! ===============================================================================
SUBROUTINE RRHDF5_RECORD(ird, datfil, rstfil, outfil,                             &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                        & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                        &
               gamma_ad, xzn, yzn, zzn, xznr, xznl, yznr, yznl, zznr, zznl,       &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                    &
               gridlx, gridly, gridlz, nx, ny, nz, nsdim,                         &
               nend, irstrt, nrstrt, trstrt, nout, tout, itstp, iars, irk, niter, &
               igeomx, igeomy, igeomz, eqs, tmax, cfl, dtini, dt, timej, nstep,   &
               small, smlrho, smallp, smallu,                                     &
               mpi_coordx, mpi_coordy, mpi_coordz, physic,CT,mnx6,mny6,mnz6)
! ===============================================================================


        USE hdf5
        USE math_constants
        USE physic_constants
        USE parallel
        USE tipos
        USE workarrays

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif


        implicit none      


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif
      
#       include "types.h"
     
#ifdef PARALELO
        integer ierror
#endif


        !------------------------------------------------------- Local Variables
      
        integer :: error, curn, dec, i, j, k, r_filer, r_filed, curR_R, nx_rec, ny_rec, nz_rec
        
        character :: c1, c2, c3, c4, c5, gc1, gc2, gc3, gc4, gc5, citers(10), mc1, mc2, mc3, &
                     mc4, mc5, mc6, mc7, mc8, mc9
        
        character(256) :: path2,fich, path_out

        
        !------------------------------------------------------- Input Variables
      
        character*(*), intent(in) :: datfil
        character*(*), intent(in) :: rstfil, outfil
        character(4), intent(in) :: basenm

        integer:: ird, plane,  mpi_coordx,   &
                  mpi_coordy, mpi_coordz

        integer(is) :: nx, ny, nz, mnx1, mny1, mnz1, mnx5, mny5, mnz5,mnx6,mny6,mnz6

      
        integer(is) :: nrst, nout1,  nout, nsdim, nend, irstrt, nrstrt,       &
                       igeomx, igeomy, igeomz, nstep, bndmnz, bndmxz,         &
                       bndmny, bndmxy, bndmxx, bndmnx, eqs,trstrt

        real(rs) :: tout1, trst, timej,                     &
                    gamma_ad, gridlx, gridly, gridlz

        real(rs) :: xznl(mnx1), xzn(nx), xznr(nx), &
                    yznl(mny1), yzn(ny), yznr(ny), &
                    zznl(mnz1), zzn(nz), zznr(nz)


        character(4)::limiter

        integer(is) :: itstp, iars, irk, niter, H_openr, H_opend

        real(rs) :: tmax, tout, cfl, dtini, small, smlrho, smallp, smallu, dt
        

        type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)
        type(ConsTransp) :: CT(-4:mnx6, -4:mny6, -4:mnz6)

      
        real(rs) :: dbufa(nx,ny),dbuft(nx,ny,mpi_dims(3))

        INTEGER(HID_T) :: file_id

        
        !     dataspaces
        integer(HID_T) :: d3space, d3filespace, d3memspace
        integer(HID_T) :: d2space, d2filespace, d2memspace
        integer(HID_T) :: d2spacea, d2filespacea, d2memspacea
        integer(HID_T) :: dxspace, dxfilespace, dxmemspace
        integer(HID_T) :: dyspace, dyfilespace, dymemspace
        integer(HID_T) :: dzspace, dzfilespace, dzmemspace
        integer(HID_T) :: dblspace, intspace

        !     groups
        integer(HID_T) :: group_id

        !     datasets
        integer(HID_T) :: d3set, d2seta, d2set, dxset, dyset, dzset, dblset, intset

        !     dimensions
        INTEGER(HSIZE_T), DIMENSION(3) :: d3dims, d3totdims
        INTEGER(HSIZE_T), DIMENSION(2) :: d2dimsa, d2totdimsa
        INTEGER(HSIZE_T), DIMENSION(3) :: d2dims, d2totdims
        INTEGER(HSIZE_T), DIMENSION(1) :: dxdims, dxtotdims
        INTEGER(HSIZE_T), DIMENSION(1) :: dydims, dytotdims
        INTEGER(HSIZE_T), DIMENSION(1) :: dzdims, dztotdims
        INTEGER(HSIZE_T), DIMENSION(1) :: dbldims = (/ 1 /), intdims = (/ 1 /)
    
        !     offsets
        INTEGER(HSIZE_T), DIMENSION(3) :: d3offset, d2offset
        INTEGER(HSIZE_T), DIMENSION(2) :: d2offseta
        INTEGER(HSIZE_T), DIMENSION(1) :: dxoffset, dyoffset, dzoffset
        
        !     attributes and properties
        integer(HID_T) :: data_transfer, file_access

        !     physical grid
        integer :: imin, imax, jmin, jmax, kmin, kmax

        !     integer and double data
        integer, dimension(1) :: intdata
        doubleprecision, dimension(1) :: dbldata

        !     MPI stuff
        integer :: hdcomm, hdinfo
  
        !     Buffers 3D
        real(rs), allocatable :: physicbuff(:,:,:)


        !-------------------------------------------------------------------------- BEGIN

      
        !     initialize dimensions
        d3dims = (/nx,ny,nz/)
        d2dims = (/nx,ny,1/)
        dxdims = (/nx/)
        dydims = (/ny/)
        dzdims = (/nz/)
      
    
        !     compute physical domain indices
        
        IF (ird .EQ. 0) THEN
           
           plane = 0
           imin = 1
           imax = nx
           jmin = 1
           jmax = ny
           kmin = 1
           kmax = nz
           allocate (physicbuff(imin:imax,jmin:jmax,kmin:kmax))
           
        ENDIF
      
       
        !     Create filename strings
        
#ifdef PARALELO
        
        hdcomm = MPI_COMM_WORLD
        hdinfo = MPI_INFO_NULL
        call MPI_INFO_CREATE(hdinfo, ierror)
      
#endif

        !     Compute offsets and total sizes

#ifdef PARALELO

        dxoffset(1)  = mpi_coordx * nx
        dxtotdims(1) = mpi_dims(1)* nx
     
        dyoffset(1)  = mpi_coordy * ny
        dytotdims(1) = mpi_dims(2)* ny
      
        dzoffset(1)  = mpi_coordz * nz
        dztotdims(1) = mpi_dims(3)* nz
           
        d3offset(1)  = dxoffset(1)
        d3offset(2)  = dyoffset(1)
        d3offset(3)  = dzoffset(1)
     
        d2offset(1)  = dxoffset(1)
        d2offset(2)  = dyoffset(1)
        d2offset(3)  = mpi_coordz

        d3totdims(1) = dxtotdims(1)
        d3totdims(2) = dytotdims(1)
        d3totdims(3) = dztotdims(1)

        d2totdims(1) =  mpi_dims(1) * nx
        d2totdims(2) =  mpi_dims(2) * ny
        d2totdims(3) =  mpi_dims(3)

      
#else
      
        dxtotdims  = dxdims
        dytotdims  = dydims
        dztotdims  = dzdims
        d3totdims  = d3dims

#endif /* offsets and dimensions*/

        !     create total dataspaces
        
        IF (ird .EQ. 0) THEN
           call h5screate_simple_f(3, d3totdims, d3space, error)     
        ENDIF

        call h5screate_simple_f(3, d2totdims, d2space, error)    
        call h5screate_simple_f(1, dxtotdims, dxspace, error)
      
        IF (nsdim .GT. 1) THEN
           call h5screate_simple_f(1, dytotdims, dyspace, error)
           call h5screate_simple_f(1, dztotdims, dzspace, error)
        endif


#ifdef PARALELO
        
        !     Create memory spaces for MPI collective output
        
        IF (ird .EQ. 0) THEN
           call h5screate_simple_f(3, d3dims, d3memspace, error)     
           call h5screate_simple_f(3, d2dims, d2memspace, error)
        ENDIF

        call h5screate_simple_f(1, dxdims, dxmemspace, error)
        
        IF (nsdim .GT. 1) THEN
           call h5screate_simple_f(1, dydims, dymemspace, error)
           call h5screate_simple_f(1, dzdims, dzmemspace, error)
        ENDIF

#endif /* buffers */

        !     create dataspaces for header
        call h5screate_simple_f(1, intdims, intspace, error)
        call h5screate_simple_f(1, dbldims, dblspace, error)

        !     create transfer properties
#ifdef PARALELO
        call h5pcreate_f(H5P_DATASET_XFER_F, data_transfer, error)
        call h5pset_dxpl_mpio_f(data_transfer, H5FD_MPIO_COLLECTIVE_F, error)
        call h5pcreate_f(H5P_FILE_ACCESS_F, file_access, error)
        call h5pset_fapl_mpio_f(file_access, hdcomm, hdinfo, error)
#endif

#ifdef PARALELO
        call MPI_INFO_FREE(hdinfo, ierror)
#endif


        !     create filename

        IF (ird .EQ. 0) THEN
           curn=R_filer
        ELSE
           curn=R_filed
        ENDIF
      
        c1=achar(48+(curn-mod(curn,10000))/10000)
        curn=mod(curn,10000)
        c2=achar(48+(curn-mod(curn,1000))/1000)
        curn=mod(curn,1000)
        c3=achar(48+(curn-mod(curn,100))/100)
        curn=mod(curn,100)
        c4=achar(48+(curn-mod(curn,10))/10)
        curn=mod(curn,10)
        c5=achar(48+curn)


        IF (ird .EQ. 0) THEN
           
           !     create new file if necessary
           IF (H_openr .EQ. 0) THEN
              IF (cartrank .EQ. 0) THEN
                 call h5fcreate_f(TRIM(outRoot)//rstfil// &
                 '-'//c1//c2//c3//c4//c5//'.h5', &
                 H5F_ACC_TRUNC_F, file_id, error &
                 )
              ENDIF
    
              curR_R = 1
         
              H_openr = 1
              IF (cartrank .EQ. 0) THEN
                 call h5fclose_f(file_id, error)
              ENDIF
           ENDIF
      
           IF (cartrank .EQ. 0) THEN
              call h5fopen_f(TRIM(outRoot)//rstfil// &
              '-'//c1//c2//c3//c4//c5//'.h5', &
              H5F_ACC_RDWR_F, file_id, error &
              )
           ENDIF

        ENDIF

        !     add new group for the current iteration
     
        curn = curR_R
     
        gc1=achar(48+(curn-mod(curn,10000))/10000)
        curn=mod(curn,10000)
        gc2=achar(48+(curn-mod(curn,1000))/1000)
        curn=mod(curn,1000)
        gc3=achar(48+(curn-mod(curn,100))/100)
        curn=mod(curn,100)
        gc4=achar(48+(curn-mod(curn,10))/10)
        curn=mod(curn,10)
        gc5=achar(48+curn)

        IF (cartrank .EQ. 0) &
        call h5gcreate_f(file_id, gc1//gc2//gc3//gc4//gc5, group_id, error)


        IF (cartrank .EQ. 0) THEN
           
           IF (ird .EQ. 0) THEN

              ! nend
              call h5dcreate_f(group_id,'nend', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nend, intdims, error )
              call h5dclose_f(intset, error)


              ! irstrt
              call h5dcreate_f(group_id,'irstrt', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,irstrt, intdims, error )
              call h5dclose_f(intset, error)


              ! nrstrt
              call h5dcreate_f(group_id,'nrstrt', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nrstrt, intdims, error )
              call h5dclose_f(intset, error) 


              !nout
              call h5dcreate_f(group_id,'nout', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nout, intdims, error )
              call h5dclose_f(intset, error) 


              ! itstp
              call h5dcreate_f(group_id,'itstp', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,itstp, intdims, error )
              call h5dclose_f(intset, error)


              ! iars
              call h5dcreate_f(group_id,'iars', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,iars, intdims, error )
              call h5dclose_f(intset, error)


              ! irk
              call h5dcreate_f(group_id,'irk', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,irk, intdims, error )
              call h5dclose_f(intset, error)


              ! niter
              call h5dcreate_f(group_id,'niter', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,niter, intdims, error )
              call h5dclose_f(intset, error)
         
            
              ! nsdim
              call h5dcreate_f(group_id,'nsdim', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nsdim, intdims, error )
              call h5dclose_f(intset, error)

         
              nx_rec = nx * mpi_dims(1)
              ny_rec = ny * mpi_dims(2)
              nz_rec = nz * mpi_dims(3)

         
              !nx
              call h5dcreate_f(group_id,'nx', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nx_rec,intdims, error )
              call h5dclose_f(intset, error)
                       
              !ny
              call h5dcreate_f(group_id,'ny', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,ny_rec,intdims, error )
              call h5dclose_f(intset, error)

              !nz
              call h5dcreate_f(group_id,'nz', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nz_rec,intdims, error )
              call h5dclose_f(intset, error)


              ! igeomx 
              call h5dcreate_f(group_id,'igeomx', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,igeomx,intdims, error )
              call h5dclose_f(intset, error)

         
              ! igeomy       
              call h5dcreate_f(group_id,'igeomy', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,igeomy,intdims, error )
              call h5dclose_f(intset, error)


              !igeomz
              call h5dcreate_f(group_id,'igeomz', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,igeomz,intdims, error )
              call h5dclose_f(intset, error)


              !eqs
              call h5dcreate_f(group_id,'eqs', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,eqs,intdims, error )
              call h5dclose_f(intset, error)


              ! tmax
              call h5dcreate_f(group_id, 'tmax', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,tmax,dbldims, error )
              call h5dclose_f(dblset, error) 


              !trstrt
              call h5dcreate_f(group_id,'trstrt', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,trstrt,intdims, error )
              call h5dclose_f(intset, error)


              !nstep
              call h5dcreate_f(group_id,'nstep', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,nstep,intdims, error )
              call h5dclose_f(intset, error)


              !tout
              call h5dcreate_f(group_id, 'tout', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,tout,dbldims, error )
              call h5dclose_f(dblset, error)
         
         
              ! cfl
              call h5dcreate_f(group_id, 'cfl', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,cfl,dbldims, error )
              call h5dclose_f(dblset, error)


              ! dtini
              call h5dcreate_f(group_id, 'dtini', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,dtini,dbldims, error )
              call h5dclose_f(dblset, error)


              ! dt
              call h5dcreate_f(group_id, 'dt', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,dt,dbldims, error )
              call h5dclose_f(dblset, error)


              ! timej
              call h5dcreate_f(group_id, 'timej', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,timej,dbldims, error )
              call h5dclose_f(dblset, error)


              ! small
              call h5dcreate_f(group_id, 'small', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,small,dbldims, error )
              call h5dclose_f(dblset, error)


              ! smlrho
              call h5dcreate_f(group_id, 'smlrho', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,smlrho,dbldims, error )
              call h5dclose_f(dblset, error)


              ! smallp
              call h5dcreate_f(group_id, 'smallp', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,smallp,dbldims, error )
              call h5dclose_f(dblset, error)


              ! smallu
              call h5dcreate_f(group_id, 'smallu', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,smallu,dbldims, error )
              call h5dclose_f(dblset, error)

         
              ! gridlx
              call h5dcreate_f(group_id, 'gridlx', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gridlx,dbldims, error )
              call h5dclose_f(dblset, error) 

              ! gridly
              call h5dcreate_f(group_id, 'gridly', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gridly,dbldims, error )
              call h5dclose_f(dblset, error) 
       
              ! gridlz
              call h5dcreate_f(group_id, 'gridlz', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gridlz,dbldims, error )
              call h5dclose_f(dblset, error) 


              ! bndmnx
              call h5dcreate_f(group_id,'bndmnx', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmnx,intdims, error )
              call h5dclose_f(intset, error)
       

              ! bndmxx
              call h5dcreate_f(group_id,'bndmxx', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmxx,intdims, error  )
              call h5dclose_f(intset, error)
        

              ! bndmny
              call h5dcreate_f(group_id,'bndmny', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmny,intdims, error )
              call h5dclose_f(intset, error)
       

              ! bndmxy
              call h5dcreate_f(group_id,'bndmxy', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmxy,intdims, error )
              call h5dclose_f(intset, error)
        

              ! bndmnz
              call h5dcreate_f(group_id,'bndmnz', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmnz,intdims, error )
              call h5dclose_f(intset, error)
       

              ! bndmxz
              call h5dcreate_f(group_id,'bndmxz', H5T_NATIVE_INTEGER, intspace, intset, error)
              call h5dwrite_f(intset, H5T_NATIVE_INTEGER,bndmxz,intdims, error )
              call h5dclose_f(intset, error)
       

              ! gamma_ad
              call h5dcreate_f(group_id, 'gamma_ad', H5T_NATIVE_DOUBLE, dblspace, dblset, error)
              call h5dwrite_f(dblset, H5T_NATIVE_DOUBLE,gamma_ad,dbldims, error )
              call h5dclose_f(dblset, error) 
   
           ENDIF
       
           !     close group
           call h5gclose_f(group_id, error)

           !     close file         
           call h5fclose_f(file_id, error)
        
        ENDIF
 

        IF (ird .EQ. 0) THEN
         
#ifdef PARALELO
           
           call h5fopen_f(TRIM(outRoot)//rstfil// &
           '-'//c1//c2//c3//c4//c5//'.h5',        &
           H5F_ACC_RDWR_F, file_id, error         &
           ,access_prp = file_access )
        
#else
        
           call h5fopen_f(TRIM(outRoot)//rstfil// &
           '-'//c1//c2//c3//c4//c5//'.h5',        &
           H5F_ACC_RDWR_F, file_id, error )
        
#endif
        
        ENDIF
      
        call h5gOPEN_f(file_id, gc1//gc2//gc3//gc4//gc5, group_id, error)

        IF (cartrank .EQ. 0) THEN
           write(6,*) 'REC_HDF, 2: ',imin, imax, jmin, jmax, kmin, kmax 
        ENDIF


        IF (ird .EQ. 0) THEN
         
#ifdef PARALELO

           !   ----------
           !   X
           !   ----------         

           call h5dcreate_f(group_id, 'xzn', H5T_NATIVE_DOUBLE, dxspace, dxset, error)           
           call h5dget_space_f(dxset, dxfilespace, error)
           call h5sselect_hyperslab_f(dxfilespace, H5S_SELECT_SET_F, dxoffset, dxdims, error)
           call h5dwrite_f(dxset, H5T_NATIVE_DOUBLE, xzn(imin:imax), dxdims, error &
                           ,file_space_id = dxfilespace, mem_space_id = dxmemspace &
                           ,xfer_prp = data_transfer )

                           
           !   ----------              
           !   Z
           !   ----------         
                           
           call h5dcreate_f(group_id, 'zzn', H5T_NATIVE_DOUBLE, dzspace, dzset, error)
           call h5dget_space_f(dzset, dzfilespace, error)
           call h5sselect_hyperslab_f(dzfilespace, H5S_SELECT_SET_F, dzoffset, dzdims, error)
           call h5dwrite_f(dzset, H5T_NATIVE_DOUBLE, zzn(kmin:kmax), dzdims, error &
                           ,file_space_id = dzfilespace, mem_space_id = dzmemspace &
                           ,xfer_prp = data_transfer )


           !   ----------                                    
           !   Y
           !   ----------         
                
           call h5dcreate_f(group_id, 'yzn', H5T_NATIVE_DOUBLE, dyspace, dyset, error)
           call h5dget_space_f(dyset, dyfilespace, error)
           call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)
           call h5dwrite_f(dyset, H5T_NATIVE_DOUBLE, yzn(jmin:jmax), dydims, error &
                           ,file_space_id = dyfilespace, mem_space_id = dymemspace &
                           ,xfer_prp = data_transfer )

#else
           call h5dwrite_f(dxset, H5T_NATIVE_DOUBLE, xzn(imin:imax), dxdims, error )
           call h5dwrite_f(dyset, H5T_NATIVE_DOUBLE, yzn(jmin:jmax), dydims, error )
           call h5dwrite_f(dzset, H5T_NATIVE_DOUBLE, zzn(kmin:kmax), dzdims, error )
#endif
         
           call h5dclose_f(dxset, error)
           call h5dclose_f(dyset, error)
           call h5dclose_f(dzset, error)

#ifdef PARALELO
           
           call h5sclose_f(dxfilespace, error)
           call h5sclose_f(dyfilespace, error)
           call h5sclose_f(dzfilespace, error)
           
#endif


         
           IF (cartrank .EQ. 0) THEN
              write(6,*) 'REC_HDF, 3: ', imin, imax, kmin,kmax, jmin, jmax  
           ENDIF
         
           !   densty
           call h5dcreate_f(group_id, 'densty', H5T_NATIVE_DOUBLE, d3space, d3set, error)

           
#ifdef PARALELO
           
           call h5dget_space_f(d3set, d3filespace, error)
           call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
           
           physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%densty
           call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
           
#ifdef PARALELO
           
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
           
#endif
           )
           call h5dclose_f(d3set, error)
           
#ifdef PARALELO
           call h5sclose_f(d3filespace, error)
#endif


           !   denstye
           call h5dcreate_f(group_id, 'denstye', H5T_NATIVE_DOUBLE, d3space, d3set, error)

           
#ifdef PARALELO
           
           call h5dget_space_f(d3set, d3filespace, error)
           call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
           
#endif
           
           physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%denstye
           call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
      
#ifdef PARALELO
           
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
           
#endif           
           )
           call h5dclose_f(d3set, error)
           
#ifdef PARALELO
           call h5sclose_f(d3filespace, error)
#endif

      
           !     pres
           call h5dcreate_f(group_id, 'pres', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
           
           call h5dget_space_f(d3set, d3filespace, error)
           call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
           
           physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%pres
           call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
           
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
          )
          call h5dclose_f(d3set, error)
          
#ifdef PARALELO
          call h5sclose_f(d3filespace, error)
#endif


          !     velx
          call h5dcreate_f(group_id, 'velx', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
          
          call h5dget_space_f(d3set, d3filespace, error)
          call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
          
          physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%velx
          call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
          
#ifdef PARALELO
          ,file_space_id = d3filespace, mem_space_id = d3memspace &
          ,xfer_prp = data_transfer &
#endif
           )
           call h5dclose_f(d3set, error)
           
#ifdef PARALELO
           call h5sclose_f(d3filespace, error)
#endif

      
           IF (nsdim .GT. 1) THEN
         
              !     vely
              call h5dcreate_f(group_id, 'vely', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
              
              call h5dget_space_f(d3set, d3filespace, error)
              call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
              
              physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%vely
              call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
              
#ifdef PARALELO
              ,file_space_id = d3filespace, mem_space_id = d3memspace &
              ,xfer_prp = data_transfer &
#endif
              )
              call h5dclose_f(d3set, error)
           
#ifdef PARALELO
              call h5sclose_f(d3filespace, error)
#endif

           ENDIF

           IF (nsdim .GT. 2) THEN
         
              !     velz
              call h5dcreate_f(group_id, 'velz', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
           
              call h5dget_space_f(d3set, d3filespace, error)
              call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
           
#endif
              physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%velz
              call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
           
#ifdef PARALELO
           
              ,file_space_id = d3filespace, mem_space_id = d3memspace &
              ,xfer_prp = data_transfer &
           
#endif
              )
              call h5dclose_f(d3set, error)
           
#ifdef PARALELO
              call h5sclose_f(d3filespace, error)
#endif
              
           ENDIF
      
           !     tracer
           call h5dcreate_f(group_id, 'tracer', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
        
           call h5dget_space_f(d3set, d3filespace, error)
           call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
        
           physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%tracer
           call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
        
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
           )
           call h5dclose_f(d3set, error)
        
#ifdef PARALELO
           call h5sclose_f(d3filespace, error)
#endif


           !     bxc
           call h5dcreate_f(group_id, 'bxc', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
        
           call h5dget_space_f(d3set, d3filespace, error)
           call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
        
           physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%bxc
           call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
        
#ifdef PARALELO
           ,file_space_id = d3filespace, mem_space_id = d3memspace &
           ,xfer_prp = data_transfer &
#endif
           )
           call h5dclose_f(d3set, error)
        
#ifdef PARALELO
           call h5sclose_f(d3filespace, error)
#endif

      
           IF (nsdim .GT. 1) THEN
         
              !     byc
              call h5dcreate_f(group_id, 'byc', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
           
              call h5dget_space_f(d3set, d3filespace, error)
              call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
              physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%byc
              call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
           
#ifdef PARALELO
              ,file_space_id = d3filespace, mem_space_id = d3memspace &
              ,xfer_prp = data_transfer &
#endif
              )
              call h5dclose_f(d3set, error)
           
#ifdef PARALELO
              call h5sclose_f(d3filespace, error)
#endif

           ENDIF

           IF (nsdim .GT. 2) THEN
         
              !     bzc
              call h5dcreate_f(group_id, 'bzc', H5T_NATIVE_DOUBLE, d3space, d3set, error)

#ifdef PARALELO
           
              call h5dget_space_f(d3set, d3filespace, error)
              call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
           
#endif
              physicbuff(imin:imax, jmin:jmax, kmin:kmax) =  physic(imin:imax, jmin:jmax, kmin:kmax)%bzc
              call h5dwrite_f(d3set, H5T_NATIVE_DOUBLE, physicbuff(imin:imax, jmin:jmax, kmin:kmax), d3dims, error &
           
#ifdef PARALELO
              ,file_space_id = d3filespace, mem_space_id = d3memspace &
              ,xfer_prp = data_transfer &
#endif
              )
              call h5dclose_f(d3set, error)
           
#ifdef PARALELO
              call h5sclose_f(d3filespace, error)
#endif
           
           ENDIF   
      
           deallocate (physicbuff)
      
        ENDIF


        !     close group
        call h5gclose_f(group_id, error)

#ifdef PARALELO
        
        !     close transfer properties
        call h5pclose_f(data_transfer, error)

        !     close file properties
        call h5pclose_f(file_access, error)
        
#endif


#ifdef PARALELO
        
        !     close memory spaces for MPI collective output
        
        IF (nsdim .GT. 1) THEN
           
           call h5sclose_f(dymemspace, error)
           call h5sclose_f(dzmemspace, error)
         
        ENDIF
      
        call h5sclose_f(dxmemspace, error)
      
        IF (ird .EQ. 0) THEN
           call h5sclose_f(d3memspace, error)
        ENDIF
      
#endif /* PARALELO */

        !     close total memory spaces

        IF (nsdim .GT. 1) THEN
           call h5sclose_f(dyspace, error)
           call h5sclose_f(dzspace, error)
        ENDIF
        
        call h5sclose_f(dxspace, error)
        
        IF (ird .EQ. 0) THEN
           call h5sclose_f(d3space, error)
        ENDIF 

        !     create memory spaces for header
        call h5sclose_f(intspace, error)
        call h5sclose_f(dblspace, error)
     
        ! WE WILL USE ONLY ONE GROUP FOR THE RESTART. 

        IF (ird /= 0) THEN
           curR_R = curR_R + 1
        ENDIF

        !     close file         
        call h5fclose_f(file_id, error)

        IF (ird .EQ. 0) THEN        
           r_filer = r_filer + 1
           H_openr = 0         
        ENDIF
        
        RETURN

        !-------------------------------------------------------------------------- END

! ===============================================================================
END SUBROUTINE rrhdf5_record
! ===============================================================================
      

! ===============================================================================      
SUBROUTINE rrhdf5_finalize
! ===============================================================================

      USE hdf5      
      USE math_constants
      USE physic_constants
      USE parallel
      USE tipos
      USE workarrays

# if defined(PARALELO) && defined(MN)
        USE mpi
# endif

      IMPLICIT NONE


# if defined(PARALELO) && !defined(MN)
        include "mpif.h"
# endif
      
# include "types.h"

#ifdef PARALELO
      integer ierror
#endif

      integer :: error

      call h5close_f(error)
      
      RETURN

! ===============================================================================      
END SUBROUTINE rrhdf5_finalize
! ===============================================================================
      
#endif


