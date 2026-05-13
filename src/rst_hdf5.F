

!*********************************
!     3D-RMHD CODE: L�STREGO
!*********************************



      
! ===============================================================================
SUBROUTINE RST_HDF5(ird, datfil, rstfil, outfil,                                  &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                        &
               H_openr, H_opend, curR_R, r_filer, r_filed,                        &
               gamma_ad, xzn, yzn, zzn, xznr, xznl, yznr, yznl, zznr, zznl,       &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                    &
               gridlx, gridly, gridlz, nx, ny, nz, nsdim,                         &
               nend, irstrt, nrstrt, trstrt, nout, tout, itstp, iars, irk, niter, &
               igeomx, igeomy, igeomz, eqs, tmax, cfl, dtini, dt, timej, nstep,   &
               small, smlrho, smallp, smallu,                                     &
               mpi_coordx, mpi_coordy, mpi_coordz, physic)
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


        integer(is) :: rst_file, rst_rec

        
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

        integer(is) :: nx, ny, nz, mnx1, mny1, mnz1, mnx5, mny5, mnz5


        integer(is) :: nrst, nout1,  nout, nsdim, nend, irstrt, nrstrt,       &
                       igeomx, igeomy, igeomz, nstep, bndmnz, bndmxz,         &
                       bndmny, bndmxy, bndmxx, bndmnx, eqs,trstrt, irstrt2

        real(rs) :: tout1, trst, timej,                     &
                    gamma_ad, gridlx, gridly, gridlz

        real(rs) :: xznl(mnx1), xzn(nx), xznr(nx), &
                    yznl(mny1), yzn(ny), yznr(ny), &
                    zznl(mnz1), zzn(nz), zznr(nz)


        character(4)::limiter

        integer(is) :: itstp, iars, irk, niter, H_openr, H_opend

        real(rs) :: tmax, tout, cfl, dtini, small, smlrho, smallp, smallu, dt, tmax2, cfl2, tout2


        type(Tphysic):: physic(-4:mnx5,-4:mny5,-4:mnz5)


        !real(rs) :: dbufa(nx,ny),dbuft(nx,ny,mpi_dims(3))

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
        !real(rs), allocatable :: physicbuff(:,:,:)

 
        !     buffers  
        real(rs) :: dbufx(mnx1), &
                    dbufy(mny1), dbufz(mnz1)


	real(rs), allocatable :: dbuf(:,:,:)

        !-------------------------------------------------------------------------- BEGIN

	rst_file=irstrt
	rst_rec=1

        !     initialize dimensions
        d3dims = (/nx,ny,nz/)
        d2dims = (/nx,ny,1/)
        dxdims = (/nx/)
        dydims = (/ny/)
        dzdims = (/nz/)
      
        !     compute physical domain indices

        plane = 0
        imin = 1
        imax = nx
        jmin = 1
        jmax = ny
        kmin = 1
        kmax = nz

           allocate (dbuf(imin:imax,jmin:jmax,kmin:kmax))


 
        !     Create filename strings

#ifdef PARALELO

        hdcomm = CART_WORLD
        hdinfo = MPI_INFO_NULL
        call MPI_INFO_CREATE(hdinfo, ierror)

!     CPU coordinates string
        
        curn=mpi_coords(1)
        mc1=achar(48+(curn-mod(curn,100))/100)
        curn=mod(curn,100)
        mc2=achar(48+(curn-mod(curn,10))/10)
        curn=mod(curn,10)
        mc3=achar(48+curn)
      
        curn=mpi_coords(2)
        mc4=achar(48+(curn-mod(curn,100))/100)
        curn=mod(curn,100)
        mc5=achar(48+(curn-mod(curn,10))/10)
        curn=mod(curn,10)
        mc6=achar(48+curn)
      
        curn=mpi_coords(3)
        mc7=achar(48+(curn-mod(curn,100))/100)
        curn=mod(curn,100)
        mc8=achar(48+(curn-mod(curn,10))/10)
        curn=mod(curn,10)
        mc9=achar(48+curn)
      
#endif

        
!     compute offsets and total sizes
#ifdef PARALELO


         dxoffset(1)  = mpi_coords(1) * nx      
         dxtotdims(1) = mpi_dims(1) * nx
         
         dyoffset(1)  = mpi_coords(2) * ny       
         dytotdims(1) = mpi_dims(2) * ny
         
         dzoffset(1)  = mpi_coords(3) * nz     
         dztotdims(1) = mpi_dims(3) * nz 


         d3offset(1)  = dxoffset(1)
         d3offset(2)  = dyoffset(1)
         d3offset(3)  = dzoffset(1)
         
         d3totdims(1) = dxtotdims(1)
         d3totdims(2) = dytotdims(1)
         d3totdims(3) = dztotdims(1)

         
         d2totdims(1) =  mpi_dims(1) * nx
         d2totdims(2) =  mpi_dims(2) * ny
         d2totdims(3) =  mpi_dims(3)

#else

         dxtotdims = dxdims
         dytotdims = dydims
         dztotdims = dzdims
         d3totdims = d3dims

#endif /* offsets and dimensions*/


      
!     create dataspaces for header
      call h5screate_simple_f(1, intdims, intspace, error)
      call h5screate_simple_f(1, dbldims, dblspace, error)


!     create total dataspaces
      call h5screate_simple_f(3, d3totdims, d3space, error)
      call h5screate_simple_f(1, dxtotdims, dxspace, error)
      
      if (nsdim.gt.1) then
         call h5screate_simple_f(1, dytotdims, dyspace, error)
         call h5screate_simple_f(1, dztotdims, dzspace, error)
      endif

#ifdef PARALELO
!     create memory spaces for MPI collective output
      call h5screate_simple_f(3, d3dims, d3memspace, error)
      call h5screate_simple_f(1, dxdims, dxmemspace, error)
      
      if (nsdim.gt.1) then
         call h5screate_simple_f(1, dydims, dymemspace, error)
         call h5screate_simple_f(1, dzdims, dzmemspace, error)
      endif

#endif /* buffers */

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


!     Create filename
      R_filer = rst_file
      curn=R_filer
      c1=achar(48+(curn-mod(curn,10000))/10000)
      curn=mod(curn,10000)
      c2=achar(48+(curn-mod(curn,1000))/1000)
      curn=mod(curn,1000)
      c3=achar(48+(curn-mod(curn,100))/100)
      curn=mod(curn,100)
      c4=achar(48+(curn-mod(curn,10))/10)
      curn=mod(curn,10)
      c5=achar(48+curn)
      

      if (cartrank .EQ. 0) then
         print*,'[rst_hdf5] Restarting from ',&
         TRIM(outRoot)//rstfil// &
         '-'//c1//c2//c3//c4//c5//'.h5'
         print*,'[rst_hdf5] Record ',rst_rec
      endif

      call h5fopen_f(TRIM(outRoot)//rstfil// &
      '-'//c1//c2//c3//c4//c5//'.h5', &
      H5F_ACC_RDWR_F, file_id, error &
      
#ifdef PARALELO
      ,access_prp = file_access &
#endif
      )

      
!     open group for the desired iteration
      curR_R = rst_rec
      
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

      call h5gopen_f(file_id, gc1//gc2//gc3//gc4//gc5, group_id, error)


      ! nend
      call h5dopen_f(group_id, 'nend', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nend, intdims, error)
      call h5dclose_f(intset, error)


      ! irstrt
      call h5dopen_f(group_id, 'irstrt', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, irstrt2, intdims, error)
      call h5dclose_f(intset, error) 
      

      ! nrstrt
      call h5dopen_f(group_id, 'nrstrt', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nrstrt, intdims, error)
      call h5dclose_f(intset, error) 

      
      ! nout
      call h5dopen_f(group_id, 'nout', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nout, intdims, error)
      call h5dclose_f(intset, error)


      ! itstp
      call h5dopen_f(group_id, 'itstp', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, itstp, intdims, error)
      call h5dclose_f(intset, error)       
      

      ! iars
      call h5dopen_f(group_id, 'iars', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, iars, intdims, error)
      call h5dclose_f(intset, error)  
      

      ! irk
      call h5dopen_f(group_id, 'irk', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, irk, intdims, error)
      call h5dclose_f(intset, error)

      
      ! niter
      call h5dopen_f(group_id, 'niter', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, niter, intdims, error)
      call h5dclose_f(intset, error)  
          

      ! nsdim
      call h5dopen_f(group_id, 'nsdim', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nsdim, intdims, error)
      call h5dclose_f(intset, error)    


      ! nx
      call h5dopen_f(group_id, 'nx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nx_rec, intdims, error)
      call h5dclose_f(intset, error) 

                
      ! ny
      call h5dopen_f(group_id, 'ny', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, ny_rec, intdims, error)
      call h5dclose_f(intset, error)  

      
      ! nz
      call h5dopen_f(group_id, 'nz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nz_rec, intdims, error)
      call h5dclose_f(intset, error)  


              nx = nx_rec / mpi_dims(1)
              ny = ny_rec / mpi_dims(2)
              nz = nz_rec / mpi_dims(3)
      
      ! igeomx
      call h5dopen_f(group_id, 'igeomx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, igeomx, intdims, error)
      call h5dclose_f(intset, error)   

      ! igeomy
      call h5dopen_f(group_id, 'igeomy', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, igeomy, intdims, error)
      call h5dclose_f(intset, error)

      ! igeomz
      call h5dopen_f(group_id, 'igeomz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, igeomz, intdims, error)
      call h5dclose_f(intset, error)  

              
      ! eqs
      call h5dopen_f(group_id, 'eqs', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, eqs, intdims, error)
      call h5dclose_f(intset, error)  


      ! tmax
      call h5dopen_f(group_id, 'tmax', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, tmax2, dbldims, error)
      call h5dclose_f(dblset, error)


      ! trstrt
      call h5dopen_f(group_id, 'trstrt', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, trstrt, intdims, error)
      call h5dclose_f(intset, error) 


      ! nstep
      call h5dopen_f(group_id, 'nstep', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, nstep, intdims, error)
      call h5dclose_f(intset, error) 


      ! tout
      call h5dopen_f(group_id, 'tout', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, tout2, dbldims, error)
      call h5dclose_f(dblset, error)


      ! cfl
      call h5dopen_f(group_id, 'cfl', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, cfl2, dbldims, error)
      call h5dclose_f(dblset, error)
      

      ! dtini
      call h5dopen_f(group_id, 'dtini', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, dtini, dbldims, error)
      call h5dclose_f(dblset, error)


      ! dt
      call h5dopen_f(group_id, 'dt', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, dt, dbldims, error)
      call h5dclose_f(dblset, error)


      ! timej
      call h5dopen_f(group_id, 'timej', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, timej, dbldims, error)
      call h5dclose_f(dblset, error)



      ! small
      call h5dopen_f(group_id, 'small', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, small, dbldims, error)
      call h5dclose_f(dblset, error)


      ! smlrho
      call h5dopen_f(group_id, 'smlrho', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, smlrho, dbldims, error)
      call h5dclose_f(dblset, error)      


      ! smallp
      call h5dopen_f(group_id, 'smallp', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, smallp, dbldims, error)
      call h5dclose_f(dblset, error)


      ! smallu
      call h5dopen_f(group_id, 'smallu', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, smallu, dbldims, error)
      call h5dclose_f(dblset, error)

      
      ! gridlx
      call h5dopen_f(group_id, 'gridlx', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gridlx, dbldims, error)
      call h5dclose_f(dblset, error)     
      
      ! gridly
      call h5dopen_f(group_id, 'gridly', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gridly, dbldims, error)
      call h5dclose_f(dblset, error)     
      
      ! gridlz
      call h5dopen_f(group_id, 'gridlz', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gridlz, dbldims, error)
      call h5dclose_f(dblset, error)  

      
      ! bndmnx
      call h5dopen_f(group_id, 'bndmnx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmnx, intdims, error)
      call h5dclose_f(intset, error)  
      
      ! bndmxx
      call h5dopen_f(group_id, 'bndmxx', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmxx, intdims, error)
      call h5dclose_f(intset, error)  
      
      ! bndmny
      call h5dopen_f(group_id, 'bndmny', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmny, intdims, error)
      call h5dclose_f(intset, error)  
      
      ! bndmxy
      call h5dopen_f(group_id, 'bndmxy', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmxy, intdims, error)
      call h5dclose_f(intset, error)  

      ! bndmnz
      call h5dopen_f(group_id, 'bndmnz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmnz, intdims, error)
      call h5dclose_f(intset, error)  
      
      ! bndmxz
      call h5dopen_f(group_id, 'bndmxz', intset, error) 
      call h5dread_f(intset, H5T_NATIVE_INTEGER, bndmxz, intdims, error)
      call h5dclose_f(intset, error)  
      
   
      ! gamma_ad
      call h5dopen_f(group_id, 'gamma_ad', dblset, error) 
      call h5dread_f(dblset, H5T_NATIVE_DOUBLE, gamma_ad, dbldims, error)
      call h5dclose_f(dblset, error) 


           !   ----------
           !   X
           !   ----------
      
      call h5dopen_f(group_id, 'xzn', dxset, error)
      
#ifdef PARALELO
      call h5dget_space_f(dxset, dxfilespace, error)
      call h5sselect_hyperslab_f(dxfilespace, H5S_SELECT_SET_F, dxoffset, dxdims, error)
#endif


     
#ifdef PARALELO
      call h5dread_f(dxset, H5T_NATIVE_DOUBLE, dbufx(imin:imax), dxdims, error &
      ,file_space_id = dxfilespace, mem_space_id = dxmemspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(dxset, H5T_NATIVE_DOUBLE, dbufx(imin:imax), dxdims, error)
#endif

      
      call h5dclose_f(dxset, error)
#ifdef PARALELO
      call h5sclose_f(dxfilespace, error)
#endif

      do i = imin,imax
         xzn(i) = dbufx(i)
      enddo


      if (nsdim.gt.2) then
         

           !   ----------
           !   Z
           !   ----------

         
         call h5dopen_f(group_id, 'zzn', dzset, error)

#ifdef PARALELO
         call h5dget_space_f(dzset, dzfilespace, error)
         call h5sselect_hyperslab_f(dzfilespace, H5S_SELECT_SET_F, dzoffset, dzdims, error)
#endif
         
#ifdef PARALELO
         call h5dread_f(dzset, H5T_NATIVE_DOUBLE, dbufz(kmin:kmax), dzdims, error &
         ,file_space_id = dzfilespace, mem_space_id = dzmemspace &
         ,xfer_prp = data_transfer)
#else
         call h5dread_f(dzset, H5T_NATIVE_DOUBLE, dbufz(kmin:kmax), dzdims, error)
#endif
         call h5dclose_f(dzset, error)
#ifdef PARALELO
         call h5sclose_f(dzfilespace, error)
#endif
      endif

      do k = kmin,kmax
         zzn(k) = dbufz(k)
      enddo


      if (nsdim.gt.1) then

         
           !   ----------
           !   Y
           !   ----------
         
         call h5dopen_f(group_id, 'yzn', dyset, error)

#ifdef PARALELO
         call h5dget_space_f(dyset, dyfilespace, error)
         call h5sselect_hyperslab_f(dyfilespace, H5S_SELECT_SET_F, dyoffset, dydims, error)        
#endif


#ifdef PARALELO
         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy(jmin:jmax), dydims, error &
         ,file_space_id = dyfilespace, mem_space_id = dymemspace &
         ,xfer_prp = data_transfer)
#else
         call h5dread_f(dyset, H5T_NATIVE_DOUBLE, dbufy(jmin:jmax), dydims, error)
#endif        
         call h5dclose_f(dyset, error)
#ifdef PARALELO
         call h5sclose_f(dyfilespace, error)
#endif

         
      do j = jmin,jmax
         yzn(j) = dbufy(j)
      enddo

	endif


!     Densty
      
      call h5dopen_f(group_id, 'densty', d3set, error)
      
#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)      

      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%densty = dbuf(i,j,k)
            enddo
         enddo
      enddo



!     denstye
      call h5dopen_f(group_id, 'denstye', d3set, error)
#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)      

      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%denstye = dbuf(i,j,k)
            enddo
         enddo
      enddo



!     pres
      call h5dopen_f(group_id, 'pres', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif

      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%pres = dbuf(i,j,k)
            enddo
         enddo
      enddo


!     velx
      call h5dopen_f(group_id, 'velx', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%velx = dbuf(i,j,k)
            enddo
         enddo
      enddo



      if (nsdim.gt.1) then
         
!     vely
         
         call h5dopen_f(group_id, 'vely', d3set, error)

#ifdef PARALELO
         call h5dget_space_f(d3set, d3filespace, error)
         call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
         
#endif
#ifdef PARALELO
         call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
         ,file_space_id = d3filespace, mem_space_id = d3memspace &
         ,xfer_prp = data_transfer)
#else
         call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
         call h5dclose_f(d3set, error)
#ifdef PARALELO
         call h5sclose_f(d3filespace, error)
#endif
      endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%vely = dbuf(i,j,k)
            enddo
         enddo
      enddo


      if (nsdim.gt.2) then
         
!     velz
         call h5dopen_f(group_id, 'velz', d3set, error)

#ifdef PARALELO
         call h5dget_space_f(d3set, d3filespace, error)
         call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
         
#endif
#ifdef PARALELO
         call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
         ,file_space_id = d3filespace, mem_space_id = d3memspace &
         ,xfer_prp = data_transfer)
#else
         call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
         call h5dclose_f(d3set, error)
#ifdef PARALELO
         call h5sclose_f(d3filespace, error)
#endif
      endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%velz = dbuf(i,j,k)
            enddo
         enddo
      enddo


!     tracer
      call h5dopen_f(group_id, 'tracer', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%tracer = dbuf(i,j,k)
            enddo
         enddo
      enddo



!     Bxc
      call h5dopen_f(group_id, 'bxc', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%bxc = dbuf(i,j,k)
            enddo
         enddo
      enddo

  
!     gravy
      call h5dopen_f(group_id, 'byc', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%byc = dbuf(i,j,k)
            enddo
         enddo
      enddo


!     bzc
      call h5dopen_f(group_id, 'bzc', d3set, error)

#ifdef PARALELO
      call h5dget_space_f(d3set, d3filespace, error)
      call h5sselect_hyperslab_f(d3filespace, H5S_SELECT_SET_F, d3offset, d3dims, error)
      
#endif
#ifdef PARALELO
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error &
      ,file_space_id = d3filespace, mem_space_id = d3memspace &
      ,xfer_prp = data_transfer)
#else
      call h5dread_f(d3set, H5T_NATIVE_DOUBLE, dbuf(imin:imax,jmin:jmax,kmin:kmax), d3dims, error)
#endif
      call h5dclose_f(d3set, error)
#ifdef PARALELO
      call h5sclose_f(d3filespace, error)
#endif
      do k = kmin,kmax
         do j = jmin,jmax
            do i = imin,imax
               physic(i, j, k)%bzc = dbuf(i,j,k)
            enddo
         enddo
      enddo


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
      if (nsdim.gt.1) then
         call h5sclose_f(dymemspace, error)
      endif
      call h5sclose_f(dxmemspace, error)
      call h5sclose_f(d3memspace, error)
#endif /* PARALELO */

!     close total memory spaces
      if (nsdim.gt.1) then
         call h5sclose_f(dyspace, error)
      endif
      call h5sclose_f(dxspace, error)
      call h5sclose_f(d3space, error)

!     create memory spaces for header
      call h5sclose_f(intspace, error)
      call h5sclose_f(dblspace, error)

!     WE WILL USE ONLY ONE GROUP FOR THE TIME BEING. PERUCHO      
!     curR_R = curR_R + 1

!     close file         
      call h5fclose_f(file_id, error)
      
!     prepare for writing a new file
!     if (curR_R.gt.RECORDER_R) then
      if (cartrank == 0) then 
         print*, &
         "Closing file "//TRIM(outRoot)//rstfil//"-" &
         //c1//c2//c3//c4//c5//".h5"
      endif
      R_filer = R_filer + 1
!     R_filed = R_filed + 1
      H_openr = 0
      H_opend = 0
!     else
!     H_openr = 1
!     endif



            
      RETURN


END SUBROUTINE RST_HDF5
