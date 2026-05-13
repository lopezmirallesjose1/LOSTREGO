!*********************************
!     3D-RMHD CODE: LOSTREGO
!*********************************

! NAME:
!       MODULOS

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
!       Define the modules to be included in the main program and other subroutines.

! ACADEMIC REFERENCE:


      
# include "types.h"   
   
!===================================MODULE======================================
! NAME                :: math_constants
! F90 SPECIFICATION   ::
! DESCRIPTION         :: define mathematical constants.
! SEE ALSO            ::
!===============================================================================

MODULE math_constants
        real(rs), parameter :: pi = 3.14159265358979323844
END MODULE math_constants


!===================================MODULE======================================
! NAME                :: physic_constants
! F90 SPECIFICATION   ::
! DESCRIPTION         :: define and print physical constants.
! SEE ALSO            ::
!===============================================================================

MODULE physic_constants

        ! ===================
          CONTAINS
        ! ===================

! ===================================================================================================      
        SUBROUTINE SHOW_PHYS_CONST(pc, me_cgs, c, kb, mp_cgs, KB1, KB2, meff)
! ===================================================================================================

                !------------------------------------------ Input Variables
                ! ...
                !------------------------------------------ Output Variables

	        real(rs) :: pc, me_cgs, c, kb, mp_cgs, KB1, KB2, meff               
 
                !------------------------------------------ Local Variables
        
                ! ...

                !------------------------------------------------------------------------- BEGIN

                pc = 1 / 3.24254215D-19   
                me_cgs = 9.1095D-28
                c = 2.99792458D10
                kb = 1.38066D-16
                mp_cgs = 1836.11 * me_cgs
                KB1 = 7.5D20
                KB2 = 4.4D-10
                meff= 1.0/9.1095D-28 * 2.0/1836.11

	        
                !------------------------------------------------------------------------- END

! ===================================================================================================
        END SUBROUTINE SHOW_PHYS_CONST
! ===================================================================================================      

END MODULE physic_constants

      
!===================================MODULE======================================
! NAME                :: tipos
! F90 SPECIFICATION   ::
! DESCRIPTION         :: define 'type' structures.
! SEE ALSO            ::
!===============================================================================

MODULE tipos
    
        ! X-Boundary coordinates.
        type TboundPosX
            real(rs) :: xzn
            real(rs) :: xznl
            real(rs) :: xznr
            real(rs) :: dx
        end type TboundPosX

        ! Y-Boundary coordinates.
        type TboundPosY
            real(rs) :: yzn
            real(rs) :: yznl
            real(rs) :: yznr
            real(rs) :: dx
        end type TboundPosY

        ! Z-Boundary coordinates.
        type TboundPosZ
            real(rs) :: zzn
            real(rs) :: zznl
            real(rs) :: zznr
            real(rs) :: dx
        end type TboundPosZ

        ! Boundary physical variables.
        type TboundPhysic
            real(rs) :: velx, vely, velz
            real(rs) :: bxc, byc, bzc
            real(rs) :: densty, denstye 
            real(rs) :: eps
            real(rs) :: pres
            real(rs) :: tracer
        end type TboundPhysic



        type TboundField
            real(rs) :: velx, vely, velz
            real(rs) :: bxc, byc, bzc
        end type TboundField



	! Boundary fluxes.
        type TboundCTxy
            real(rs) :: fbxy, rflxy
        end type TboundCTxy

        type TboundCTxz
            real(rs) :: fbxz, rflxz
        end type TboundCTxz

        type TboundCTyx
            real(rs) :: fbyx, rflxx
        end type TboundCTyx

        type TboundCTyz
            real(rs) :: fbyz, rflxz
        end type TboundCTyz

        type TboundCTzx
            real(rs) :: fbzx, rflxx
        end type TboundCTzx

        type TboundCTzy
            real(rs) :: fbzy, rflxy
        end type TboundCTzy


        ! Physical variables (.....RESTRT).
        type TPlanePhysic
            real(rs) :: velx, vely, velz
            real(rs) :: bxc, byc, bzc
            real(rs) :: densty, denstye 
            real(rs) :: eps
            real(rs) :: pres
            real(rs) :: tracer
        end type TPlanePhysic


        
        ! Physical and relativistic variables.
        type Tphysic
            real(rs) :: velx, vely, velz
            real(rs) :: ovelx, ovely, ovelz 
            real(rs) :: momenx, momeny, momenz
            real(rs) :: omomenx, omomeny, omomenz
            real(rs) :: densty, denstye
            real(rs) :: odensty, odenstye
            real(rs) :: eps
            real(rs) :: oeps
            real(rs) :: pres
            real(rs) :: opres
            real(rs) :: denstw, denstew
            real(rs) :: odenstw, odenstew
            real(rs) :: energy
            real(rs) :: oenergy
            real(rs) :: tracer, tracer2
            real(rs) :: otracer
            real(rs) :: bxc, byc, bzc, bxc2, byc2, bzc2 
            real(rs) :: obxc, obyc, obzc
            real(rs) :: gravx, gravy, gravz
	    real(rs) :: bxzn, byzn, bzzn, bxct, byct, bzct
        end type Tphysic

	type ConsTransp
	    real(rs) :: bx, by, bz, obx, oby, obz
	    real(rs) :: fbyx,fbzx,fbxy,fbzy,fbxz,fbyz
	    real(rs) :: omega_x, omega_y, omega_z
	    real(rs) :: rflxx, rflxy, rflxz
	end type ConsTransp

        type TstepInfoReal
        real(rs) :: dt
        real(rs) :: sound
        real(rs) :: v
        end type TstepInfoReal
        
        type TstepInfoInt
        integer(is) :: ic
        integer(is) :: jc
        integer(is) :: kc
        end type TstepInfoInt

        ! ===================
        ! BEGIN :: MEMORY OPTIMIZATION
        ! ===================
        
        ! One-dimensional working arrays (wa) 
        type Twaspace

           ! The definition of "Working Arrays (wa)" depends on the particular script 
           ! where they are explicitly defined. The objective of this formulation 
           ! is to improve the ordering of the definitions attending future modifications.

           !====================
           ! 1D - vector
           !====================
           
        
           real(rs), pointer :: wa1init100(:), wa1init101(:), wa1init102(:), &
                                wa1init103(:), wa1init104(:), wa1init105(:), &
                                wa1init106(:), wa1init107(:), wa1init108(:), &
				wa1init109(:)

           real(rs), pointer :: wa1rest100(:), wa1rest101(:), wa1rest102(:), &
                                wa1rest103(:), wa1rest104(:), wa1rest105(:), &
                                wa1rest106(:), wa1rest107(:), wa1rest108(:), &
			        wa1rest109(:)

           real(rs), pointer :: wa1get100(:), wa1get101(:), wa1get102(:), &
                                wa1get103(:), wa1get104(:), wa1get105(:), &
                                wa1get106(:), wa1get107(:), wa1get108(:), &
                                wa1get109(:)       

                                
           !====================
           ! -4 : mn5 vector
           !====================

           real(rs), pointer :: wa5plm100(:)
                
           real(rs), pointer :: wa5eos100(:), wa5eos101(:), wa5eos102(:), &
                                wa5eos103(:), wa5eos104(:), wa5eos105(:), &
                                wa5eos106(:), wa5eos107(:), wa5eos108(:), &
                                wa5eos109(:), wa5eos110(:)


           real(rs), pointer :: wa5bis100(:), wa5bis101(:), wa5bis102(:), &
                                wa5bis103(:), wa5bis104(:), wa5bis105(:), &
                                wa5bis106(:), wa5bis107(:), wa5bis108(:), &
                                wa5bis109(:)


           real(rs), pointer :: wa5get100(:), wa5get101(:), wa5get102(:), &
                                wa5get103(:), wa5get104(:), wa5get105(:), &
                                wa5get106(:), wa5get107(:), wa5get108(:), &
                                wa5get109(:), wa5get110(:), wa5get111(:), &
                                wa5get112(:), wa5get113(:), wa5get114(:), &
                                wa5get115(:), wa5get116(:), wa5get117(:), &
                                wa5get118(:), wa5get119(:), wa5get120(:), &
                                wa5get121(:), wa5get122(:), wa5get123(:)

           real(rs), pointer :: wa5hlle100(:), wa5hlle101(:), wa5hlle102(:), &
                                wa5hlle103(:), wa5hlle104(:), wa5hlle105(:), &
                                wa5hlle106(:), wa5hlle107(:), wa5hlle108(:), &
                                wa5hlle109(:), wa5hlle110(:), wa5hlle111(:), &
                                wa5hlle112(:), wa5hlle113(:), wa5hlle114(:), &
                                wa5hlle115(:), wa5hlle116(:), wa5hlle117(:), &
                                wa5hlle118(:), wa5hlle119(:), wa5hlle120(:), &
                                wa5hlle121(:), wa5hlle122(:), wa5hlle123(:), &
                                wa5hlle124(:), wa5hlle125(:), wa5hlle126(:), &
                                wa5hlle127(:), wa5hlle128(:), wa5hlle129(:), &
                                wa5hlle130(:), wa5hlle131(:), wa5hlle132(:), &
                                wa5hlle133(:), wa5hlle134(:), wa5hlle135(:), &
                                wa5hlle136(:), wa5hlle137(:), wa5hlle138(:), &
                                wa5hlle139(:), wa5hlle140(:), wa5hlle141(:), &
                                wa5hlle142(:)

           real(rs), pointer :: wa5hlld100(:), wa5hlld101(:), wa5hlld102(:), &
                                wa5hlld103(:), wa5hlld104(:), wa5hlld105(:), &
                                wa5hlld106(:), wa5hlld107(:), wa5hlld108(:), &
                                wa5hlld109(:), wa5hlld110(:), wa5hlld111(:), &
                                wa5hlld112(:), wa5hlld113(:), wa5hlld114(:), &
                                wa5hlld115(:), wa5hlld116(:), wa5hlld117(:), &
                                wa5hlld118(:), wa5hlld119(:), wa5hlld120(:), &
                                wa5hlld121(:), wa5hlld122(:), wa5hlld123(:), &
                                wa5hlld124(:), wa5hlld125(:), wa5hlld126(:), &
                                wa5hlld127(:), wa5hlld128(:), wa5hlld129(:), &
                                wa5hlld130(:), wa5hlld131(:), wa5hlld132(:)
	   integer(is), pointer :: 	       wa5hlld133(:)

                                
           real(rs), pointer :: wa5mhydrow100(:), wa5mhydrow101(:), wa5mhydrow102(:), &
                                wa5mhydrow103(:), wa5mhydrow104(:), wa5mhydrow105(:), &
                                wa5mhydrow106(:), wa5mhydrow107(:), wa5mhydrow108(:), &
                                wa5mhydrow109(:), wa5mhydrow110(:), wa5mhydrow111(:), &
                                wa5mhydrow112(:), wa5mhydrow113(:), wa5mhydrow114(:), &
                                wa5mhydrow115(:), wa5mhydrow116(:), wa5mhydrow117(:), &
                                wa5mhydrow118(:), wa5mhydrow119(:), wa5mhydrow120(:), &
                                wa5mhydrow121(:), wa5mhydrow122(:), wa5mhydrow123(:)
                                
	   integer(is), pointer :: wa5mhydrow124(:)

	   real(rs), pointer :: wa5mhydrow125(:) 
          	

           real(rs), pointer :: wa5sweeps100(:), wa5sweeps101(:), wa5sweeps102(:), &
                                wa5sweeps103(:), wa5sweeps104(:), wa5sweeps105(:), &
                                wa5sweeps106(:), wa5sweeps107(:), wa5sweeps108(:), &
                                wa5sweeps109(:), wa5sweeps110(:), wa5sweeps111(:), &
                                wa5sweeps112(:), wa5sweeps113(:), wa5sweeps114(:), &
                                wa5sweeps115(:), wa5sweeps116(:), wa5sweeps117(:), &
                                wa5sweeps118(:), wa5sweeps119(:), wa5sweeps120(:)


           real(rs), pointer :: wa5tstep100(:), wa5tstep101(:), wa5tstep102(:), &
                                wa5tstep103(:), wa5tstep104(:), wa5tstep105(:), &
                                wa5tstep106(:), wa5tstep107(:), wa5tstep108(:), &
                                wa5tstep109(:), wa5tstep110(:), wa5tstep111(:), &
                                wa5tstep112(:), wa5tstep113(:), wa5tstep114(:), &
                                wa5tstep115(:), wa5tstep116(:), wa5tstep117(:), &
                                wa5tstep118(:), wa5tstep119(:), wa5tstep120(:), &
                                wa5tstep121(:), wa5tstep122(:), wa5tstep123(:), &
                                wa5tstep124(:), wa5tstep125(:), wa5tstep126(:), &
                                wa5tstep127(:), wa5tstep128(:), wa5tstep129(:), &
                                wa5tstep130(:), wa5tstep131(:), wa5tstep132(:), &
                                wa5tstep133(:), wa5tstep134(:), wa5tstep135(:), &
                                wa5tstep136(:), wa5tstep137(:), wa5tstep138(:)

 	   integer(is), pointer :: wa5tstep139(:)	
	   real(rs), pointer    :: wa5tstep140(:), wa5tstep141(:)
                               
           !====================
           ! -4 : mn6 vector
           !====================


           real(rs), pointer :: wa6mhydrow100(:), wa6mhydrow101(:), wa6mhydrow102(:), &
                                wa6mhydrow103(:), wa6mhydrow104(:), wa6mhydrow105(:), &
                                wa6mhydrow106(:), wa6mhydrow107(:), wa6mhydrow108(:), &
                                wa6mhydrow109(:), wa6mhydrow110(:), wa6mhydrow111(:), &
                                wa6mhydrow112(:), wa6mhydrow113(:), wa6mhydrow114(:), &
                                wa6mhydrow115(:), wa6mhydrow116(:), wa6mhydrow117(:), &
                                wa6mhydrow118(:), wa6mhydrow119(:), wa6mhydrow120(:), &
                                wa6mhydrow121(:), wa6mhydrow122(:), wa6mhydrow123(:), &
                                wa6mhydrow124(:), wa6mhydrow125(:), wa6mhydrow126(:), &
                                wa6mhydrow127(:), wa6mhydrow128(:), wa6mhydrow129(:), &
                                wa6mhydrow130(:), wa6mhydrow131(:), wa6mhydrow132(:), &
                                wa6mhydrow133(:), wa6mhydrow134(:), wa6mhydrow135(:), &
                                wa6mhydrow136(:), wa6mhydrow137(:), wa6mhydrow138(:), &
                                wa6mhydrow139(:), wa6mhydrow140(:), wa6mhydrow141(:), &
                                wa6mhydrow142(:), wa6mhydrow143(:), wa6mhydrow144(:), &
                                wa6mhydrow145(:), wa6mhydrow146(:), wa6mhydrow147(:), & 
                                wa6mhydrow148(:), wa6mhydrow149(:), wa6mhydrow150(:), &
                                wa6mhydrow151(:), wa6mhydrow152(:), wa6mhydrow153(:), &
                                wa6mhydrow154(:), wa6mhydrow155(:), wa6mhydrow156(:), &
                                wa6mhydrow157(:), wa6mhydrow158(:)

           real(rs), pointer :: wa6sweeps100(:), wa6sweeps101(:), wa6sweeps102(:), &
                                wa6sweeps103(:), wa6sweeps104(:), wa6sweeps105(:), &
                                wa6sweeps106(:), wa6sweeps107(:), wa6sweeps108(:)

           real(rs), pointer :: wa6tstep101(:), wa6tstep102(:),                 &
                                wa6tstep103(:), wa6tstep104(:), wa6tstep105(:), &
                                wa6tstep106(:), wa6tstep107(:), wa6tstep108(:), &
                                wa6tstep109(:), wa6tstep110(:), wa6tstep111(:), &
                                wa6tstep112(:), wa6tstep113(:), wa6tstep114(:), &
                                wa6tstep115(:), wa6tstep116(:), wa6tstep117(:), &
                                wa6tstep118(:), wa6tstep119(:), wa6tstep120(:), &
                                wa6tstep121(:), wa6tstep122(:), wa6tstep123(:), &
                                wa6tstep124(:), wa6tstep125(:), wa6tstep126(:), &
                                wa6tstep127(:), wa6tstep128(:), wa6tstep129(:), &
                                wa6tstep130(:), wa6tstep131(:), wa6tstep132(:), &
                                wa6tstep133(:), wa6tstep134(:), wa6tstep135(:), &
                                wa6tstep136(:), wa6tstep137(:), wa6tstep138(:), &
                                wa6tstep139(:), wa6tstep140(:), wa6tstep141(:), &
                                wa6tstep142(:), wa6tstep143(:), wa6tstep144(:), &
                                wa6tstep145(:), wa6tstep146(:), wa6tstep147(:), & 
                                wa6tstep148(:), wa6tstep149(:), wa6tstep150(:), &
                                wa6tstep151(:), wa6tstep152(:)                             

           real(rs), pointer :: wa6init100(:), wa6init101(:)

           real(rs), pointer :: wa6rest100(:), wa6rest101(:)

           real(rs), pointer :: wa6MP100(:), wa6MP101(:), wa6MP102(:), &
                                wa6MP103(:), wa6MP104(:)

           
           !====================
           ! 3D - vector
           !====================

           
           !mnx1,mny1,mnz1 arrays for tracer evolution
           real(rs), pointer :: wa3sweeps100(:,:,:)

	   real(rs), pointer :: wa3tstep100(:,:,:)

           !0:mnx1, 0:mny1, 0:mnz1 arrays for CT-magnetic fluxes
           real(rs), pointer :: wa3sweeps101(:,:,:), wa3sweeps102(:,:,:), wa3sweeps103(:,:,:),     &
                                wa3sweeps104(:,:,:), wa3sweeps105(:,:,:), wa3sweeps106(:,:,:),     &
                                wa3sweeps107(:,:,:), wa3sweeps108(:,:,:), wa3sweeps109(:,:,:)

           
        end type Twaspace

        ! ===================
        ! END :: MEMORY OPTIMIZATION
        ! ===================
        
END MODULE tipos


!===================================MODULE======================================
! NAME                :: memoria
! F90 SPECIFICATION   ::
! DESCRIPTION         :: establish memory attributes.
! SEE ALSO            ::
!===============================================================================
      
MODULE memoria

	USE tipos

        integer(is), parameter :: maxVars = 200

        integer(is) :: nbvars
        data           nbvars  /0/

        real(rs) :: maxMem, usedMem
        data        maxMem, usedMem  /0.0,0.0/

        type Tvar
            character(32) :: vname
            integer(is)   :: esize
            integer(is)   :: bsize
        end type Tvar

        type(Tvar) :: vtable(maxVars)

        ! TboundPhysic arrays used in intercambiaBoundPhysic
        type(TboundPhysic), pointer :: boundsXs1(:,:,:), boundsYs1(:,:,:),boundsZs1(:,:,:) 
        type(TboundPhysic), pointer :: boundsXr1(:,:,:), boundsYr1(:,:,:),boundsZr1(:,:,:)
        type(TboundPhysic), pointer :: boundsXs2(:,:,:), boundsYs2(:,:,:),boundsZs2(:,:,:) 
        type(TboundPhysic), pointer :: boundsXr2(:,:,:), boundsYr2(:,:,:),boundsZr2(:,:,:)

        type(TboundCTxy), pointer :: CTboundsXs1xy(:,:,:)
        type(TboundCTxy), pointer :: CTboundsXr1xy(:,:,:)
        type(TboundCTxy), pointer :: CTboundsXs2xy(:,:,:)
        type(TboundCTxy), pointer :: CTboundsXr2xy(:,:,:)

        type(TboundCTyz), pointer :: CTboundsYs1yz(:,:,:)
        type(TboundCTyz), pointer :: CTboundsYr1yz(:,:,:)
        type(TboundCTyz), pointer :: CTboundsYs2yz(:,:,:)
        type(TboundCTyz), pointer :: CTboundsYr2yz(:,:,:)

        type(TboundCTzx), pointer :: CTboundsZs1zx(:,:,:)
        type(TboundCTzx), pointer :: CTboundsZr1zx(:,:,:)
        type(TboundCTzx), pointer :: CTboundsZs2zx(:,:,:)
        type(TboundCTzx), pointer :: CTboundsZr2zx(:,:,:)

        type(TboundCTxz), pointer :: CTTboundsXs1xz(:,:,:)                                    
        type(TboundCTxz), pointer :: CTTboundsXr1xz(:,:,:)                                       
        type(TboundCTxz), pointer :: CTTboundsXs2xz(:,:,:)                                       
        type(TboundCTxz), pointer :: CTTboundsXr2xz(:,:,:)                                       

        type(TboundCTyx), pointer :: CTTboundsYs1yx(:,:,:)                   
        type(TboundCTyx), pointer :: CTTboundsYr1yx(:,:,:)                   
        type(TboundCTyx), pointer :: CTTboundsYs2yx(:,:,:)                   
        type(TboundCTyx), pointer :: CTTboundsYr2yx(:,:,:)                   

        type(TboundCTzy), pointer :: CTTboundsZs1zy(:,:,:)
        type(TboundCTzy), pointer :: CTTboundsZr1zy(:,:,:)
        type(TboundCTzy), pointer :: CTTboundsZs2zy(:,:,:)
        type(TboundCTzy), pointer :: CTTboundsZr2zy(:,:,:)

        ! ===================
          CONTAINS
        ! ===================      

! ===================================================================================================
        SUBROUTINE MEMPUSH (nele, elesize, var)
! ===================================================================================================
          
                implicit none
        
                !------------------------------------------ Input Variables

                integer(is)  :: nele, elesize
                character(*) :: var

                !------------------------------------------ Output Variables

                ! ...

                !------------------------------------------ Local Variables

                character(16) :: message1

                !------------------------------------------------------------------------- BEGIN
        
                usedMem = usedMem + REAL(nele) * REAL(elesize)
                maxMem = MAX (maxMem, usedMem)
                
                IF (nbvars .LT. maxVars) THEN
                   
                   nbvars = nbvars + 1
                   vtable(nbvars)%vname = var
                   vtable(nbvars)%esize = nele
                   vtable(nbvars)%bsize = elesize
                   
                ELSE
                   
                   message1='Increase maxVars'
                   call ERROR ('MEMPUSH', message1, __FILE__, __LINE__)
                   
                ENDIF
      
                !--------------------------------------------------------------------------- END

! ===================================================================================================
        END SUBROUTINE MEMPUSH
! ===================================================================================================

        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------

! ===================================================================================================      
        SUBROUTINE MEMPOP (var)
! ===================================================================================================
      
                implicit none
        
                !------------------------------------------ Input Variables
        
                character(*) :: var

                !------------------------------------------ Output Variables

                ! ...
                
                !------------------------------------------ Local Variables
                
                integer(is) :: ii
		character(27) :: message1

                !------------------------------------------------------------------------- BEGIN

                IF (nbvars .GT. 0) THEN
                   
                   ii= 1

                   ! LGT, LLT uses ASCII ordering.
                   DO WHILE ((LGT(vtable(ii)%vname,var) .OR. LLT(vtable(ii)%vname,var)) .AND. &
                              ii .LE. nbvars)
                      ii= ii+1
                   ENDDO
                   
                   IF (ii .LE. nbvars) THEN
                      
                      usedMem = usedMem - REAL(vtable(ii)%esize) * REAL(vtable(ii)%bsize)

                      ! We fill the empty space in the table with the last element.
                      
                      IF (nbvars .GT. 1) THEN
                         vtable(ii) = vtable(nbvars)
                      ENDIF
                      
                      nbvars = nbvars - 1
                      
                   ELSE
                      
                      message1='Variable not found in table'
                      call ERROR ('MEMPOP', message1, __FILE__, __LINE__)
                   
                   ENDIF
                   
                ENDIF
                
                !--------------------------------------------------------------------------- END
 
! ===================================================================================================               
        END SUBROUTINE MEMPOP
! ===================================================================================================
      
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------

! ===================================================================================================      
        SUBROUTINE SHOWMEM ( )
! ===================================================================================================
      
                implicit none

                !------------------------------------------ Input Variables

                ! ...

                !------------------------------------------ Output Variables

                ! ...

                !------------------------------------------ Local Variables
                
                integer(is)   :: ii
                integer(is*2) :: kk

                !------------------------------------------------------------------------- BEGIN
                
                IF (nbvars .GT. 0) THEN
                   
                   write(*,*) '------------ Memory report -------------'
                   write(*,*) ' Number of allocated variables: ', nbvars
                   
                   DO ii= 1, nbvars
                      
                      kk = (vtable(ii)%esize*vtable(ii)%bsize)
                      write(*,*) ii, vtable(ii)%esize, vtable(ii)%bsize, kk, &
                                 TRIM(vtable(ii)%vname)
                      
                   ENDDO
                   
                   write(*,*) '----------------------------------------'
                   
                ENDIF
      
                !--------------------------------------------------------------------------- END

! ===================================================================================================      
        END SUBROUTINE SHOWMEM
! ===================================================================================================
      
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------

! ===================================================================================================      
        SUBROUTINE GETBOUNDSMEMORYSPACE (nx, ny, nz)
! ===================================================================================================

                implicit none

                !------------------------------------------ Input Variables

                integer(is) :: nx, ny, nz

                !------------------------------------------ Output Variables

                ! ...
                
                !------------------------------------------ Local Variables

                integer(is) :: ierr

                
                !------------------------------------------------------------------------- BEGIN
 
                allocate (boundsXr1(5,ny,nz), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__) 
                ENDIF

                allocate (boundsYr1(nx,5,nz), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

 
                allocate (boundsZr1(nx,ny,5), stat= ierr)
 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                
                allocate (boundsXr2(5,ny,nz), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                
                allocate (boundsYr2(nx,5,nz), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                
                allocate (boundsZr2(nx,ny,5), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF
 

                allocate (boundsXs1(5,ny,nz), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                
                allocate (boundsYs1(nx,5,nz), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                
                allocate (boundsZs1(nx,ny,5), stat= ierr)
 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                
                allocate (boundsXs2(5,ny,nz), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF 


                allocate (boundsYs2(nx,5,nz), stat= ierr)

                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                
                allocate (boundsZs2(nx,ny,5), stat= ierr)
 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__ ) 
                ENDIF

                !--------------------------------------------------------------------------- END

! ===================================================================================================                
        END SUBROUTINE GETBOUNDSMEMORYSPACE
! ===================================================================================================

        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
    
! ===================================================================================================      
        SUBROUTINE GETBOUNDSMEMORYSPACECT (nx, ny, nz)
! ===================================================================================================

                implicit none

                !------------------------------------------ Input Variables

                integer(is) :: nx, ny, nz

                !------------------------------------------ Output Variables

                ! ...
                
                !------------------------------------------ Local Variables

                integer(is) :: ierr

                
                !------------------------------------------------------------------------- BEGIN
 

                allocate (CTboundsXr1xy(2,ny+1,nz), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__) 
                ENDIF


                allocate (CTboundsYr1yz(nx,2,nz+1), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

 

                allocate (CTboundsZr1zx(nx+1,ny,2), stat= ierr) 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTboundsXr2xy(2,ny+1,nz), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTboundsYr2yz(nx,2,nz+1), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTboundsZr2zx(nx+1,ny,2), stat= ierr) 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF
 


                allocate (CTboundsXs1xy(2,ny+1,nz), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTboundsYs1yz(nx,2,nz+1), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTboundsZs1zx(nx+1,ny,2), stat= ierr) 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTboundsXs2xy(2,ny+1,nz), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF 



                allocate (CTboundsYs2yz(nx,2,nz+1), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTboundsZs2zx(nx+1,ny,2), stat= ierr) 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__ ) 
                ENDIF

                !--------------------------------------------------------------------------- END

! ===================================================================================================                
        END SUBROUTINE GETBOUNDSMEMORYSPACECT
! ===================================================================================================


        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
    
! ===================================================================================================      
        SUBROUTINE GETBOUNDSMEMORYSPACECTT (nx, ny, nz)
! ===================================================================================================

                implicit none

                !------------------------------------------ Input Variables

                integer(is) :: nx, ny, nz

                !------------------------------------------ Output Variables

                ! ...
                
                !------------------------------------------ Local Variables

                integer(is) :: ierr

                
                !------------------------------------------------------------------------- BEGIN
 

                allocate (CTTboundsXr1xz(2,ny,nz+1), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__) 
                ENDIF


                allocate (CTTboundsYr1yx(nx+1,2,nz), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

 

                allocate (CTTboundsZr1zy(nx,ny+1,2), stat= ierr) 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTTboundsXr2xz(2,ny,nz+1), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTTboundsYr2yx(nx+1,2,nz), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTTboundsZr2zy(nx,ny+1,2), stat= ierr) 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF
 


                allocate (CTTboundsXs1xz(2,ny,nz+1), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTTboundsYs1yx(nx+1,2,nz), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTTboundsZs1zy(nx,ny+1,2), stat= ierr) 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTTboundsXs2xz(2,ny,nz+1), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF 



                allocate (CTTboundsYs2yx(nx+1,2,nz), stat= ierr)
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__  ) 
                ENDIF

                

                allocate (CTTboundsZs2zy(nx,ny+1,2), stat= ierr) 
                IF (ierr .NE. 0) THEN
                   call ERROR ('GETBOUNDSMEMORYSPACE',        &
                               'Check allocation of bounds',    &
                               __FILE__, __LINE__ ) 
                ENDIF

                !--------------------------------------------------------------------------- END

! ===================================================================================================                
        END SUBROUTINE GETBOUNDSMEMORYSPACECTT
! ===================================================================================================
      
      
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------

     
! ===================================================================================================
        SUBROUTINE FREEBOUNDSMEMORYSPACE ()
! ===================================================================================================
        
                implicit none

                !------------------------------------------------------------------------- BEGIN

                deallocate (boundsXs1)              
                deallocate (boundsYs1)
                deallocate (boundsZs1)
                deallocate (boundsXs2)
                deallocate (boundsYs2)
                deallocate (boundsZs2)
                deallocate (boundsXr1)
                deallocate (boundsYr1)
                deallocate (boundsZr1)
                deallocate (boundsXr2)
                deallocate (boundsYr2)
                deallocate (boundsZr2)

                !------------------------------------------------------------------------- END

! ===================================================================================================
        END SUBROUTINE FREEBOUNDSMEMORYSPACE
! ===================================================================================================

        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
      
! ===================================================================================================
        SUBROUTINE FREEBOUNDSMEMORYSPACECT ()
! ===================================================================================================
        
                implicit none

                !------------------------------------------------------------------------- BEGIN

                deallocate (CTboundsXs1xy)
                deallocate (CTboundsYs1yz)
                deallocate (CTboundsZs1zx)
                deallocate (CTboundsXs2xy)
                deallocate (CTboundsYs2yz)
                deallocate (CTboundsZs2zx)
                deallocate (CTboundsXr1xy)
                deallocate (CTboundsYr1yz)
                deallocate (CTboundsZr1zx)
                deallocate (CTboundsXr2xy)
                deallocate (CTboundsYr2yz)
                deallocate (CTboundsZr2zx)

                !------------------------------------------------------------------------- END

! ===================================================================================================
        END SUBROUTINE FREEBOUNDSMEMORYSPACECT
! ===================================================================================================

        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
      
! ===================================================================================================
        SUBROUTINE FREEBOUNDSMEMORYSPACECTT ()
! ===================================================================================================
        
                implicit none

                !------------------------------------------------------------------------- BEGIN

                deallocate (CTTboundsXs1xz)
                deallocate (CTTboundsYs1yx)
                deallocate (CTTboundsZs1zy)
                deallocate (CTTboundsXs2xz)
                deallocate (CTTboundsYs2yx)
                deallocate (CTTboundsZs2zy)
                deallocate (CTTboundsXr1xz)
                deallocate (CTTboundsYr1yx)
                deallocate (CTTboundsZr1zy)
                deallocate (CTTboundsXr2xz)
                deallocate (CTTboundsYr2yx)
                deallocate (CTTboundsZr2zy)

                !------------------------------------------------------------------------- END

! ===================================================================================================
        END SUBROUTINE FREEBOUNDSMEMORYSPACECTT
! ===================================================================================================
      
      
END MODULE memoria

      
!===================================MODULE======================================
! NAME                :: parallel
! F90 SPECIFICATION   ::
! DESCRIPTION         :: Declare 'parallel' attributes for PARALELO.F
! SEE ALSO            ::
!===============================================================================

MODULE parallel

        integer(is)              :: sizeTotal(3), utrace, traza
        integer(is), allocatable :: iniTramo(:,:), finTramo(:,:), sizeTramo(:,:)

# ifdef OPENMP
        integer, parameter :: MAXTHREADS = 128
# else
        integer, parameter :: MAXTHREADS = 1
# endif

        ! number of OpenMP threads
        integer :: omp_th, curIter
      
        ! coords of OpenMP threads
        integer, dimension(:), allocatable :: ocurx, ocury, ocurz
      
        ! number of OpenMP threads per dimension
        integer :: onumx, onumy, onumz

        ! MPI dims, coords and periods
        integer :: mpi_dims(3), mpi_coords(3), cartrank, CART_WORLD
        logical :: periods(3)
        
        ! neighbor ranks
        integer :: xleft, xright, yleft, yright, zleft, zright

END MODULE parallel

      
!===================================MODULE======================================
! NAME                :: workarrays
! F90 SPECIFICATION   ::
! DESCRIPTION         :: Declare and get memory for workarrays. 
! SEE ALSO            ::
!===============================================================================

MODULE workarrays
      
        USE tipos
        USE parallel

        ! Command line arguments
        character(len=256) :: outRoot
        character(len=256) :: path_wkd

        ! Save Twaspace as wa(MAXTHREADS)
        type(Twaspace), save :: wa(MAXTHREADS)

        integer(is) :: nele, maxele

        ! TstepInfoReal & TstepInfoInt vectors (size 0:nbproc-1)
        type(TstepInfoReal), pointer :: wass100(:)
        type(TstepInfoInt),  pointer :: wass101(:)


        ! ===================
          CONTAINS
        ! ===================

! ===================================================================================================
        SUBROUTINE GETWAMEMORYSPACE (mn6, mn5, mn1, mnx5, mny5, mnz5, nx, ny, nz,   &
                                     mpi_dims1, mpi_dims2, mpi_dims3)
! ===================================================================================================
                                     
                implicit none

                !------------------------------------------ Input Variables

                integer(is):: mn6, mn5, mn1, mnx5, mny5, mnz5, nx, ny,  &
                              nz, mpi_dims1, mpi_dims2, mpi_dims3


                !------------------------------------------ Output Variables

                ! ...

                !------------------------------------------ Local Variables

                integer(is):: tid, nth, ierr

# ifdef OPENMP
                integer(is):: omp_get_max_threads
# endif

                !------------------------------------------------------------------------- BEGIN

# ifdef OPENMP

                nth = omp_get_max_threads()
                IF (nth .GT. MAXTHREADS) THEN
                   call ERROR ('GETWAMEMMORYSPACE',       &
                               'nth > MAXTHREADS', __FILE__, __LINE__) 
                ENDIF
      
# else
                nth = 1
# endif

                DO tid = 1, nth
                   
                   !====================
                   ! 1D - vector
                   !====================

                   allocate (wa(tid)%wa1init100(-4:mn5), wa(tid)%wa1init101(-4:mn5), wa(tid)%wa1init102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1init103(-4:mn5), wa(tid)%wa1init104(-4:mn5), wa(tid)%wa1init105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1init106(-4:mn5), wa(tid)%wa1init107(-4:mn5), wa(tid)%wa1init108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1init109(-4:mn5), stat=ierr)

                   allocate (wa(tid)%wa1rest100(-4:mn5), wa(tid)%wa1rest101(-4:mn5), wa(tid)%wa1rest102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1rest103(-4:mn5), wa(tid)%wa1rest104(-4:mn5), wa(tid)%wa1rest105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1rest106(-4:mn5), wa(tid)%wa1rest107(-4:mn5), wa(tid)%wa1rest108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1rest109(-4:mn5), stat=ierr)
 

                   allocate (wa(tid)%wa1get100(-4:mn5), wa(tid)%wa1get101(-4:mn5), wa(tid)%wa1get102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1get103(-4:mn5), wa(tid)%wa1get104(-4:mn5), wa(tid)%wa1get105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1get106(-4:mn5), wa(tid)%wa1get107(-4:mn5), wa(tid)%wa1get108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa1get109(-4:mn5), stat=ierr) 


                   IF (ierr .NE. 0) THEN
                      call ERROR ('GETWAMEMORYSPACE',    &
                                  'Allocation of wa1',   &
                                  __FILE__, __LINE__) 
                   ENDIF


                   !====================
                   ! -4 : mn5 vector
                   !====================

                   allocate (wa(tid)%wa5plm100(-4:mn5), stat=ierr)
                
                   allocate (wa(tid)%wa5eos100(-4:mn5), wa(tid)%wa5eos101(-4:mn5), wa(tid)%wa5eos102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5eos103(-4:mn5), wa(tid)%wa5eos104(-4:mn5), wa(tid)%wa5eos105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5eos106(-4:mn5), wa(tid)%wa5eos107(-4:mn5), wa(tid)%wa5eos108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5eos109(-4:mn5), wa(tid)%wa5eos110(-4:mn5), stat=ierr)


                   allocate (wa(tid)%wa5bis100(-4:mn5), wa(tid)%wa5bis101(-4:mn5), wa(tid)%wa5bis102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5bis103(-4:mn5), wa(tid)%wa5bis104(-4:mn5), wa(tid)%wa5bis105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5bis106(-4:mn5), wa(tid)%wa5bis107(-4:mn5), wa(tid)%wa5bis108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5bis109(-4:mn5), stat=ierr)



                   allocate (wa(tid)%wa5get100(-4:mn5), wa(tid)%wa5get101(-4:mn5), wa(tid)%wa5get102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5get103(-4:mn5), wa(tid)%wa5get104(-4:mn5), wa(tid)%wa5get105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5get106(-4:mn5), wa(tid)%wa5get107(-4:mn5), wa(tid)%wa5get108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5get109(-4:mn5), wa(tid)%wa5get110(-4:mn5), wa(tid)%wa5get111(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5get112(-4:mn5), wa(tid)%wa5get113(-4:mn5), wa(tid)%wa5get114(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5get115(-4:mn5), wa(tid)%wa5get116(-4:mn5), wa(tid)%wa5get117(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5get118(-4:mn5), wa(tid)%wa5get119(-4:mn5), wa(tid)%wa5get120(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5get121(-4:mn5), wa(tid)%wa5get122(-4:mn5), wa(tid)%wa5get123(-4:mn5), stat=ierr)

                   allocate (wa(tid)%wa5hlle100(-4:mn5), wa(tid)%wa5hlle101(-4:mn5), wa(tid)%wa5hlle102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle103(-4:mn5), wa(tid)%wa5hlle104(-4:mn5), wa(tid)%wa5hlle105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle106(-4:mn5), wa(tid)%wa5hlle107(-4:mn5), wa(tid)%wa5hlle108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle109(-4:mn5), wa(tid)%wa5hlle110(-4:mn5), wa(tid)%wa5hlle111(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle112(-4:mn5), wa(tid)%wa5hlle113(-4:mn5), wa(tid)%wa5hlle114(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle115(-4:mn5), wa(tid)%wa5hlle116(-4:mn5), wa(tid)%wa5hlle117(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle118(-4:mn5), wa(tid)%wa5hlle119(-4:mn5), wa(tid)%wa5hlle120(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle121(-4:mn5), wa(tid)%wa5hlle122(-4:mn5), wa(tid)%wa5hlle123(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle124(-4:mn5), wa(tid)%wa5hlle125(-4:mn5), wa(tid)%wa5hlle126(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle127(-4:mn5), wa(tid)%wa5hlle128(-4:mn5), wa(tid)%wa5hlle129(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle130(-4:mn5), wa(tid)%wa5hlle131(-4:mn5), wa(tid)%wa5hlle132(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle133(-4:mn5), wa(tid)%wa5hlle134(-4:mn5), wa(tid)%wa5hlle135(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle136(-4:mn5), wa(tid)%wa5hlle137(-4:mn5), wa(tid)%wa5hlle138(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle139(-4:mn5), wa(tid)%wa5hlle140(-4:mn5), wa(tid)%wa5hlle141(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlle142(-4:mn5), stat=ierr)

                   allocate (wa(tid)%wa5hlld100(-4:mn5), wa(tid)%wa5hlld101(-4:mn5), wa(tid)%wa5hlld102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld103(-4:mn5), wa(tid)%wa5hlld104(-4:mn5), wa(tid)%wa5hlld105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld106(-4:mn5), wa(tid)%wa5hlld107(-4:mn5), wa(tid)%wa5hlld108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld109(-4:mn5), wa(tid)%wa5hlld110(-4:mn5), wa(tid)%wa5hlld111(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld112(-4:mn5), wa(tid)%wa5hlld113(-4:mn5), wa(tid)%wa5hlld114(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld115(-4:mn5), wa(tid)%wa5hlld116(-4:mn5), wa(tid)%wa5hlld117(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld118(-4:mn5), wa(tid)%wa5hlld119(-4:mn5), wa(tid)%wa5hlld120(-4:mn5), stat=ierr)                   
                   allocate (wa(tid)%wa5hlld121(-4:mn5), wa(tid)%wa5hlld122(-4:mn5), wa(tid)%wa5hlld123(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld124(-4:mn5), wa(tid)%wa5hlld125(-4:mn5), wa(tid)%wa5hlld126(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld127(-4:mn5), wa(tid)%wa5hlld128(-4:mn5), wa(tid)%wa5hlld129(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5hlld130(-4:mn5), wa(tid)%wa5hlld131(-4:mn5), wa(tid)%wa5hlld132(-4:mn5), stat=ierr) 
		   allocate (wa(tid)%wa5hlld133(-4:mn5), stat=ierr)

                   
                   allocate (wa(tid)%wa5mhydrow100(-4:mn5), wa(tid)%wa5mhydrow101(-4:mn5), wa(tid)%wa5mhydrow102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5mhydrow103(-4:mn5), wa(tid)%wa5mhydrow104(-4:mn5), wa(tid)%wa5mhydrow105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5mhydrow106(-4:mn5), wa(tid)%wa5mhydrow107(-4:mn5), wa(tid)%wa5mhydrow108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5mhydrow109(-4:mn5), wa(tid)%wa5mhydrow110(-4:mn5), wa(tid)%wa5mhydrow111(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5mhydrow112(-4:mn5), wa(tid)%wa5mhydrow113(-4:mn5), wa(tid)%wa5mhydrow114(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5mhydrow115(-4:mn5), wa(tid)%wa5mhydrow116(-4:mn5), wa(tid)%wa5mhydrow117(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5mhydrow118(-4:mn5), wa(tid)%wa5mhydrow119(-4:mn5), wa(tid)%wa5mhydrow120(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5mhydrow121(-4:mn5), wa(tid)%wa5mhydrow122(-4:mn5), wa(tid)%wa5mhydrow123(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5mhydrow124(-4:mn5), wa(tid)%wa5mhydrow125(-4:mn5), stat=ierr)
                   

                   allocate (wa(tid)%wa5sweeps100(-4:mn5), wa(tid)%wa5sweeps101(-4:mn5), wa(tid)%wa5sweeps102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5sweeps103(-4:mn5), wa(tid)%wa5sweeps104(-4:mn5), wa(tid)%wa5sweeps105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5sweeps106(-4:mn5), wa(tid)%wa5sweeps107(-4:mn5), wa(tid)%wa5sweeps108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5sweeps109(-4:mn5), wa(tid)%wa5sweeps110(-4:mn5), wa(tid)%wa5sweeps111(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5sweeps112(-4:mn5), wa(tid)%wa5sweeps113(-4:mn5), wa(tid)%wa5sweeps114(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5sweeps115(-4:mn5), wa(tid)%wa5sweeps116(-4:mn5), wa(tid)%wa5sweeps117(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5sweeps118(-4:mn5), wa(tid)%wa5sweeps119(-4:mn5), wa(tid)%wa5sweeps120(-4:mn6), stat=ierr)
                   

                   allocate (wa(tid)%wa5tstep100(-4:mn5), wa(tid)%wa5tstep101(-4:mn5), wa(tid)%wa5tstep102(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep103(-4:mn5), wa(tid)%wa5tstep104(-4:mn5), wa(tid)%wa5tstep105(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep106(-4:mn5), wa(tid)%wa5tstep107(-4:mn5), wa(tid)%wa5tstep108(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep109(-4:mn5), wa(tid)%wa5tstep110(-4:mn5), wa(tid)%wa5tstep111(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep112(-4:mn5), wa(tid)%wa5tstep113(-4:mn5), wa(tid)%wa5tstep114(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep115(-4:mn5), wa(tid)%wa5tstep116(-4:mn5), wa(tid)%wa5tstep117(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep118(-4:mn5), wa(tid)%wa5tstep119(-4:mn5), wa(tid)%wa5tstep120(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep121(-4:mn5), wa(tid)%wa5tstep122(-4:mn5), wa(tid)%wa5tstep123(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep124(-4:mn5), wa(tid)%wa5tstep125(-4:mn5), wa(tid)%wa5tstep126(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep127(-4:mn5), wa(tid)%wa5tstep128(-4:mn5), wa(tid)%wa5tstep129(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep130(-4:mn5), wa(tid)%wa5tstep131(-4:mn5), wa(tid)%wa5tstep132(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep133(-4:mn5), wa(tid)%wa5tstep134(-4:mn5), wa(tid)%wa5tstep135(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep136(-4:mn5), wa(tid)%wa5tstep137(-4:mn5), wa(tid)%wa5tstep138(-4:mn5), stat=ierr)
                   allocate (wa(tid)%wa5tstep139(-4:mn5), wa(tid)%wa5tstep140(-4:mn5), wa(tid)%wa5tstep141(-4:mn5), stat=ierr)
           
                   IF (ierr .NE. 0) THEN
                      call ERROR ('GETWAMEMORYSPACE',     &
                                  'Allocation of wa5',    &
                                  __FILE__, __LINE__) 
                   ENDIF


                   !====================
                   ! -4 : mn6 vector
                   !====================


                   allocate (wa(tid)%wa6mhydrow100(-4:mn6), wa(tid)%wa6mhydrow101(-4:mn6), wa(tid)%wa6mhydrow102(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow103(-4:mn6), wa(tid)%wa6mhydrow104(-4:mn6), wa(tid)%wa6mhydrow105(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow106(-4:mn6), wa(tid)%wa6mhydrow107(-4:mn6), wa(tid)%wa6mhydrow108(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow109(-4:mn6), wa(tid)%wa6mhydrow110(-4:mn6), wa(tid)%wa6mhydrow111(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow112(-4:mn6), wa(tid)%wa6mhydrow113(-4:mn6), wa(tid)%wa6mhydrow114(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow115(-4:mn6), wa(tid)%wa6mhydrow116(-4:mn6), wa(tid)%wa6mhydrow117(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow118(-4:mn6), wa(tid)%wa6mhydrow119(-4:mn6), wa(tid)%wa6mhydrow120(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow121(-4:mn6), wa(tid)%wa6mhydrow122(-4:mn6), wa(tid)%wa6mhydrow123(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow124(-4:mn6), wa(tid)%wa6mhydrow125(-4:mn6), wa(tid)%wa6mhydrow126(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow127(-4:mn6), wa(tid)%wa6mhydrow128(-4:mn6), wa(tid)%wa6mhydrow129(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow130(-4:mn6), wa(tid)%wa6mhydrow131(-4:mn6), wa(tid)%wa6mhydrow132(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow133(-4:mn6), wa(tid)%wa6mhydrow134(-4:mn6), wa(tid)%wa6mhydrow135(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow136(-4:mn6), wa(tid)%wa6mhydrow137(-4:mn6), wa(tid)%wa6mhydrow138(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow139(-4:mn6), wa(tid)%wa6mhydrow140(-4:mn6), wa(tid)%wa6mhydrow141(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow142(-4:mn6), wa(tid)%wa6mhydrow143(-4:mn6), wa(tid)%wa6mhydrow144(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow145(-4:mn6), wa(tid)%wa6mhydrow146(-4:mn6), wa(tid)%wa6mhydrow147(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow148(-4:mn6), wa(tid)%wa6mhydrow149(-4:mn6), wa(tid)%wa6mhydrow150(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow151(-4:mn6), wa(tid)%wa6mhydrow152(-4:mn6), wa(tid)%wa6mhydrow153(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow154(-4:mn6), wa(tid)%wa6mhydrow155(-4:mn6), wa(tid)%wa6mhydrow156(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6mhydrow157(-4:mn6), wa(tid)%wa6mhydrow158(-4:mn6), stat=ierr)

                   allocate (wa(tid)%wa6sweeps100(-4:mn6), wa(tid)%wa6sweeps101(-4:mn6), wa(tid)%wa6sweeps102(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6sweeps103(-4:mn6), wa(tid)%wa6sweeps104(-4:mn6), wa(tid)%wa6sweeps105(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6sweeps106(-4:mn6), wa(tid)%wa6sweeps107(-4:mn6), wa(tid)%wa6sweeps108(-4:mn6), stat=ierr)

                   allocate (wa(tid)%wa6MP100(-4:mn6), wa(tid)%wa6MP101(-4:mn6), wa(tid)%wa6MP102(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6MP103(-4:mn6), wa(tid)%wa6MP104(-4:mn6), stat=ierr)

                   allocate (wa(tid)%wa6tstep101(-4:mn6), wa(tid)%wa6tstep102(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep103(-4:mn6), wa(tid)%wa6tstep104(-4:mn6), wa(tid)%wa6tstep105(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep106(-4:mn6), wa(tid)%wa6tstep107(-4:mn6), wa(tid)%wa6tstep108(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep109(-4:mn6), wa(tid)%wa6tstep110(-4:mn6), wa(tid)%wa6tstep111(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep112(-4:mn6), wa(tid)%wa6tstep113(-4:mn6), wa(tid)%wa6tstep114(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep115(-4:mn6), wa(tid)%wa6tstep116(-4:mn6), wa(tid)%wa6tstep117(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep118(-4:mn6), wa(tid)%wa6tstep119(-4:mn6), wa(tid)%wa6tstep120(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep121(-4:mn6), wa(tid)%wa6tstep122(-4:mn6), wa(tid)%wa6tstep123(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep124(-4:mn6), wa(tid)%wa6tstep125(-4:mn6), wa(tid)%wa6tstep126(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep127(-4:mn6), wa(tid)%wa6tstep128(-4:mn6), wa(tid)%wa6tstep129(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep130(-4:mn6), wa(tid)%wa6tstep131(-4:mn6), wa(tid)%wa6tstep132(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep133(-4:mn6), wa(tid)%wa6tstep134(-4:mn6), wa(tid)%wa6tstep135(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep136(-4:mn6), wa(tid)%wa6tstep137(-4:mn6), wa(tid)%wa6tstep138(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep139(-4:mn6), wa(tid)%wa6tstep140(-4:mn6), wa(tid)%wa6tstep141(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep142(-4:mn6), wa(tid)%wa6tstep143(-4:mn6), wa(tid)%wa6tstep144(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep145(-4:mn6), wa(tid)%wa6tstep146(-4:mn6), wa(tid)%wa6tstep147(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep148(-4:mn6), wa(tid)%wa6tstep149(-4:mn6), wa(tid)%wa6tstep150(-4:mn6), stat=ierr)
                   allocate (wa(tid)%wa6tstep151(-4:mn6), wa(tid)%wa6tstep152(-4:mn6), stat=ierr)                           

                   allocate (wa(tid)%wa6init100(-4:mn6), wa(tid)%wa6init101(-4:mn6), stat=ierr)

                   allocate (wa(tid)%wa6rest100(-4:mn6), wa(tid)%wa6rest101(-4:mn6), stat=ierr)


                   IF (ierr .NE. 0) THEN
                      call ERROR ('GETWAMEMORYSPACE',    &
                                  'Allocation of wa6',   &
                                  __FILE__, __LINE__) 
                   ENDIF
         
                ENDDO

        
                !====================
                ! 3D - vector
                !====================

                DO tid = 1, nth

                   allocate (wa(tid)%wa3sweeps100(mnx5, mny5, mnz5), stat = ierr)

                   allocate (wa(tid)%wa3sweeps101(-4:mn6, -4:mn6, -4:mn6), stat = ierr)
                   allocate (wa(tid)%wa3sweeps102(-4:mn6, -4:mn6, -4:mn6), stat = ierr)
                   allocate (wa(tid)%wa3sweeps103(-4:mn6, -4:mn6, -4:mn6), stat = ierr)
                   allocate (wa(tid)%wa3sweeps104(-4:mn6, -4:mn6, -4:mn6), stat = ierr)
                   allocate (wa(tid)%wa3sweeps105(-4:mn6, -4:mn6, -4:mn6), stat = ierr)
                   allocate (wa(tid)%wa3sweeps106(-4:mn6, -4:mn6, -4:mn6), stat = ierr)
        
                   allocate (wa(tid)%wa3sweeps107(-4:mn6, -4:mn6, -4:mn6), stat = ierr)
                   allocate (wa(tid)%wa3sweeps108(-4:mn6, -4:mn6, -4:mn6), stat = ierr)
                   allocate (wa(tid)%wa3sweeps109(-4:mn6, -4:mn6, -4:mn6), stat = ierr)

                   allocate(wa(tid)%wa3tstep100(mnx5,mny5,mnz5), stat=ierr)
        

                   IF (ierr .NE. 0) THEN
                      call ERROR ('GETWAMEMORYSPACE', 'Allocation of wa3', &
                                  __FILE__, __LINE__)
                   ENDIF    
        

		ENDDO	


                allocate (wass100(0:mpi_dims1*mpi_dims2*mpi_dims3-1),         &
                          wass101(0:mpi_dims1*mpi_dims2*mpi_dims3-1), stat= ierr)


                IF (ierr .NE. 0) THEN
                    call ERROR ('GETWAMEMORYSPACE', 'Allocation of wass', &
                                 __FILE__, __LINE__)
                ENDIF
 
                !--------------------------------------------------------------------------- END
        
! ===================================================================================================
        END SUBROUTINE GETWAMEMORYSPACE
! ===================================================================================================


        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------


! ===================================================================================================
        SUBROUTINE FREEWAMEMORYSPACE()
! ===================================================================================================

                implicit none

                !------------------------------------------ Input Variables

                ! ...
                
                !------------------------------------------ Output Variables

                ! ...
                
                !------------------------------------------ Local Variables

                integer(is):: tid, nth, ierr

# ifdef OPENMP
                integer(is):: omp_get_max_threads
# endif

                !------------------------------------------------------------------------- BEGIN
                    
# ifdef OPENMP
                nth = omp_get_max_threads()
# else
                nth = 1
# endif

                    
                !====================
                ! 3D - vector
                !====================    

   		DO tid=1, nth
		                    
                   deallocate (wa(tid)%wa3sweeps100, wa(tid)%wa3sweeps101, wa(tid)%wa3sweeps102, &
                               wa(tid)%wa3sweeps103, wa(tid)%wa3sweeps104, wa(tid)%wa3sweeps105, &
                               wa(tid)%wa3sweeps106, wa(tid)%wa3sweeps107, wa(tid)%wa3sweeps108, &
                               wa(tid)%wa3sweeps109)

	 	   deallocate (wa(tid)%wa3tstep100)    

		ENDDO
                            
                DO tid= 1, nth


                   !====================
                   ! 1D - vector
                   !====================

                   deallocate (wa(tid)%wa1init100, wa(tid)%wa1init101, wa(tid)%wa1init102)
                   deallocate (wa(tid)%wa1init103, wa(tid)%wa1init104, wa(tid)%wa1init105)
                   deallocate (wa(tid)%wa1init106, wa(tid)%wa1init107, wa(tid)%wa1init108)
		   deallocate (wa(tid)%wa1init109)

                   deallocate (wa(tid)%wa1rest100, wa(tid)%wa1rest101, wa(tid)%wa1rest102)
                   deallocate (wa(tid)%wa1rest103, wa(tid)%wa1rest104, wa(tid)%wa1rest105)
                   deallocate (wa(tid)%wa1rest106, wa(tid)%wa1rest107, wa(tid)%wa1rest108)
                   deallocate (wa(tid)%wa1rest109)


                   deallocate (wa(tid)%wa1get100, wa(tid)%wa1get101, wa(tid)%wa1get102)
                   deallocate (wa(tid)%wa1get103, wa(tid)%wa1get104, wa(tid)%wa1get105)
                   deallocate (wa(tid)%wa1get106, wa(tid)%wa1get107, wa(tid)%wa1get108)
                   deallocate (wa(tid)%wa1get109) 

                   !====================
                   ! -4 : mn5 vector
                   !====================

                   deallocate (wa(tid)%wa5plm100)
                
                   deallocate (wa(tid)%wa5eos100, wa(tid)%wa5eos101, wa(tid)%wa5eos102)
                   deallocate (wa(tid)%wa5eos103, wa(tid)%wa5eos104, wa(tid)%wa5eos105)
                   deallocate (wa(tid)%wa5eos106, wa(tid)%wa5eos107, wa(tid)%wa5eos108)
                   deallocate (wa(tid)%wa5eos109, wa(tid)%wa5eos110)

                   deallocate (wa(tid)%wa5bis100, wa(tid)%wa5bis101, wa(tid)%wa5bis102)
                   deallocate (wa(tid)%wa5bis103, wa(tid)%wa5bis104, wa(tid)%wa5bis105)
                   deallocate (wa(tid)%wa5bis106, wa(tid)%wa5bis107, wa(tid)%wa5bis108)
                   deallocate (wa(tid)%wa5bis109)

                   deallocate (wa(tid)%wa5get100, wa(tid)%wa5get101, wa(tid)%wa5get102)
                   deallocate (wa(tid)%wa5get103, wa(tid)%wa5get104, wa(tid)%wa5get105)
                   deallocate (wa(tid)%wa5get106, wa(tid)%wa5get107, wa(tid)%wa5get108)
                   deallocate (wa(tid)%wa5get109, wa(tid)%wa5get110, wa(tid)%wa5get111)
                   deallocate (wa(tid)%wa5get112, wa(tid)%wa5get113, wa(tid)%wa5get114)
                   deallocate (wa(tid)%wa5get115, wa(tid)%wa5get116, wa(tid)%wa5get117)
                   deallocate (wa(tid)%wa5get118, wa(tid)%wa5get119, wa(tid)%wa5get120)
                   deallocate (wa(tid)%wa5get121, wa(tid)%wa5get122, wa(tid)%wa5get123)

                   deallocate (wa(tid)%wa5hlle100, wa(tid)%wa5hlle101, wa(tid)%wa5hlle102)
                   deallocate (wa(tid)%wa5hlle103, wa(tid)%wa5hlle104, wa(tid)%wa5hlle105)
                   deallocate (wa(tid)%wa5hlle106, wa(tid)%wa5hlle107, wa(tid)%wa5hlle108)
                   deallocate (wa(tid)%wa5hlle109, wa(tid)%wa5hlle110, wa(tid)%wa5hlle111)
                   deallocate (wa(tid)%wa5hlle112, wa(tid)%wa5hlle113, wa(tid)%wa5hlle114)
                   deallocate (wa(tid)%wa5hlle115, wa(tid)%wa5hlle116, wa(tid)%wa5hlle117)
                   deallocate (wa(tid)%wa5hlle118, wa(tid)%wa5hlle119, wa(tid)%wa5hlle120)
                   deallocate (wa(tid)%wa5hlle121, wa(tid)%wa5hlle122, wa(tid)%wa5hlle123)
                   deallocate (wa(tid)%wa5hlle124, wa(tid)%wa5hlle125, wa(tid)%wa5hlle126)
                   deallocate (wa(tid)%wa5hlle127, wa(tid)%wa5hlle128, wa(tid)%wa5hlle129)
                   deallocate (wa(tid)%wa5hlle130, wa(tid)%wa5hlle131, wa(tid)%wa5hlle132)
                   deallocate (wa(tid)%wa5hlle133, wa(tid)%wa5hlle134, wa(tid)%wa5hlle135)
                   deallocate (wa(tid)%wa5hlle136, wa(tid)%wa5hlle137, wa(tid)%wa5hlle138)
                   deallocate (wa(tid)%wa5hlle139, wa(tid)%wa5hlle140, wa(tid)%wa5hlle141)
                   deallocate (wa(tid)%wa5hlle142)


                   deallocate (wa(tid)%wa5hlld100, wa(tid)%wa5hlld101, wa(tid)%wa5hlld102)
                   deallocate (wa(tid)%wa5hlld103, wa(tid)%wa5hlld104, wa(tid)%wa5hlld105)
                   deallocate (wa(tid)%wa5hlld106, wa(tid)%wa5hlld107, wa(tid)%wa5hlld108)
                   deallocate (wa(tid)%wa5hlld109, wa(tid)%wa5hlld110, wa(tid)%wa5hlld111)
                   deallocate (wa(tid)%wa5hlld112, wa(tid)%wa5hlld113, wa(tid)%wa5hlld114)
                   deallocate (wa(tid)%wa5hlld115, wa(tid)%wa5hlld116, wa(tid)%wa5hlld117)
                   deallocate (wa(tid)%wa5hlld118, wa(tid)%wa5hlld119, wa(tid)%wa5hlld120)
                   deallocate (wa(tid)%wa5hlld121, wa(tid)%wa5hlld122, wa(tid)%wa5hlld123)
                   deallocate (wa(tid)%wa5hlld124, wa(tid)%wa5hlld125, wa(tid)%wa5hlld126)
                   deallocate (wa(tid)%wa5hlld127, wa(tid)%wa5hlld128, wa(tid)%wa5hlld129)
                   deallocate (wa(tid)%wa5hlld130, wa(tid)%wa5hlld131, wa(tid)%wa5hlld132)
	           deallocate (wa(tid)%wa5hlld133)		


                   deallocate (wa(tid)%wa5mhydrow100, wa(tid)%wa5mhydrow101, wa(tid)%wa5mhydrow102)
                   deallocate (wa(tid)%wa5mhydrow103, wa(tid)%wa5mhydrow104, wa(tid)%wa5mhydrow105)
                   deallocate (wa(tid)%wa5mhydrow106, wa(tid)%wa5mhydrow107, wa(tid)%wa5mhydrow108)
                   deallocate (wa(tid)%wa5mhydrow109, wa(tid)%wa5mhydrow110, wa(tid)%wa5mhydrow111)
                   deallocate (wa(tid)%wa5mhydrow112, wa(tid)%wa5mhydrow113, wa(tid)%wa5mhydrow114)
                   deallocate (wa(tid)%wa5mhydrow115, wa(tid)%wa5mhydrow116, wa(tid)%wa5mhydrow117)
                   deallocate (wa(tid)%wa5mhydrow118, wa(tid)%wa5mhydrow119, wa(tid)%wa5mhydrow120)
                   deallocate (wa(tid)%wa5mhydrow121, wa(tid)%wa5mhydrow122, wa(tid)%wa5mhydrow123)
                   deallocate (wa(tid)%wa5mhydrow124, wa(tid)%wa5mhydrow125)
                   

                   deallocate (wa(tid)%wa5sweeps100, wa(tid)%wa5sweeps101, wa(tid)%wa5sweeps102)
                   deallocate (wa(tid)%wa5sweeps103, wa(tid)%wa5sweeps104, wa(tid)%wa5sweeps105)
                   deallocate (wa(tid)%wa5sweeps106, wa(tid)%wa5sweeps107, wa(tid)%wa5sweeps108)
                   deallocate (wa(tid)%wa5sweeps109, wa(tid)%wa5sweeps110, wa(tid)%wa5sweeps111)
                   deallocate (wa(tid)%wa5sweeps112, wa(tid)%wa5sweeps113, wa(tid)%wa5sweeps114)
                   deallocate (wa(tid)%wa5sweeps115, wa(tid)%wa5sweeps116, wa(tid)%wa5sweeps117)
                   deallocate (wa(tid)%wa5sweeps118, wa(tid)%wa5sweeps119, wa(tid)%wa5sweeps120)
                   

                   deallocate (wa(tid)%wa5tstep100, wa(tid)%wa5tstep101, wa(tid)%wa5tstep102)
                   deallocate (wa(tid)%wa5tstep103, wa(tid)%wa5tstep104, wa(tid)%wa5tstep105)
                   deallocate (wa(tid)%wa5tstep106, wa(tid)%wa5tstep107, wa(tid)%wa5tstep108)
                   deallocate (wa(tid)%wa5tstep109, wa(tid)%wa5tstep110, wa(tid)%wa5tstep111)
                   deallocate (wa(tid)%wa5tstep112, wa(tid)%wa5tstep113, wa(tid)%wa5tstep114)
                   deallocate (wa(tid)%wa5tstep115, wa(tid)%wa5tstep116, wa(tid)%wa5tstep117)
                   deallocate (wa(tid)%wa5tstep118, wa(tid)%wa5tstep119, wa(tid)%wa5tstep120)
                   deallocate (wa(tid)%wa5tstep121, wa(tid)%wa5tstep122, wa(tid)%wa5tstep123)
                   deallocate (wa(tid)%wa5tstep124, wa(tid)%wa5tstep125, wa(tid)%wa5tstep126)
                   deallocate (wa(tid)%wa5tstep127, wa(tid)%wa5tstep128, wa(tid)%wa5tstep129)
                   deallocate (wa(tid)%wa5tstep130, wa(tid)%wa5tstep131, wa(tid)%wa5tstep132)
                   deallocate (wa(tid)%wa5tstep133, wa(tid)%wa5tstep134, wa(tid)%wa5tstep135)
                   deallocate (wa(tid)%wa5tstep136, wa(tid)%wa5tstep137, wa(tid)%wa5tstep138)
                   deallocate (wa(tid)%wa5tstep139, wa(tid)%wa5tstep140, wa(tid)%wa5tstep141)
           

                   !====================
                   ! -4 : mn6 vector
                   !====================


                   deallocate (wa(tid)%wa6mhydrow100, wa(tid)%wa6mhydrow101, wa(tid)%wa6mhydrow102)
                   deallocate (wa(tid)%wa6mhydrow103, wa(tid)%wa6mhydrow104, wa(tid)%wa6mhydrow105)
                   deallocate (wa(tid)%wa6mhydrow106, wa(tid)%wa6mhydrow107, wa(tid)%wa6mhydrow108)
                   deallocate (wa(tid)%wa6mhydrow109, wa(tid)%wa6mhydrow110, wa(tid)%wa6mhydrow111)
                   deallocate (wa(tid)%wa6mhydrow112, wa(tid)%wa6mhydrow113, wa(tid)%wa6mhydrow114)
                   deallocate (wa(tid)%wa6mhydrow115, wa(tid)%wa6mhydrow116, wa(tid)%wa6mhydrow117)
                   deallocate (wa(tid)%wa6mhydrow118, wa(tid)%wa6mhydrow119, wa(tid)%wa6mhydrow120)
                   deallocate (wa(tid)%wa6mhydrow121, wa(tid)%wa6mhydrow122, wa(tid)%wa6mhydrow123)
                   deallocate (wa(tid)%wa6mhydrow124, wa(tid)%wa6mhydrow125, wa(tid)%wa6mhydrow126)
                   deallocate (wa(tid)%wa6mhydrow127, wa(tid)%wa6mhydrow128, wa(tid)%wa6mhydrow129)
                   deallocate (wa(tid)%wa6mhydrow130, wa(tid)%wa6mhydrow131, wa(tid)%wa6mhydrow132)
                   deallocate (wa(tid)%wa6mhydrow133, wa(tid)%wa6mhydrow134, wa(tid)%wa6mhydrow135)
                   deallocate (wa(tid)%wa6mhydrow136, wa(tid)%wa6mhydrow137, wa(tid)%wa6mhydrow138)
                   deallocate (wa(tid)%wa6mhydrow139, wa(tid)%wa6mhydrow140, wa(tid)%wa6mhydrow141)
                   deallocate (wa(tid)%wa6mhydrow142, wa(tid)%wa6mhydrow143, wa(tid)%wa6mhydrow144)
                   deallocate (wa(tid)%wa6mhydrow145, wa(tid)%wa6mhydrow146, wa(tid)%wa6mhydrow147)
                   deallocate (wa(tid)%wa6mhydrow148, wa(tid)%wa6mhydrow149, wa(tid)%wa6mhydrow150)
                   deallocate (wa(tid)%wa6mhydrow151, wa(tid)%wa6mhydrow152, wa(tid)%wa6mhydrow153)
                   deallocate (wa(tid)%wa6mhydrow154, wa(tid)%wa6mhydrow155, wa(tid)%wa6mhydrow156)
                   deallocate (wa(tid)%wa6mhydrow157, wa(tid)%wa6mhydrow158)

                   deallocate (wa(tid)%wa6sweeps100, wa(tid)%wa6sweeps101, wa(tid)%wa6sweeps102)
                   deallocate (wa(tid)%wa6sweeps103, wa(tid)%wa6sweeps104, wa(tid)%wa6sweeps105)
                   deallocate (wa(tid)%wa6sweeps106, wa(tid)%wa6sweeps107, wa(tid)%wa6sweeps108)

                   deallocate (wa(tid)%wa6MP100, wa(tid)%wa6MP101, wa(tid)%wa6MP102)
                   deallocate (wa(tid)%wa6MP103, wa(tid)%wa6MP104)

                   deallocate (wa(tid)%wa6tstep101, wa(tid)%wa6tstep102)
                   deallocate (wa(tid)%wa6tstep103, wa(tid)%wa6tstep104, wa(tid)%wa6tstep105)
                   deallocate (wa(tid)%wa6tstep106, wa(tid)%wa6tstep107, wa(tid)%wa6tstep108)
                   deallocate (wa(tid)%wa6tstep109, wa(tid)%wa6tstep110, wa(tid)%wa6tstep111)
                   deallocate (wa(tid)%wa6tstep112, wa(tid)%wa6tstep113, wa(tid)%wa6tstep114)
                   deallocate (wa(tid)%wa6tstep115, wa(tid)%wa6tstep116, wa(tid)%wa6tstep117)
                   deallocate (wa(tid)%wa6tstep118, wa(tid)%wa6tstep119, wa(tid)%wa6tstep120)
                   deallocate (wa(tid)%wa6tstep121, wa(tid)%wa6tstep122, wa(tid)%wa6tstep123)
                   deallocate (wa(tid)%wa6tstep124, wa(tid)%wa6tstep125, wa(tid)%wa6tstep126)
                   deallocate (wa(tid)%wa6tstep127, wa(tid)%wa6tstep128, wa(tid)%wa6tstep129)
                   deallocate (wa(tid)%wa6tstep130, wa(tid)%wa6tstep131, wa(tid)%wa6tstep132)
                   deallocate (wa(tid)%wa6tstep133, wa(tid)%wa6tstep134, wa(tid)%wa6tstep135)
                   deallocate (wa(tid)%wa6tstep136, wa(tid)%wa6tstep137, wa(tid)%wa6tstep138)
                   deallocate (wa(tid)%wa6tstep139, wa(tid)%wa6tstep140, wa(tid)%wa6tstep141)
                   deallocate (wa(tid)%wa6tstep142, wa(tid)%wa6tstep143, wa(tid)%wa6tstep144)
                   deallocate (wa(tid)%wa6tstep145, wa(tid)%wa6tstep146, wa(tid)%wa6tstep147)
                   deallocate (wa(tid)%wa6tstep148, wa(tid)%wa6tstep149, wa(tid)%wa6tstep150)
                   deallocate (wa(tid)%wa6tstep151, wa(tid)%wa6tstep152)                           

                   deallocate (wa(tid)%wa6init100, wa(tid)%wa6init101)

                   deallocate (wa(tid)%wa6rest100, wa(tid)%wa6rest101)


                   IF (ierr .NE. 0) THEN
                      call ERROR ('GETWAMEMORYSPACE',    &
                                  'Allocation of wa6',   &
                                  __FILE__, __LINE__) 
                   ENDIF
         
     
               ENDDO

               !--------------------------------------------------------------------------- END

! ===================================================================================================
        END SUBROUTINE FREEWAMEMORYSPACE
! ===================================================================================================

END MODULE workarrays



