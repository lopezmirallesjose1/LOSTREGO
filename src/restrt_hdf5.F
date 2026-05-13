! NAME:
!       RESTRT

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

! DESCRIPTION:
!       Reads and writes files for restarting jobs, using HDF5 directives.

! ACADEMIC REFERENCE:


! LIST OF MODIFICATIONS:
!       17 April 2019. Subroutine's core.


# ifdef HDF
      
      subroutine restrt( basenm, ird, mn6, mn5, mn1, mnx1, mnx5, mny1,&
      mny5, mnz1, mnz5, freqpl, dtini,                    & 
      H_openr, H_opend, rst_file, rst_rec, rst_filed,    &
      datfil, rstfil, outfil, tdlfil, suffixd, suffix,   &
      nrst, nout1,                                       & 
      nx, nxo, ny0, ny, nyold, nyh, nz, nzo, nyt,        &
      nxtot, nytot, nztot,                               &
      nsdim, igodu, igeomx, igeomy, igeomz, nstep,       &
      bndmxz, bndmny, bndmxy, bndmnz, bndmxx, bndmnx,    &
      small, epsiln, omg1, omg2, boundPos,               &
      dt, dk1, dk2, psh1, psh2, temb, machb, csb2,       &
      nn, mm1, mm2, tem0, tem1,                          &
      oang, ipr, ms, nmod, A0, fra, frs, pha,            &
      dely, delz, tout1, trst, gamma_ad,                 &
      gridlx, gridly, gridlz, pb, rhob, rhoeb, velb,     &
      timej, grav, rhoa, pa, dpady, uwi, vwi, vk,        &
      radius, cdenst, shift, fg, rm, brt, rorb,          &
      t0, t1, vptbd, rptbd, pptbd, alpn, betn, gamn,     &
      rbn, m_0, R_b, c, pc, rho_0a, q0, me, mp,          &
      xznl, xzn, xznr, yznl, yzn, yznr, zznl, zzn, zznr, &
      delx, iextrax, iextraz, gridlx2, gridlz2,          &
      nx2, nz2, ixsym, izsym, ncl, ic, jc, kc, mclouds,  &
      physic, xCPU, yCPU, zCPU )

                   
           
      USE constantes
      USE parallel
      USE tipos
      USE workarrays
      
#if defined(PARALELO) && defined(MN) 
      USE mpi
#endif
      IMPLICIT NONE
#if defined(PARALELO) &&        !defined(MN)
      include "mpif.h"
#endif
#      include "types.h"

!     ny0 is the size of the homogeneous grid
!     ny is the total size of the local grid
!     nyh is the global size of the homogeneous grid
!     nyt is the global size of the grid
!-------------------------------------------------------Input Variables
      
      
      character(4), intent(in) :: basenm

      integer(is):: ird, mnx1, mnx5, mny1, mny5, mnz1, mnz5, mn6, mn5, mn1, freqpl,   &
      iextrax, iextraz, nx2, nz2, ixsym, izsym, ny0, igodu, &
      H_openr, H_opend, r_filer, r_filed, nxo, nzo, nxtot, nztot, nytot


      integer(is):: ic(2000), jc(2000), kc(2000)  

      real(rs) :: dtini, delx, small, gridlx2, gridlz2,        &
      epsiln, omg1, omg2, alpn, betn, gamn,        &
      rbn, m_0, R_b, c, pc, shift, rho_0a, q0,     &
      radius, cdenst, fg, rm, me, mp, me_cgs, kb, mp_cgs 

      integer :: xCPU, yCPU, zCPU

      type(TboundPosY) :: boundPos(-4:5)

      character*(*), intent(in) :: datfil
      character*(*), intent(in) :: rstfil, outfil, tdlfil
      character*(*), intent(in) :: suffixd
      character*(*), intent(in) :: suffix

!------------------------------------------------------Output Variables

      integer(is) :: nrst, nout1, nx, ny, nyold, nyh, nz, nyt,        &
      nsdim, igeomx, igeomy, igeomz, nstep,            &
      bndmxz, bndmny, bndmxy, bndmnz, bndmxx, bndmnx,  &
      ipr, ms, nmod, mclouds, brt, ncl

      real(rs) :: dt, dk1, dk2, psh1, psh2, temb, machb, csb2, nn,    &
      mm1, mm2, dely, delz, tout1, trst,                  &
      gamma_ad, gridlx, gridly, gridlz, pb, rhob, rhoeb,  & 
      velb, timej, tem0, tem1, &
      uwi, vwi, vk, rorb, oang, t0, t1, vptbd, rptbd, pptbd,     &
      A0, fra(25), frs(25), pha(25)

      real(rs) :: grav(-4:mny5), rhoa(-4:mny5),pa(-4:mny5), dpady(-4:mny5)

      real(rs) :: xznl(mnx1), xzn(mnx1), xznr(mnx1), &
      yznl(mny1), yzn(mny1), yznr(mny1), &
      zznl(mnz1), zzn(mnz1), zznr(mnz1)

      type(Tphysic)::  physic(-4:mnx5,-4:mny5,-4:mnz5)

!-------------------------------------------------------Local Variables
      character(256) :: fich

      character(4), external :: int2ch4
      integer(is) :: rst_file, rst_rec, rst_filed, curR_R

      integer(is) :: tid, i, j, k, l, mpi_coordx, mpi_coordy, mpi_coordz, ll, plane

      integer(is) :: jj, ierr, nx3, nz3, &
                    iic, jjc, kkc, idum, nrad, ci, cj, ck

      integer(is) :: dxoff, dyoff, dzoff

      real(rs):: scr1, scrc1, scr2, scrc2, scr3, scrc3, scr4, scrc4, & 
      gamm1, c1, c2, c3, ar, fr, dfdr, erro

      real(rs) :: tt, theta1, theta2, yz, orho, WL, xc, yc, zc, ac, bc, cc,  &
      distc, rhocl, delyh, rhot, rangl, sangl, radj, &
      delc, dist, ran1, gridlxo, z1, z2, rr, rstar
  
      real(rs) :: mom0(2), parmax0(4), parmin0(2), yextra(ny0/4), &
      paext(ny0/4), dpadyext(ny0/4), gravext(ny0/4),  &
      xznl0(mnx1-nx2), zznl0(mnx1-nx2)

      real(rs):: yzncl(2000)

      real(rs):: tem_min, tem_max, hb, gamb, epsb, etars, radi, & 
      rc, teme, he, game, pe, epse, rhoe, rhoa3 


      real(rs), pointer :: rho_sy(:),xe_s(:),tem_s(:),   & 
      eps_sy(:),p_s(:),gam_s(:),    &
      h_s(:),cs_s(:),              &
      dpdrho_s(:),dpdeps_s(:)  
      
      real(rs) :: rad(-4:mn5), rho2(-4:mn5), vel2(-4:mn5) 
      real(rs), allocatable :: scrch1(:), scrch2(:), scrch3(:), scrch4(:)
      real(rs), pointer :: rhocol0(:,:), profs0(:,:,:)

      type(TplanePhysic), pointer :: xplane0(:,:), yplane0(:,:,:), zplane0(:,:)

!--------------------------------------------------------------------------BEGIN
      
      tid = 1

      r_filed = rst_filed
      r_filer = rst_file
      plane  = 0  
      curR_R = 0
      H_openr = 0
      H_opend = 0
      
      allocate (scrch1(-4:mn6), stat= ierr)
      allocate (scrch2(-4:mn6), stat= ierr)
      allocate (scrch3(-4:mn6), stat= ierr)
      allocate (scrch4(-4:mn6), stat= ierr)


 ! -4:mn5 vectors
      rho_sy=> wa(1)%wa201
      xe_s  => wa(1)%wa202
      tem_s => wa(1)%wa203
      eps_sy=> wa(1)%wa204
      p_s   => wa(1)%wa205
      gam_s => wa(1)%wa206
      h_s   => wa(1)%wa207
      cs_s  => wa(1)%wa208
      dpdrho_s => wa(1)%wa209
      dpdeps_s => wa(1)%wa210


      IF (ird == 0) THEN
         
! Write Restart

         mpi_coordx = mpi_coords(1)
         mpi_coordy = mpi_coords(2)
         mpi_coordz = mpi_coords(3)
         if (cartrank == 0) then
            write(6,*) 'RESTRT, 1: ', cartrank, mpi_coordx, mpi_coordy, mpi_coordz 
            !call flush(6)
         endif
         
         call rrhdf5_record( ird, plane, datfil, rstfil, outfil,        &
         basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
         H_openr, H_opend, curR_R, r_filer, r_filed,                    &
         gamma_ad, xzn, yzn, zzn, timej,                                &
         bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
         gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
         nstep, igeomx, igeomy, igeomz,                                 &
         mpi_coordx, mpi_coordy, mpi_coordz,                            &
         pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &
         psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
         radius, cdenst, shift, fg, rm, brt, me, mp,                    &
         m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
         ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
         t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

         rst_file = r_filer
        
         nrst = 0
         trst = 0.0

!     call filnam( 'r', basenm, datfil, rstfil, outfil, tdlfil, suffixd, suffix )
         
      ELSE IF (ird == 2) THEN
         mpi_coordx = mpi_dims(1)
         mpi_coordy = mpi_dims(2)
         mpi_coordz = mpi_dims(3)
! Write Data

         if (basenm(1:1) == 'H') then

            if (mpi_dims(2) /= 1) then

! Write XY plane information

               mpi_coordx = mpi_dims(1)
               mpi_coordy = mpi_dims(2)
               mpi_coordz = mpi_dims(3)
               plane = 1
               call rrhdf5_record(ird, plane, datfil, rstfil, outfil,        &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                    &
               gamma_ad, xzn, yzn, zzn, timej,                                &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
               gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
               nstep, igeomx, igeomy, igeomz,                                 &
               mpi_coordx, mpi_coordy, mpi_coordz,                            &
               pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &
               psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
               radius, cdenst, shift, fg, rm, brt, me, mp,                    &
               m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
               ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
               t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

! Write YZ plane information
               plane = 2
               call rrhdf5_record(ird, plane, datfil, rstfil, outfil,        &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                    &
               gamma_ad, xzn, yzn, zzn, timej,                                &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
               gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
               nstep, igeomx, igeomy, igeomz,                                 &
               mpi_coordx, mpi_coordy, mpi_coordz,                            &
               pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &
               psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
               radius, cdenst, shift, fg, rm, brt, me, mp,                    &
               m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
               ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
               t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

! Write XZ plane information
               plane = 3
               call rrhdf5_record(ird, plane, datfil, rstfil, outfil,        &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                    &
               gamma_ad, xzn, yzn, zzn, timej,                                &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
               gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
               nstep, igeomx, igeomy, igeomz,                                 &
               mpi_coordx, mpi_coordy, mpi_coordz,                            &
               pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &
               psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
               radius, cdenst, shift, fg, rm, brt, me, mp,                    &
               m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
               ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
               t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

               rst_filed = r_filed
            else
! Sequential case
               call rrhdf5_record(ird, plane, datfil, rstfil, outfil,        &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                    &
               gamma_ad, xzn, yzn, zzn, timej,                                &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
               gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
               nstep, igeomx, igeomy, igeomz,                                 &
               mpi_coordx, mpi_coordy, mpi_coordz,                            &
               pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &
               psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
               radius, cdenst, shift, fg, rm, brt, me, mp,                    &
               m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
               ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
               t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

               rst_filed = r_filed
               

            endif
         else

            

! Write XY plane information
            plane = 1
            call rrhdf5_record( ird, plane, datfil, rstfil, outfil,        &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                    &
               gamma_ad, xzn, yzn, zzn, timej,                                &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
               gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
               nstep, igeomx, igeomy, igeomz,                                 &
               mpi_coordx, mpi_coordy, mpi_coordz,                            &
               pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &
               psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
               radius, cdenst, shift, fg, rm, brt, me, mp,                    &
               m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
               ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
               t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

! Write YZ plane information
            plane = 2
            call rrhdf5_record(ird, plane, datfil, rstfil, outfil,        &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                    &
               gamma_ad, xzn, yzn, zzn, timej,                                &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
               gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
               nstep, igeomx, igeomy, igeomz,                                 &
               mpi_coordx, mpi_coordy, mpi_coordz,                            &
               pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &
               psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
               radius, cdenst, shift, fg, rm, brt, me, mp,                    &
               m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
               ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
               t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

! Write XZ plane information
            plane = 3
            call rrhdf5_record(ird, plane, datfil, rstfil, outfil,        &
               basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    & 
               H_openr, H_opend, curR_R, r_filer, r_filed,                    &
               gamma_ad, xzn, yzn, zzn, timej,                                &
               bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
               gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
               nstep, igeomx, igeomy, igeomz,                                 &
               mpi_coordx, mpi_coordy, mpi_coordz,                            &
               pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &
               psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
               radius, cdenst, shift, fg, rm, brt, me, mp,                    &
               m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
               ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
               t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)

            rst_filed = r_filed

         endif
        
!     call filnam('d', basenm, datfil, rstfil, outfil, tdlfil, suffixd, suffix)

         nout1 = 0
         tout1 = 0.0

      ELSE IF (ird == 1) THEN
         call rst_hdf5( rst_file, rst_rec, rstfil,                           &
              basenm, mnx1, mnx5, mny1, mny5, mnz1, mnz5,                    &
              H_openr, H_opend, r_filer, r_filed,                            &
              gamma_ad, xzn, yzn, zzn, timej, dt,                            &
              bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,                &
              gridlx, gridly, gridlz, nx, ny, ny0, nyold, nz, nsdim,         &
              nxo, nzo, nstep, igeomx, igeomy, igeomz,                       &
              mpi_coordx, mpi_coordy, mpi_coordz,                            &
              pb, rhob, rhoeb, velb, temb, machb, csb2, dk1, dk2,            &      
              psh1, psh2, mm1, mm2, nn, tem0, tem1, alpn, betn, gamn, rbn,   &
              radius, cdenst, shift, fg, rm, brt, me, mp,                    &
              m_0, R_b, c, pc, rho_0a, q0, oang, ipr, ms, nmod, A0,          &
              ncl, ic, jc, kc, mclouds, uwi, vwi, vk, rorb,                  &
              t0, t1, vptbd, rptbd, pptbd, grav, rhoa, pa, physic)


              rst_file = r_filer

              dt = dtini

              dxoff  = mpi_coords(1) * nx      
              dyoff  = mpi_coords(2) * ny       
              dzoff  = mpi_coords(3) * nz     


!              print*,'before set up',mpi_coords(1),mpi_coords(2),mpi_coords(3),dxoff,dzoff,nxtot,nxo,nztot,nzo

              if (  dyoff < nyold .and. dxoff >= (nxtot-nxo)/2 .and. dxoff < (nxtot+nxo)/2  &
                  .and. dzoff >= (nztot-nzo)/2 .and. dzoff < (nztot+nzo)/2 ) then  

!              print*,'setting up 1',mpi_coords(1),mpi_coords(2),mpi_coords(3),nx,ny,nz 
!              print*,'setting up 2',mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%pres,physic(nx,ny,nz)%pres
!              print*,'setting up 3',mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%densty,physic(nx,ny,nz)%densty
!              print*,'setting up 4',mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(1,1,1)%denstye,physic(nx,ny,nz)%denstye


              do k = 1, nz
                 do j= 1, ny
                    do i= 1, nx
                       rho_sy(i) = physic(i,j,k)%densty
                       xe_s(i)   = physic(i,j,k)%denstye/physic(i,j,k)%densty
                       tem_s(i)  = physic(i,j,k)%pres/(xe_s(i)*physic(i,j,k)%densty/me +       &
                                    (1.0-xe_s(i))*physic(i,j,k)%densty/mp)
                    enddo
                    
!                    print*,'before synge',mpi_coords(1),mpi_coords(2),mpi_coords(3),physic(nx,ny,nz)%pres

                    call eos_synge( 1, tid, 1, nx, me, mp, rho_sy, xe_s, tem_s, eps_sy,   &
                                    p_s, gam_s, h_s, cs_s, dpdrho_s, dpdeps_s)


!                    print*,'after synge',mpi_coords(1),mpi_coords(2),mpi_coords(3)


                    do i = 1, nx
                        scrch1(i) =    1.0 / SQRT( 1.0                           &
                                      - physic(i,j,k)%velx * physic(i,j,k)%velx &
                                      - physic(i,j,k)%vely * physic(i,j,k)%vely &
                                      - physic(i,j,k)%velz * physic(i,j,k)%velz)

                        physic(i,j,k)%pres  = p_s(i)
                        physic(i,j,k)%opres = physic(i,j,k)%pres

                        scrch2(i)    = h_s(i)

                        physic(i,j,k)%denstw  = physic(i,j,k)%densty * scrch1(i)
                        physic(i,j,k)%denstew = physic(i,j,k)%denstye * scrch1(i)
                        physic(i,j,k)%momenx  = physic(i,j,k)%denstw * scrch2(i) &
                                                * scrch1(i) * physic(i,j,k)%velx
                        physic(i,j,k)%momeny  = physic(i,j,k)%denstw * scrch2(i) &
                                                * scrch1(i) * physic(i,j,k)%vely
                        physic(i,j,k)%momenz  = physic(i,j,k)%denstw * scrch2(i) &
                                                * scrch1(i) * physic(i,j,k)%velz
                        physic(i,j,k)%energy  = physic(i,j,k)%denstw * scrch2(i) &
                                                * scrch1(i) - physic(i,j,k)%pres   &
                                               - physic(i,j,k)%denstw
                    enddo
                 enddo
              enddo 

              endif

      ELSE 
         call ERROR( 'restrt', 'Wrong ird parameter', __FILE__, __LINE__ )
      ENDIF

      RETURN
!----------------------------------------------------------------------------END
      END subroutine restrt
#endif
