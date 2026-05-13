
!*********************************
!     3D-RMHD CODE: LÓSTREGO
!*********************************
      
! NAME:
!       EOS (Equation of State)

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
!       Calculate the Equation of State (EoS) of an ideal gas.

! ACADEMIC REFERENCE:
!       Relativistic Hydrodynamics. L.Rezzolla & O. Zanotti (Chap. 2.4.6)



! ===============================================================================      
SUBROUTINE EOS (input, m, n, g, eps, rho, p, dpdrho, dpdeps, h, cs)
! ===============================================================================
      
        ! The Equation of State (EoS) for an ideal gas can be written as:
      
        ! p   :: pressure
        ! rho :: density
        ! eps :: internal energy
        ! g   :: adiabatic index

        ! h :: specific enthalpy

        ! cs :: speed of sound

      
        ! input :: (1) pressure (p), density (rho) --> 
        !              enthalpy (h), internal energy (eps), speed of sound (cs).
        !          (2) internal energy (eps), density (rho) -->
        !              enthalpy (h), pressure (p), speed of sound (cs)

      
        implicit none
#include "types.h"
      
        !------------------------------------------------------- Input Variables

        integer(is):: input, m, n

        real(rs)   :: g

        ! size -4:mn6
        real(rs)   :: p(-4:*), rho(-4:*)

        !------------------------------------------------------- Output Variables

        !size -4:mn6
        real(rs)   :: eps(-4:*), dpdrho(-4:*), dpdeps(-4:*), h(-4:*), cs(-4:*)

        !------------------------------------------------------- Local Variables
      
        integer(is):: i

        real(rs)   :: gam1

        !------------------------------------------------------------------------- BEGIN      

        gam1 = g - 1.D0
	
        ! The variables of state are calculated at indices m to n.
      
        DO i= m, n
        
           IF (input .EQ. 1) THEN
              eps(i) = p(i)/rho(i)/gam1
           ELSE IF (input .EQ. 2) THEN
              p(i) = gam1 * rho(i)*eps(i)
           ENDIF

           dpdrho(i) = gam1*eps(i)
           dpdeps(i) = gam1*rho(i)
           h(i)      = 1.D0 + eps(i) + p(i)/rho(i)
           cs(i)     = DSQRT((dpdrho(i) + p(i)*dpdeps(i)/rho(i)/rho(i))/h(i))
         
        ENDDO

	RETURN          

        !-------------------------------------------------------------------------END

! ===============================================================================        
END SUBROUTINE EOS
! ===============================================================================

      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      
! NAME:
!       EOS_SYNGE (Synge Equation of State)

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
!       Calculate the Equation of State of a multicomponent Synge fluid (non-degenarate 
!       relativistic fluid).
     
! ACADEMIC REFERENCE:
      
!       "An upwind numerical scheme for relativistic hydrodynamics with a 
!        general equation of state." S.A.E.G.Falle and S.S.Komissarov. MNRAS. 
!        278, 586-602. 1996.

!       "A multidimensional relativistic hydrodynamic code with a general 
!        equation of state." E. Choi and P.J. Wiita. AJSS, 191:113-123. 2010

!       "A numerical simulation of the evolution and fate of a Fanaroff-Riley type 
!        I jet. The case of 3C 31." M. Perucho and J.M. Marti. MNRAS. 382, 526-542.2007

!        See also: "Relativistic Hydrodynamics." L.Rezzolla & O. Zanotti (Chap. 2.4.2)
!        See also: "Linear and non-linear Kelvin-Helmhotz instabilities in relativistic 
!        flows." M.Perucho. 2005

      

! ===============================================================================      
SUBROUTINE EOS_SYNGE (input, tid, m, n, me, mp, rho, xe, tem, eps,   & 
                      p, gam, h, cs, dpdrho, dpdeps)
! ===============================================================================

                      
        ! For a multi-component relativistic gas (with i species), we adopt the
	! following nomenclature:

        ! z_i :: z_i = m_i/Tem
        ! h   :: specific enthalpy
        ! m_i :: mass of the i-component
        ! rho :: density
        ! tem :: temperature
        ! eps :: specific internal energy

        ! p   :: total pressure (sum of the partial i-pressures).
        ! n_i :: number density of the i-component.

        ! meff :: effective mass
        ! Xe   :: lepton fraction mass :: Xe=me/(me+mp)
        ! Xp   :: proton fraction mass :: Xp=mp/(me+mp)=Xe-1

        ! gam :: adiabatic factor

        ! cs :: sound speed
    
        ! input :: (1) density (rho), mass fraction (xe), temperature (tem) --> 
        !              int. energy (eps), pressure (p)...
        !          (2) density (rho), mass fraction (xe), int. energy (eps) --> 
        !              temperature (tem), pressure (p)...
        !          (3) mass fraction (xe), temperature (tem)                --> 
        !              int. energy (eps), pressure (p)...
	!          (4) density (rho), mass fraction (xe), pressure (p)      -->
	!              int. energy (eps), temperature (tem)...


        
        ! units :: c=1, K_b=1              
                      

        USE workarrays 
        implicit none

# include "types.h"

        !------------------------------------------------------- Input Variables 
      
        integer(is):: input, m, n, tid

        real(rs)   ::  me, mp

        ! size -4:mn6
        real(rs)   :: tem(-4:*), rho(-4:*), xe(-4:*)


        !------------------------------------------------------ Output Variables 
      
        ! size -4:mn6
        real(rs)   :: dpdrho(-4:*), dpdeps(-4:*), h(-4:*), cs(-4:*), gam(-4:*), &
                      p(-4:*), eps(-4:*)


        !------------------------------------------------------- Local Variables

        integer(is):: i, it, itmax

        real(rs)   ::  errlim , errmax

        real(rs), pointer   ::  ze(:), zp(:), meff1(:),              &
                                ge(:), dge(:), gp(:), dgp(:),        &
                                f(:), df(:), oldtem(:), relerr(:)             
                             

        ! Limits for Newton-Raphson iterations.                         
        errlim = 1.D-7
        itmax  = 100	

        ze     => wa(tid)%wa5eos100
        zp     => wa(tid)%wa5eos101
        meff1  => wa(tid)%wa5eos102             
        ge     => wa(tid)%wa5eos103
        dge    => wa(tid)%wa5eos104
        gp     => wa(tid)%wa5eos105 
        dgp    => wa(tid)%wa5eos106       
        f      => wa(tid)%wa5eos107
        df     => wa(tid)%wa5eos108
        oldtem => wa(tid)%wa5eos109 
        relerr => wa(tid)%wa5eos110

       
        !------------------------------------------------------------------------- BEGIN

        ! (1) Density (rho), mass fraction (xe) and temperature (tem)

      
        IF (input .EQ. 1)  THEN
           
           DO i = m, n     

		  
              ze(i) = me / tem(i)
              zp(i) = mp / tem(i)


              ! Inverse effective mass term.
              meff1(i) = (1.D0-xe(i))/mp + xe(i)/me

              dpdrho(i) = tem(i)*meff1(i)
              p(i) = rho(i)*dpdrho(i)
              ! Lepton specific enthalpy.
              ge(i)=5.D0/2.D0/ze(i)+DSQRT(9.D0/4.D0/ze(i)**2+1.D0)

              ! Proton specific enthalpy.
              gp(i)=5.D0/2.D0/zp(i)+DSQRT(9.D0/4.D0/zp(i)**2+1.D0)

              ! Derivatives of lepton/proton specific enthalpy (to calculate gamma factor).
              dge(i)=-5.D0/2.D0/ze(i)**2-9.D0/4.D0/ze(i)**3/ &
                      DSQRT(9.D0/4.D0/ze(i)**2+1.D0)
              dgp(i)=-5.D0/2.D0/zp(i)**2-9.D0/4.D0/zp(i)**3/ &
                      DSQRT(9.D0/4.D0/zp(i)**2+1.D0)     


           ENDDO
      
           DO i = m, n 

              h(i)= (1.D0-xe(i))*gp(i) + xe(i)*ge(i)
              eps(i) = h(i) - dpdrho(i) - 1.D0
              gam(i) = 1.D0 - 1.D0/(1.D0+((1.D0-xe(i))*dgp(i)*zp(i)**2/mp + & 
                       xe(i)*dge(i)*ze(i)**2/me)/meff1(i))
              dpdeps(i) = rho(i)*(gam(i)-1.D0)
              cs(i) = DSQRT(gam(i)*p(i) / (h(i)*rho(i)))

           ENDDO

	   RETURN

        ! ----------------

        ! (2) Density (rho), mass fraction (xe) and specific internal energy (eps)

        ELSE IF (input .EQ. 2)  THEN

           DO i = m, n
              meff1(i) = (1.D0-xe(i))/mp + xe(i)/me
           ENDDO

           ! The temperature (tem) is used as the initial seed for Newton Raphson.
           DO it = 1, itmax      

              DO i = m, n
                 
                 ze(i) = me / tem(i)
                 zp(i) = mp / tem(i)

                 dpdrho(i) = tem(i)*meff1(i)

                 ! Lepton specific enthalpy.
                 ge(i)=5.D0/2.D0/ze(i)+DSQRT(9.D0/4.D0/ze(i)**2+1.D0)

                 ! Proton specific enthalpy.
                 gp(i)=5.D0/2.D0/zp(i)+DSQRT(9.D0/4.D0/zp(i)**2+1.D0)

                 dge(i)=-5.D0/2.D0/ze(i)**2-9.D0/4.D0/ze(i)**3/  &
                         DSQRT(9.D0/4.D0/ze(i)**2+1.D0)
                 dgp(i)=-5.D0/2.D0/zp(i)**2-9.D0/4.D0/zp(i)**3/  &
                         DSQRT(9.D0/4.D0/zp(i)**2+1.D0)
                 
              ENDDO
  
              ! =============================
              !NEWTON-RAPHSON ITERATIONS
              ! =============================


              DO i = m, n
                 
                 ! f(i) = h(tem)-p(tem)/rho-1-eps
                 f(i) = (1.D0-xe(i))*gp(i) + xe(i)*ge(i) - dpdrho(i) - 1.D0 - eps(i)

                 ! df(i) = dh(tem)-dp(tem)/rho (the internal energy does NOT depend on tem)
                 df(i) = -(1.D0-xe(i))*dgp(i)*zp(i)**2/mp -   &
                           xe(i)*dge(i)*ze(i)**2/me - meff1(i)

		 !df(i) = (1.D0-xe(i))*dgp(i)+xe(i)*dge(i)-meff1(i)

                 ! (i+1)-ROOT selection. 
		 ! -------------------------------------------

                 oldtem(i) = tem(i)
                 tem(i) = oldtem(i) - f(i)/df(i)
                 relerr(i) = ABS((tem(i)-oldtem(i))/oldtem(i))
                 
              ENDDO


              ! Readjust the temperature when negative values are found.
              DO i = m, n
                 IF (tem(i) .LE. 0.D0) THEN
                    
                    tem(i) = oldtem(i)/2.D0
                    relerr(i) = ABS((tem(i)-oldtem(i))/oldtem(i))
                    
                 ENDIF
              ENDDO 

              errmax = 0.D0  
              DO i = m, n
                 IF (relerr(i) .GT. errmax) THEN
                    errmax = relerr(i)
                 ENDIF
              ENDDO

              ! Stop Newton-Raphson if the error limit is achieved before the
              ! maximum number of iterations.
              IF (errmax .LE. errlim) exit

           ENDDO

	   ! =============================
           ! END NEWTON-RAPHSON ITERATIONS
           ! =============================


           DO i = m, n
              
              p(i) = rho(i)*dpdrho(i)
              h(i) = 1.D0 + eps(i) + dpdrho(i)
              gam(i) = 1.D0 - 1.D0/(1.D0+((1.D0-xe(i))*dgp(i)*zp(i)**2/mp +    &
                       xe(i)*dge(i)*ze(i)**2/me)/meff1(i))
              dpdeps(i) = rho(i)*(gam(i)-1.D0)
              cs(i) = DSQRT(gam(i)*p(i)/(h(i)*rho(i)))
              
           ENDDO

	   RETURN

        ! ----------------

        ! (3) Mass fraction (xe) and Temperature (tem)

        ELSE IF (input .EQ. 3) THEN

           DO i=m, n

              ze(i) = me/tem(i)
              zp(i) = mp/tem(i)
              meff1(i) = (1.D0-xe(i))/mp + xe(i)/me

              ge(i)=5.D0/2.D0/ze(i)+DSQRT(9.D0/4.D0/ze(i)**2+1.D0)
              gp(i)=5.D0/2.D0/zp(i)+DSQRT(9.D0/4.D0/zp(i)**2+1.D0)

              dge(i)=-5.D0/2.D0/ze(i)**2-9.D0/4.D0/ze(i)**3/   &
                      DSQRT(9.D0/4.D0/ze(i)**2+1.D0)
              dgp(i)=-5.D0/2.D0/zp(i)**2-9.D0/4.D0/zp(i)**3/   &
                      DSQRT(9.D0/4.D0/zp(i)**2+1.D0)

              gam(i) = 1.D0 - 1.D0/(1.D0+( (1.D0-xe(i))*dgp(i)*zp(i)**2/mp + & 
                       xe(i)*dge(i)*ze(i)**2/me)/meff1(i))

           ENDDO

	   RETURN


	ENDIF


        !------------------------------------------------------------------------- END

! ===============================================================================        
END SUBROUTINE EOS_SYNGE
! ===============================================================================

