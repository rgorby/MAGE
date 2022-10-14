module glsolution
    use gltypes
    use glutils
    implicit none


    !---------------------------------------------
    ! Routines to Calculate the Gibson-Low CME Model
    ! B, j, v and density, pressure, temperature
    ! port of IDL module GIBLOW in SSW version as of July, 2022
    ! and Anna Malanushenko (HAO UCAR) modifications 
    ! comments and structure modified by AJM
    !---------------------------------------------
    contains

        !> Calculate Gibson-Low Solution Derivatives 
        !> for Density
        !>
        subroutine calcDerivs(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution

            real(rp), dimension(:,:), allocatable :: st, ct, sph, cph, mu
            real(rp), dimension(:,:), allocatable :: sTbig
            real(rp), dimension(:,:), allocatable :: cTbig
            real(rp), dimension(:,:), allocatable :: sPbig
            real(rp), dimension(:,:), allocatable :: cPbig
            real(rp), dimension(:, :), allocatable :: gfun, gfunot
            real(rp), dimension(:, :), allocatable :: dSdT, dgfundr
            real(rp), dimension(:, :), allocatable :: dSdR, dbPdt, dbPdR, d2gfundr2
            real(rp), dimension(:, :), allocatable :: d2SdR2, d2SdT2, dbTdR, dbRdT
            real(rp), dimension(:, :), allocatable :: dTdmu, dTdR, dRdmu, dPidmu, dmudlam
            real(rp), dimension(:, :), allocatable :: dthetadmu, dPidR, dRdlam, d2SdTdR

            real(rp) :: cs, ss 

            allocate (st(State%js:State%je, State%ks:State%ke))
            allocate (ct(State%js:State%je, State%ks:State%ke))
            allocate (mu(State%js:State%je, State%ks:State%ke))
            allocate (sph(State%js:State%je, State%ks:State%ke))
            allocate (cph(State%js:State%je, State%ks:State%ke))
            allocate (sTbig(State%js:State%je, State%ks:State%ke))
            allocate (cTbig(State%js:State%je, State%ks:State%ke))
            allocate (sPbig(State%js:State%je, State%ks:State%ke))
            allocate (cPbig(State%js:State%je, State%ks:State%ke))
            allocate (gfunot(State%js:State%je, State%ks:State%ke))
            allocate (gfun(State%js:State%je, State%ks:State%ke))
            allocate (dSdT(State%js:State%je, State%ks:State%ke))
            allocate (dgfundr(State%js:State%je, State%ks:State%ke))
            allocate (dSdR(State%js:State%je, State%ks:State%ke))
            allocate (dbPdt(State%js:State%je, State%ks:State%ke))
            allocate (dbPdR(State%js:State%je, State%ks:State%ke))
            allocate (d2gfundr2(State%js:State%je, State%ks:State%ke))
            allocate (d2SdT2(State%js:State%je, State%ks:State%ke))
            allocate (dbTdR(State%js:State%je, State%ks:State%ke))
            allocate (dbRdT(State%js:State%je, State%ks:State%ke))
            allocate (dTdmu(State%js:State%je, State%ks:State%ke))
            allocate (dTdR(State%js:State%je, State%ks:State%ke))
            allocate (dRdmu(State%js:State%je, State%ks:State%ke))
            allocate (dPidmu(State%js:State%je, State%ks:State%ke))
            allocate (dthetadmu(State%js:State%je, State%ks:State%ke))
            allocate (dPidR(State%js:State%je, State%ks:State%ke))
            allocate (dRdlam(State%js:State%je, State%ks:State%ke))
            allocate (dmudlam(State%js:State%je, State%ks:State%ke))
            allocate (d2SdTdR(State%js:State%je, State%ks:State%ke))
            
            ! ;  See eq A19
            ! ; 
            ! ;  first need to calculate tderiv = dT/drlam, where T = Pi + Bstrength (squared)
            ! ;  in the Rcap,Thcap system
            ! ;
            st = sin(State%thpB)
            ct = cos(State%thpB)
            sph = sin(State%phpB)
            cph = cos(State%phpB)

            sTbig = sin(State%thcap)
            cTbig = cos(State%thcap)
            sPbig = sin(State%phcap)
            cPbig = cos(State%phcap)

            call zero2tiny2d(st) 
            call zero2tiny2d(ct) 
            call zero2tiny2d(sph)
            call zero2tiny2d(cph)
            
            call zero2tiny2d(sTbig) 
            call zero2tiny2d(cTbig)
            call zero2tiny2d(sPbig)
            call zero2tiny2d(cPbig)

            mu = State%rcap*sTbig
            gfunot = sin(Model%alnot*Model%rbub)/Model%alnot/Model%rbub - cos(Model%alnot*Model%rbub)
            
            dSdT = State%stream*2*cTbig/sTbig
            dgfundr = cos(Model%alnot*State%rcap)/State%rcap - sin(Model%alnot*State%rcap)/State%rcap/State%rcap/Model%alnot + &
                    Model%alnot*sin(Model%alnot*State%rcap)
            dSdR = (4.*Model%ao*pi/Model%alnot/Model%alnot)*(dgfundr*Model%rbub**2./gfunot - 2*State%rcap)*sTbig**2.
            
            d2SdTdR = dSdR*2.*cTbig/sTbig
            d2gfundr2 = (-State%rcap*State%rcap*Model%alnot*sin(Model%alnot*State%rcap) - State%rcap*cos(Model%alnot*State%rcap) + &
                        Model%alnot*State%rcap**3.*Model%alnot*cos(Model%alnot*State%rcap) - &
                        cos(Model%alnot*State%rcap)*State%rcap + 2.*sin(Model%alnot*State%rcap)/Model%alnot)/State%rcap**3.
            d2SdR2 = (4.*Model%ao*pi/Model%alnot/Model%alnot)*(d2gfundr2*Model%rbub**2/gfunot - 2.)*sTbig**2.
            d2SdT2 = State%stream*2.*(cTbig**2. - sTbig**2.)/(sTbig**2)
            
            dRdmu = sTbig
            dthetadmu = cTbig/State%rcap
            dPidR = Model%ao*dSdR
            dPidmu = Model%ao*(dRdmu*dSdR + dthetadmu*dSdT)

            dTdR = -4*(dSdT**2)/(State%rcap**5)/(sTbig**2)
            dTdR = dTdR + 2*dSdT*d2SdTdR/(State%rcap**4)/(sTbig**2)
            dTdR = dTdR - 2*(dSdR**2)/(State%rcap**3)/(sTbig**2)
            dTdR = dTdR + 2*dSdR*d2SdR2/(State%rcap**2)/(sTbig**2)

            dTdR = dTdR - 2*Model%alnot*Model%alnot*(State%stream**2)/(State%rcap**3)/(sTbig**2)
            dTdR = dTdR + 2*Model%alnot*Model%alnot*State%stream*dSdR/(State%rcap**2)/(sTbig**2)

            dTdR = dTdR/(8.*pi)
            dTdR = dTdR + dPidR

            dTdmu = -2*((dSdT/State%rcap)**2)/(mu**3)
            dTdmu = dTdmu + 2*(dSdT/State%rcap)*(dRdmu*(-dSdT/(State%rcap**2) + &
                            d2SdTdR/State%rcap) + dthetadmu*d2SdT2/State%rcap)/(mu**2)
            dTdmu = dTdmu - 2*(dSdR**2)/(mu**3)
            dTdmu = dTdmu + 2*dSdR*(dRdmu*d2SdR2 + dthetadmu*d2SdTdR)/(mu**2)
            dTdmu = dTdmu - 2*Model%alnot*Model%alnot*(State%stream**2)/(mu**3)
            dTdmu = dTdmu + 2*Model%alnot*Model%alnot*State%stream*(dRdmu*dSdR + dthetadmu*dSdT)/(mu**2)
            dTdmu = dTdmu/8./pi
            dTdmu = dTdmu + dPidmu

            dRdlam = (State%rlam - Model%xo*st*cph)/State%rcap
            ! ;
            ! ;  remember we are in a rotated mu frame
            ! ;
            cs = cos(Model%sigma)
            ss = sin(Model%sigma)
            dmudlam = (State%rlam*(st**2*(cph**2 + sph**2*cs**2) + &
                        ct**2*ss**2 - 2*st*ct*sph*ss*cs) - &
                        Model%xo*st*cph)/mu
            State%tderivR = dTdR*dRdlam
            State%tderivmu = dTdmu*dmudlam
            State%tderiv = State%tderivR + State%tderivmu
            if (State%outside_count .gt. 0) then
                where (Solution%inside_mask(State%ri, :, :) .lt. 1.) 
                    State%tderiv = Model%outScale**2.*State%tderivout
                end where
            end if
        end subroutine calcDerivs

        !> Calculate Gibson-Low Solution Density
        !> 
        subroutine calcDensity(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution

            real(rp), dimension(:,:), allocatable :: d1, d2, d3, d4

            allocate(d1(State%js:State%je, State%ks:State%ke), d2(State%js:State%je, State%ks:State%ke))
            allocate(d3(State%js:State%je, State%ks:State%ke), d4(State%js:State%je, State%ks:State%ke))

            call calcDerivs(Model, State, Solution)

            ! ;  the force F could include the effects of the self-similar
            ! ;  acceleration (alpha ne 0) - otherwise its just gravity.
            ! ; Note that to properly include alpha ne 0, however, phiss needs
            ! ; to be defined as a differential equation
            ! ;
            State%F = GMm/(State%rsquig**2)/(Rsolar_cgs**2.) + Model%alpha*State%rsquig*Rsolar_cgs*Mp_cgs

            d1 = -((State%rlam/State%rsquig)**2.)*(1. - ((State%rlam/State%rsquig)**2.))*State%tderiv
            d2 = 2.*(State%rlam/State%rsquig)*(Model%apar/(State%rsquig**2.))*State%glpi
            d3 = (State%rlam/State%rsquig)*(Model%apar/(State%rsquig**2))* &
                (1. - 2.*((State%rlam/State%rsquig)**2.))*State%blittlerlamb**2./4./pi
            d4 = ((State%rlam/State%rsquig)**2.)*((Model%apar**2)/(State%rsquig**2.) + &
                                    2*(Model%apar/State%rsquig))* &
                                    (State%blittlethlamb**2. + State%blittlephlamb**2.)/4./pi/State%rlam

            ! Equilibrium Density
            Solution%dens(State%ri,:,:)  = (d1 + d2 + d3 + d4) / (State%F * Rsolar_cgs)
            deallocate (d1, d2, d3, d4)

            ! ;  put in background density (atmospheric constants)
            ! ;  first radial power law HE
            ! ; see Gibson et al 1999 - eq 3
            ! ; https://ui.adsabs.harvard.edu/link_gateway/1999JGR...104.9691G/doi:10.1029/98JA02681 
            ! ;
            State%densin = (Model%aa + Model%cc + Model%ee)*State%rsquig**(-3.)
            State%densout = (Model%aa*State%rsquig**(-Model%bb) + &
                                Model%cc*State%rsquig**(-Model%dd) + &
                                Model%ee*State%rsquig**(-Model%ff)) 

            if (State%inside_count .gt. 0) then 
                where (Solution%inside_mask(State%ri, :, :) .gt. 0.1) 
                    State%densback = State%densin
                end where
            end if
            if (State%outside_count .gt. 0) then 
                where (Solution%inside_mask(State%ri, :, :) .lt. 1.) 
                    State%densback = State%densout
                end where
            end if
            State%DensbackHEonly = State%densback
            ! ; 
            ! ; now add total pressure continuity part inside bubble
            ! ;
            if (State%inside_count .gt. 0) then 
                where (Solution%inside_mask(State%ri, :, :) .gt. 0.1) 
                    State%densback = State%densback + State%dbackin*Model%phiss**3.
                end where
            end if
            ! Add background to solution density
            Solution%dens(State%ri,:,:) = Solution%dens(State%ri,:,:) + State%Densback
            ! ;
            ! ;  put dens in ss coords
            ! ;
            Solution%dens(State%ri,:,:) = Solution%dens(State%ri,:,:)/Model%phiss**3.
            State%Densback = State%Densback/Model%phiss**3.

        end subroutine calcDensity

        !> Calculate Gibson-Low Solution Pressure
        !>
        !>
        subroutine calcPressure(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution

            ! ;
            ! ; calculate pressure from Gibson and Low eq. A18
            ! ;
            Solution%Pres(State%ri, :, :) = (State%rlam**2/State%rsquig**2.) & 
                                            * (1. - (State%rlam**2/State%rsquig**2.)) &
                                            * (State%blittlerlamb**2.)/(8.*pi) &
                                            + (State%rlam**2./State%rsquig**2.)*State%glpi
            ! ; now we have to match total pressure Pmag+Pgas at the bubble interface
            ! ;  note that the INNER bubble solution is defined so that Pmag and Pgas are zero there
            ! ;  so, we calculate for each rlam, at the theta of the bubble boundary, the OUTER
            ! ;  solution total pressure, and add it to the INNER gas pressure 
            ! ;  (for all theta within the bubble at that rlam).  Then we must also
            ! ;  add a similar density increment to the INNER density, because the OUTER solution
            ! ; is potential and so has density and pressure in HE balance 
            ! ;
            State%presback = 0.
            State%presout = 0.

            if (State%inside_count .gt. 0) then
                if (Model%isLoud .and. Model%isDebug) write(*,*) 'calling calcInside'
                call calcInside(Model, State, Solution)
                where (Solution%inside_mask(State%ri, :, :) .gt. 0.1) State%presback = State%Pbackin*Model%phiss**4.
            end if
            ! ;  now we put in additional hydrostatic background pressure to keep
            ! ;  things positive
            State%presout = (Model%aa/(Model%bb + 1.))*State%rsquig**(-Model%bb - 1.) &
                            + (Model%cc/(Model%dd + 1.))*State%rsquig**(-Model%dd - 1.) &
                            + (Model%ee/(Model%ff + 1.))*State%rsquig**(-Model%ff - 1.)

            State%presout = State%presout*GMm/Rsolar_cgs
            ! ;
            ! ; see Gibson et al 1998, eqs 3 and 4 - note we are assuming alpha (He) = 0
            ! ;
            State%presin = 0.
            State%presin_0 = 0.
            ! ; presin and densin only matter for CME solution
            ! ; and are put in as a means of stopping the background density
            ! ; from blowing up as the structure self-similarly expands
            ! ; since density scales as phiss^-3 and pressure as phiss^-4
            ! ; and rsquig = rpB/phiss
            ! ; however, note this means that for rsquig less than 1, the
            ! ; background density won't fall off as defined by aa,bb,cc,dd,ee,f
            ! ; but rather by an r^3 power law.  Thus, if fitting a magnetostatic
            ! ; solution it is best to set phiss=1
            ! ;
            State%presin = (1./4.)*(Model%aa + Model%cc + Model%ee)*State%rsquig**(-4.)
            State%presin_0 = -(1./4.)*(Model%aa + Model%cc + Model%ee)
            State%presin_0 = State%presin_0 + (Model%aa/(Model%bb + 1.)) + &
                             (Model%cc/(Model%dd + 1.)) + &
                             (Model%ee/(Model%ff + 1.))
            State%presin = (State%presin + State%presin_0)*GMm/Rsolar_cgs
            ! ; thus presin and presout will match smoothly at rsquig = 1
            ! ;  
            where (State%rsquig .lt. 1.) State%presback = State%presback + State%presin
            where (State%rsquig .ge. 1.) State%presback = State%presback + State%presout

            Solution%Pres(State%ri, :, :) = Solution%Pres(State%ri, :, :) + State%presback
            
            !; self-similar transform
            Solution%Pres(State%ri, :, :) = Solution%Pres(State%ri, :, :)/Model%phiss**4.
        end subroutine calcPressure

        !> Calculate Outside (of bubble) Gibson-Low Solution
        !>
        !>
        subroutine calcOutside(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution

            real(rp) :: r1
            real(rp), dimension(:, :), allocatable :: al1, al2, psi0, f, g, h
            real(rp), dimension(:, :), allocatable :: st, mu, sph, cph
            real(rp), dimension(:, :), allocatable :: psi1, i, psi3, dpsi0dr, dpsi0dth
            real(rp), dimension(:, :), allocatable :: d2psi0dth2, dpsi0dr2, dpsi0drdth, dfdr
            real(rp), dimension(:, :), allocatable :: dfdmu, dfdrdmu, dfdr2, dgdr, dgdmu
            real(rp), dimension(:, :), allocatable :: dgdrdmu, dgdr2, dhdr, dhdmu, d2hdmu2
            real(rp), dimension(:, :), allocatable :: dhdrdmu, dhdr2, dpsi1dr, dpsi1dr2
            real(rp), dimension(:, :), allocatable :: dpsi1dmu, d2psi1dmu2, dpsi1drdmu, dpsi1dth
            real(rp), dimension(:, :), allocatable :: dpsi1drdth, d2psi1dth2, didr, didmu
            real(rp), dimension(:, :), allocatable :: didrdmu, didr2, dpsi3dr, dpsi3dmu
            real(rp), dimension(:, :), allocatable :: d2psi3dmu2, dpsi3drdmu, dpsi3dr2, dpsi3dth
            real(rp), dimension(:, :), allocatable :: dpsi3drdth, d2psi3dth2, dAdr, dAdth
            real(rp), dimension(:, :), allocatable :: dAdrdth, dAdr2, dAdth2, brlambout1
            real(rp), dimension(:, :), allocatable :: bthlambout1, bphlambout1, dbth1dr, dbr1dth
            real(rp), dimension(:, :), allocatable :: jrlambout1, jthlambout1, jphlambout1
            real(rp), dimension(:, :), allocatable :: blittleX1, blittleY1
            real(rp), dimension(:, :), allocatable :: blittleZ1, jlittleX1, jlittleY1, jlittleZ1
            real(rp), dimension(:, :), allocatable :: blittleY2, blittleZ, jlittleY2, jlittleZ
            real(rp), dimension(:, :), allocatable :: blittleX, blittleY, jlittleX, jlittleY
            real(rp), dimension(:, :), allocatable :: streal, ctreal, spreal, cpreal, dbr1dr
            real(rp), dimension(:, :), allocatable :: dbrdr
            real(rp), dimension(:, :), allocatable :: yt2


            allocate (st(State%js:State%je, State%ks:State%ke), mu(State%js:State%je, State%ks:State%ke), sph(State%js:State%je, State%ks:State%ke),  cph(State%js:State%je, State%ks:State%ke))
            allocate (al1(State%js:State%je, State%ks:State%ke), al2(State%js:State%je, State%ks:State%ke), psi0(State%js:State%je, State%ks:State%ke), f(State%js:State%je, State%ks:State%ke), g(State%js:State%je, State%ks:State%ke), h(State%js:State%je, State%ks:State%ke))
            allocate (psi1(State%js:State%je, State%ks:State%ke), i(State%js:State%je, State%ks:State%ke), psi3(State%js:State%je, State%ks:State%ke), dpsi0dr(State%js:State%je, State%ks:State%ke), dpsi0dth(State%js:State%je, State%ks:State%ke))
            allocate (d2psi0dth2(State%js:State%je, State%ks:State%ke), dpsi0dr2(State%js:State%je, State%ks:State%ke), dpsi0drdth(State%js:State%je, State%ks:State%ke), dfdr(State%js:State%je, State%ks:State%ke))
            allocate (dfdmu(State%js:State%je, State%ks:State%ke), dfdrdmu(State%js:State%je, State%ks:State%ke), dfdr2(State%js:State%je, State%ks:State%ke), dgdr(State%js:State%je, State%ks:State%ke), dgdmu(State%js:State%je, State%ks:State%ke))
            allocate (dgdrdmu(State%js:State%je, State%ks:State%ke), dgdr2(State%js:State%je, State%ks:State%ke), dhdr(State%js:State%je, State%ks:State%ke), dhdmu(State%js:State%je, State%ks:State%ke), d2hdmu2(State%js:State%je, State%ks:State%ke))
            allocate (dhdrdmu(State%js:State%je, State%ks:State%ke), dhdr2(State%js:State%je, State%ks:State%ke), dpsi1dr(State%js:State%je, State%ks:State%ke), dpsi1dr2(State%js:State%je, State%ks:State%ke))
            allocate (dpsi1dmu(State%js:State%je, State%ks:State%ke), d2psi1dmu2(State%js:State%je, State%ks:State%ke), dpsi1drdmu(State%js:State%je, State%ks:State%ke), dpsi1dth(State%js:State%je, State%ks:State%ke))
            allocate (dpsi1drdth(State%js:State%je, State%ks:State%ke), d2psi1dth2(State%js:State%je, State%ks:State%ke), didr(State%js:State%je, State%ks:State%ke), didmu(State%js:State%je, State%ks:State%ke))
            allocate (didrdmu(State%js:State%je, State%ks:State%ke), didr2(State%js:State%je, State%ks:State%ke), dpsi3dr(State%js:State%je, State%ks:State%ke), dpsi3dmu(State%js:State%je, State%ks:State%ke))
            allocate (d2psi3dmu2(State%js:State%je, State%ks:State%ke), dpsi3drdmu(State%js:State%je, State%ks:State%ke), dpsi3dr2(State%js:State%je, State%ks:State%ke), dpsi3dth(State%js:State%je, State%ks:State%ke))
            allocate (dpsi3drdth(State%js:State%je, State%ks:State%ke), d2psi3dth2(State%js:State%je, State%ks:State%ke), dAdr(State%js:State%je, State%ks:State%ke), dAdth(State%js:State%je, State%ks:State%ke))
            allocate (dAdrdth(State%js:State%je, State%ks:State%ke), dAdr2(State%js:State%je, State%ks:State%ke), dAdth2(State%js:State%je, State%ks:State%ke), brlambout1(State%js:State%je, State%ks:State%ke))
            allocate (bthlambout1(State%js:State%je, State%ks:State%ke), bphlambout1(State%js:State%je, State%ks:State%ke), dbth1dr(State%js:State%je, State%ks:State%ke), dbr1dth(State%js:State%je, State%ks:State%ke))
            allocate (jrlambout1(State%js:State%je, State%ks:State%ke), jthlambout1(State%js:State%je, State%ks:State%ke), jphlambout1(State%js:State%je, State%ks:State%ke))
            allocate (blittleX1(State%js:State%je, State%ks:State%ke), blittleY1(State%js:State%je, State%ks:State%ke))
            allocate (blittleZ1(State%js:State%je, State%ks:State%ke), jlittleX1(State%js:State%je, State%ks:State%ke), jlittleY1(State%js:State%je, State%ks:State%ke), jlittleZ1(State%js:State%je, State%ks:State%ke))
            allocate (blittleY2(State%js:State%je, State%ks:State%ke), blittleZ(State%js:State%je, State%ks:State%ke), jlittleY2(State%js:State%je, State%ks:State%ke), jlittleZ(State%js:State%je, State%ks:State%ke))
            allocate (blittleX(State%js:State%je, State%ks:State%ke), blittleY(State%js:State%je, State%ks:State%ke), jlittleX(State%js:State%je, State%ks:State%ke), jlittleY(State%js:State%je, State%ks:State%ke))
            allocate (streal(State%js:State%je, State%ks:State%ke), ctreal(State%js:State%je, State%ks:State%ke), spreal(State%js:State%je, State%ks:State%ke), cpreal(State%js:State%je, State%ks:State%ke), dbr1dr(State%js:State%je, State%ks:State%ke))
            allocate (dbrdr(State%js:State%je, State%ks:State%ke))
            allocate (yt2(State%js:State%je, State%ks:State%ke))

            ! ;
            ! ;  Calculate the external (to bubble) field, Appendix B1
            ! ;
            ! ;  we are going to rotate the outside field so that its axis
            ! ;  of symmetry lies along the vector r1 defined by offset
            ! ;  xo -- basically a rotation about the y axis 
            ! ;  we do this by transforming through cartesian coords
            ! ;  note we have now changed naming conventions to be
            ! ;  consistent with right-hand coordinate system
            ! ;
            ! ;

            State%xtr = State%rpb*sin(State%thpb)*cos(State%phpb)
            State%ytr = State%rpb*sin(State%thpb)*sin(State%phpb)
            State%ztr = State%rpb*cos(State%thpb)

            !; now transform to spherical coords (radius is unchanged)
            State%rout = State%rpb

            ! ;
            ! ; the tilde coords are the new x and z with up defined as above
            ! ; rotated 90 about z, and then 90 about y
            ! ; from the physical coordinate system
            ! ; (so that the r1 vector points up and the outside current sheet 
            ! ;  is in the y-z plane)
            ! ;  this is because the inner solution has its axis of symmetry naturally
            ! ;  at the equator, and the outer at the pole
            ! ;
            ! ; first rotate about z
            State%xtilde = -State%ytr
            yt2 = State%xtr

            !; now about xtilde
            State%ytilde = -State%ztr
            State%ztilde = yt2

            State%rat = (State%ztilde/State%rout)
            where (State%rat .gt. 1) State%rat = 1
            where (State%rat .lt. -1) State%rat = -1

            State%thout = acos(State%rat)
            State%phout = atan2(State%ytilde, State%xtilde)
    
            ! ;
            ! ;
            ! ; make sure phout is betweeen 0. and 2.*!dpi
            ! ;
            where (State%phout .lt. 0.) State%phout = 2.*pi + State%phout
            where (State%phout .gt. 2.*pi) State%phout = State%phout - 2.*pi
            
            st = sin(State%thout)
            mu = cos(State%thout)
            sph = sin(State%phout)
            cph = cos(State%phout)

            call zero2tiny2d(st)
            call zero2tiny2d(mu)
            call zero2tiny2d(sph)
            call zero2tiny2d(cph)
    
            r1 = sqrt(Model%xo*Model%xo)

            ! ;
            ! ;    break down bits of Stream function
            ! ;
            al1 = r1*r1 - Model%rbub*Model%rbub
            al2 = 2.*r1*r1 - Model%rbub*Model%rbub
            ! ;
            ! ; eq. B2
            psi0 = mu
            ! ;
            ! ; eq. B3
            f = r1*(State%rout*State%rout + al1) - State%rout*mu*al2
            g = al1*al1 + r1*r1*State%rout*State%rout - 2.*State%rout*r1*al1*mu
            h = 1./(sqrt(g))
            psi1 = (-1/Model%rbub)*f*h
            ! ;
            ! ; eq. B4
            i = State%rout*State%rout + r1*r1 - 2.*State%rout*r1*mu
            psi3 = (1./Model%rbub)*sqrt(i)
            ! ;
            ! ; eq. B5
            State%streamout = psi0 + psi1 + psi3

            ! ;
            ! ;  now take derivatives for magnetic fields, etc.
            ! ;
            dpsi0dr = 0.
            dpsi0dth = -st
            d2psi0dth2 = -mu
            dpsi0dr2 = 0.
            dpsi0drdth = 0.

            dfdr = 2.*r1*State%rout - mu*al2
            dfdmu = -State%rout*al2
            dfdrdmu = -al2
            dfdr2 = 2.*r1

            dgdr = 2.*State%rout*r1*r1 - 2.*r1*al1*mu
            dgdmu = -2.*State%rout*r1*al1
            dgdrdmu = -2.*r1*al1
            dgdr2 = 2.*r1*r1

            dhdr = -dgdr/2./(g**1.5)
            dhdmu = -dgdmu/2./(g**1.5)
            d2hdmu2 = 3.*dgdmu/4./(g**2.5)
            dhdrdmu = -dgdrdmu/2./(g**1.5) + 3.*dgdmu*dgdr/4./(g**2.5)
            dhdr2 = -dgdr2/2/(g**1.5) + 3*dgdr*dgdr/4./(g**2.5)

            dpsi1dr = (-1./Model%rbub)*(dfdr*h + f*dhdr)
            dpsi1dr2 = (-1./Model%rbub)*(2.*dfdr*dhdr + dfdr2*h + f*dhdr2)

            dpsi1dmu = (-1/Model%rbub)*(dfdmu*h + f*dhdmu)
            d2psi1dmu2 = (-1./Model%rbub)*(2.*dfdmu*dhdmu + f*d2hdmu2)
            dpsi1drdmu = (-1./Model%rbub)*(dfdmu*dhdr + dfdr*dhdmu + dfdrdmu*h + f*dhdrdmu)

            dpsi1dth = -st*dpsi1dmu
            dpsi1drdth = -st*dpsi1drdmu
            d2psi1dth2 = -mu*dpsi1dmu + (st**2.)*d2psi1dmu2

            didr = 2*State%rout - 2*r1*mu
            didmu = -2*State%rout*r1
            didrdmu = -2*r1
            didr2 = 2.

            dpsi3dr = (1./Model%rbub)*didr/2./sqrt(i)
            dpsi3dmu = (1./Model%rbub)*didmu/2./sqrt(i)
            d2psi3dmu2 = -(1./Model%rbub)*didmu/4./(i**(1.5))
            dpsi3drdmu = (1./Model%rbub)*(didrdmu/2./sqrt(i) - didmu*didr/4./(i**1.5))
            dpsi3dr2 = (1./Model%rbub)*(didr2/2./sqrt(i) - (didr)**2./4./(i**1.5))

            dpsi3dth = -st*dpsi3dmu
            dpsi3drdth = -st*dpsi3drdmu
            d2psi3dth2 = -mu*dpsi3dmu + (st**2.)*d2psi3dmu2

            dAdr = dpsi0dr + dpsi1dr + dpsi3dr
            dAdth = dpsi0dth + dpsi1dth + dpsi3dth
            dAdrdth = dpsi0drdth + dpsi1drdth + dpsi3drdth
            dAdr2 = dpsi0dr2 + dpsi1dr2 + dpsi3dr2
            dAdth2 = d2psi0dth2 + d2psi1dth2 + d2psi3dth2
            
            ! ;
            ! ;  now calculate field
            ! ;
            ! ;  eq B1
            
            brlambout1 = dAdth/(State%rout*State%rout*st)
            bthlambout1 = -dAdr/(State%rout*st)
            bphlambout1 = 0.

            ! where (brlambout1 /= brlambout1) brlambout1 = 0.
            ! where (bthlambout1 /= bthlambout1) bthlambout1 = 0.
            ! where (bphlambout1 /= bphlambout1) bphlambout1 = 0.

            !; and calculate currents

            dbth1dr = (dAdr/State%rout - dAdr2)/(State%rout*st)
            dbr1dth = ((-mu/st)*dAdth + dAdth2)/(State%rout*State%rout*st)

            jrlambout1 = 0.
            jthlambout1 = 0.
            jphlambout1 = (1./State%rout)*(bthlambout1 + State%rout*dbth1dr - dbr1dth)

            ! ;
            ! ; But now we need to get blittle and derivatives of blittle with respect to
            ! ; the physical coord r, which we need for density and pressure  (eq. A18-A19)
            ! ;  We do so by going through the cartesian coordinates.
            ! ;

            blittleX1 = brlambout1*st*cph + bthlambout1*mu*cph
            blittleY1 = brlambout1*st*sph + bthlambout1*mu*sph
            blittleZ1 = brlambout1*mu - bthlambout1*st
            jlittleX1 = jrlambout1*st*cph + jthlambout1*mu*cph
            jlittleY1 = jrlambout1*st*sph + jthlambout1*mu*sph
            jlittleZ1 = jrlambout1*mu - jthlambout1*st

            jphlambout1 = 0.
            ! ;
            ! ; but now we have to translate to the xtilde-rotated frame for x and z
            ! ;
            blittleY2 = blittleZ1
            blittleZ = -blittleY1
            jlittleY2 = jlittleZ1
            jlittleZ = -jlittleY1
            jlittleX1 = 0.
            jlittleZ1 = 0.
            ! ;
            ! ; and now about the z-rotated frame for x and y
            ! ;
            blittleX = blittleY2
            blittleY = -blittleX1
            jlittleX = jlittleY2
            jlittleY = -jlittleX1
            jlittleX1 = 0.
            jlittleZ1 = 0.
            ! ;
            ! ; and finally back to real physical spherical coords
            ! ;
            streal = sin(State%thpb)
            ctreal = cos(State%thpb)
            spreal = sin(State%phpb)
            cpreal = cos(State%phpb)

            call zero2tiny2d(streal)
            call zero2tiny2d(ctreal)
            call zero2tiny2d(spreal)
            call zero2tiny2d(cpreal)

            State%brlambout = blittleY*streal*spreal + blittleX*streal*cpreal + blittleZ*ctreal
            State%bthlambout = blittleY*ctreal*spreal + blittleX*ctreal*cpreal - blittleZ*streal
            State%bphlambout = -blittleX*spreal + blittleY*cpreal
            State%jrlambout = jlittleY*streal*spreal + jlittleX*streal*cpreal + jlittleZ*ctreal
            State%jthlambout = jlittleY*ctreal*spreal + jlittleX*ctreal*cpreal - jlittleZ*streal
            State%jphlambout = -jlittleX*spreal + jlittleY*cpreal

            jlittleX = 0.
            jlittleY = 0.
            jlittleZ = 0.
            ! ;  now calculate tderivout (needed for density  A19)
            ! ;
            dbr1dr = (-2.*dAdth/State%rout + dAdrdth)/(State%rout*State%rout*st)
            State%tderivout = (brlambout1*dbr1dr + bthlambout1*dbth1dr)/4./pi

            ! ;
            ! ;   and finally flip the sign about the equator
            ! ;
            ! ;	south = where(phout gt 0.  and phout lt double(!dpi))
            where ((State%phout .gt. pi) .and. (State%phout .lt. 2.*pi))
                State%brlambout = -1*State%brlambout
                State%bthlambout = -1*State%bthlambout
                State%bphlambout = -1*State%bphlambout
            end where

            ! Take care of NaNs
            ! where (State%brlambout /= State%brlambout) State%brlambout = 0.
            ! where (State%bthlambout /= State%bthlambout) State%bthlambout = 0.
            ! where (State%bphlambout /= State%bphlambout) State%bphlambout = 0.
            ! where (State%tderivout /= State%tderivout) State%tderivout = 0.
            ! where (State%streamout /= State%streamout) State%streamout = 0.

            deallocate (al1, al2, psi0, f, g, h, st, mu, sph, cph)
            deallocate (psi1, i, psi3, dpsi0dr, dpsi0dth)
            deallocate (d2psi0dth2, dpsi0dr2, dpsi0drdth, dfdr)
            deallocate (dfdmu, dfdrdmu, dfdr2, dgdr, dgdmu)
            deallocate (dgdrdmu, dgdr2, dhdr, dhdmu, d2hdmu2)
            deallocate (dhdrdmu, dhdr2, dpsi1dr, dpsi1dr2)
            deallocate (dpsi1dmu, d2psi1dmu2, dpsi1drdmu, dpsi1dth)
            deallocate (dpsi1drdth, d2psi1dth2, didr, didmu)
            deallocate (didrdmu, didr2, dpsi3dr, dpsi3dmu)
            deallocate (d2psi3dmu2, dpsi3drdmu, dpsi3dr2, dpsi3dth)
            deallocate (dpsi3drdth, d2psi3dth2, dAdr, dAdth)
            deallocate (dAdrdth, dAdr2, dAdth2, brlambout1)
            deallocate (bthlambout1, bphlambout1, dbth1dr, dbr1dth)
            deallocate (jrlambout1, jthlambout1, jphlambout1)
            deallocate (blittleX1, blittleY1)
            deallocate (blittleZ1, jlittleX1, jlittleY1, jlittleZ1)
            deallocate (blittleY2, blittleZ, jlittleY2, jlittleZ)
            deallocate (blittleX, blittleY, jlittleX, jlittleY)
            deallocate (streal, ctreal, spreal, cpreal, dbr1dr)
            deallocate (dbrdr)
            deallocate (yt2)
        end subroutine calcOutside

        !> Calculate Gibson-Low Fields
        !> 
        !> Calls calcOutside to determine outside Stream function
        subroutine calcFields(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            real(rp), dimension(:, :), allocatable :: st, ct, sph, cph, gfun, gfunot
            real(rp), dimension(:, :), allocatable :: dSdT, dgfundr
            real(rp), dimension(:, :), allocatable :: dSdR, Q
            real(rp), dimension(:, :), allocatable :: dbPdt, dbPdR, d2gfundr2
            real(rp), dimension(:, :), allocatable :: d2SdR2, d2SdT2, dbTdR, dbRdT
            real(rp), dimension(:, :), allocatable :: ctreal, streal, cpreal, spreal
            real(rp), dimension(:, :), allocatable :: blittleR, blittleT, blittleP
            real(rp), dimension(:, :), allocatable :: blittleX, blittleY1, blittleZ1
            real(rp), dimension(:, :), allocatable :: blittley, blittlez
            real(rp), dimension(:, :), allocatable :: jlittleR, jlittleT, jlittleP
            real(rp), dimension(:, :), allocatable :: jlittleX, jlittleY1, jlittleZ1
            real(rp), dimension(:, :), allocatable :: jlittley, jlittlez
            real(rp), dimension(:, :), allocatable :: jlittlerlamb, jlittlethlamb, jlittlephlamb
            real(rp) :: dlamdr
 
            allocate (st(State%js:State%je, State%ks:State%ke), ct(State%js:State%je, State%ks:State%ke))
            allocate (sph(State%js:State%je, State%ks:State%ke), cph(State%js:State%je, State%ks:State%ke))
            allocate (gfun(State%js:State%je, State%ks:State%ke), gfunot(State%js:State%je, State%ks:State%ke))
            allocate (dSdT(State%js:State%je, State%ks:State%ke), dgfundr(State%js:State%je, State%ks:State%ke))
            allocate (dSdR(State%js:State%je, State%ks:State%ke), blittleR(State%js:State%je, State%ks:State%ke))
            allocate (blittleT(State%js:State%je, State%ks:State%ke), Q(State%js:State%je, State%ks:State%ke))
            allocate (blittleP(State%js:State%je, State%ks:State%ke))
            allocate (dbPdt(State%js:State%je, State%ks:State%ke), dbPdR(State%js:State%je, State%ks:State%ke))
            allocate (d2gfundr2(State%js:State%je, State%ks:State%ke))
            allocate (d2SdR2(State%js:State%je, State%ks:State%ke), d2SdT2(State%js:State%je, State%ks:State%ke))
            allocate (dbTdR(State%js:State%je, State%ks:State%ke), dbRdT(State%js:State%je, State%ks:State%ke))
            allocate (jlittleR(State%js:State%je, State%ks:State%ke))
            allocate (jlittleT(State%js:State%je, State%ks:State%ke), jlittleP(State%js:State%je, State%ks:State%ke))
            allocate (blittleX(State%js:State%je, State%ks:State%ke), blittleY1(State%js:State%je, State%ks:State%ke))
            allocate (blittleZ1(State%js:State%je, State%ks:State%ke), jlittleX(State%js:State%je, State%ks:State%ke))
            allocate (jlittleY1(State%js:State%je, State%ks:State%ke), jlittleZ1(State%js:State%je, State%ks:State%ke))
            allocate (ctreal(State%js:State%je, State%ks:State%ke), streal(State%js:State%je, State%ks:State%ke))
            allocate (cpreal(State%js:State%je, State%ks:State%ke), spreal(State%js:State%je, State%ks:State%ke))
            allocate (blittley(State%js:State%je, State%ks:State%ke), blittlez(State%js:State%je, State%ks:State%ke))
            allocate (jlittley(State%js:State%je, State%ks:State%ke), jlittlez(State%js:State%je, State%ks:State%ke))

            st = 0.0; ct = 0.0; sph = 0.0; cph = 0.0; gfun = 0.0; gfunot = 0.0
            dSdT = 0.0; dgfundr = 0.0
            dSdR = 0.0; blittleR = 0.0; blittleT = 0.0; Q = 0.0; blittleP = 0.0
            dbPdt = 0.0; dbPdR = 0.0; d2gfundr2 = 0.0
            d2SdR2 = 0.0; d2SdT2 = 0.0; dbTdR = 0.0; dbRdT = 0.0; jlittleR = 0.0
            jlittleT = 0.0; jlittleP = 0.0; blittleX = 0.0; blittleY1 = 0.0
            blittleZ1 = 0.0; jlittleX = 0.0; jlittleY1 = 0.0; jlittleZ1 = 0.0
            ctreal = 0.0; streal = 0.0; cpreal = 0.0; spreal = 0.0
            blittley = 0.0; blittlez = 0.0; jlittley = 0.0; jlittlez = 0.0

            st = sin(State%thcap)
            ct = cos(State%thcap)
            sph = sin(State%phcap)
            cph = cos(State%phcap)

            call zero2tiny2d(st) ! where (abs(st) .lt. tiny) st=tiny*sign(st)
            call zero2tiny2d(ct) ! where (abs(ct) .lt. tiny) ct=tiny*sign(ct)
            call zero2tiny2d(sph) ! where (abs(sph) .lt. tiny) sph=tiny*sign(sph)
            call zero2tiny2d(cph) ! where (abs(cph) .lt. tiny) cph=tiny*sign(cph)

            ! eq. B8
            gfun = sin(Model%alnot*State%rcap)/(Model%alnot*State%rcap) - cos(Model%alnot*State%rcap)
            gfunot = sin(Model%alnot*Model%rbub)/(Model%alnot*Model%rbub) - cos(Model%alnot*Model%rbub)

            ! Stream function A in the region Rcap lt rbub
            ! eq. B7
            where (Solution%inside_mask(State%ri, :, :) .gt. 0.1) 
                State%stream = (4.*Model%ao*pi/Model%alnot**2.)*(gfun*Model%rbub**2./gfunot - State%rcap**2.)*st**2.
            end where
            ! in the region Rcap ge rbub
            !  call routine that calculates Stream function, etc
            !  APPENDIX B1

            if (State%outside_count .gt. 0) then
                if (Model%isLoud .and. Model%isDebug) write(*,*) 'giblow_fieldcalc: calling calcOutside'
                call calcOutside(Model, State, Solution)
                ! Set Outside Stream function
                where (Solution%inside_mask(State%ri, :, :) .lt. 1.) State%stream = Model%outScale*State%streamout
            end if
            
            ! ;
            ! ;  define pressure
            ! ;
            State%glpi = Model%ao*State%stream + Model%Pio
            ! Set Outside Pressure
            if (State%outside_count .gt. 0) where (Solution%inside_mask(State%ri, :, :) .lt. 1.) State%glpi = Model%Pio
            
            ! ;
            ! ; now calculate little b (blittle), that goes with this Stream function
            ! ; and pressure in these bubble coordinates
            ! ; To do this, we need to use the derivatives of Stream with respect
            ! ; to the cap coords.
            ! ;
            dSdT = State%stream*2.*ct/st
            dgfundr = cos(Model%alnot*State%rcap)/State%rcap &
                        - sin(Model%alnot*State%rcap)/(State%rcap**2.*Model%alnot) &
                        + Model%alnot*sin(Model%alnot*State%rcap)
            dSdR = (4.*Model%ao*pi/Model%alnot**2.)*(dgfundr*Model%rbub**2./gfunot - 2.*State%rcap)*st**2.
            
            !;  eq. B6
            blittleR = dSdT/(State%rcap**2.*st)
            blittleT = -dSdR/(State%rcap*st)
            ! ;  ADD THE PHI FIELD!!!!!
            ! ; 
            ! ; NOTE! negative sign on BlittleP is to make a left-handed wind
            ! ;  about the apple core
            Q = Model%alnot*State%stream
            blittleP = -Q/(State%rcap*st)
            ! TODO: Check NaN
            !where (blittleR /= blittleR) blittleR = 0.
            !where (blittleT /= blittleT) blittleT = 0.
            !where (blittleP /= blittleP) blittleP = 0.
            ! ;
            ! ; also calculate currents by taking curl
            ! ;
            dbPdt = -(1./State%rcap/st)*Model%alnot*dSdT
            dbPdt = dbPdt + (ct/State%rcap/st/st)*Q

            dbPdR = -(1./State%rcap/st)*Model%alnot*dSdR
            dbPdR = dbPdR + (1./State%rcap/State%rcap/st)*Q

            d2gfundr2 = 2.*gfun/State%rcap/State%rcap - Model%alnot*Model%alnot*gfun
            
            d2SdR2 = ((((Model%rbub**2)/gfunot)*d2gfundr2 - 2.)/ &
                        (((Model%rbub**2)/gfunot)*gfun - State%rcap**2.))*State%stream
            d2SdT2 = (4.*((ct/st)**2) - 2./st/st)*State%stream
            
            dbTdR = (1./State%rcap/State%rcap/st)*dSdR - (1./State%rcap/st)*d2SdR2
            dbRdT = -(ct/State%rcap/State%rcap/st/st)*dSdT + (1/State%rcap/State%rcap/st)*d2SdT2

            jlittleR = (1./State%rcap/st)*(ct*blittleP + st*dbPdt)
            jlittleT = -(1./State%rcap)*(blittleP + State%rcap*dbPdR)
            jlittleP = (1./State%rcap)*(blittleT + State%rcap*dbTdR - dbRdT)


            ! ;
            ! ; But now we need to transform to blittle coordinates with respect to
            ! ; the physical coord r, which we need for density and pressure (eq. A19)
            ! ; (note the derivative calculation for this interior field happens in giblow - tderiv)
            ! ;  We do so by going through the cartesian coordinates.
            ! ;  note this has also now been rewritten in right-hand coordinates
            ! ;
            blittleX = blittleR*st*cph + blittleT*ct*cph - blittleP*sph
            blittleY1 = blittleR*st*sph + blittleT*ct*sph + blittleP*cph
            blittleZ1 = blittleR*ct - blittleT*st
            
            jlittleX = jlittleR*st*cph + jlittleT*ct*cph - jlittleP*sph
            jlittleY1 = jlittleR*st*sph + jlittleT*ct*sph + jlittleP*cph
            jlittleZ1 = jlittleR*ct - jlittleT*st

            ! ;
            ! ; but now we have to translate to the sigma rotated (about x axis) frame for y and z
            ! ;  note clockwise rotation-- moves back into physical coordinates
            ! ;
            blittleY = cos(Model%sigma)*blittleY1 + sin(Model%sigma)*blittleZ1
            blittleZ = cos(Model%sigma)*blittleZ1 - sin(Model%sigma)*blittleY1
            jlittleY = cos(Model%sigma)*jlittleY1 + sin(Model%sigma)*jlittleZ1
            jlittleZ = cos(Model%sigma)*jlittleZ1 - sin(Model%sigma)*jlittleY1

            ctreal = cos(State%thpb)
            streal = sin(State%thpb)
            cpreal = cos(State%phpb)
            spreal = sin(State%phpb)

            call zero2tiny2d(ctreal)
            call zero2tiny2d(streal)
            call zero2tiny2d(cpreal)
            call zero2tiny2d(spreal)

            ! ;
            ! ;  these are the "little b" magnetic field components- the little b r
            ! ;  is needed in calculating pressure and density in the pB coords
            ! ;  see Appendix A, esp. eqs. A18-A19
            State%blittlerlamb = blittleY*streal*spreal + blittleX*streal*cpreal + blittleZ*ctreal
            State%blittlethlamb = blittleY*ctreal*spreal + blittleX*ctreal*cpreal - blittleZ*streal
            State%blittlephlamb = -blittleX*spreal + blittleY*cpreal
            State%jlittlerlamb = jlittleY*streal*spreal + jlittleX*streal*cpreal + jlittleZ*ctreal
            State%jlittlethlamb = jlittleY*ctreal*spreal + jlittleX*ctreal*cpreal - jlittleZ*streal
            State%jlittlephlamb = -jlittleX*spreal + jlittleY*cpreal

            !;   we need to put in the proper values for the outside field
            if (State%outside_count .ne. 0) then
                where (Solution%inside_mask(State%ri, :, :) .lt. 1.)
                    State%blittlerlamb = Model%outScale*State%brlambout
                    State%blittlethlamb = Model%outScale*State%bthlambout
                    State%blittlephlamb = Model%outScale*State%bphlambout
                    State%jlittlerlamb = Model%outScale*State%jrlambout
                    State%jlittlethlamb = Model%outScale*State%jthlambout
                    State%jlittlephlamb = Model%outScale*State%jphlambout
                end where
            end if

            ! ;  and finally, here are the magnetic field coordinates in the rpB,thetapB
            ! ;  phipB coordinate systems !  we have to also get rid of selfsim stuff
            ! ;
            ! ;  eqs.  A2-A4
            
            ! ;  dlamdr is one, because rlam = rsquig - apar.
            dlamdr = 1.

            Solution%b(State%ri, :, :, XDIR) = State%blittlerlamb*(State%rlam/State%rsquig)**2./Model%phiss**2.
            Solution%b(State%ri, :, :, YDIR) = State%blittlethlamb*(State%rlam/State%rsquig)*dlamdr/Model%phiss**2.
            Solution%b(State%ri, :, :, ZDIR) = State%blittlephlamb*(State%rlam/State%rsquig)*dlamdr/Model%phiss**2.
            Solution%j(State%ri, :, :, XDIR) = State%jlittlerlamb*(State%rlam/State%rsquig)**2./Model%phiss**2.
            Solution%j(State%ri, :, :, YDIR) = State%jlittlethlamb*(State%rlam/State%rsquig)*dlamdr/Model%phiss**2.
            Solution%j(State%ri, :, :, ZDIR) = State%jlittlephlamb*(State%rlam/State%rsquig)*dlamdr/Model%phiss**2.
            
            deallocate (sph, cph, gfun, gfunot)
            deallocate (dSdT, dgfundr)
            deallocate (dSdR, blittleR, blittleT, Q, blittleP)
            deallocate (dbPdt, dbPdR, d2gfundr2)
            deallocate (d2SdR2, d2SdT2, dbTdR, dbRdT, jlittleR)
            deallocate (jlittleT, jlittleP, blittleX, blittleY1)
            deallocate (blittleZ1, jlittleX, jlittleY1, jlittleZ1)
            deallocate (ctreal, streal, cpreal, spreal)
            deallocate (jlittley, jlittlez, blittley, blittlez)       
        end subroutine calcFields

        !> determines values of density, pressure, and magnetic pressure
        !> along the bubble boundary -- thus, for a given value of r (physical, pre-ss radial coord)
        !> quantities will be calculated using the OUTSIDE solution that need to be added to the INSIDE
        !> gas pressure and density so that total pressure is constant across the bubble boundary
        !> NOTE, inside solution Stream function (and so Pgas and Pmag) is zero at surface
        !> also, outside solution is potential so density and pressure will be
        !> in radial  HEbalance
        subroutine calcInside(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            
            real(rp) :: r1
            real(rp), dimension(:, :), allocatable :: st, mu
            real(rp), dimension(:, :), allocatable :: al1, al2, psi0, f, g, h, psi1
            real(rp), dimension(:, :), allocatable :: i, psi3, dpsi0dr, dpsi0dth, dpsi0dr2
            real(rp), dimension(:, :), allocatable :: dpsi0drdth, dpsi0dth2, dfdr, dfdmu
            real(rp), dimension(:, :), allocatable :: dfdrdmu, dfdr2, dfdmu2, dgdr, dgdmu
            real(rp), dimension(:, :), allocatable :: dgdrdmu, dgdr2, dgdmu2, dhdr, dhdmu
            real(rp), dimension(:, :), allocatable :: dhdrdmu, dhdr2, dhdmu2, dpsi1dr
            real(rp), dimension(:, :), allocatable :: dpsi1dr2, dpsi1dmu, dpsi1drdmu, dpsi1dmu2
            real(rp), dimension(:, :), allocatable :: dpsi1dth, dpsi1drdth, dpsi1dth2, didr
            real(rp), dimension(:, :), allocatable :: didmu, didrdmu, didr2, didmu2
            real(rp), dimension(:, :), allocatable :: dpsi3dr, dpsi3dmu, dpsi3drdmu, dpsi3dr2
            real(rp), dimension(:, :), allocatable :: dpsi3dmu2, dpsi3dth, dpsi3drdth, dpsi3dth2
            real(rp), dimension(:, :), allocatable :: dAdr, dAdth, dAdrdth, dAdr2, dAdth2
            real(rp), dimension(:, :), allocatable :: dbrdr, dbthdr, dbrdth, dbthdth
            real(rp), dimension(:, :), allocatable :: rnolam, Bpresin
            real(rp), dimension(:, :), allocatable :: dldr, lr, ufunc, dufdr
            real(rp), dimension(:, :), allocatable :: dthdr, dstuffdr, dstuffdth, dpdr
            real(rp) :: dlamdr

            allocate (al1(State%js:State%je, State%ks:State%ke), st(State%js:State%je, State%ks:State%ke), mu(State%js:State%je, State%ks:State%ke))
            allocate (al2(State%js:State%je, State%ks:State%ke), psi0(State%js:State%je, State%ks:State%ke), f(State%js:State%je, State%ks:State%ke), g(State%js:State%je, State%ks:State%ke), h(State%js:State%je, State%ks:State%ke), psi1(State%js:State%je, State%ks:State%ke))
            allocate (i(State%js:State%je, State%ks:State%ke), psi3(State%js:State%je, State%ks:State%ke), dpsi0dr(State%js:State%je, State%ks:State%ke), dpsi0dth(State%js:State%je, State%ks:State%ke), dpsi0dr2(State%js:State%je, State%ks:State%ke))
            allocate (dpsi0drdth(State%js:State%je, State%ks:State%ke), dpsi0dth2(State%js:State%je, State%ks:State%ke), dfdr(State%js:State%je, State%ks:State%ke), dfdmu(State%js:State%je, State%ks:State%ke))
            allocate (dfdrdmu(State%js:State%je, State%ks:State%ke), dfdr2(State%js:State%je, State%ks:State%ke), dfdmu2(State%js:State%je, State%ks:State%ke), dgdr(State%js:State%je, State%ks:State%ke), dgdmu(State%js:State%je, State%ks:State%ke))
            allocate (dgdrdmu(State%js:State%je, State%ks:State%ke), dgdr2(State%js:State%je, State%ks:State%ke), dgdmu2(State%js:State%je, State%ks:State%ke), dhdr(State%js:State%je, State%ks:State%ke), dhdmu(State%js:State%je, State%ks:State%ke))
            allocate (dhdrdmu(State%js:State%je, State%ks:State%ke), dhdr2(State%js:State%je, State%ks:State%ke), dhdmu2(State%js:State%je, State%ks:State%ke), dpsi1dr(State%js:State%je, State%ks:State%ke))
            allocate (dpsi1dr2(State%js:State%je, State%ks:State%ke), dpsi1dmu(State%js:State%je, State%ks:State%ke), dpsi1drdmu(State%js:State%je, State%ks:State%ke), dpsi1dmu2(State%js:State%je, State%ks:State%ke))
            allocate (dpsi1dth(State%js:State%je, State%ks:State%ke), dpsi1drdth(State%js:State%je, State%ks:State%ke), dpsi1dth2(State%js:State%je, State%ks:State%ke), didr(State%js:State%je, State%ks:State%ke))
            allocate (didmu(State%js:State%je, State%ks:State%ke), didrdmu(State%js:State%je, State%ks:State%ke), didr2(State%js:State%je, State%ks:State%ke), didmu2(State%js:State%je, State%ks:State%ke))
            allocate (dpsi3dr(State%js:State%je, State%ks:State%ke), dpsi3dmu(State%js:State%je, State%ks:State%ke), dpsi3drdmu(State%js:State%je, State%ks:State%ke), dpsi3dr2(State%js:State%je, State%ks:State%ke))
            allocate (dpsi3dmu2(State%js:State%je, State%ks:State%ke), dpsi3dth(State%js:State%je, State%ks:State%ke), dpsi3drdth(State%js:State%je, State%ks:State%ke), dpsi3dth2(State%js:State%je, State%ks:State%ke))
            allocate (dAdr(State%js:State%je, State%ks:State%ke), dAdth(State%js:State%je, State%ks:State%ke), dAdrdth(State%js:State%je, State%ks:State%ke), dAdr2(State%js:State%je, State%ks:State%ke), dAdth2(State%js:State%je, State%ks:State%ke))
            allocate (dbrdr(State%js:State%je, State%ks:State%ke), dbthdr(State%js:State%je, State%ks:State%ke), dbrdth(State%js:State%je, State%ks:State%ke), dbthdth(State%js:State%je, State%ks:State%ke))
            allocate (rnolam(State%js:State%je, State%ks:State%ke), Bpresin(State%js:State%je, State%ks:State%ke))
            allocate (dldr(State%js:State%je, State%ks:State%ke), lr(State%js:State%je, State%ks:State%ke), ufunc(State%js:State%je, State%ks:State%ke), dufdr(State%js:State%je, State%ks:State%ke))
            allocate (dthdr(State%js:State%je, State%ks:State%ke), dstuffdr(State%js:State%je, State%ks:State%ke), dstuffdth(State%js:State%je, State%ks:State%ke), dpdr(State%js:State%je, State%ks:State%ke))

            al1 = 0.0; st = 0.0; mu = 0.0
            al2 = 0.0; psi0 = 0.0; f = 0.0; g = 0.0; h = 0.0; psi1 = 0.0
            i = 0.0; psi3 = 0.0; dpsi0dr = 0.0; dpsi0dth = 0.0; dpsi0dr2 = 0.0
            dpsi0drdth = 0.0; dpsi0dth2 = 0.0; dfdr = 0.0; dfdmu = 0.0
            dfdrdmu = 0.0; dfdr2 = 0.0; dfdmu2 = 0.0; dgdr = 0.0; dgdmu = 0.0
            dgdrdmu = 0.0; dgdr2 = 0.0; dgdmu2 = 0.0; dhdr = 0.0; dhdmu = 0.0
            dhdrdmu = 0.0; dhdr2 = 0.0; dhdmu2 = 0.0; dpsi1dr = 0.0
            dpsi1dr2 = 0.0; dpsi1dmu = 0.0; dpsi1drdmu = 0.0; dpsi1dmu2 = 0.0
            dpsi1dth = 0.0; dpsi1drdth = 0.0; dpsi1dth2 = 0.0; didr = 0.0
            didmu = 0.0; didrdmu = 0.0; didr2 = 0.0; didmu2 = 0.0
            dpsi3dr = 0.0; dpsi3dmu = 0.0; dpsi3drdmu = 0.0; dpsi3dr2 = 0.0
            dpsi3dmu2 = 0.0; dpsi3dth = 0.0; dpsi3drdth = 0.0; dpsi3dth2 = 0.0
            dAdr = 0.0; dAdth = 0.0; dAdrdth = 0.0; dAdr2 = 0.0; dAdth2 = 0.0
            dbrdr = 0.0; dbthdr = 0.0; dbrdth = 0.0; dbthdth = 0.0
            rnolam = 0.0; Bpresin = 0.0
            dldr = 0.0; lr = 0.0; ufunc = 0.0; dufdr = 0.0
            dthdr = 0.0; dstuffdr = 0.0; dstuffdth = 0.0; dpdr = 0.0

            r1 = sqrt(Model%xo*Model%xo)

            ! Stretched r
            State%rout = State%rlam

            ! ; calculate theta along the boundary, and rotate so that the
            ! ; bubble is at the equator
            State%rat = ((r1*r1 + State%rout*State%rout - Model%rbub*Model%rbub)/(2*State%rout*r1))
            where (State%rat .gt. 1.) State%rat = 1
            where (State%rat .lt. -1.) State%rat = -1
            State%thout = acos(State%rat)

            st = sin(State%thout)
            mu = cos(State%thout)

            call zero2tiny2d(st) ! where (abs(ct) .lt. tiny) ct=tiny*sign(ct)
            call zero2tiny2d(mu) ! where (abs(st) .lt. tiny) st=tiny*sign(st)

            al1 = r1*r1 - Model%rbub*Model%rbub
            al2 = 2*r1*r1 - Model%rbub*Model%rbub
           
            ! ;
            ! ;  eq. B2
            psi0 = mu

            ! ;
            ! ;  eq. B3
            f = r1*(State%rout*State%rout + al1) - State%rout*mu*al2
            g = al1*al1 + r1*r1*State%rout*State%rout - 2.*State%rout*r1*al1*mu
            h = 1./(sqrt(g))
            psi1 = (-1./Model%rbub)*f*h

            ! ;
            ! ;  eq. B4
            i = State%rout*State%rout + r1*r1 - 2.*State%rout*r1*mu
            psi3 = (1./Model%rbub)*sqrt(i)

            ! ;
            ! ;  eq. B5
            State%streamout = psi0 + psi1 + psi3

            ! ;
            ! ;  now take derivatives for magnetic fields, etc.
            ! ;
            dpsi0dr = 0.
            dpsi0dth = -st
            dpsi0dr2 = 0.
            dpsi0drdth = 0.
            dpsi0dth2 = -mu

            dfdr = 2.*r1*State%rout - mu*al2
            dfdmu = -State%rout*al2
            dfdrdmu = -al2
            dfdr2 = 2.*r1
            dfdmu2 = 0.

            dgdr = 2.*State%rout*r1*r1 - 2.*r1*al1*mu
            dgdmu = -2.*State%rout*r1*al1
            dgdrdmu = -2.*r1*al1
            dgdr2 = 2.*r1*r1
            dgdmu2 = 0.

            dhdr = -dgdr/2./(g**1.5)
            dhdmu = -dgdmu/2./(g**1.5)
            dhdrdmu = -dgdrdmu/2./(g**1.5) + 3.*dgdmu*dgdr/4./(g**2.5)
            dhdr2 = -dgdr2/2./(g**1.5) + 3.*dgdr*dgdr/4./(g**2.5)
            dhdmu2 = 1.5*(dgdmu**2)/2./(g**2.5)

            dpsi1dr = (-1./Model%rbub)*(dfdr*h + f*dhdr)
            dpsi1dr2 = (-1./Model%rbub)*(2.*dfdr*dhdr + dfdr2*h + f*dhdr2)

            dpsi1dmu = (-1./Model%rbub)*(dfdmu*h + f*dhdmu)
            dpsi1drdmu = (-1./Model%rbub)*(dfdmu*dhdr + dfdr*dhdmu + dfdrdmu*h + f*dhdrdmu)
            dpsi1dmu2 = (-1./Model%rbub)*(dfdmu2*h + 2.*dfdmu*dhdmu + f*dhdmu2)

            dpsi1dth = -st*dpsi1dmu
            dpsi1drdth = -st*dpsi1drdmu
            dpsi1dth2 = -mu*dpsi1dmu + st**2*dpsi1dmu2

            didr = 2.*State%rout - 2.*r1*mu
            didmu = -2.*State%rout*r1
            didrdmu = -2.*r1
            didr2 = 2.
            didmu2 = 0.

            dpsi3dr = (1./Model%rbub)*didr/2./sqrt(i)
            dpsi3dmu = (1./Model%rbub)*didmu/2./sqrt(i)
            dpsi3drdmu = (1./Model%rbub)*(didrdmu/2./sqrt(i) - didmu*didr/4./(i**1.5))
            dpsi3dr2 = (1./Model%rbub)*(didr2/2./sqrt(i) - (didr)**2./4./(i**1.5))
            dpsi3dmu2 = (1./Model%rbub)*didmu2/2./sqrt(i) - (1./Model%rbub)*(didmu**2.)/4./(i**1.5)

            dpsi3dth = -st*dpsi3dmu
            dpsi3drdth = -st*dpsi3drdmu
            dpsi3dth2 = -mu*dpsi3dmu + st**2.*dpsi3dmu2

            dAdr = dpsi0dr + dpsi1dr + dpsi3dr
            dAdth = dpsi0dth + dpsi1dth + dpsi3dth
            dAdrdth = dpsi0drdth + dpsi1drdth + dpsi3drdth
            dAdr2 = dpsi0dr2 + dpsi1dr2 + dpsi3dr2
            dAdth2 = dpsi0dth2 + dpsi1dth2 + dpsi3dth2
            ! ;
            ! ;  now calculate magnetic field
            ! ;
            ! ;  eq. B1
            State%brlambout = dAdth/(State%rout**2*st)
            State%bthlambout = -dAdr/(State%rout*st)
            State%bphlambout = 0.
            
            ! ;
            ! ;  now calculate field derivatives
            ! ;
            dbrdr = (-2.*dAdth/State%rout + dAdrdth)/(State%rout**2.*st)
            dbthdr = (dAdr/State%rout - dAdr2)/State%rout/st

            dbrdth = -mu*dAdth/(State%rout**2.*st**2)
            dbrdth = dbrdth + dAdth2/(State%rout**2.*st)
            dbthdth = mu*dAdr/(State%rout*st**2)
            dbthdth = dbthdth - dAdrdth/(State%rout*st)

            ! Back to rsquig
            rnolam = State%rout - Model%apar
            State%presin = (State%rout**2./rnolam**2.)*(1. - (State%rout**2./rnolam**2.))*(State%brlambout**2.)/8./pi
            State%presin = State%presin/Model%phiss**4.
            ! ;  here are the magnetic pressure in the rpB,thetapB
            ! ;  phipB coordinate systems !  we have to also get rid of selfsim stuff
            ! ;
            ! ;  dlamdr is one, because rlam = rsquig + apar.
            dlamdr = 1.
            State%Bpresin = (( State%brlambout*(State%rout**2./rnolam**2.)/Model%phiss**2.)**2. &
                            + (State%bthlambout*(State%rout/rnolam)*dlamdr/Model%phiss**2.)**2.)/8./pi
            State%Pbackin = (State%presin + State%Bpresin)*Model%outScale**2.

            ! ;
            ! ;  now calculate the background density
            ! ;
            dldr = -Model%apar/(rnolam*rnolam)
            lr = (Model%apar + rnolam)/rnolam

            ufunc = (r1*r1 + State%rout*State%rout - Model%rbub*Model%rbub)/(2*State%rout*r1)
            dufdr = (1./r1) - ufunc/State%rout

            dthdr = -(1./sqrt(1. - ufunc*ufunc))*dufdr
            dstuffdr = -(lr*dldr*State%brlambout*State%brlambout + lr*lr*State%brlambout*dbrdr + &
                            .5*dldr*State%bthlambout*State%bthlambout + lr*State%bthlambout*dbthdr)/(4.*pi)
            dstuffdth = -(lr*lr*State%brlambout*dbrdth + lr*State%bthlambout*dbthdth)/(4.*pi)
            dpdr = dstuffdr + dstuffdth*dthdr

            F = GMm/((State%rout - Model%apar)**2)/(Rsolar_cgs**2) + Model%alpha*(State%rout - Model%apar)*Rsolar_cgs*Mp_cgs
            State%Dbackin = dpdr/(F*Rsolar_cgs*Model%phiss**3.)*Model%outScale**2.

            deallocate (al1)
            deallocate (al2, psi0, f, g, h, psi1)
            deallocate (i, psi3, dpsi0dr, dpsi0dth, dpsi0dr2)
            deallocate (dpsi0drdth, dpsi0dth2, dfdr, dfdmu)
            deallocate (dfdrdmu, dfdr2, dfdmu2, dgdr, dgdmu)
            deallocate (dgdrdmu, dgdr2, dgdmu2, dhdr, dhdmu)
            deallocate (dhdrdmu, dhdr2, dhdmu2, dpsi1dr)
            deallocate (dpsi1dr2, dpsi1dmu, dpsi1drdmu, dpsi1dmu2)
            deallocate (dpsi1dth, dpsi1drdth, dpsi1dth2, didr)
            deallocate (didmu, didrdmu, didr2, didmu2)
            deallocate (dpsi3dr, dpsi3dmu, dpsi3drdmu, dpsi3dr2)
            deallocate (dpsi3dmu2, dpsi3dth, dpsi3drdth, dpsi3dth2)
            deallocate (dAdr, dAdth, dAdrdth, dAdr2, dAdth2)
            deallocate (dbrdr, dbthdr, dbrdth, dbthdth)
            deallocate (rnolam, Bpresin)
            deallocate (dldr, lr, ufunc, dufdr)
            deallocate (dthdr, dstuffdr, dstuffdth, dpdr)
        end subroutine calcInside

        !> Calculate final Gibson-Low Solution scaler values based on total dens, temp and pressure
        !>
        !>
        subroutine finishSolution(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            real(rp), dimension(:, :), allocatable :: dphidt

            allocate(dphidt(State%js:State%je, State%ks:State%ke))
            ! ; also velocity (in km/sec)
            ! ; note this is radial

            dphidt = sqrt((Model%eta*Model%phiss - 2.*Model%alpha)/Model%phiss)
            Solution%v(State%ri, :, :, XDIR)= (State%rpb*Rsolar*dphidt)/Model%phiss
    
            where (Solution%dens(State%ri, :, :) .gt. 0) 
                Solution%Temp(State%ri, :, :) = Solution%Pres(State%ri, :, :) / &
                                                2./kboltz/Solution%dens(State%ri, :, :)
            end where

            if (Model%bonly .gt. 0) then
                Solution%Temp(State%ri, :, :) = Solution%Temp(State%ri, :, :)*0. + Model%Isothermal
                Solution%dens(State%ri, :, :) = State%DensbackHEonly
                Solution%Pres(State%ri, :, :) = gas_R*Solution%Temp(State%ri, :, :)*Mp_cgs*Solution%dens(State%ri, :, :)
                if (State%inside_count .gt. 0) then
                    where (Solution%inside_mask(State%ri, :, :) .gt. 0.1) 
                        Solution%dens(State%ri, :, :) = Model%bonly*Solution%dens(State%ri, :, :)
                    end where
                end if
                State%cavinside = .false.
                where (State%rcap .lt. (Model%c2bonly*Model%rbub)) State%cavinside = .true.
                if (any(State%cavinside)) where (State%cavinside) Solution%dens(State%ri, :, :) = Model%c1bonly* Solution%dens(State%ri, :, :)
            end if
    
            ! Zero-out values outside of bubble for fields and velocity
            ! set scalars to default value
            if (Model%bubbleonly .eq. 1) then
                where (Solution%inside_mask(State%ri, :, :) .lt. 1.0)
                    Solution%dens(State%ri, :, :) = 1.d-5
                    Solution%temp(State%ri, :, :)= 1.d-5
                    Solution%pres(State%ri, :, :) = 1.d-5
                    Solution%b(State%ri, :, :, 1) = 0.
                    Solution%b(State%ri, :, :, 2) = 0.
                    Solution%b(State%ri, :, :, 3) = 0.
                    Solution%j(State%ri, :, :, 1) = 0.
                    Solution%j(State%ri, :, :, 2) = 0.
                    Solution%j(State%ri, :, :, 3) = 0.
                    Solution%v(State%ri, :, :, 1) = 0.
                    Solution%v(State%ri, :, :, 2) = 0.
                    Solution%v(State%ri, :, :, 3) = 0.
                end where
            end if

            where (Solution%dens(State%ri, :, :)  /= Solution%dens(State%ri, :, :) )
                Solution%dens(State%ri, :, :) = 0.
            end where
            where (Solution%temp(State%ri, :, :)  /= Solution%temp(State%ri, :, :) )
                Solution%temp(State%ri, :, :)= 0.
            end where
            where (Solution%pres(State%ri, :, :)  /= Solution%pres(State%ri, :, :) )
                Solution%pres(State%ri, :, :) = 0.
            end where
        end subroutine finishSolution

        !> Precondition model magnetic field based on CME geometry 
        !> this normalizes the bubble field by the user input Bmax and 
        !> the geometry max bmag value (field monotonically increases, 
        !> peaks at t0, depending on geometry parameters then exp decay
        !> over time.)
        !> 
        function calcModelBmax(Model, State, Solution, nt) result(ao)
            type(glModel_T), intent(inout) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            integer, intent(in) :: nt
            integer :: i
            real(rp) :: ao, bmagmax, zeta, time, val
            real(rp), allocatable, dimension(:,:) :: bmag, st, ct, sph, cph
            
            bmagmax = 0.0
            if(Model%isLoud) write(*,"(1X,A24,2X,1F)") "velmult scaled:", Model%velmult
            State%r = Model%frontheight
            ! do while bmag is monotonically increasing, save max, exit if bmagmax decreases from previous
            if(Model%isLoud) write(*,"(1X,A24)") "Find bmagmax:"
            if(Model%isLoud) write(*,"(1X,A24)") "------------------------"
            do i=1, nt
                zeta = Model%x0 - Model%r0 - Model%apar + (nt - i - 1.)*2.*Model%r0/nt
                if(Model%isLoud) write(*,"(1X,A14,2X,1F,1I)") "zeta, nt:", zeta, i
                if (zeta > 0.01) then
                    Model%time = (Model%frontheight/zeta - 1.)/Model%s_eta
                    call generateGLSolution(Solution, Model, State)
                    st = sin(State%thpb)                    
                    ct = cos(State%thpb)
                    sph = sin(State%phpb)
                    cph = cos(State%phpb)
                    bmag = sqrt((Solution%b(1,:,:,XDIR)*st*cph + Solution%b(1,:,:,YDIR)*ct*cph - Solution%b(1,:,:,ZDIR)*sph)**2. + &
                                (Solution%b(1,:,:,XDIR)*st*sph + Solution%b(1,:,:,YDIR)*ct*sph + Solution%b(1,:,:,ZDIR)*cph)**2. + &
                                (Solution%b(1,:,:,XDIR)*ct - Solution%b(1,:,:,YDIR)*st)**2.)
                    val = maxval(bmag)
                    if (bmagmax < val) then 
                        bmagmax = val
                    else 
                        if(Model%isLoud) write(*,"(1X,A36,2X,2F)") "bmagmax (1) < val (2), exit loop:", bmagmax, val
                        exit
                    end if
                end if
            end do
            if(Model%isLoud) write(*,"(1X,A24,2X,1F,1I,1F)") "bmagmax, nt, time:", bmagmax, i, Model%time
            if(Model%isLoud) write(*,"(1X,A24)") "------------------------"
            ao = Model%bmax / bmagmax 
        end function calcModelBmax

        !> For a previously initialized take Model
        !> and calculate the Gibson-Low solution and populate
        !> Solution to be used in other boundary condition code
        !> or ingested into LOS code
        subroutine generateGLSolution(Solution, Model, State)
            class(glSolution_T), intent(inout)  :: Solution
            class(glModel_T), intent(inout) :: Model
            class(glState_T), intent(inout)  :: State
            integer :: i, rdim

            rdim = size(State%r)
            if(Model%isLoud .and. (Model%isDebug)) write(*,*) "Generating GL Solution"

            ! time assumed in seconds
            Model%phiss = sqrt(Model%eta)*Model%time + 1.0 
            if (Model%isLoud .and. (Model%isDebug)) then
                write(*,"(1X,A14,2X,F)") "eta: ", Model%eta 
                write(*,"(1X,A14,2X,F)") "time: ", Model%time 
                write(*,"(1X,A14,2X,F)") "phiss: ", Model%phiss 
            end if
            
            do i = 1, rdim
                State%rpb = State%r(i)
                State%ri = i

                ! Calculate cap coordinate transformation
                call calcCoords(Model, State)
                ! Determine the bubble mask
                where (State%rcap .lt. Model%rbub) Solution%inside_mask(State%ri,:,:) = 1.              
                State%inside_count = sum(Solution%inside_mask(State%ri,:,:))
                State%outside_count = sum(1. - Solution%inside_mask(State%ri,:,:))   
                ! Calculate Outside Field solution
                ! Outside Field Pressure and Dens contributions
                ! and Inside fields
                call calcFields(Model, State, Solution)
                ! Calculate total and background pressures
                ! inside bubble
                call calcPressure(Model, State, Solution)
                ! Calculate background and total density 
                ! inside bubble
                call calcDensity(Model, State, Solution)
                call finishSolution(Model, State, Solution)
            end do

            if(Model%isLoud .and. Model%isDebug) write(*,*) "GL Solution Complete"
            if(Model%isLoud .and. Model%isDebug) write(*,*) "inside_count: ", sum(Solution%inside_mask)
        end subroutine generateGLSolution 
end module