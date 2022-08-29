module solution
    use gltypes
    use glutils
    implicit none

    contains
        !---------------
        ! Solution
        !---------------
    
        !> Calculate Gibson-Low Solution Derivatives 
        !>
        !>
        subroutine calcDerivs(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution

            real(rp), dimension(:,:), allocatable :: sph
            real(rp), dimension(:,:), allocatable :: cph
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

            integer :: dim1, dim2, cs, ss

            dim1 = State%Nj
            dim2 = State%Nk
            allocate (sph(dim1, dim2))
            allocate (cph(dim1, dim2))
            allocate (sTbig(dim1, dim2))
            allocate (cTbig(dim1, dim2))
            allocate (sPbig(dim1, dim2))
            allocate (cPbig(dim1, dim2))
            allocate (gfunot(dim1, dim2))
            allocate (gfun(dim1, dim2))
            allocate (dSdT(dim1, dim2))
            allocate (dgfundr(dim1, dim2))
            allocate (dSdR(dim1, dim2))
            allocate (dbPdt(dim1, dim2))
            allocate (dbPdR(dim1, dim2))
            allocate (d2gfundr2(dim1, dim2))
            allocate (d2SdT2(dim1, dim2))
            allocate (dbTdR(dim1, dim2))
            allocate (dbRdT(dim1, dim2))
            allocate (dTdmu(dim1, dim2))
            allocate (dTdR(dim1, dim2))
            allocate (dRdmu(dim1, dim2))
            allocate (dPidmu(dim1, dim2))
            allocate (dthetadmu(dim1, dim2))
            allocate (dPidR(dim1, dim2))
            allocate (dRdlam(dim1, dim2))
            allocate (dmudlam(dim1, dim2))
            allocate (d2SdTdR(dim1, dim2))
            

            ! ;  See eq A19
            ! ; 
            ! ;  first need to calculate tderiv = dT/drlam, where T = Pi + Bstrength (squared)
            ! ;  in the Rcap,Thcap system
            ! ;
            State%st = sin(State%thpB)
            State%ct = cos(State%thpB)
            sph = sin(State%phpB)
            cph = cos(State%phpB)

            sTbig = sin(State%thcap)
            cTbig = cos(State%thcap)
            sPbig = sin(State%phcap)
            cPbig = cos(State%phcap)

            call zero2tiny2d(State%st) ! where (abs(State%st) .lt. tiny) State%st=tiny*sign(State%st)
            call zero2tiny2d(State%ct) ! where (abs(State%ct) .lt. tiny) State%ct=tiny*sign(State%ct)
            call zero2tiny2d(sTbig) ! where (abs(sTbig) .lt. tiny) sTbig=tiny*sign(sTbig)
            call zero2tiny2d(cTbig) ! where (abs(cTbig) .lt. tiny) cTbig=tiny*sign(cTbig)

            State%mu = State%rcap*sTbig
            gfunot = sin(Model%alnot*Model%rbub)/Model%alnot/Model%rbub - cos(Model%alnot*Model%rbub)
            dSdT = State%stream*2*cTbig/sTbig
            dgfundr = cos(Model%alnot*State%rcap)/State%rcap - sin(Model%alnot*State%rcap)/State%rcap/State%rcap/Model%alnot + &
                    Model%alnot*sin(Model%alnot*State%rcap)
            dSdR = (4.*Model%ao*pi/Model%alnot/Model%alnot)*(dgfundr*Model%rbub**2/gfunot - 2*State%rcap)*sTbig**2
            d2SdTdR = dSdR*2.*cTbig/sTbig
            d2gfundr2 = (-State%rcap*State%rcap*Model%alnot*sin(Model%alnot*State%rcap) - State%rcap*cos(Model%alnot*State%rcap) + &
                        Model%alnot*State%rcap**3.*Model%alnot*cos(Model%alnot*State%rcap) - &
                        cos(Model%alnot*State%rcap)*State%rcap + 2*sin(Model%alnot*State%rcap)/Model%alnot)/State%rcap**3.
            d2SdR2 = (4.*Model%ao*pi/Model%alnot/Model%alnot)*(d2gfundr2*Model%rbub**2/gfunot - 2)*sTbig**2
            d2SdT2 = State%stream*2*(cTbig**2 - sTbig**2)/(sTbig**2)
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

            dTdR = dTdR/8./pi
            dTdR = dTdR + dPidR

            dTdmu = -2*((dSdT/State%rcap)**2)/(State%mu**3)
            dTdmu = dTdmu + 2*(dSdT/State%rcap)*(dRdmu*(-dSdT/(State%rcap**2) + &
                                                d2SdTdR/State%rcap) + dthetadmu*d2SdT2/State%rcap)/(State%mu**2)
            dTdmu = dTdmu - 2*(dSdR**2)/(State%mu**3)
            dTdmu = dTdmu + 2*dSdR*(dRdmu*d2SdR2 + dthetadmu*d2SdTdR)/(State%mu**2)
            dTdmu = dTdmu - 2*Model%alnot*Model%alnot*(State%stream**2)/(State%mu**3)
            dTdmu = dTdmu + 2*Model%alnot*Model%alnot*State%stream*(dRdmu*dSdR + dthetadmu*dSdT)/(State%mu**2)
            dTdmu = dTdmu/8./pi
            dTdmu = dTdmu + dPidmu

            dRdlam = (State%rlam - Model%xo*State%st*cph)/State%rcap

            cs = cos(Model%sigma)
            ss = sin(Model%sigma)
            dmudlam = (State%rlam*(State%st**2*(cph**2 + sph**2*cs**2) + &
                        State%ct**2*ss**2 - 2*State%st*State%ct*sph*ss*cs) - &
                        Model%xo*State%st*cph)/State%mu
            State%tderivR = dTdR*dRdlam
            State%tderivmu = dTdmu*dmudlam
            State%tderiv = State%tderivR + State%tderivmu
            if (State%outside_count .gt. 0) where (Solution%inside_mask(State%ri, :, :) .lt. 1.) State%tderiv = Model%outScale**2.*State%tderivout

        end subroutine calcDerivs

        !> Calculate Gibson-Low Solution Density
        !> 
        subroutine calcDensity(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            real(rp), dimension(:, :), allocatable :: d1, d2, d3, d4
            
            ! ;  the force F could include the effects of the self-similar
            ! ;  acceleration (alpha ne 0) - otherwise its just gravity.
            ! ; Note that to properly include alpha ne 0, however, phiss needs
            ! ; to be defined as a differential equation
            ! ;
            State%F = GMm/(State%rsquig**2)/(Rsun**2) + Model%alpha*State%rsquig*Rsun*mprot

            d1 = -((State%rlam/State%rsquig)**2)*(1 - ((State%rlam/State%rsquig)**2))*State%tderiv
            d2 = 2*(State%rlam/State%rsquig)*(Model%apar/(State%rsquig**2))*State%glpi
            d3 = (State%rlam/State%rsquig)*(Model%apar/(State%rsquig**2))* &
                (1 - 2*((State%rlam/State%rsquig)**2))*State%blittlerlamb**2/4./pi
            d4 = ((State%rlam/State%rsquig)**2)*((Model%apar**2)/(State%rsquig**2) + &
                                    2*(Model%apar/State%rsquig))*(State%blittlethlamb**2 + State%blittlephlamb**2)/4./pi/State%rlam
            Solution%dens(State%ri,:,:)  = d1 + d2 + d3 + d4

            Solution%dens(State%ri,:,:)  = Solution%dens(State%ri,:,:) / State%F / Rsun
            

            if ((State%rcap(1,1) .gt. Model%rbub) .and. (Model%isLoud))  then
                write(*,"(1X,A20,2X,E13.6)") "Max d1: ", maxval(d1)
                write(*,"(1X,A20,2X,E13.6)") "Min d1: ",  minval(d1)
                write(*,"(1X,A20,2X,E13.6)") "Max d2: ",  maxval(d2)
                write(*,"(1X,A20,2X,E13.6)") "Min d2: ", minval(d2)
                write(*,"(1X,A20,2X,E13.6)") "Max d3: ",  maxval(d3)
                write(*,"(1X,A20,2X,E13.6)") "Min d3: ", minval(d3)
                write(*,"(1X,A20,2X,E13.6)") "Max d4: ",  maxval(d4)
                write(*,"(1X,A20,2X,E13.6)") "Min d4: ", minval(d4)
                write(*,"(1X,A20,2X,E13.6)") "Max densnoback: ",  maxval(Solution%dens(State%ri,:,:))
                write(*,"(1X,A20,2X,E13.6)") "Min densnoback: ", minval(Solution%dens(State%ri,:,:))
                write(*,"(1X,A20,2X,E13.6)") "Max F: ",  maxval(State%F)
                write(*,"(1X,A20,2X,E13.6)") "Min F: ", minval(State%F)
            end if
            ! ;  put in background density
            ! ;  first radial power law HE
            ! ; see Gibson et al 1998, eq 3 
            ! ;
            State%densin = (Model%aa + Model%cc + Model%ee)*State%rsquig**(-3.)
            State%densout = (Model%aa*State%rsquig**(-Model%bb) + Model%cc*State%rsquig**(-Model%dd) + Model%ee*State%rsquig**(-Model%ff))
            State%densback = 0.0

            if (State%inside_count .gt. 0) where (Solution%inside_mask(State%ri, :, :) .gt. 0.) State%densback = State%densin
            if (State%outside_count .gt. 0) where (Solution%inside_mask(State%ri, :, :) .lt. 1.) State%densback = State%densout
            State%DensbackHEonly = State%densback
            ! ; 
            ! ; now add total pressure continuity part
            ! ;
            if (State%inside_count .gt. 0) where (Solution%inside_mask(State%ri, :, :) .gt. 0.) State%densback = State%densback + State%dbackin*(Model%phiss**3)
            Solution%dens(State%ri,:,:) = Solution%dens(State%ri,:,:) + State%Densback
            ! ;
            ! ;  put dens in ss coords
            ! ;
            Solution%dens(State%ri,:,:) = Solution%dens(State%ri,:,:)/(Model%phiss**3)
            State%Densback = State%Densback/(Model%phiss**3)
            if ((State%rcap(1,1) .gt. Model%rbub) .and. (Model%isLoud)) then
                write(*,*) "r: ", State%rcap(1,1)
                write(*,"(1X,A20,2X,E13.6)") "max densin: ", maxval(State%densin)
                write(*,"(1X,A20,2X,E13.6)") "min densin: ", minval(State%densin)
                write(*,"(1X,A20,2X,E13.6)") "max densout: ", maxval(State%densout)
                write(*,"(1X,A20,2X,E13.6)") "min densout: ", minval(State%densout)
                write(*,"(1X,A20,2X,E13.6)") "max densback: ", maxval(State%densback)
                write(*,"(1X,A20,2X,E13.6)") "min densback: ", minval(State%densback)
                write(*,"(1X,A20,2X,E13.6)") "max dens: ", maxval(Solution%dens(State%ri,:,:))
                write(*,"(1X,A20,2X,E13.6)") "min dens: ", minval(Solution%dens(State%ri,:,:))
            end if
        end subroutine calcDensity

        !> Calculate Gibson-Low Solution Pressure
        !>
        !>
        subroutine calcPresssure(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            ! ;
            ! ; calculate pressure from Gibson and Low eq. A18
            ! ;
            Solution%Pres(State%ri, :, :) = (State%rlam**2/State%rsquig**2) & 
                                            * (1 - (State%rlam**2/State%rsquig**2)) &
                                            * (State%blittlerlamb**2)/(8*pi) &
                                            + (State%rlam**2/State%rsquig**2)*State%glpi
            ! ; now we have to match total pressure Pmag+Pgas at the bubble interface
            ! ;  note that the INNER bubble solution is defined so that Pmag and Pgas are zero there
            ! ;  so, we calculate for each rlam, at the theta of the bubble boundary, the OUTER
            ! ;  solution total pressure, and add it to the INNER gas pressure 
            ! ;  (for all theta within the bubble at that rlam).  Then we must also
            ! ;  add a similar density increment to the INNER density, because the OUTER solution
            ! ; is potential and so has density and pressure in HE balance 
            ! ;
            State%presback = Solution%Pres(State%ri, :, :) * 0.

            if (State%inside_count .gt. 0) then
                if (Model%isDebug) write(*,*) 'calling calcInside'
                call calcInside(Model, State, Solution)
                where (Solution%inside_mask(State%ri, :, :) .eq. 1) State%presback = State%Pbackin*(Model%phiss**4)
            end if
            ! ;  now we put in additional hydrostatic background pressure to keep
            ! ;  things positive
            State%presout = (Model%aa/(Model%bb + 1.))*State%rsquig**(-Model%bb - 1.) &
                            + (Model%cc/(Model%dd + 1.))*State%rsquig**(-Model%dd - 1.) &
                            + (Model%ee/(Model%ff + 1.))*State%rsquig**(-Model%ff - 1.)

            State%presout = State%presout*GMm/Rsun
            ! ;
            ! ; see Gibson et al 1998, eqs 3 and 4 - note we are assuming alpha (He) = 0
            ! ;
            
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
            State%presin_0 = State%presin_0 + (Model%aa/(Model%bb + 1.)) &
                        + (Model%cc/(Model%dd + 1.)) &
                        + (Model%ee/(Model%ff + 1.))
            State%presin = State%presin + State%presin_0
            State%presin = State%presin*GMm/Rsun
            ! ; thus presin and presout will match smoothly at rsquig = 1
            ! ;  
            where (State%rsquig .lt. 1) State%presback = State%presback + State%presin
            where (State%rsquig .ge. 1) State%presback = State%presback + State%presout

            Solution%Pres(State%ri, :, :) = Solution%Pres(State%ri, :, :) + State%presback
            
            !; self-similar transform
            Solution%Pres(State%ri, :, :) = Solution%Pres(State%ri, :, :)/Model%phiss**4

        end subroutine calcPresssure

        !> Calculate Outside (of bubble) Gibson-Low Solution
        !>
        !>
        subroutine calcOutside(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution

            real(rp) :: r1
            integer :: dim1, dim2
            real(rp), dimension(:, :), allocatable :: al1, al2, psi0, f, g, h
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
            real(rp), dimension(:, :), allocatable :: sph, cph, blittleX1, blittleY1
            real(rp), dimension(:, :), allocatable :: blittleZ1, jlittleX1, jlittleY1, jlittleZ1
            real(rp), dimension(:, :), allocatable :: blittleY2, blittleZ, jlittleY2, jlittleZ
            real(rp), dimension(:, :), allocatable :: blittleX, blittleY, jlittleX, jlittleY
            real(rp), dimension(:, :), allocatable :: streal, ctreal, spreal, cpreal, dbr1dr
            real(rp), dimension(:, :), allocatable :: dbrdr
            real(rp), dimension(:, :), allocatable :: yt2

            dim1 = State%Nj
            dim2 = State%Nk
            allocate (al1(dim1, dim2), al2(dim1, dim2), psi0(dim1, dim2), f(dim1, dim2), g(dim1, dim2), h(dim1, dim2))
            allocate (psi1(dim1, dim2), i(dim1, dim2), psi3(dim1, dim2), dpsi0dr(dim1, dim2), dpsi0dth(dim1, dim2))
            allocate (d2psi0dth2(dim1, dim2), dpsi0dr2(dim1, dim2), dpsi0drdth(dim1, dim2), dfdr(dim1, dim2))
            allocate (dfdmu(dim1, dim2), dfdrdmu(dim1, dim2), dfdr2(dim1, dim2), dgdr(dim1, dim2), dgdmu(dim1, dim2))
            allocate (dgdrdmu(dim1, dim2), dgdr2(dim1, dim2), dhdr(dim1, dim2), dhdmu(dim1, dim2), d2hdmu2(dim1, dim2))
            allocate (dhdrdmu(dim1, dim2), dhdr2(dim1, dim2), dpsi1dr(dim1, dim2), dpsi1dr2(dim1, dim2))
            allocate (dpsi1dmu(dim1, dim2), d2psi1dmu2(dim1, dim2), dpsi1drdmu(dim1, dim2), dpsi1dth(dim1, dim2))
            allocate (dpsi1drdth(dim1, dim2), d2psi1dth2(dim1, dim2), didr(dim1, dim2), didmu(dim1, dim2))
            allocate (didrdmu(dim1, dim2), didr2(dim1, dim2), dpsi3dr(dim1, dim2), dpsi3dmu(dim1, dim2))
            allocate (d2psi3dmu2(dim1, dim2), dpsi3drdmu(dim1, dim2), dpsi3dr2(dim1, dim2), dpsi3dth(dim1, dim2))
            allocate (dpsi3drdth(dim1, dim2), d2psi3dth2(dim1, dim2), dAdr(dim1, dim2), dAdth(dim1, dim2))
            allocate (dAdrdth(dim1, dim2), dAdr2(dim1, dim2), dAdth2(dim1, dim2), brlambout1(dim1, dim2))
            allocate (bthlambout1(dim1, dim2), bphlambout1(dim1, dim2), dbth1dr(dim1, dim2), dbr1dth(dim1, dim2))
            allocate (jrlambout1(dim1, dim2), jthlambout1(dim1, dim2), jphlambout1(dim1, dim2))
            allocate (sph(dim1, dim2), cph(dim1, dim2), blittleX1(dim1, dim2), blittleY1(dim1, dim2))
            allocate (blittleZ1(dim1, dim2), jlittleX1(dim1, dim2), jlittleY1(dim1, dim2), jlittleZ1(dim1, dim2))
            allocate (blittleY2(dim1, dim2), blittleZ(dim1, dim2), jlittleY2(dim1, dim2), jlittleZ(dim1, dim2))
            allocate (blittleX(dim1, dim2), blittleY(dim1, dim2), jlittleX(dim1, dim2), jlittleY(dim1, dim2))
            allocate (streal(dim1, dim2), ctreal(dim1, dim2), spreal(dim1, dim2), cpreal(dim1, dim2), dbr1dr(dim1, dim2))
            allocate (dbrdr(dim1, dim2))
            allocate (yt2(dim1, dim2))

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
            State%rout = State%rlam

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

            r1 = sqrt(Model%xo*Model%xo)

            State%st = sin(State%thout)
            State%mu = cos(State%thout)

            call zero2tiny2d(State%st)
            call zero2tiny2d(State%mu)

            ! ;
            ! ;    break down bits of stream function
            ! ;
            al1 = State%rout*State%rout - Model%rbub*Model%rbub
            al2 = 2*State%rout*State%rout - Model%rbub*Model%rbub
            ! ;
            ! ; eq. B2
            psi0 = State%mu
            ! ;
            ! ; eq. B3
            f = r1*(State%rout*State%rout + al1) - State%rout*State%mu*al2
            g = al1*al1 + r1*r1*State%rout*State%rout - 2*State%rout*r1*al1*State%mu
            h = 1/(sqrt(g))
            psi1 = (-1/Model%rbub)*f*h
            ! ;
            ! ; eq. B4
            i = State%rout*State%rout + r1*r1 - 2*State%rout*r1*State%mu
            psi3 = (1/Model%rbub)*sqrt(i)
            ! ;
            ! ; eq. B5
            State%streamout = psi0 + psi1 + psi3
            ! ;
            ! ;  now take derivatives for magnetic fields, etc.
            ! ;
            dpsi0dr = 0.
            dpsi0dth = -State%st
            d2psi0dth2 = -State%mu
            dpsi0dr2 = 0.
            dpsi0drdth = 0.

            dfdr = 2*r1*State%rout - State%mu*al2
            dfdmu = -State%rout*al2
            dfdrdmu = -al2
            dfdr2 = 2*r1

            dgdr = 2*State%rout*r1*r1 - 2*r1*al1*State%mu
            dgdmu = -2*State%rout*r1*al1
            dgdrdmu = -2*r1*al1
            dgdr2 = 2*r1*r1

            dhdr = -dgdr/2./(g**1.5)
            dhdmu = -dgdmu/2./(g**1.5)
            d2hdmu2 = 3.*dgdmu/4./(g**2.5)
            dhdrdmu = -dgdrdmu/2./(g**1.5) + 3.*dgdmu*dgdr/4./(g**2.5)
            dhdr2 = -dgdr2/2/(g**1.5) + 3*dgdr*dgdr/4./(g**2.5)

            dpsi1dr = (-1/Model%rbub)*(dfdr*h + f*dhdr)
            dpsi1dr2 = (-1/Model%rbub)*(2.*dfdr*dhdr + dfdr2*h + f*dhdr2)

            dpsi1dmu = (-1/Model%rbub)*(dfdmu*h + f*dhdmu)
            d2psi1dmu2 = (-1/Model%rbub)*(2.*dfdmu*dhdmu + f*d2hdmu2)
            dpsi1drdmu = (-1/Model%rbub)*(dfdmu*dhdr + dfdr*dhdmu + dfdrdmu*h + f*dhdrdmu)

            dpsi1dth = -State%st*dpsi1dmu
            dpsi1drdth = -State%st*dpsi1drdmu
            d2psi1dth2 = -State%mu*dpsi1dmu + (State%st**2)*d2psi1dmu2

            didr = 2*State%rout - 2*r1*State%mu
            didmu = -2*State%rout*r1
            didrdmu = -2*r1
            didr2 = 2.

            dpsi3dr = (1/Model%rbub)*didr/2./sqrt(i)
            dpsi3dmu = (1/Model%rbub)*didmu/2./sqrt(i)
            d2psi3dmu2 = -(1/Model%rbub)*didmu/4./(i**(1.5))
            dpsi3drdmu = (1/Model%rbub)*(didrdmu/2./sqrt(i) - didmu*didr/4./(i**1.5))
            dpsi3dr2 = (1/Model%rbub)*(didr2/2./sqrt(i) - (didr)**2./4./(i**1.5))

            dpsi3dth = -State%st*dpsi3dmu
            dpsi3drdth = -State%st*dpsi3drdmu
            d2psi3dth2 = -State%mu*dpsi3dmu + (State%st**2)*d2psi3dmu2

            dAdr = dpsi0dr + dpsi1dr + dpsi3dr
            dAdth = dpsi0dth + dpsi1dth + dpsi3dth
            dAdrdth = dpsi0drdth + dpsi1drdth + dpsi3drdth
            dAdr2 = dpsi0dr2 + dpsi1dr2 + dpsi3dr2
            dAdth2 = d2psi0dth2 + d2psi1dth2 + d2psi3dth2
            ! ;
            ! ;  now calculate field
            ! ;
            ! ;  eq B1
            brlambout1 = dAdth/State%rout/State%rout/State%st
            bthlambout1 = -dAdr/State%rout/State%st
            bphlambout1 = 0.
            
            !; and calculate currents

            dbth1dr = (dAdr/State%rout - dAdr2)/State%rout/State%st
            dbr1dth = ((-State%mu/State%st)*dAdth + dAdth2)/State%rout/State%rout/State%st

            jrlambout1 = 0.
            jthlambout1 = 0.
            jphlambout1 = (1./State%rout)*(bthlambout1 + State%rout*dbth1dr - dbr1dth)

            ! ;
            ! ; But now we need to get blittle and derivatives of blittle with respect to
            ! ; the physical coord r, which we need for density and pressure  (eq. A18-A19)
            ! ;  We do so by going through the cartesian coordinates.
            ! ;

            sph = sin(State%phout)
            cph = cos(State%phout)
            call zero2tiny2d(sph)
            call zero2tiny2d(cph)

            blittleX1 = brlambout1*State%st*cph + bthlambout1*State%mu*cph
            blittleY1 = brlambout1*State%st*sph + bthlambout1*State%mu*sph
            blittleZ1 = brlambout1*State%mu - bthlambout1*State%st
            jlittleX1 = jrlambout1*State%st*cph + jthlambout1*State%mu*cph
            jlittleY1 = jrlambout1*State%st*sph + jthlambout1*State%mu*sph
            jlittleZ1 = jrlambout1*State%mu - jthlambout1*State%st

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
            ! ; and finally back to spherical coords
            ! ;
            streal = sin(State%thpb)
            ctreal = cos(State%thpb)
            spreal = sin(State%phpb)
            cpreal = cos(State%phpb)

            State%brlambout = blittleY*streal*spreal + blittleX*streal*cpreal + blittleZ*ctreal
            State%bthlambout = blittleY*ctreal*spreal + blittleX*ctreal*cpreal - blittleZ*streal
            State%bphlambout = -blittleX*spreal + blittleY*cpreal
            State%jrlambout = jlittleY*streal*spreal + jlittleX*streal*cpreal + jlittleZ*ctreal
            State%jthlambout = jlittleY*ctreal*spreal + jlittleX*ctreal*cpreal - jlittleZ*streal
            State%jphlambout = -jlittleX*spreal + jlittleY*cpreal

            jlittleX = 0.
            jlittleY = 0.
            jlittleZ = 0.
            ! ;  now calculated tderivout (needed for density  A19)
            ! ;
            dbr1dr = (-2.*dAdth/State%rout + dAdrdth)/State%rout/State%rout/State%st
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

            ! TODO: is this nan check necessary?
            where (State%brlambout /= State%brlambout)
                State%brlambout = 0.
                State%bthlambout = 0.
                State%bphlambout = 0.
                State%tderivout = 0.
                State%streamout = 0.
            end where

            deallocate (al1, al2, psi0, f, g, h)
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
            deallocate (sph, cph, blittleX1, blittleY1)
            deallocate (blittleZ1, jlittleX1, jlittleY1, jlittleZ1)
            deallocate (blittleY2, blittleZ, jlittleY2, jlittleZ)
            deallocate (blittleX, blittleY, jlittleX, jlittleY)
            deallocate (streal, ctreal, spreal, cpreal, dbr1dr)
            deallocate (dbrdr)
            deallocate (yt2)
        end subroutine calcOutside

        !> Calculate Gibson-Low Fields
        !> 
        !> Calls calcOutside after 
        subroutine calcFields(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            real(rp), dimension(:, :), allocatable :: sph, cph, gfun, gfunot
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
            real(rp) :: val
            character(len=10) :: arr_name
            integer :: dim1, dim2

            dim1 = State%Nj
            dim2 = State%Nk

            allocate (sph(dim1, dim2), cph(dim1, dim2), gfun(dim1, dim2), gfunot(dim1, dim2))
            allocate (dSdT(dim1, dim2), dgfundr(dim1, dim2))
            allocate (dSdR(dim1, dim2), blittleR(dim1, dim2), blittleT(dim1, dim2), Q(dim1, dim2), blittleP(dim1, dim2))
            allocate (dbPdt(dim1, dim2), dbPdR(dim1, dim2), d2gfundr2(dim1, dim2))
            allocate (d2SdR2(dim1, dim2), d2SdT2(dim1, dim2), dbTdR(dim1, dim2), dbRdT(dim1, dim2), jlittleR(dim1, dim2))
            allocate (jlittleT(dim1, dim2), jlittleP(dim1, dim2), blittleX(dim1, dim2), blittleY1(dim1, dim2))
            allocate (blittleZ1(dim1, dim2), jlittleX(dim1, dim2), jlittleY1(dim1, dim2), jlittleZ1(dim1, dim2))
            allocate (ctreal(dim1, dim2), streal(dim1, dim2), cpreal(dim1, dim2), spreal(dim1, dim2))
            allocate (blittley(dim1, dim2), blittlez(dim1, dim2), jlittley(dim1, dim2), jlittlez(dim1, dim2))

            ! sph = 0.0; cph = 0.0; gfun = 0.0; gfunot = 0.0
            ! dSdT = 0.0; dgfundr = 0.0
            ! dSdR = 0.0; blittleR = 0.0; blittleT = 0.0; Q = 0.0; blittleP = 0.0
            ! dbPdt = 0.0; dbPdR = 0.0; d2gfundr2 = 0.0
            ! d2SdR2 = 0.0; d2SdT2 = 0.0; dbTdR = 0.0; dbRdT = 0.0; jlittleR = 0.0
            ! jlittleT = 0.0; jlittleP = 0.0; blittleX = 0.0; blittleY1 = 0.0
            ! blittleZ1 = 0.0; jlittleX = 0.0; jlittleY1 = 0.0; jlittleZ1 = 0.0
            ! ctreal = 0.0; streal = 0.0; cpreal = 0.0; spreal = 0.0
            ! State%jlittlerlamb = 0.0; State%jlittlethlamb = 0.0; State%jlittlephlamb = 0.0
            ! blittley = 0.0; blittlez = 0.0; jlittley = 0.0; jlittlez = 0.0

            !debug
            !write(555,* ) "[EP] in giblow_fieldcalc dim1, dim2:", dim1, dim2

            State%st = sin(State%thcap)
            State%ct = cos(State%thcap)
            sph = sin(State%phcap)
            cph = cos(State%phcap)

            !debug
            !write(555,*) '[EP] in giblow_model: State%st,State%ct,sph,cph before zero2tiny2d'
            !write(555,*) dim1, dim2
            !write(555,*) minval(State%st), maxval(State%st)
            !write(555,*) minval(State%ct), maxval(State%ct)
            !write(555,*) minval(sph), maxval(sph)
            !write(555,*) minval(cph), maxval(cph)
            call zero2tiny2d(State%st) ! where (abs(State%st) .lt. tiny) State%st=tiny*sign(State%st)
            call zero2tiny2d(State%ct) ! where (abs(State%ct) .lt. tiny) State%ct=tiny*sign(State%ct)
            call zero2tiny2d(sph) ! where (abs(sph) .lt. tiny) sph=tiny*sign(sph)
            call zero2tiny2d(cph) ! where (abs(cph) .lt. tiny) cph=tiny*sign(cph)

            !write(555,*) '[EP] in giblow_model: State%st,State%ct,sph,cph after zero2tiny2d'
            !write(555,*) minval(State%st), maxval(State%st)
            !write(555,*) minval(State%ct), maxval(State%ct)
            !write(555,*) minval(sph), maxval(sph)
            !write(555,*) minval(cph), maxval(cph)

            ! eq. B8
            gfun = sin(Model%alnot*State%rcap)/Model%alnot/State%rcap - cos(Model%alnot*State%rcap)
            gfunot = sin(Model%alnot*Model%rbub)/Model%alnot/Model%rbub - cos(Model%alnot*Model%rbub)

            ! eq. B7
            State%stream = (4.*Model%ao*pi/Model%alnot/Model%alnot)*(gfun*Model%rbub**2/gfunot - State%rcap**2)*State%st**2

            ! if (Model%isDebug) then
            !     write(*,*) "gfun: ", gfun
            !     write(*,*) "gfunot: ", gfunot
            !     write(*,*) "stream: ", stream
            ! end if
            
            ! in the region Rcap ge rbub
            !  call routine that calculates stream function, etc
            !  APPENDIX B1

            if (State%outside_count .ne. 0) then
                if (Model%isDebug) write(*,*) 'giblow_fieldcalc: calling calcOutside'
                call calcOutside(Model, State, Solution)

                where (Solution%inside_mask(State%ri, :, :) .lt. 1.) State%stream = Model%outScale*State%streamout
            end if
            
            ! ;
            ! ;  define pressure
            ! ;
            State%glpi = Model%ao*State%stream + Model%Pio
            if (State%outside_count .ne. 0) where (Solution%inside_mask(State%ri, :, :) .lt. 1.) State%glpi = Model%Pio
            
            ! ;
            ! ; now calculate little b (blittle), that goes with this Stream function
            ! ; and pressure in these bubble coordinates
            ! ; To do this, we need to use the derivatives of Stream with respect
            ! ; to the cap coords.
            ! ;
            dSdT = State%stream*2*State%ct/State%st
            dgfundr = cos(Model%alnot*State%rcap)/State%rcap &
                        - sin(Model%alnot*State%rcap)/State%rcap/State%rcap/Model%alnot &
                        + Model%alnot*sin(Model%alnot*State%rcap)
            dSdR = (4.*Model%ao*pi/Model%alnot/Model%alnot)*(dgfundr*Model%rbub**2/gfunot - 2.*State%rcap)*State%st**2
            
            !;  eq. B6
            blittleR = (1/State%rcap/State%st)*((1/State%rcap)*dSdT)
            blittleT = -(1/State%rcap/State%st)*dSdR
            
            ! ;  ADD THE PHI FIELD!!!!!
            ! ; 
            ! ; NOTE! negative sign on BlittleP is to make a left-handed wind
            ! ;  about the apple core
            Q = Model%alnot*State%stream
            blittleP = -(1/State%rcap/State%st)*Q

            ! ;
            ! ; also calculate currents by taking curl
            ! ;
            dbPdt = -(1./State%rcap/State%st)*Model%alnot*dSdT
            dbPdt = dbPdt + (State%ct/State%rcap/State%st/State%st)*Q
            dbPdR = -(1./State%rcap/State%st)*Model%alnot*dSdR
            dbPdR = dbPdR + (1./State%rcap/State%rcap/State%st)*Q
            d2gfundr2 = 2.*gfun/State%rcap/State%rcap - Model%alnot*Model%alnot*gfun
            d2SdR2 = ((((Model%rbub**2)/gfunot)*d2gfundr2 - 2.)/(((Model%rbub**2)/gfunot)*gfun - State%rcap*State%rcap))*State%stream
            d2SdT2 = (4.*((State%ct/State%st)**2) - 2./State%st/State%st)*State%stream
            dbTdR = (1/State%rcap/State%rcap/State%st)*dSdR - (1/State%rcap/State%st)*d2SdR2
            dbRdT = -(State%ct/State%rcap/State%rcap/State%st/State%st)*dSdT + (1/State%rcap/State%rcap/State%st)*d2SdT2

            !debug
            !write(555,*) '[EP] in giblow_model: Thcap, Stream, Model%alnot, Model%rbub'
            !write(555,*) minval(Thcap), maxval(Thcap)
            !write(555,*) minval(Stream), maxval(Stream)
            !write(555,*) Model%alnot
            !write(555,*) Model%rbub

            !write(555,*) '[EP] in giblow_model: State%rcap, State%st, dSdT, dSdR, Q'
            !write(555,*) minval(State%rcap), maxval(State%rcap)
            !write(555,*) minval(State%st), maxval(State%st)
            !write(555,*) minval(dSdT), maxval(dSdT)
            !write(555,*) minval(dSdR), maxval(dSdR)
            !write(555,*) minval(Q), maxval(Q)

            jlittleR = (1./State%rcap/State%st)*(State%ct*blittleP + State%st*dbPdt)
            jlittleT = -(1./State%rcap)*(blittleP + State%rcap*dbPdR)
            jlittleP = (1./State%rcap)*(blittleT + State%rcap*dbTdR - dbRdT)
            ! ;
            ! ; But now we need to transform to blittle coordinates with respect to
            ! ; the physical coord r, which we need for density and pressure (eq. A19)
            ! ; (note the derivative calculation for this interior field happens in giblow - tderiv)
            ! ;  We do so by going through the cartesian coordinates.
            ! ;  note this has also now been rewritten in right-hand coordinates
            ! ;
            blittleX = blittleR*State%st*cph + blittleT*State%ct*cph - blittleP*sph
            blittleY1 = blittleR*State%st*sph + blittleT*State%ct*sph + blittleP*cph
            blittleZ1 = blittleR*State%ct - blittleT*State%st
            jlittleX = jlittleR*State%st*cph + jlittleT*State%ct*cph - jlittleP*sph
            jlittleY1 = jlittleR*State%st*sph + jlittleT*State%ct*sph + jlittleP*cph
            jlittleZ1 = jlittleR*State%ct - jlittleT*State%st

            !debug
            !write(555,*) '[EP] in giblow_model: blittleR, blittleT, blittleP'
            !write(555,*) minval(blittleR), maxval(blittleR)
            !write(555,*) minval(blittleT), maxval(blittleT)
            !write(555,*) minval(blittleP), maxval(blittleP)
            ! ;
            ! ; but now we have to translate to the sigma rotated (about x axis) frame for y and z
            ! ;  note clockwise rotation-- moves back into physical coordinates
            ! ;
            associate(sigma=>Model%sigma)

            blittleY = cos(sigma)*blittleY1 + sin(sigma)*blittleZ1
            blittleZ = cos(sigma)*blittleZ1 - sin(sigma)*blittleY1
            jlittleY = cos(sigma)*jlittleY1 + sin(sigma)*jlittleZ1
            jlittleZ = cos(sigma)*jlittleZ1 - sin(sigma)*jlittleY1

            end associate

            !debug
            !write(555,*) '[EP] in giblow_model: blittleX, blittleY, blittleZ'
            !write(555,*) minval(blittleX), maxval(blittleX)
            !write(555,*) minval(blittleY), maxval(blittleY)
            !write(555,*) minval(blittleZ), maxval(blittleZ)

            ctreal = cos(State%thpb)
            streal = sin(State%thpb)
            cpreal = cos(State%phpb)
            spreal = sin(State%phpb)
            
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

            dlamdr = 1
            !debug
            !write(555,*) '[EP] in giblow_model: State%blittlerlamb, State%rlam, State%rsquig, phiss'
            !write(555,*) minval(State%blittlerlamb), maxval(State%blittlerlamb)
            !write(555,*) minval(State%rlam), maxval(State%rlam)
            !write(555,*) minval(State%rsquig), maxval(State%rsquig)
            !write(555,*) phiss

            ! ;  and finally, here are the magnetic field coordinates in the rpB,thetapB
            ! ;  phipB coordinate systems !  we have to also get rid of selfsim stuff
            ! ;
            ! ;  eqs.  A2-A4
            
            ! ;  dlamdr is one, because rlam = rsquig - apar.
            Solution%b(State%ri, :, :, 1) = (State%blittlerlamb*(State%rlam/State%rsquig)**2)/(Model%phiss**2)
            Solution%b(State%ri, :, :, 2) = State%blittlethlamb*(State%rlam/State%rsquig)*dlamdr/(Model%phiss**2)
            Solution%b(State%ri, :, :, 3) = State%blittlephlamb*(State%rlam/State%rsquig)*dlamdr/(Model%phiss**2)
            Solution%j(State%ri, :, :, 1) = (State%jlittlerlamb*(State%rlam/State%rsquig)**2)/(Model%phiss**2)
            Solution%j(State%ri, :, :, 2) = State%jlittlethlamb*(State%rlam/State%rsquig)*dlamdr/(Model%phiss**2)
            Solution%j(State%ri, :, :, 3) = State%jlittlephlamb*(State%rlam/State%rsquig)*dlamdr/(Model%phiss**2)

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

        !> Calculate Inside (of bubble) Gibson-Low Solution
        !>
        !>
        subroutine calcInside(Model, State, Solution)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            ! real(rp), dimension(:, :), allocatable :: r1, State%rout, State%thout, State%mu, al1
            
            real(rp) :: r1
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
            real(rp), dimension(:, :), allocatable :: rnolam, dlamdr, Bpresin
            real(rp), dimension(:, :), allocatable :: dldr, lr, ufunc, dufdr
            real(rp), dimension(:, :), allocatable :: dthdr, dstuffdr, dstuffdth, dpdr
            integer :: dim1, dim2

            dim1 = State%Nj
            dim2 = State%Nk

            allocate (al1(dim1, dim2))
            allocate (al2(dim1, dim2), psi0(dim1, dim2), f(dim1, dim2), g(dim1, dim2), h(dim1, dim2), psi1(dim1, dim2))
            allocate (i(dim1, dim2), psi3(dim1, dim2), dpsi0dr(dim1, dim2), dpsi0dth(dim1, dim2), dpsi0dr2(dim1, dim2))
            allocate (dpsi0drdth(dim1, dim2), dpsi0dth2(dim1, dim2), dfdr(dim1, dim2), dfdmu(dim1, dim2))
            allocate (dfdrdmu(dim1, dim2), dfdr2(dim1, dim2), dfdmu2(dim1, dim2), dgdr(dim1, dim2), dgdmu(dim1, dim2))
            allocate (dgdrdmu(dim1, dim2), dgdr2(dim1, dim2), dgdmu2(dim1, dim2), dhdr(dim1, dim2), dhdmu(dim1, dim2))
            allocate (dhdrdmu(dim1, dim2), dhdr2(dim1, dim2), dhdmu2(dim1, dim2), dpsi1dr(dim1, dim2))
            allocate (dpsi1dr2(dim1, dim2), dpsi1dmu(dim1, dim2), dpsi1drdmu(dim1, dim2), dpsi1dmu2(dim1, dim2))
            allocate (dpsi1dth(dim1, dim2), dpsi1drdth(dim1, dim2), dpsi1dth2(dim1, dim2), didr(dim1, dim2))
            allocate (didmu(dim1, dim2), didrdmu(dim1, dim2), didr2(dim1, dim2), didmu2(dim1, dim2))
            allocate (dpsi3dr(dim1, dim2), dpsi3dmu(dim1, dim2), dpsi3drdmu(dim1, dim2), dpsi3dr2(dim1, dim2))
            allocate (dpsi3dmu2(dim1, dim2), dpsi3dth(dim1, dim2), dpsi3drdth(dim1, dim2), dpsi3dth2(dim1, dim2))
            allocate (dAdr(dim1, dim2), dAdth(dim1, dim2), dAdrdth(dim1, dim2), dAdr2(dim1, dim2), dAdth2(dim1, dim2))
            allocate (dbrdr(dim1, dim2), dbthdr(dim1, dim2), dbrdth(dim1, dim2), dbthdth(dim1, dim2))
            allocate (rnolam(dim1, dim2), dlamdr(dim1, dim2), Bpresin(dim1, dim2))
            allocate (dldr(dim1, dim2), lr(dim1, dim2), ufunc(dim1, dim2), dufdr(dim1, dim2))
            allocate (dthdr(dim1, dim2), dstuffdr(dim1, dim2), dstuffdth(dim1, dim2), dpdr(dim1, dim2))

            ! al1 = 0.0
            ! al2 = 0.0; psi0 = 0.0; f = 0.0; g = 0.0; h = 0.0; psi1 = 0.0
            ! i = 0.0; psi3 = 0.0; dpsi0dr = 0.0; dpsi0dth = 0.0; dpsi0dr2 = 0.0
            ! dpsi0drdth = 0.0; dpsi0dth2 = 0.0; dfdr = 0.0; dfdmu = 0.0
            ! dfdrdmu = 0.0; dfdr2 = 0.0; dfdmu2 = 0.0; dgdr = 0.0; dgdmu = 0.0
            ! dgdrdmu = 0.0; dgdr2 = 0.0; dgdmu2 = 0.0; dhdr = 0.0; dhdmu = 0.0
            ! dhdrdmu = 0.0; dhdr2 = 0.0; dhdmu2 = 0.0; dpsi1dr = 0.0
            ! dpsi1dr2 = 0.0; dpsi1dmu = 0.0; dpsi1drdmu = 0.0; dpsi1dmu2 = 0.0
            ! dpsi1dth = 0.0; dpsi1drdth = 0.0; dpsi1dth2 = 0.0; didr = 0.0
            ! didmu = 0.0; didrdmu = 0.0; didr2 = 0.0; didmu2 = 0.0
            ! dpsi3dr = 0.0; dpsi3dmu = 0.0; dpsi3drdmu = 0.0; dpsi3dr2 = 0.0
            ! dpsi3dmu2 = 0.0; dpsi3dth = 0.0; dpsi3drdth = 0.0; dpsi3dth2 = 0.0
            ! dAdr = 0.0; dAdth = 0.0; dAdrdth = 0.0; dAdr2 = 0.0; dAdth2 = 0.0
            ! dbrdr = 0.0; dbthdr = 0.0; dbrdth = 0.0; dbthdth = 0.0
            ! rnolam = 0.0; dlamdr = 0.0; Bpresin = 0.0
            ! dldr = 0.0; lr = 0.0; ufunc = 0.0; dufdr = 0.0
            ! dthdr = 0.0; dstuffdr = 0.0; dstuffdth = 0.0; dpdr = 0.0
            ! rat = 0.0

            r1 = Model%xo

            ! where(inside) State%rout = rpb
            ! State%rout = rpb
            State%rout = State%rlam

            ! ; calculte theta along the boundary, and rotate so that the
            ! ; bubble is at the equator
            State%rat = ((r1*r1 + State%rout*State%rout - Model%rbub*Model%rbub)/(2*State%rout*r1))
            where (State%rat .gt. 1.) State%rat = 1
            where (State%rat .lt. -1.) State%rat = -1
            State%thout = acos(State%rat)

            State%mu = cos(State%thout)
            State%st = sin(State%thout)

            call zero2tiny2d(State%mu) ! where (abs(State%st) .lt. tiny) State%st=tiny*sign(State%st)
            call zero2tiny2d(State%ct) ! where (abs(State%ct) .lt. tiny) State%ct=tiny*sign(State%ct)

            al1 = r1*r1 - Model%rbub*Model%rbub
            al2 = 2*r1*r1 - Model%rbub*Model%rbub
           
            ! ;
            ! ;  eq. B2
            psi0 = State%mu

            ! ;
            ! ;  eq. B3
            f = r1*(State%rout*State%rout + al1) - State%rout*State%mu*al2
            g = al1*al1 + r1*r1*State%rout*State%rout - 2.*State%rout*r1*al1*State%mu
            h = 1./(sqrt(g))
            psi1 = (-1./Model%rbub)*f*h

            ! ;
            ! ;  eq. B4
            i = State%rout*State%rout + r1*r1 - 2.*State%rout*r1*State%mu
            psi3 = (1/Model%rbub)*sqrt(i)

            ! ;
            ! ;  eq. B5
            State%streamout = psi0 + psi1 + psi3

            ! ;
            ! ;  now take derivatives for magnetic fields, etc.
            ! ;
            dpsi0dr = 0.
            dpsi0dth = -State%st
            dpsi0dr2 = 0.
            dpsi0drdth = 0.
            dpsi0dth2 = -1.*State%mu

            dfdr = 2.*r1*State%rout - State%mu*al2
            dfdmu = -State%rout*al2
            dfdrdmu = -al2
            dfdr2 = 2.*r1
            dfdmu2 = 0.

            dgdr = 2.*State%rout*r1*r1 - 2.*r1*al1*State%mu
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

            dpsi1dth = -State%st*dpsi1dmu
            dpsi1drdth = -State%st*dpsi1drdmu
            dpsi1dth2 = -State%mu*dpsi1dmu + State%st**2*dpsi1dmu2

            didr = 2.*State%rout - 2.*r1*State%mu
            didmu = -2.*State%rout*r1
            didrdmu = -2.*r1
            didr2 = 2.
            didmu2 = 0.

            dpsi3dr = (1./Model%rbub)*didr/2./sqrt(i)
            dpsi3dmu = (1./Model%rbub)*didmu/2./sqrt(i)
            dpsi3drdmu = (1./Model%rbub)*(didrdmu/2./sqrt(i) - didmu*didr/4./(i**1.5))
            dpsi3dr2 = (1./Model%rbub)*(didr2/2./sqrt(i) - (didr)**2./4./(i**1.5))
            dpsi3dmu2 = (1./Model%rbub)*didmu2/2./sqrt(i) - (1./Model%rbub)*(didmu**2.)/4./(i**1.5)

            dpsi3dth = -State%st*dpsi3dmu
            dpsi3drdth = -State%st*dpsi3drdmu
            dpsi3dth2 = -State%mu*dpsi3dmu + State%st**2.*dpsi3dmu2

            dAdr = dpsi0dr + dpsi1dr + dpsi3dr
            dAdth = dpsi0dth + dpsi1dth + dpsi3dth
            dAdrdth = dpsi0drdth + dpsi1drdth + dpsi3drdth
            dAdr2 = dpsi0dr2 + dpsi1dr2 + dpsi3dr2
            dAdth2 = dpsi0dth2 + dpsi1dth2 + dpsi3dth2
            ! ;
            ! ;  now calculate field
            ! ;
            ! ;  eq. B1
            State%brlambout = dAdth/State%rout**2/State%st
            State%bthlambout = -dAdr/State%rout/State%st
            State%bphlambout = 0.
            ! ;
            ! ;  now calculate field derivatives
            ! ;
            dbrdr = (-2.*dAdth/State%rout + dAdrdth)/State%rout/State%rout/State%st
            dbthdr = (dAdr/State%rout - dAdr2)/State%rout/State%st

            dbrdth = -State%mu*dAdth/State%rout/State%rout/(State%st)**2
            dbrdth = dbrdth + dAdth2/State%rout/State%rout/State%st
            dbthdth = State%mu*dAdr/State%rout/(State%st)**2
            dbthdth = dbthdth - dAdrdth/State%rout/State%st

            rnolam = State%rout - Model%apar
            State%presin = (State%rout**2/rnolam**2)*(1 - (State%rout**2/rnolam**2))*(State%brlambout**2)/8./pi
            State%presin = State%presin/(Model%phiss**4)
            ! ;  here are the magnetic pressure in the rpB,thetapB
            ! ;  phipB coordinate systems !  we have to also get rid of selfsim stuff
            ! ;
            ! ;  dlamdr is one, because rlam = rsquig + apar.
            dlamdr = 1.
            State%Bpresin = (( State%brlambout*(State%rout**2/rnolam**2)/(Model%phiss**2))**2 &
                            + (State%bthlambout*(State%rout/rnolam)*dlamdr/(Model%phiss**2))**2)/8./pi
            State%Pbackin = State%presin +  State%Bpresin
            State%Pbackin =  State%Pbackin*Model%outScale*Model%outScale
            ! ;
            ! ;  now calculate the background density
            ! ;
            dldr = -Model%apar/rnolam/rnolam
            lr = (Model%apar + rnolam)/rnolam

            ufunc = (r1*r1 + State%rout*State%rout - Model%rbub*Model%rbub)/(2*State%rout*r1)
            dufdr = (1./r1) - ufunc/State%rout

            dthdr = -(1./sqrt(1. - ufunc*ufunc))*dufdr
            dstuffdr = -(lr*dldr*State%brlambout*State%brlambout + lr*lr*State%brlambout*dbrdr + &
                            .5*dldr*State%bthlambout*State%bthlambout + lr*State%bthlambout*dbthdr)/4./pi
            dstuffdth = -(lr*lr*State%brlambout*dbrdth + lr*State%bthlambout*dbthdth)/4./pi
            dpdr = dstuffdr + dstuffdth*dthdr

            F = GMm/((State%rout - Model%apar)**2)/(Rsun**2) + Model%alpha*(State%rout - Model%apar)*Rsun*mprot
            State%Dbackin = dpdr/F/Rsun/(Model%phiss**3)
            State%Dbackin = State%Dbackin*Model%outScale*Model%outScale

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
            deallocate (rnolam, dlamdr, Bpresin)
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

            allocate(dphidt(State%Nj, State%Nk))

            ! ; also velocity (in km/sec)
            ! ; note this is radial

            dphidt = sqrt((Model%eta*Model%phiss - 2.*Model%alpha)/Model%phiss)
            Solution%v(State%ri, :, :, XDIR)= (State%rpb*6.96d5*dphidt)/Model%phiss
    
            where (Solution%dens(State%ri, :, :) .gt. 0) Solution%Temp(State%ri, :, :) = Solution%Pres(State%ri, :, :) / 2./kboltz/Solution%dens(State%ri, :, :)
            ! if ((gl_verbose.eq.1).and.(maxval(temp).eq.0)) then
            !     print *, 'warning: all temp values are zeroes, the bubble is likely outside of the computational domain'
            ! end if
            ! ;Temp = 6.02d-9*Pres/dens/1.6696d-24
    
            if (Model%bonly .gt. 0) then
                Solution%Temp(State%ri, :, :) = Solution%Temp(State%ri, :, :)*0. + Model%Isothermal
                Solution%dens(State%ri, :, :) = State%DensbackHEonly
                Solution%Pres(State%ri, :, :) = gas_R*Solution%Temp(State%ri, :, :)*mprot*Solution%dens(State%ri, :, :)
                if (State%inside_count .gt. 0) then
                    where (Solution%inside_mask(State%ri, :, :) .eq. 1) 
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
                where (Solution%inside_mask(State%ri, :, :) .eq. 0)
                    Solution%dens(State%ri, :, :) = 1d-5
                    Solution%temp(State%ri, :, :)= 1d-5
                    Solution%pres(State%ri, :, :) = 1d-5
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
        end subroutine finishSolution

        !> For a initialized Model, State, Solution
        !>
        !>
        subroutine generateGLSolution(Model, State, Solution)
            type(glModel_T), intent(inout) :: Model
            type(glState_T), intent(inout)  :: State
            type(glSolution_T), intent(inout)  :: Solution
            integer :: i, rdim

            rdim = size(State%r)
            if(Model%isDebug) write(*,*) "Generating GL Solution"

            ! time assumed in seconds
            Model%phiss = sqrt(Model%eta)*Model%time + 1.0 

            do i = 1, rdim
                State%rpb = State%r(i)
                State%ri = i

                call calcCoords(Model, State)

                where (State%rcap .lt. Model%rbub) Solution%inside_mask(State%ri,:,:) = 1.              
                State%inside_count = sum(Solution%inside_mask(State%ri,:,:))
                State%outside_count = sum(1. - Solution%inside_mask(State%ri,:,:))   
    
                call calcFields(Model, State, Solution)
                call calcPresssure(Model, State, Solution)
                call calcDerivs(Model, State, Solution)
                call calcDensity(Model, State, Solution)
                call finishSolution(Model, State, Solution)

                ! if (Model%isDebug) then
                !     write(*,"(1X,A20,2X,E13.6)") "current r: ", State%r(i)
                !     write(*,"(1X,A20,2X,E13.6)") "bubble r: ",  Model%rbub
                !     write(*,"(1X,A20,2X,E13.6)") "Max rlam: ",  maxval(State%rlam)
                !     write(*,"(1X,A20,2X,E13.6)") "Min rlam: ",  minval(State%rlam)
                !     write(*,"(1X,A20,2X,E13.6)") "Max rcap: ",  maxval(State%rcap)
                !     write(*,"(1X,A20,2X,E13.6)") "Min rcap: ",  minval(State%rcap)
                !     write(*,"(1X,A20,2X,E13.6)") "Max rsquig: ",  maxval(State%rsquig)
                !     write(*,"(1X,A20,2X,E13.6)") "Min rsquig: ",  minval(State%rsquig)
                !     write(*,"(1X,A20,2X,E13.6)") "Max densback: ",  maxval(State%densback)
                !     write(*,"(1X,A20,2X,E13.6)") "Min densback: ",  minval(State%densback)
                !     write(*,"(1X,A20,2X,E13.6)") "Max densin: ",  maxval(State%densin)
                !     write(*,"(1X,A20,2X,E13.6)") "Min densin: ",  minval(State%densin)
                !     write(*,"(1X,A20,2X,E13.6)") "Max densout: ",  maxval(State%densout)
                !     write(*,"(1X,A20,2X,E13.6)") "Min densout: ",  minval(State%densout)
                !     write(*,"(1X,A20,2X,E13.6)") "Max solDens: ",  maxval(Solution%dens(i,:,:))
                !     write(*,"(1X,A20,2X,E13.6)") "Min solDens: ",  minval(Solution%dens(i,:,:))
                !     write(*,"(1X,A20,2X,F13.1)") "inside_count: ", State%inside_count
                !     write(*,"(1X,A20,2X,F13.1)") "outside_count: ", State%outside_count
                ! end if
            end do
            if (Model%isLoud) then
                write(*,"(1X,A20,2X,E13.6)") "max br phi=0.0: ", maxval(Solution%b(:,:,1,XDIR))
                write(*,"(1X,A20,2X,E13.6)") "min br phi=0.0: ", minval(Solution%b(:,:,1,XDIR))
                write(*,"(1X,A20,2X,E13.6)") "max btheta phi=0.0: ", maxval(Solution%b(:,:,1,YDIR))
                write(*,"(1X,A20,2X,E13.6)") "min btheta phi=0.0: ", minval(Solution%b(:,:,1,YDIR))
                write(*,"(1X,A20,2X,E13.6)") "max bphi phi=0.0: ", maxval(Solution%b(:,:,1,ZDIR))
                write(*,"(1X,A20,2X,E13.6)") "min bphi phi=0.0: ", minval(Solution%b(:,:,1,ZDIR))
                write(*,"(1X,A20,2X,E13.6)") "max dens phi=0.0: ", maxval(Solution%dens(:,:,1))
                write(*,"(1X,A20,2X,E13.6)") "min dens phi=0.0: ", minval(Solution%dens(:,:,1))
                write(*,*) "inside bubble"
                write(*,"(1X,A20,2X,E13.6)") "max dens: ", maxval(Solution%dens*merge(1,0,Solution%inside_mask > 0.0))
                write(*,"(1X,A20,2X,E13.6)") "min dens: ", minval(Solution%dens*merge(1,0,Solution%inside_mask > 0.0))
                write(*,*) "outside bubble"
                write(*,"(1X,A20,2X,E13.6)") "max dens: ", maxval(Solution%dens*merge(1,0,Solution%inside_mask < 1.0))
                write(*,"(1X,A20,2X,E13.6)") "min dens: ", minval(Solution%dens*merge(1,0,Solution%inside_mask < 1.0))
            end if
            if(Model%isDebug) write(*,*) "GL Solution Complete"
        end subroutine generateGLSolution 
end module