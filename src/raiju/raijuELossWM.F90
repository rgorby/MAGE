module raijuELossWM
    use ioh5
    use files
    
    use raijudefs
    use raijutypes
    use raijuspecieshelper, only : spcIdx

    implicit none

    integer, parameter, private :: MAXIOVAR = 50

    contains

    subroutine initEWM(eWM, configFname, xmlInp, shGrid)
        type(eLossWM_T), intent(inout) :: eWM
        character(len=strLen) :: configFname
        type(XML_Input_T), intent(in) :: xmlInp
        type(ShellGrid_T), intent(in) :: shGrid

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        type(IOVar_T) :: tauVAR

        ! Determine if we are expected to provide info to output file
        call xmlInp%Set_Val(eWM%doOutput,"/Kaiju/RAIJU/losses/doOutput",eWM%doOutput)

        ! We need Kp from wind file, try to get that first
        call xmlInp%Set_Val(eWM%KpTS%wID,"/Kaiju/Gamera/wind/tsfile","NONE")
        call eWM%KpTS%initTS("Kp",doLoudO=.false.)

        call ClearIO(IOVars)

        if(.not. ioExist(configfname, "Tau", "waveModel")) then
            write(*,*) "This config file not structured for RAIJU, use genRAIJU.py. Good day."
            stop
        endif

        !Chorus wave
        call AddInVar(IOVars,"Kp")  
        call AddInVar(IOVars,"MLT") 
        call AddInVar(IOVars,"L") 
        call AddInVar(IOVars,"Ek") 
        call AddInVar(IOVars,"Tau")
        call ReadVars(IOVars,.false.,configFname, "waveModel")

        eWM%Nkp  = IOVars(FindIO(IOVars, "Kp"))%N
        eWM%Nmlt = IOVars(FindIO(IOVars, "MLT"))%N
        eWM%Nl   = IOVars(FindIO(IOVars, "L"))%N
        eWM%Ne   = IOVars(FindIO(IOVars, "Ek"))%N

        tauVAR = IOVars(FindIO(IOVars, "Tau"))

        ! Do some dimension checks
        if (tauVAR%Nr /= 4) then
            write(*,*) "tauDim:",tauVAR%Nr
            write(*,*) 'Currently only support tau model files in the form tau(Kp,MLT,L,Ek)'
            stop
        endif

        if (eWM%Nkp /=  tauVAR%dims(1) .or. eWM%Nmlt /=  tauVAR%dims(2) .or. eWM%Nl /=  tauVAR%dims(3) .or. eWM%Ne /=  tauVAR%dims(4)) then
            write(*,*) "tauDims:",tauVAR%dims,"Nk:",eWM%Nkp,"Nm:",eWM%Nmlt,"Nl:",eWM%Nl,"Ne:",eWM%Ne
            write(*,*) 'Dimensions of tau arrays are not compatible'
            stop
        endif

        ! If still here we're good to allocate and store data

        allocate(eWM%Kp1D    (eWM%Nkp ))
        allocate(eWM%MLT1D   (eWM%Nmlt))
        allocate(eWM%L1D     (eWM%Nl  ))
        allocate(eWM%Energy1D(eWM%Ne  ))
        allocate(eWM%Tau4D(eWM%Nkp,eWM%Nmlt,eWM%Nl,eWM%Ne))

        call IOArray1DFill(IOVars,"Kp" , eWM%Kp1D    )
        call IOArray1DFill(IOVars,"MLT", eWM%MLT1D   )
        call IOArray1DFill(IOVars,"L"  , eWM%L1D     )
        call IOArray1DFill(IOVars,"Ek" , eWM%Energy1D)
        call IOArray4DFill(IOVars,"Tau", eWM%Tau4D   )

        call ClearIO(IOVars)

        !Array order check: array is in acsending order
           !Chorus
        if(eWM%Kp1D(1) > eWM%Kp1D(eWM%Nkp)) then
            write(*,*) "Kp: ",eWM%Kp1D
            write(*,*) "reorder wave model so Kp is in ascending order"
            stop
        end if

        if(eWM%MLT1D(1) > eWM%MLT1D(eWM%Nmlt)) then
            write(*,*) "MLT: ",eWM%MLT1D
            write(*,*) "reorder wave model so MLT is in ascending order"
            stop
        end if

        if(eWM%L1D(1) > eWM%L1D(eWM%Nl)) then
            write(*,*) "L: ",eWM%L1D
            write(*,*) "reorder wave model so L shell is in ascending order"
            stop
        end if

        if(eWM%Energy1D(1) > eWM%Energy1D(eWM%Ne)) then
            write(*,*) "Ek: ",eWM%Energy1D
            write(*,*) "reorder wave model so Ek is in ascending order"
            stop
        end if

        if (eWM%doOutput) then
            allocate(eWM%wPS    (shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg))
            allocate(eWM%wHISS  (shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg))
            allocate(eWM%wCHORUS(shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg))
        endif
         

    end subroutine initEWM


    function calcELossRate_WM(Model, Grid, State, i, j, k) result(lossRate2)
        type(raijuModel_T) , intent(inout) :: Model
        type(raijuGrid_T)  , intent(in) :: Grid
        type(raijuState_T) , intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp), dimension(2) :: lossRate2

        integer :: psphIdx, eleIdx
        real(rp) :: NpsphPnt
            !! Density [#/cc] of plasmasphere at point i,j
        real(rp) :: L, MLT, E, Kp
            !! L shell and MLT of given point, channel energy in MeV, current Kp
        real(rp) :: wLBlend, wNBlend
            !! L-weighting of blending between IMAG and PS. 0=PS
            !! Density-weighting between Chorus and Hiss
        real(rp) :: wPS, wHISS, wCHORUS
            !! Weights for plasma sheet, plasmasphere hiss, and chorus models

        associate(eWM=>Model%eLossWM)

            lossRate2 = 0.0
                !! (1) = value, (2) = type

            L = sqrt(State%xyzMin(i,j,1)**2 + State%xyzMin(i,j,2)**2)  ! [Re]
            MLT = atan2(State%xyzMin(i,j,2),State%xyzMin(i,j,1))/pi*12.D0+12.D0
            E = abs(Grid%alamc(k) * State%bvol(i,j)**(-2./3.)) * 1.0E-6  ! [MeV]
            Kp = eWM%KpTS%evalAt(State%t)

            psphIdx = spcIdx(Grid, F_PSPH)
            NpsphPnt = State%Den(i,j,psphIdx+1)  ! Add 1 cause we're grabbing from density, which has bulk as first element
            
            ! Calculate blending
            wNBlend = dlog(NpsphPnt/eWM%NpsphLow) / dlog(eWM%NpsphHigh/eWM%NpsphLow)
            call ClampValue(wNBlend, 0.0_rp, 1.0_rp)
                !! 1 => Psphere Hiss, 0 => Other
            wLBlend = RampDown(L, eWM%ChorusLMax, eWM%PsheetLMin - eWM%ChorusLMax)
                !! 1 => Chorus, 0 => PS
            
            ! If psphere density is high enough, we always apply Hiss regardless of L
                ! This means Hiss model needs to appropriately handle any L value its given
            ! If psphere density is low, we choose between Chorus and plasma sheet strong-scattering based on L
            ! Also note: the weighting is not k-dependent, but if its not too costly to do for each k then ¯\_(ツ)_/¯
            
            ! Calculate weights
            wHISS = wNBlend
            wCHORUS = (1 - wNBlend)*wLBlend
            wPS = (1 - wNBlend)*(1-wLBlend)

            ! Make weighting variable, idk if this will work well
            lossRate2(2) = lossRate2(2) + int(wHISS*10)
            lossRate2(2) = lossRate2(2) + int(wCHORUS*100)
            lossRate2(2) = lossRate2(2) + int(wPS*1000)

            ! Calculate loss rates and accumulate
            if (wHiss > TINY) then
                lossRate2(1) = lossRate2(1) + wHISS*calcHissRate(MLT, L, E, Kp)
            endif

            ! Save output info if expected
            ! All electron k's will use this function, only first k should write something
            eleIdx = Grid%k2spc(k)
            if (eWM%doOutput .and. k .eq. Grid%spc(eleIdx)%kStart) then
                eWM%wPS(i,j)     = wPS
                eWM%wHISS(i,j)   = wHISS
                eWM%wCHORUS(i,j) = wCHORUS
            endif

        end associate

    end function calcELossRate_WM

    function calcHissRate(MLTin, Lin, Ein, Kpin) result(rateK)
        ! Empirical lifetime against plasmaspheric hiss pitch angle diffusion, based on Orlova et al. 2015JA021878.
        ! Improvements relative to 2014GL060100: 1. Hiss wave intensity distribution model is based on new data 
        ! (O14 was based on single-component E field in CRRES data. O16 used Spasojevic+2015 model based on EMFISIS B data on VAP); 
        ! 2. Wave spectrum is assumed differently (O14 assume Gaussian spectrum based on CRRES data).
        ! Electron lifetime tau(L,E,MLT,Kp) = tau_av(L,E)/g(MLT)/h(Kp), 
        ! where 1.5<L<5.5, E=log10(Ek[MeV]) for 1 KeV < Ek < 10 MeV.
        ! log10(tau_av(L,E)) = a1+a2*L+a3*E+...+a20*E^3, when E >= f(L).
        ! f(L) = 0.1328*L^2-2.1463*L+3.7857.
        ! g(MLT) = 10^g0(MLT)/G0
        ! h(Kp) = 10^h0(Kp)/H0
        !   G0 = int_0^24(10^g0(MLT))dMLT / 24 = 782.3.
        !   g0(MLT) = b2*MLT^2 + b1*MLT + b0
        !   H0 = 1315.
        !   h0(Kp) = c2*Kp^2 + c1*Kp + c0
        real(rp), intent(in) :: MLTin, Lin
            !! L in Re
        real(rp), intent(in) :: Ein, Kpin
            !! E in MeV

        real(rp) :: tau_av, tau, rateK
            !! tau in s, rate in 1/s
        real(rp) :: MLT, L, E, Kp
            !! Parames we clamp if needed
        real(rp) :: L2, L3, L4, fL, E2, E3, E4, E5, LE
            !! Polynomials
        real(rp), dimension(20) :: le_pol
            !! Container for polynomials
        real(rp) :: b0, b1, b2, G0, g0_MLT, g_MLT, c0, c1, c2, H0, h0_Kp, h_Kp
            !! Coefficients for g_MLT and h_Kp
        real(rp), dimension(20), parameter :: a1_20 = [77.323, -92.641, -55.754, 44.497, 48.981, 8.9067, -10.704, &
                                                         -15.711, -3.3326, 1.5189, 1.294, 2.2546, 0.31889, -0.85916, & 
                                                         -0.22182, 0.034318, 0.097248, -0.12192, -0.062765, 0.0063218]
            !! Coefficients for polynomials

        rateK = 0.0

        MLT = MLTin
        Kp = Kpin
        L = Lin
        call ClampValue(L,1.5_rp,5.5_rp)
        L2 = L*L
        fL = 0.1328*L2 - 2.1463*L + 3.7857
        E = log10(Ein)
        call ClampValue(E,max(-3.0_rp,fL),1.0_rp) ! 1 keV to 1 Mev

        E2 = E*E
        E3 = E2*E
        E4 = E3*E
        E5 = E4*E
        L3 = L2*L
        L4 = L3*L
        LE = L*E
        le_pol = (/1.D0,L,E,L2,LE,E2,L3,L2*E,L*E2,E3,L4,L3*E,L2*E2,L*E3,E4,L*E4,L2*E3,L4*E,L2*L3,E5/)

        b0 = 2.080
        b1 = 0.1773
        b2 = -0.007338
        G0 = 782.3
        g0_MLT = b2*MLT*MLT + b1*MLT + b0
        g_MLT = 10**g0_MLT/G0
        c0 = 2.598
        c1 = 0.2321
        c2 = -0.01414
        H0 = 1315.0
        Kp = min(Kp,5.0_rp) ! 0<Kp<5, 1.2% data in Kp bin of 4.3-7.6
        h0_Kp = c2*Kp*Kp + c1*Kp + c0
        h_Kp = 10**h0_Kp/H0

        tau_av = 10.0**(dot_product(a1_20,le_pol))*86400.D0 ! seconds
        tau = tau_av/g_MLT/h_Kp
        rateK = 1.0/tau

        !write(*,*) Ein, L, fL, E, tau, rateK, '--'
        !stop

    end function calcHissRate

    subroutine eWMOutput(Model, Grid, State, gStr, doGhostsO)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: gStr
        logical, optional, intent(in) :: doGhostsO

        character(len=strLen) :: gStrWM
        integer :: is, ie, js, je, ks, ke
        logical :: doGhosts
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        if (doGhosts) then
            is = Grid%shGrid%isg
            ie = Grid%shGrid%ieg
            js = Grid%shGrid%jsg
            je = Grid%shGrid%jeg
        else
            is = Grid%shGrid%is
            ie = Grid%shGrid%ie
            js = Grid%shGrid%js
            je = Grid%shGrid%je
        endif

        write(gStrWM, "(A, A)") trim(gStr), "/eWM"
        !Reset IO chain
        call ClearIO(IOVars)

        ! Fat output for precip info
        call AddOutVar(IOVars, "wPS"    , Model%eLossWM%wPS    (is:ie,js:je), dStr="0-1 weighting for plasma sheet strong scattering losses")
        call AddOutVar(IOVars, "wHISS"  , Model%eLossWM%wHISS  (is:ie,js:je), dStr="0-1 weighting for plasmasphere hiss losses")
        call AddOutVar(IOVars, "wCHORUS", Model%eLossWM%wCHORUS(is:ie,js:je), dStr="0-1 weighting for inner mag chorus wave losses")
        call WriteVars(IOVars,.true.,Model%raijuH5, gStrWM)

    end subroutine eWMOutput

end module raijuELossWM