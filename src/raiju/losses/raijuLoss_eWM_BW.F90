module raijuLoss_eWM_BW
    ! Bao-Wang electron Wave Model - based loss and precipitation model

    use kdefs, only : mec2, vc_cgs
    use raijudefs
    use raijutypes
    use raijuSpeciesHelper
    use raijuLoss_SS

    use arrayutil

    implicit none

    integer, parameter, private :: MAXIOVAR = 50

    type, extends(baseRaijuLoss_T) :: raiLoss_eWM_BW_T
        logical :: reqsGood = .false.

        type(raiLoss_SS_T) :: loss_SS
            !! Our personal instance of the strong scattering implementation to apply where needed

        ! -- Model -- !
        real(rp) :: NpsphHigh  = 100.0  ! [#/cc]
        real(rp) :: NpsphLow   = 10.0   ! [#/cc]
        real(rp) :: ChorusLMax = 7.0  ! [Re]
        real(rp) :: PsheetLMin = 8.0  ! [Re]
        real(rp) :: ChorusEMin = 1.1  ! [keV]
        
        type(TimeSeries_T) :: KpTS
            !! Kp data from wind file

        ! -- Grid -- !
        ! Chorus info
        integer :: Nkp, Nmlt, Nl, Ne
            !! Number of bins for Kp, MLT, L shell, and Energy
        real(rp), allocatable, dimension(:) :: Kp1D
            !! 1D array of Kp dimension for Tau4D
        real(rp), allocatable, dimension(:) :: MLT1D
            !! 1D array of MLT dimension for Tau4D
        real(rp), allocatable, dimension(:) :: L1D
            !! 1D array of L shell dimension for Tau4D [Re]
        real(rp), allocatable, dimension(:) :: Energy1D
            !! 1D array of energy dimension for Tau4D [MeV]
        real(rp), allocatable, dimension(:,:,:,:) :: Tau4D
            !! Tau(Kp, MLT, L, E) table electron scattering table [seconds]

        ! -- State -- !
        real(rp), allocatable, dimension(:,:) :: wPS
        real(rp), allocatable, dimension(:,:) :: wHISS
        real(rp), allocatable, dimension(:,:) :: wCHORUS

        real(rp), allocatable, dimension(:,:,:) :: tauTotal
            !! (i,j,k) Final taus considering all processes and weights
        

        contains
        
        procedure :: doInit => eWM_BW_LossInit
        procedure :: doUpdate => eWM_BW_DoUpdate
        procedure :: isValidSpc => eWM_BW_LossValidSpc
        procedure :: calcTau    => eWM_BW_CalcTau
        procedure :: doOutput => eWM_BW_DoOutput

    end type raiLoss_eWM_BW_T

    contains

    subroutine eWM_BW_LossInit(this, Model, Grid, xmlInp)
        class(raiLoss_eWM_BW_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(XML_Input_T) , intent(in) :: xmlInp
        
        type(raijuSpecies_T) :: spc_ele
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        type(IOVar_T) :: tauVAR

        ! Only valid for Earth
        if (trim(toUpper(Model%planet%name)) .ne. "EARTH") then
            write(*,*)"WARNING in raijuLoss_eWM_BW: Not simulating Earth, idk what waves are like elswhere"
            write(*,*)"No electron wave model happening"
            return
        else
            this%reqsGood = .true.
        endif
        this%isPrecip = .true.

        ! If we are still here then we are gonna be used, set everything up

        spc_ele = Grid%spc(spcIdx(Grid, F_HOTE))

        ! Init our instance of RAIJU strong scattering model
        call this%loss_SS%doInit(Model, Grid, xmlInp)

        ! We need Kp from wind file, try to get that first
        call xmlInp%Set_Val(this%KpTS%wID,"/Kaiju/Gamera/wind/tsfile","NONE")
        call this%KpTS%initTS("Kp",doLoudO=.false.)

        call ClearIO(IOVars)

        if(.not. ioExist(Model%configFName, "Tau", "waveModel")) then
            write(*,*) "This config file not structured for RAIJU, use genRAIJU.py. Good day."
            stop
        endif

        !Chorus wave
        call AddInVar(IOVars,"Kp")  
        call AddInVar(IOVars,"MLT") 
        call AddInVar(IOVars,"L") 
        call AddInVar(IOVars,"Ek") 
        call AddInVar(IOVars,"Tau")
        call ReadVars(IOVars,.false.,Model%configFname, "waveModel")

        this%Nkp  = IOVars(FindIO(IOVars, "Kp" ))%N
        this%Nmlt = IOVars(FindIO(IOVars, "MLT"))%N
        this%Nl   = IOVars(FindIO(IOVars, "L"  ))%N
        this%Ne   = IOVars(FindIO(IOVars, "Ek" ))%N

        tauVAR = IOVars(FindIO(IOVars, "Tau"))

        ! Do some dimension checks
        if (tauVAR%Nr /= 4) then
            write(*,*) "tauDim:",tauVAR%Nr
            write(*,*) 'Currently only support tau model files in the form tau(Kp,MLT,L,Ek)'
            stop
        endif

        if (this%Nkp /=  tauVAR%dims(1) .or. this%Nmlt /=  tauVAR%dims(2) .or. this%Nl /=  tauVAR%dims(3) .or. this%Ne /=  tauVAR%dims(4)) then
            write(*,*) "tauDims:",tauVAR%dims,"Nk:",this%Nkp,"Nm:",this%Nmlt,"Nl:",this%Nl,"Ne:",this%Ne
            write(*,*) 'Dimensions of tau arrays are not compatible'
            stop
        endif

        ! If still here we're good to allocate and store data

        allocate(this%Kp1D    (this%Nkp ))
        allocate(this%MLT1D   (this%Nmlt))
        allocate(this%L1D     (this%Nl  ))
        allocate(this%Energy1D(this%Ne  ))
        allocate(this%Tau4D(this%Nkp,this%Nmlt,this%Nl,this%Ne))

        call IOArray1DFill(IOVars,"Kp" , this%Kp1D    )
        call IOArray1DFill(IOVars,"MLT", this%MLT1D   )
        call IOArray1DFill(IOVars,"L"  , this%L1D     )
        call IOArray1DFill(IOVars,"Ek" , this%Energy1D)
        call IOArray4DFill(IOVars,"Tau", this%Tau4D   )

        call ClearIO(IOVars)

        !Array order check: array is in acsending order
           !Chorus
        if(this%Kp1D(1) > this%Kp1D(this%Nkp)) then
            write(*,*) "Kp: ",this%Kp1D
            write(*,*) "reorder wave model so Kp is in ascending order"
            stop
        end if

        if(this%MLT1D(1) > this%MLT1D(this%Nmlt)) then
            write(*,*) "MLT: ",this%MLT1D
            write(*,*) "reorder wave model so MLT is in ascending order"
            stop
        end if

        if(this%L1D(1) > this%L1D(this%Nl)) then
            write(*,*) "L: ",this%L1D
            write(*,*) "reorder wave model so L shell is in ascending order"
            stop
        end if

        if(this%Energy1D(1) > this%Energy1D(this%Ne)) then
            write(*,*) "Ek: ",this%Energy1D
            write(*,*) "reorder wave model so Ek is in ascending order"
            stop
        end if

        
        allocate(this%wPS    (Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg))
        allocate(this%wHISS  (Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg))
        allocate(this%wCHORUS(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg))

        allocate(this%tauTotal(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, spc_ele%kStart:spc_ele%kEnd))

        
    end subroutine eWM_BW_LossInit


    subroutine eWM_BW_DoUpdate(this, Model, Grid, State)
        !! We are called at the start of every step, once coupling info, boundary conditions, etc. have beenr resolved
        !! We calculate all of the weighting, loss taus, for all electrons, and store it here
        !! When our calcTau gets called for an i,j,k, we'll just report the right value from our State
        class(raiLoss_eWM_BW_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State

        integer :: psphIdx, eleIdx
        logical :: isGood
        integer :: i,j,k
        real(rp) :: NpsphPnt
            !! Density [#/cc] of plasmasphere at point i,j
        real(rp) :: L, MLT, E, Kp
            !! L shell and MLT of given point, channel energy in MeV, current Kp
        real(rp) :: wLBlend, wNBlend
            !! L-weighting of blending between IMAG and PS. 0=PS
            !! Density-weighting between Chorus and Hiss
        real(rp) :: tauPS, tauHiss, tauCHORUS

        ! Zero everyone out in prep for new values
        call fillArray(this%wPS     , 0.0_rp)
        call fillArray(this%wHISS   , 0.0_rp)
        call fillArray(this%wCHORUS , 0.0_rp)
        call fillArray(this%tauTotal, HUGE  )

        call this%loss_SS%doUpdate(Model, Grid, State)  ! Just in case its doUpdate actually does something in the future

        ! Things that aren't i,j,k-dependent
        psphIdx = spcIdx(Grid, F_PSPH)
        eleIdx = spcIdx(Grid, F_HOTE)
        Kp = this%KpTS%evalAt(State%t)

        associate(sh=>Grid%shGrid, spc=>Grid%spc(eleIdx))

            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,k,isGood,L,MLT,E,NpsphPnt,wNBlend,wLBlend,tauPS,tauHiss,tauCHORUS)
            do j=sh%jsg,sh%jeg
                do i=sh%isg,sh%ieg
                    isGood = State%active(i,j)
                    if (.not. isGood) then
                        cycle
                    endif

                    L = sqrt(State%xyzMincc(i,j,1)**2 + State%xyzMincc(i,j,2)**2)  ! [Re]
                    MLT = atan2(State%xyzMincc(i,j,2),State%xyzMincc(i,j,1))/pi*12.D0+12.D0

                    NpsphPnt = State%Den(psphIdx)%data(i,j)  ! [#/cc]
                    
                    ! Calculate blending
                    wNBlend = dlog(NpsphPnt/this%NpsphLow) / dlog(this%NpsphHigh/this%NpsphLow)
                    call ClampValue(wNBlend, 0.0_rp, 1.0_rp)
                        !! 1 => Psphere Hiss, 0 => Other
                    wLBlend = RampDown(L, this%ChorusLMax, this%PsheetLMin - this%ChorusLMax)
                        !! 1 => Chorus, 0 => PS
                    
                    ! If psphere density is high enough, we always apply Hiss regardless of L
                        ! This means Hiss model needs to appropriately handle any L value its given
                    ! If psphere density is low, we choose between Chorus and plasma sheet strong-scattering based on L
                    ! Also note: the weighting is not k-dependent, but if its not too costly to do for each k then ¯\_(ツ)_/¯
                    
                    ! Calculate weights
                    this%wHISS(i,j) = wNBlend
                    this%wCHORUS(i,j) = (1 - wNBlend)*wLBlend
                    this%wPS(i,j) = (1 - wNBlend)*(1-wLBlend)

                    ! Now calculate loss taus and accumulate
                
                    do k=spc%kStart,spc%kEnd        
                        
                        E = abs(Grid%alamc(k) * State%bvol_cc(i,j)**(-2./3.)) * 1.0E-6  ! [MeV]

                        if (this%wHISS(i,j) > TINY) then
                            tauHiss = calcHissTau(MLT, L, E, Kp)
                        else
                            tauHiss = HUGE
                        endif

                        if (this%wPS(i,j) > TINY) then
                            tauPS = this%loss_SS%calcTau(Model, Grid, State, i, j, k)
                        else
                            tauPS = HUGE
                        endif

                        if (this%wCHORUS(i,j) > TINY) then
                            !! Implement chorus tau calculation here
                            continue
                        else
                            tauCHORUS = HUGE
                        endif

                        this%tauTotal(i,j,k) = this%wHISS(i,j)*tauHiss + this%wPS(i,j)*tauPS  + this%wCHORUS(i,j)*tauCHORUS 
                    enddo
                enddo
            enddo

        end associate

    end subroutine eWM_BW_DoUpdate


    function eWM_BW_CalcTau(this, Model, Grid, State, i, j, k) result(tau)
        !! We are not doing any heavy calculations since we did it in DoUpdate
        !! Just give back thr right value after some checks
        class(raiLoss_eWM_BW_T), intent(in) :: this
        type(raijuModel_T ), intent(in) :: Model
        type(raijuGrid_T  ), intent(in) :: Grid
        type(raijuState_T ), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau

        associate(spc => Grid%spc(Grid%k2spc(k)))

        tau = HUGE
        if (.not. this%reqsGood) then
            return
        endif

        if (k < spc%kStart .or. k > spc%kEnd) then
            return
        endif

        ! If still here then the requester is right to be calling us, give back our tau
        tau = this%tauTotal(i,j,k)

        end associate

    end function eWM_BW_CalcTau


    function eWM_BW_LossValidSpc(this, spc) result(isValid)
        class(raiLoss_eWM_BW_T), intent(in) :: this
        type(raijuSpecies_T), intent(in) :: spc
        logical :: isValid
        
        isValid = .false.
        if (.not. this%reqsGood) return

        if ( (spc%spcType .eq. RAIJUELE) ) then
            isValid = .true.
        endif
    end function eWM_BW_LossValidSpc


    subroutine eWM_BW_DoOutput(this, Model, Grid, State, gStr, doGhostsO)
        class(raiLoss_eWM_BW_T), intent(in) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: gStr
            !! Group ("Step#X") to append our data to in a sub-group ("eWM")
        logical, intent(in), optional :: doGhostsO


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
        associate(spc=>Grid%spc(spcIdx(Grid, F_HOTE)))
        call AddOutVar(IOVars, "wPS"    , this%wPS    (is:ie,js:je), dStr="0-1 weighting for plasma sheet strong scattering losses")
        call AddOutVar(IOVars, "wHISS"  , this%wHISS  (is:ie,js:je), dStr="0-1 weighting for plasmasphere hiss losses")
        call AddOutVar(IOVars, "wCHORUS", this%wCHORUS(is:ie,js:je), dStr="0-1 weighting for inner mag chorus wave losses")
        call AddOutVar(IOVars, "tauTotal", this%tauTotal(is:ie,js:je,spc%kStart:spc%kEnd), dStr="[s] combined lifetimes of all contributing loss processes")
        call WriteVars(IOVars,.true.,Model%raijuH5, gStrWM)
        end associate
    end subroutine eWM_BW_DoOutput


    function calcHissTau(MLTin, Lin, Ein, Kpin) result(tau)
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
        !rateK = 1.0/tau

        !write(*,*) Ein, L, fL, E, tau, rateK, '--'
        !stop

    end function calcHissTau


end module raijuLoss_eWM_BW