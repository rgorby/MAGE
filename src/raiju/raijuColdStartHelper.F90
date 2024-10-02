module raijuColdStartHelper

    use raijudefs
    use raijutypes
    use imaghelper
    use earthhelper

    use raijuetautils
    use raijuloss_CX

    implicit none

    contains

    subroutine raijuGeoColdStart(Model, Grid, State, t0, dstModel)
        !! Cold start RAIJU assuming we are at Earth sometime around 21st century
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        real(rp), intent(in) :: t0
            !! Target time to pull SW values from
        real(rp), intent(in) :: dstModel
            !! Current dst of global model

        integer :: sIdx_p, sIdx_e
        real(rp) :: dstReal, dstTarget
        real(rp) :: dps_preCX, dps_postCX, dps_ele

        sIdx_p = spcIdx(Grid, F_HOTP)
        sIdx_e = spcIdx(Grid, F_HOTE)

        ! Calc our target RC dst
        write(*,*) "WARNING: Setting QTRC from raijuColdStart, idk if we should be in charge of this"
        dstReal = GetSWVal('symh', Model%tsF, t0)
        dstTarget = dstReal - dstModel

        ! Start by nuking all etas, we will set it all up ourselves
        State%eta = 0.0

        ! Init psphere
        call setRaijuInitPsphere(Model, Grid, State, Model%psphInitKp)

        ! Init hot protons
        call raiColdStart_initHOTP(Model, Grid, State, t0, dstTarget)
        dps_preCX  = spcEta2DPS(Model, Grid, State, Grid%spc(sIdx_p), State%active .eq. RAIJUACTIVE)  ! Note: final argument is making a logical mask for isGood
        ! Hit it with some charge exchange
        call raiColdStart_applyCX(Model, Grid, State, Grid%spc(sIdx_p))
        dps_postCX = spcEta2DPS(Model, Grid, State, Grid%spc(sIdx_p), State%active .eq. RAIJUACTIVE)
        

        write(*,*)"Lazy raijuGeoColdStart: not rescaling proton eta to target dst, just adding electrons"

        ! Calc moments to update pressure and density
        call EvalMoments(Grid, State)

        ! Use HOTP moments to set electrons
        call raiColdStart_initHOTE(Model, Grid, State)
        dps_ele = spcEta2DPS(Model, Grid, State, Grid%spc(sIdx_e), State%active .eq. RAIJUACTIVE)

        write(*,'(a,f7.2)') "  Target DPS-Dst       : ",dstTarget
        write(*,'(a,f7.2)') "  Hot proton pre-loss  : ",dps_preCX
        write(*,'(a,f7.2)') "            post-loss  : ",dps_postCX
        write(*,'(a,f7.2)') "  Hot electron DPS-Dst :",dps_ele
        !write(*,*)" Hot electron DPS-Dst:",dps_ele

    end subroutine raijuGeoColdStart


    subroutine raiColdStart_initHOTP(Model, Grid, State, t0, dstTarget)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        real(rp), intent(in) :: t0
            !! Target time to pull SW values from
        real(rp), intent(in) :: dstTarget

        real(rp) :: dstTarget_p
        logical :: isInTM03
        integer :: i,j,sIdx
        integer, dimension(2) :: ij_TM
        real(rp) :: vSW, dSW, dPS_emp, pPS_emp, ktPS_emp
        real(rp) :: x0_TM, y0_TM, T0_TM, Bvol0_TM, P0_ps, N0_ps
        real(rp) :: L, vm, P_rc, D_rc, kt_rc, P_final, D_final

        P_final = 0.0
        D_final = 0.0
        sIdx = spcIdx(Grid, F_HOTP)

        associate(sh=>Grid%shGrid, spc=>Grid%spc(sIdx))

        ! Initialize TM03
        call InitTM03(Model%tsF,t0)

        ! Scale target Dst down to account for electrons contributing stuff later
        dstTarget_p = dstTarget / (1.0 + 1.0/Model%tiote)
        call SetQTRC(dstTarget_p) ! This sets a global QTRC_P0 inside earthhelper.F90

        ! Get Borovsky statistical values
        vSW = abs(GetSWVal("Vx",Model%tsF,t0))
        dSW =     GetSWVal("D" ,Model%tsF,t0)
        dPS_emp  = 0.292*(dSW**0.49)
        ktPS_emp = -3.65 + 0.0190*vSW*1.0e-3 !m/s=>km/s
        ktPS_emp = max(ktPS_emp,TINY)
        pPS_emp = DkT2P(dPS_emp, ktPS_emp)

        ! Get reference TM value at -10 Re
        x0_TM = -10.0-TINY
        y0_TM = 0.0
        ! Empirical temperature
        call EvalTM03([x0_TM,y0_TM,0.0_rp],N0_ps,P0_ps,isInTM03)
        T0_TM = DP2kT(N0_ps, P0_ps)
        ! Model FTV
        ij_TM = minloc( sqrt( (State%xyzMincc(:,:,XDIR)-x0_TM)**2 + (State%xyzMincc(:,:,YDIR)**2) ) )
        Bvol0_TM = State%bvol_cc(ij_TM(IDIR), ij_TM(JDIR))
        if (.not. isInTM03) then
            write(*,*) "This should not happen w/ TM03, you should figure this out ..."
            stop
        endif

        ! Now set our initial density and pressure profile
        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg
                
                if (State%active(i,j) .eq. RAIJUINACTIVE) cycle

                L = norm2(State%xyzMincc(i,j,XDIR:YDIR))
                vm = State%bvol_cc(i,j)**(-2./3.)

                kt_rc = T0_TM*(Bvol0_TM/State%bvol_cc(i,j))**(2./3.)
                kt_rc = min(kt_rc, 4.0*T0_TM)  ! Limit cap. Not a big fan, but without cap we get stuff that's too energetic and won't go away (until FLC maybe)

                call EvalTM03_SM(State%xyzMincc(i,j,:),N0_ps,P0_ps,isInTM03)

                P_rc = P_QTRC(L)  ! From earthhelper.F90

                if (.not. isInTM03) then
                    N0_ps = dPS_emp
                    P0_ps = pPS_emp
                endif

                if (P0_ps > P_rc) then
                    P_final = P0_ps
                    D_final = PkT2Den(P0_ps, kt_rc)
                else
                    P_final = P_rc
                    D_final = PkT2Den(P_rc, kt_rc)
                endif

                ! Finally map it to HOTP etas
                call DkT2SpcEta(Model, spc, State%eta(i,j,spc%kStart:spc%kEnd), D_final, kt_rc, vm)

            enddo
        enddo

        end associate

    end subroutine raiColdStart_initHOTP


    subroutine raiColdStart_applyCX(Model, Grid, State, spc)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        type(raijuSpecies_T), intent(in) :: spc

        integer :: i,j,k
        type(raiLoss_CX_T) :: lossCX
        type(XML_Input_T) :: nullXML  ! Needed for raiLoss inits, but CX doesn't actually need it, so make a dummy one
        real(rp) :: tCX = 12*3600  ! [s] Amount of time to apply CX for
        real(rp) :: tau

        call lossCX%doInit(Model, Grid, nullXML)
        
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,tau)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                do k = spc%kStart,spc%kEnd
                    tau = lossCX%calcTau(Model, Grid, State, i, j, k)
                    State%eta(i,j,k) = State%eta(i,j,k)*exp(-tCX/tau)
                enddo
            enddo
        enddo

    end subroutine raiColdStart_applyCX


    subroutine raiColdStart_initHOTE(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: sIdx_e, sIdx_p
        integer :: i,j
        real(rp) :: kt_p, kt_e, den, vm

        sIdx_p = spcIdx(Grid, F_HOTP)
        sIdx_e = spcIdx(Grid, F_HOTE)

        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,vm,den,kt_p,kt_e)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (State%active(i,j) .eq. RAIJUINACTIVE) cycle

                vm = State%bvol_cc(i,j)**(-2./3.)
                den = State%Den(i,j,sIdx_p+1)
                kt_p = DP2kT(den, State%Press(i,j,sIdx_p+1))
                kt_e = kt_p / Model%tiote
                call DkT2SpcEta(Model, Grid%spc(sIdx_e), &
                                State%eta(i,j,Grid%spc(sIdx_e)%kStart:Grid%spc(sIdx_e)%kEnd), &
                                den, kt_e, vm)
            enddo
        enddo

    end subroutine raiColdStart_initHOTE


    function GetSWVal(vID,fID,t) result(qSW)
        character(len=*), intent(in) :: vID,fID
        real(rp), intent(in) :: t
        real(rp) :: qSW

        type(TimeSeries_T) :: tsQ
        tsQ%wID = trim(fID)
        call tsQ%initTS(trim(vID),doLoudO=.false.)
        qSW = tsQ%evalAt(t)
    end function GetSWVal

end module raijuColdStartHelper