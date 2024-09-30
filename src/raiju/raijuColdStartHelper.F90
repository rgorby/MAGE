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

        type(raiLoss_CX_T) :: lossCX
        real(rp) :: dstReal, dstTarget

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
        ! CX RC
        ! Rescale RC to target dst
        ! EvalMoments, set electrons based using Maxwellian temperature

    end subroutine raijuGeoColdStart


    subroutine raiColdStart_initHOTP(Model, Grid, State, t0, dstTarget)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        real(rp), intent(in) :: t0
            !! Target time to pull SW values from
        real(rp), intent(in) :: dstTarget

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

        call SetQTRC(dstTarget) ! This sets a global QTRC_P0 inside earthhelper.F90

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

        contains

    end subroutine raiColdStart_initHOTP


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