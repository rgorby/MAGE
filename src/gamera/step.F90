!Routines to advance the grid

module step
    use gamtypes
    use gamutils
    use bcs
    use gdefs
    use output
    use multifluid

    implicit none

    contains

    !Calls BCs for all directions
    subroutine EnforceBCs(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr ! must be inout so that BC functions can be called
        type(State_T), intent(inout) :: State

        integer :: n
        character(len=strLen) :: BCID
        logical :: ownBC

        !Loop over BCs for this grid and call them
        do n=1,Gr%NumBC
            if (allocated(Gr%externalBCs(n)%p)) then
                SELECT CASE(Gr%externalBCs(n)%p%bcDir())
                    CASE(INI)
                        ownBC = Gr%hasLowerBC(IDIR)
                    CASE(OUTI)
                        ownBC = Gr%hasUpperBC(IDIR)
                    CASE(INJ)
                        ownBC = Gr%hasLowerBC(JDIR)
                    CASE(OUTJ)
                        ownBC = Gr%hasUpperBC(JDIR)
                    CASE(INK)
                        ownBC = Gr%hasLowerBC(KDIR)
                    CASE(OUTK)
                        ownBC = Gr%hasUpperBC(KDIR)
                    CASE DEFAULT
                        write (*,*) 'Warning, BC ignored for unowned boundary'
                        ownBC = .false.
                END SELECT
                if(ownBC) then
                    write (BCID, '(A,I0)') "BC#", n
                    call Tic(BCID)
                    call Gr%externalBCs(n)%p%doBC(Model,Gr,State)
                    call Toc(BCID)
                endif
            endif
        enddo

    end subroutine

    !Calculates timestep
    function CalcDT(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State
        real(rp) :: CalcDT

        real(rp) :: dtMin,dtOld,dtijk
        integer :: i,j,k
        integer :: is,ie,js,je,ks,ke
        logical :: isDisaster,isBad
        character(len=strLen) :: eStr

        if(Model%fixedTimestep) then
            CalcDT = Model%dt
            return
        endif

        dtOld = Model%dt
 
        dtMin = HUGE
       

        is = Gr%isDT
        ie = Gr%ieDT
        js = Gr%jsDT
        je = Gr%jeDT
        ks = Gr%ksDT
        ke = Gr%keDT

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(dtijk) reduction(min:dtMin)      
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    call CellDT(Model,Gr,State,i,j,k,dtijk)
                    dtMin = min(dtijk,dtMin)

                enddo
            enddo
        enddo

        if(dtMin == HUGE) then
            write(*,*) "<No timestep was smaller than HUGE, exiting ...>"
            isDisaster = .true.
        endif

        CalcDT = dtMin

        !Whether to crash or resort to desperate measures
        isDisaster = .false. 
        isBad      = .false.

        !Check for horrible things
        if (Model%ts > 0) then
            !Check for sudden drop in dt
            if (dtOld/dtMin >= 10) then
                write(eStr,*) "<Drop in timestep by 10x (",dtOld,"=>",dtMin,"), exiting ...>"
                isBad = .true.
            endif

            !Check for too small dt
            if (dtMin <= TINY) then
                write(eStr,*) "<Timestep too small (",dtMin,"), exiting ...>"
                isDisaster = .true. !We're boned
            endif

            !Check for slower but significant timestep drop
            if ( (Model%dt0>TINY) .and. (dtMin/Model%dt0 <= Model%limDT0) ) then
                write(eStr,*) "<Timestep less than limDT0 of initial (",Model%dt0,"=>",dtMin,"), exiting ...>"
                isBad = .true.
            endif

            if ( Model%doCPR .and. (Model%dt0>TINY) .and. (dtMin/Model%dt0 <= Model%limCPR) ) then
                !Patient is dying, attempt CPR to keep patient alive
                write(eStr,*) "<Patient is dying, trying to resuscitate. Clear!>"
                isBad = .true.
            endif

            !Done testing, now is the time for action!
            if (isDisaster .or. isBad) then

                if (Model%doCPR .and. isBad .and. (.not. isDisaster)) then
                !Doing CPR, it's bad but not a disaster. Try to save the patient
                    call StateCPR(Model,Gr,State)
                else
                !Meh, let it die
                    write(*,*) trim(eStr)
                    !Call black box and implode
                    call BlackBox(Model,Gr,State,dtOld)
                    write(*,*) 'Commiting suicide in 300s ...'
                    call sleep(300) !Sleep before blowing up
                    write(*,*) 'Goodbye cruel world'
                    stop !Self destruct

                endif
            endif !isDis or isBad

        endif !Terrible things check

        !Make sure we don't overstep the end of the simulation
        if ( (Model%t+CalcDT) > Model%tFin ) then
            CalcDT = max(Model%tFin-Model%t,TINY)
        endif

        !Apply a limit of 10x the initial timestep in case it starts static
        if (Model%dt0>TINY .and. CalcDT>(10.0*Model%dt0)) then
            CalcDT = 10.0 * Model%dt0
        endif

    end function CalcDT

    subroutine BlackBox(Model,Gr,State,dt0)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp), intent(in) :: dt0

        integer :: i,j,k,Nb
        integer :: is,ie,js,je,ks,ke
        logical :: isBad

        is = Gr%isDT
        ie = Gr%ieDT
        js = Gr%jsDT
        je = Gr%jeDT
        ks = Gr%ksDT
        ke = Gr%keDT

        Nb = 0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(isBad) reduction(+:Nb)  
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    call ChkCell(Model,Gr,State,i,j,k,dt0,isBad)
                    if (isBad) Nb = Nb+1
                enddo
            enddo
        enddo
        write(*,*) 'Bad cells: ', Nb
        write(*,*) 'Writing black box ...'
        !Do last output
        call WriteBlackBox(Model,Gr,State)

    end subroutine BlackBox
    
    !Check cell for bad things
    subroutine ChkCell(Model,Gr,State,i,j,k,dt0,isBad)
        integer, intent(in) :: i,j,k
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp),intent(in) :: dt0
        logical, intent(out) :: isBad

        real(rp) :: dtijk
        real(rp), dimension(NVAR) :: pW,pCon
        integer :: s,iG,jG,kG
        character(len=strLen) :: oStr
        isBad = .false.

        call CellDT(Model,Gr,State,i,j,k,dtijk)
        isBad = ( (dt0/dtijk)>10 ) .or. (dtijk < TINY)
        if (isBad) then
            !Display information
            !$OMP CRITICAL
            iG = i+Gr%ijkShift(IDIR)
            jG = j+Gr%ijkShift(JDIR)
            kG = k+Gr%ijkShift(KDIR)

            write(*,'(A,3I5)')     '<------- Bad Cell @ ijk = ',iG,jG,kG
            write(*,'(A,3es12.2)') 'xyz       = ', Gr%xyzcc(i,j,k,:)
            
            oStr = 'Bxyz [' // trim(Model%gamOut%bID) // '] = '
            write(*,'(A,3es12.2)') trim(oStr),State%Bxyz(i,j,k,:)*Model%gamOut%bScl

            do s=0,Model%nSpc
                !Get prim variables
                pCon = State%Gas(i,j,k,:,BLK)
                call CellC2P(Model,pCon,pW)

                if (s > 0) then
                    write(*,'(A,I0)') 'Fluid ', s
                else
                    write(*,'(A)') 'Bulk '
                endif
                !Den and pressure
                oStr = '   D/P [' // trim(Model%gamOut%dID) // ',' // trim(Model%gamOut%pID) // '] = '
                write(*,'(A,2es12.2)') trim(oStr),pW(DEN)*Model%gamOut%dScl,pW(PRESSURE)*Model%gamOut%pScl
                !Velocity
                oStr = '   Vxyz [' // trim(Model%gamOut%vID) // ']    = '
                write(*,'(A,3es12.2)') trim(oStr),pW(VELX:VELZ)*Model%gamOut%vScl
            enddo
            write(*,'(A)') '------->'
            !$OMP END CRITICAL

        endif

    end subroutine ChkCell

    !Try to do CPR on the plasma state, slow down fastest speeds
    subroutine StateCPR(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State

        integer :: i,j,k
        real(rp) :: dtijk

        if (Model%doMultiF .or. Model%doResistive) then
            write(*,*) 'CPR not yet implemented for these options, bailing ...'
            stop
        endif
        
        
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,dtijk)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie
                    call CellDT(Model,Gr,State,i,j,k,dtijk)
                    if (dtijk/Model%dt0 <= Model%limCPR) then
                        !Fix this cell
                        call CellCPR(Model,Gr,State,i,j,k)
                    endif
                    
                enddo
            enddo
        enddo

    end subroutine StateCPR

    !Perform CPR on a cell, try to nudge the timestep up a bit
    subroutine CellCPR(Model,Gr,State,i,j,k)
        integer, intent(in) :: i,j,k
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State
        

        real(rp), parameter :: alpha = 0.9 !Slow-down factor
        real(rp) :: D,Vf,Va,Cs,MagB
        real(rp), dimension(NVAR) :: pCon,pW
        real(rp), dimension(NDIM) :: B
        logical :: doVa,doCs,doVf,isFloored !Which speeds to slow
    
    !TODO: Remove this replicated code
    !Get three speeds, fluid/sound/alfven
        pCon = State%Gas(i,j,k,:,BLK)
        call CellC2P(Model,pCon,pW)
        Vf = norm2(pW(VELX:VELZ))

        call CellPress2Cs(Model,pCon,Cs)

        Va = 0.0
        D = pCon(DEN)

        if (Model%doMHD) then
            B = State%Bxyz(i,j,k,:)

            if (Model%doBackground) then
                B = B + Gr%B0(i,j,k,:)
            endif

            MagB = norm2(B)
            Va = MagB/sqrt(D)
            !Boris correct Alfven speed
            if (Model%doBoris) then
                Va = Model%Ca*Va/sqrt(Model%Ca**2.0 + Va**2.0)
            endif
        endif !doMHD

    !Identify which speeds are too fast
        doVa = .false.
        doCs = .false.
        doVf = .false.
        !Decide which speed is the problem
        if ( (Vf>=Cs) .and. (Vf>=Va) ) then
            !Fast flow speed
            doVf = .true.
        else if ( (Cs>=Va) .and. (Cs>=Va) ) then
            doCs = .true.
        else if ( (Va>=Cs) .and. (Va>=Vf) ) then
            !Fast Alfven speed
            doVa = .true.
        endif

        !For now disabling doing anything about Va
        doVa = .false.
        
    !Slow down speeds by a bit
        if (doVf) then
            !Directly slow speed
            pW(VELX:VELZ) = alpha*pW(VELX:VELZ)
        endif

        if (doCs) then
            !Reduce pressure
            pW(PRESSURE) = alpha*alpha*pW(PRESSURE)
        endif

        if (doVa) then
            !Nudge up density
            pW(DEN) = pW(DEN)/(alpha*alpha)
        endif
    !Check for floored cell
        isFloored = (pW(DEN)     <dFloor+TINY) .or. &
                  & (pW(PRESSURE)<pFloor+TINY)
        if (isFloored) then
            !This cell was just floored and it's causing us problems, so take away its velocity
            pW(VELX:VELZ) = 0.0
        endif          

    !Store changed values
        call CellP2C(Model,pW,pCon)
        State%Gas(i,j,k,:,BLK) = pCon
    end subroutine CellCPR

    subroutine CellDT(Model,Gr,State,i,j,k,dtijk)
        integer, intent(in) :: i,j,k
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp),intent(out) :: dtijk


        real(rp) ::  ke,e,P, dt
        real(rp) :: Vx,Vy,Vz,rho,Bx,By,Bz
        real(rp) :: dl,MagB,MagV,Vfl,Valf,Cs,vCFL,Diff,Vdiff

        !Get three speeds, fluid/sound/alfven
        rho = State%Gas(i,j,k,DEN,BLK)
        Vx  = State%Gas(i,j,k,MOMX,BLK)/rho
        Vy  = State%Gas(i,j,k,MOMY,BLK)/rho
        Vz  = State%Gas(i,j,k,MOMZ,BLK)/rho
        MagV = sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
        

        ke = 0.5*rho*(MagV**2.0)
        e = State%Gas(i,j,k,ENERGY,BLK) - ke
        P = (Model%gamma-1)*e

        !Handle multifluid case for sound speed/max flow
        if (Model%doMultiF) then
            Cs = MultiFCs(Model,State%Gas(i,j,k,:,:))
            MagV = MultiFSpeed(Model,State%Gas(i,j,k,:,:))
        else
            Cs = sqrt(Model%gamma*P/rho)
        endif

        Vfl   = MagV
        Valf  = 0.0
        VDiff = 0.0

        if (Model%doMHD) then
            Bx = State%Bxyz(i,j,k,XDIR)
            By = State%Bxyz(i,j,k,YDIR)
            Bz = State%Bxyz(i,j,k,ZDIR)
            if (Model%doBackground) then
                Bx = Bx + Gr%B0(i,j,k,XDIR)
                By = By + Gr%B0(i,j,k,YDIR)
                Bz = Bz + Gr%B0(i,j,k,ZDIR)
            endif
            !Use B to calculate Alfven speed for CFL speed
            MagB = sqrt(Bx**2.0+By**2.0+Bz**2.0)
            Valf = MagB/sqrt(rho)
            !Boris correct Alfven speed
            if (Model%doBoris) then
                Valf = Model%Ca*Valf/sqrt(Model%Ca*Model%Ca + Valf*Valf)
            endif
            
            if(Model%doResistive) then
               ! Asume t ~ x^2/(2Diff)
               Diff = maxval((/State%Deta(i,j,k,XDIR),State%Deta(i+1,j,k,XDIR), &
                               State%Deta(i,j,k,YDIR),State%Deta(i,j+1,k,YDIR), &
                               State%Deta(i,j,k,ZDIR),State%Deta(i,j,k+1,ZDIR)/))
               Vdiff = 2.0d0*Diff/minval((/Gr%di(i,j,k),Gr%dj(i,j,k),Gr%dk(i,j,k)/))
            end if
        endif

        vCFL = Vfl + sqrt(Cs**2.0 + Valf**2.0) + Vdiff
        
        !Use min length for timestep calculation
        dl = minval((/Gr%di(i,j,k),Gr%dj(i,j,k),Gr%dk(i,j,k)/))
        dtijk = Model%CFL*dl/vCFL

    end subroutine CellDT

end module step
