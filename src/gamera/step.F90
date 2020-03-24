!Routines to advance the grid

module step
    use gamtypes
    use gamutils
    use bcs
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
        !Loop over BCs for this grid and call them
        do n=1,Gr%NumBC
            if (allocated(Gr%externalBCs(n)%p)) then
                write (BCID, '(A,I0)') "BC#", n
                call Tic(BCID)
                call Gr%externalBCs(n)%p%doBC(Model,Gr,State)
                call Toc(BCID)
            endif
        enddo

    end subroutine

    !Calculates timestep
    function CalcDT(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp) :: CalcDT

        real(rp) :: dtMin,dtOld,dtijk
        integer :: i,j,k
        integer :: is,ie,js,je,ks,ke
        logical :: isDisaster

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

        CalcDT = dtMin
        isDisaster = .false.
        if (Model%ts > 0) then
        !Check for horrible things

            !Check for sudden drop in dt
            if (dtOld/dtMin >= 10) then
                write(*,*) "<Drop in timestep by 10x, exiting ...>"
                isDisaster = .true.
            endif

            !Check for too small dt
            if (dtMin <= TINY) then
                write(*,*) "<Timestep too small, exiting ...>"
                isDisaster = .true.
            endif

            !Check for slower but significant timestep drop
            if ( (Model%dt0>TINY) .and. (Model%dt0/dtMin >=100) ) then
                write(*,*) "<Timestep less than 1% of initial, exiting ...>"
                isDisaster = .true.
            endif

            !Call black box and implode
            if (isDisaster) then
                call BlackBox(Model,Gr,State,dtOld)
                write(*,*) 'Commiting suicide in 120s ...'
                call sleep(120) !Sleep for 2 minuts before blowing up
                write(*,*) 'Goodbye cruel world'
                stop
            endif !Self-destruct

        endif !Disaster check

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
        integer :: s
        isBad = .false.

        call CellDT(Model,Gr,State,i,j,k,dtijk)
        isBad = ( (dt0/dtijk)>10 ) .or. (dtijk < TINY)
        if (isBad) then
            !Display information
            !$OMP CRITICAL
            write(*,*) '<---------------->'
            write(*,*) 'ijk = ',i,j,k
            write(*,*) 'xyz = ', Gr%xyzcc(i,j,k,:)
            write(*,*) 'Bxyz = ', State%Bxyz(i,j,k,:)
            pCon = State%Gas(i,j,k,:,BLK)
            call CellC2P(Model,pCon,pW)
            write(*,'(A,5es12.2)') 'Bulk (PRIM) = ', pW
            if (Model%doMultiF) then
                do s=1,Model%nSpc
                    pCon = State%Gas(i,j,k,:,s)
                    call CellC2P(Model,pCon,pW)
                    write(* ,'(A,I0,A,5es12.2)') "Fluid" , s, ' = ', pW
                enddo
            endif
            write(*,*) '<---------------->'
            !$OMP END CRITICAL

        endif

    end subroutine ChkCell

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
