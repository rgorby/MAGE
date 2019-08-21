!Routines related to multifluid simulations

module multifluid
    use gamtypes
    use xml_input
    use gamutils
    
    implicit none

    !Information for a given fluid species
    type FLUID_T
        character(len=strLen) :: idStr="NONE"
        real(rp) :: dVac = TINY
    end type FLUID_T

    !Array of fluid descriptors (1-nSpc)
    type(FLUID_T), dimension(:), allocatable :: Spcs

    contains

    subroutine InitMultiF(Model,xmlInp)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: xmlInp

        integer :: n
        real(rp) :: gDFloor !Global density floor
        character(len=strLen) :: sID,sTag

        !FIXME: Reading sim/dFloor here instead of reordering MPI & OMP versions
        call xmlInp%Set_Val(gDFloor,"sim/dFloor",dFloor)

        !Find number of species
        call xmlInp%Set_Val(Model%nSpc,'multifluid/nSpc',1)

        !Allocate number of descriptors
        allocate(Spcs(Model%nSpc))

        !Loop through species and get information
        do n=1,Model%nSpc
            write(sID ,'(A,I0)') "fluid" , n
            write(sTag,'(A,I0,A)') "fluid" , n,"/id"

            call xmlInp%Set_Val(Spcs(n)%idStr,sTag,sID)
            write(sTag,'(A,I0,A)') "fluid" , n,"/dVac"            
            call xmlInp%Set_Val(Spcs(n)%dVac,sTag,gDFloor)
        enddo

    end subroutine InitMultiF

    !Convert multi-species conserved variables into bulk flow quantities
    subroutine MultiF2Bulk(Model,U)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U(NVAR,BLK:Model%nSpc)

        integer :: n
        logical, dimension(Model%nSpc) :: isGood
        real(rp), dimension(NVAR) :: pW,pCon

        !Which fluids in this cell are good
        isGood = ( U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        
        do n=1,Model%nSpc
            if (.not. isGood(n)) then
                !Enforce global (not species specific floors)
                pCon = U(:,n)
                call CellC2P(Model,pCon,pW)
                pW(DEN)      = max(pW(DEN),dFloor)
                pW(PRESSURE) = max(pW(PRESSURE),pFloor)
                call CellP2C(Model,pW,pCon)
            endif
        enddo

        do n=1,NVAR
            U(n,BLK) = sum(U(n,1:Model%nSpc),mask=isGood)
        enddo
        
    end subroutine MultiF2Bulk
    
    !Convert full Gas state to bulk, ie include loop
    subroutine State2Bulk(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg
                    call MultiF2Bulk(Model,State%Gas(i,j,k,:,:))
                enddo
            enddo
        enddo

    end subroutine State2Bulk
    
    function MultiFCs(Model,U) result(MaxCs)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: U(NVAR,BLK:Model%nSpc)
        real(rp) :: MaxCs

        integer :: n
        real(rp) :: Csn
        logical, dimension(Model%nSpc) :: isGood
        !Which fluids in this cell are good
        isGood = ( U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        MaxCs = 0.0

        do n=1,Model%nSpc
            if (isGood(n)) then
                call CellPress2Cs(Model,U(:,n),Csn)
                MaxCs = max(MaxCs,Csn)
            endif
        enddo

    end function MultiFCs

end module multifluid
