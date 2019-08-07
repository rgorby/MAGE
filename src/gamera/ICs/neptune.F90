!Setup for Neptune magnetosphere (no MI coupling)

module useric
    use types
    use gamutils
    use xml_input
    use bcs
    use msphutils
    use wind
    
    implicit none

    !Various global would go here
    real(rp) :: Rho0,P0
    real(rp) :: ThD,CosThD,SinThD

    ! type for remix BC
    type, extends(baseBC_T) :: CoroInnerBC_T

        !Main electric field structures
        real(rp), allocatable, dimension(:,:,:,:) :: inEijk,inExyz

        contains

        procedure :: doInit => InitCoroInner
        procedure :: doBC => CoroInner

    end type CoroInnerBC_T

    contains

    subroutine initUser(Model,Grid,State,inpXML)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML


    !Get defaults from input deck
        call inpXML%Set_Val(Rho0 ,"prob/Rho0",1.0_rp)
        call inpXML%Set_Val(P0   ,"prob/P0"  ,0.001_rp)
        call inpXML%Set_Val(ThD  ,"prob/ThD" ,46.0_rp)

        CosThD = cos(ThD*deg2rad)
        SinThD = sin(ThD*deg2rad)

        !Set magnetosphere parameters
        call setMagsphere(Model,inpXML)

    !Set BCs
        ! deallocate default BCs
        deallocate(Grid%ExternalBCs(INI )%p)
        deallocate(Grid%ExternalBCs(OUTI)%p)
        deallocate(Grid%ExternalBCs(INJ )%p)
        deallocate(Grid%ExternalBCs(OUTJ)%p)
        deallocate(Grid%ExternalBCs(INK )%p)
        deallocate(Grid%ExternalBCs(OUTK)%p)

        !Set BCs (spherical, RPT)
        allocate(CoroInnerBC_T      :: Grid%externalBCs(INI )%p)
        allocate(WindBC_T           :: Grid%externalBCs(OUTI)%p)
        allocate(lfmInBC_T          :: Grid%externalBCs(INJ )%p)
        allocate(lfmOutBC_T         :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(OUTK)%p)
    !Set ICs
        call VectorPot2Flux(Model,Grid,State,TiltDipoleVP)
        call GasIC2State(Model,Grid,State,GasIC)

    !Set bounds
        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

        !Set MG bounds
        Grid%isMG = Grid%is
        Grid%ieMG = Grid%ie
        Grid%jsMG = Grid%js
        Grid%jeMG = Grid%je
        Grid%ksMG = Grid%ks
        Grid%keMG = Grid%ke
    !Set hack functions
        if ( (Model%Ri == Model%NumRi) .or. (Model%Ri == 1) ) then
            Model%HackE => EFix
        endif

        contains

        subroutine TiltDipoleVP(x,y,z,Ax,Ay,Az)
            
            real(rp), intent(in) :: x,y,z
            real(rp), intent(out) :: Ax,Ay,Az

            real(rp), dimension(NDIM) :: A,m,r,rhat

            m = [0,0,1]
            m = [SinThD,0.0_rp,CosThD]

            r = [x,y,z]
            rhat = r/norm2(r)

            A = M0*cross(m,rhat)/(dot_product(r,r))
            Ax = A(XDIR)
            Ay = A(YDIR)
            Az = A(ZDIR)
        end subroutine TiltDipoleVP

        subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
            real(rp), intent(in) :: x,y,z
            real(rp), intent(out) :: D,Vx,Vy,Vz,P

            D = Rho0
            P = P0
            Vx = 0.0
            Vy = 0.0
            Vz = 0.0
        end subroutine GasIC

    end subroutine initUser

!---
    ! subroutine postBCInitUser(Model,Grid,State)
    !     type(Model_T), intent(inout) :: Model
    !     type(Grid_T), intent(inout) :: Grid
    !     type(State_T), intent(inout) :: State

    !     SELECT type(pWind=>Grid%externalBCs(OUTI)%p)
    !         TYPE IS (WindBC_T)
    !             if (associated(pWind%getWind)) then
    !                 write(*,*) 'Using solar wind BC from file ...'
    !             else
    !                 write(*,*) 'Using solar wind BC from subroutine ...'
    !                 pWind%getWind => SolarWindTS
    !             endif
    !         CLASS DEFAULT
    !             write(*,*) 'Could not find Wind BC in remix IC'
    !             stop
    !     END SELECT

    ! end subroutine postBCInitUser

!---
    !Fixes electric field before application
    subroutine EFix(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State

        write(*,*) 'Not doing anything yet ...'
    end subroutine EFix
!---
    !Initialization for Coro Inner BC
    subroutine InitCoroInner(bc,Model,Grid,State,xmlInp)
        class(CoroInnerBC_T), intent(inout) :: bc
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        integer :: PsiShells
        PsiShells = 5
        !Are we on the inner (REMIX) boundary
        if (Model%Ri == 1) then

            !Create holders for coupling electric field
            allocate(bc%inExyz(1:PsiShells,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM))
            allocate(bc%inEijk(1:PsiShells+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM))
            bc%inExyz = 0.0
            bc%inEijk = 0.0
        endif
    end subroutine InitCoroInner
!---
    !Inner-I BC for inner corotation boundary
    subroutine CoroInner(bc,Model,Grid,State)
        class(CoroInnerBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        write(*,*) 'Not doing anything yet ...'
    end subroutine CoroInner

end module useric