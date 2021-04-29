!Various routines to handle plasmapuase location in CHIMP
module plasmaputils
    use kdefs
    use chmpdefs
    use chmpunits
    use ebtypes
    use xml_input
    use streamline

    implicit none

    integer,  private :: Nl=80,Nphi=72
    real(rp), private :: Lin=2.0,Lout=6.0,phi0=0.0,phi1=2.0*PI ! ranges to calculate Lpp
    real(rp), private :: dL,dphi

    contains

    !Initialize the plasmapuase parameters
    subroutine initPP(ebState,inpXML,doStaticO)
        type(ebState_T), intent(inout)   :: ebState
        type(XML_Input_T), intent(inout) :: inpXML
        logical, optional, intent(in) :: doStaticO

        logical :: doStatic

        if (present(doStaticO)) then
            doStatic = doStaticO
        else
            doStatic = .false.
        endif

        !initializing the wave model
        call inpXML%Set_Val(Nl,'pp/Nl',Nl)
        call inpXML%Set_Val(Lin,'pp/Lin',Lin)
        call inpXML%Set_Val(Lout,'pp/Lout',Lout)
        call inpXML%Set_Val(Nphi,'pp/Nphi',Nphi)
        call inpXML%Set_Val(phi0,'pp/phi0',phi0)
        call inpXML%Set_Val(phi1,'pp/phi1',phi1)

        dL = (Lout-Lin)/Nl
        dPhi = (phi1-phi0)/Nphi

        !Allocate arrays that holds PP information
        !Always do eb1
        allocate(ebState%eb1%Lpp(Nphi+1))
        ebState%eb1%Lpp = 0.0

        if (.not. doStatic) then
            allocate(ebState%eb2%Lpp(Nphi+1))
            ebState%eb2%Lpp = 0.0
        end if

    end subroutine initPP

    ! Determine the location of the PP over desired MLT range
    subroutine calcLppMLT(Model,ebState,t,Lpp)
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T), intent(inout)   :: ebState 
        real(rp), intent(in) :: t
        real(rp), dimension(Nphi), intent(inout) :: Lpp

        real(rp) :: phi
        integer :: i

        ! loop through MLT and L in equatorial plane
        !$OMP PARALLEL DO default(none) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,phi) &
        !$OMP shared (Nphi,phi0,dPhi,Lpp,t,ebState,Model)
        do i=1, Nphi+1
            phi = phi0+(i-1)*dPhi
            call findLpp(Model,ebState,t,phi,Lpp(i))
        end do

    end subroutine calcLppMLT

    ! Find the location of the PP for a given MLT
    subroutine findLpp(Model,ebState,t,phi,Lpp)
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T), intent(inout)   :: ebState 
        real(rp), intent(in) :: t,phi
        real(rp), intent(inout) :: Lpp

        real(rp), dimension(2) :: bD !flux tube avg. density to calculate the gradient
        real(rp), dimension(NDIM) :: x 
        real(rp) :: L
        integer :: j
        real(rp) :: gradD, gradPP=-100 ! pp location criteria

        !Data for tracing
        type(ebTrc_T) :: ebTrc

        ! Move out in L in equatorial plane until PP criteria is met
        do j=1,Nl
            L = Lin+j*dL
            x = [L*cos(phi),L*sin(phi),0.0_rp] ! in equatorial plane
            call SliceFL(Model,ebState,x,t,ebTrc)
            bD(2) = ebTrc%bD
            if(j==1) bD(1) = bD(2) ! start off with no gradient
            
            ! gradD = (bD(2)-bD(1))/dL 
            ! if (gradD <= gradPP) then ! Check if pp location criteria is met
            if (bD(2) <= 50.0_rp .and. L > 3.0) then !using characteristic density for now
                Lpp = L
                exit
            else
                bD(1) = bD(2) !copy current state to previous state
            end if
        end do

    end subroutine findLpp

    ! calculate distance to plasmapuase from a point xeq
    subroutine deltaLpp(ebState,xeq,t,dLpp)
        type(ebState_T), intent(in)   :: ebState 
        real(rp), dimension(NDIM), intent(in) :: xeq
        real(rp), intent(in) :: t
        real(rp), intent(inout) :: dLpp

        real(rp) :: L,phi
        real(rp) :: dt,wp,wt
        integer  :: i0,i1
        real(rp) :: Lpp1,Lpp2,Lpp

        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        L = norm2(xeq)
        phi = atan2(xeq(2),xeq(1))
        if (phi <= 0) phi = phi + 2.0*PI ! to make 0-2pi range

        !get indices of MLT that point is between
        i0 = floor(phi/dPhi)+1
        i1 = ceiling(phi/dPhi)+1

        !FIX ME: Check if doStatic case is done correctly
        if (ebState%doStatic) then
            wp = (phi-i0*dPhi)/dPhi
            Lpp = eb1%Lpp(i0)+wp*(eb1%Lpp(i1)-eb1%Lpp(i0))
        else
            !linear interpolation to phi of location in each eb slice
            wp = (phi-i0*dPhi)/dPhi
            Lpp1 = eb1%Lpp(i0)+wp*(eb1%Lpp(i1)-eb1%Lpp(i0))
            Lpp2 = eb2%Lpp(i0)+wp*(eb2%Lpp(i1)-eb2%Lpp(i0))

            !linear interpolation in time
            dt = eb2%time-eb1%time
            wt = (t-eb1%time)/dt
            Lpp = Lpp1+wt*(Lpp2-Lpp1)
            ! write(*,*) "eb1%time, eb2%time ", eb1%time,eb2%time
            ! write(*,*) "eb1Lpp(i0),eb1Lpp(i1),eb2Lpp(i0),eb2Lpp(i1)",eb1%Lpp(i0),eb1%Lpp(i1),eb2%Lpp(i0),eb2%Lpp(i1)
            ! write(*,*) "Lpp1,Lpp2,Lpp",Lpp1,Lpp2,Lpp
        end if

        dLpp = L-Lpp

        end associate
    end subroutine deltaLpp

end module plasmaputils