!Routines to handle RCM inner magnetosphere model
module rcmimag
    use volttypes
    use ioh5
    use files
    use earthhelper

    implicit none

    contains

    !Initialize RCM inner magnetosphere model
    subroutine initRCM(iXML)
        type(XML_Input_T), intent(in) :: iXML

    end subroutine initRCM

    !Evaluate eq map at a given point
    !Returns density (#/cc) and pressure (nPa)
    subroutine EvalRCM(r,phi,t,imW)
        real(rp), intent(in) :: r,phi,t
        real(rp), intent(out) :: imW(NVARIMAG)

        imW = 0.0

    end subroutine EvalRCM

end module rcmimag
