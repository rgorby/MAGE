module sifapp

    use kdefs
    use volttypes
    
    use siftypes
    use sifstarter


    implicit none

    contains

    subroutine initSifVolt(vApp, sifApp, iXML)
        type(voltApp_T), intent(in) :: vApp
        type(sifApp_T), intent(inout) :: sifApp
        type(XML_Input_T), intent(in) :: iXML

    ! Init core model
        call sifInitCore(sifApp, iXML)


    end subroutine initSifVolt


end module sifapp