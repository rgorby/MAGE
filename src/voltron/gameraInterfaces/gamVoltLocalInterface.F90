!Code for reading and writing gamera data from voltron on a single OpenMP node

module gamvoltlocalinterface
    use gamapp
    use voltapp

    implicit none

    type gamVoltLocalInterface_T
    end type gamVoltLocalInterface_T

    contains

    subroutine initializeVoltronFromGamera(gvLocal, gameraApp, voltronApp)
    end subroutine initializeVoltronFromGamer

    subroutine transferGameraToVoltron(gvLocal)
    end subroutine transferGameraToVoltron

    subroutine transferVoltronToGamera(gvLocal)
    end subroutine transferVoltronToGamera

end module gamvoltlocalinterface

