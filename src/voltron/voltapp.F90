! Collection of data and objects for the voltron middle man

module voltapp
    use mixapp
    use gamvoltlocalinterface

    implicit none

    type voltApp_T
        type(mixApp_T) :: remixApp

        type(gamVoltLocalInterface_T) :: gamVoltLocalInterface
    end type voltApp_T

    contains

end module voltapp

