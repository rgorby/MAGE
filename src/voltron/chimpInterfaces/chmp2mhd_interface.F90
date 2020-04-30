! Converting chimp data to gamera data

module chmp2mhd_interface
    use gamtypes
    use volttypes
    use gamapp
    use ebtypes
    
    implicit none

    contains

    subroutine init_chmp2Mhd(chmp2mhd, ebTrcApp, gamApp)
        type(chmp2Mhd_T), intent(inout) :: chmp2mhd
        type(ebTrcApp_T), intent(inout) :: ebTrcApp
        type(gamApp_T)  , intent(in)    :: gamApp

        associate(Gr=>gamApp%Grid)
        !Changing squish to do nodes
        allocate(chmp2mhd%xyzSquish(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,2))
        allocate(chmp2mhd%isGood   (Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1))

        chmp2mhd%xyzSquish = 0.0
        chmp2mhd%isGood = .false.
        
        end associate

    end subroutine init_chmp2Mhd

end module chmp2mhd_interface
