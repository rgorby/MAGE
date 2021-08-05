
module planethelper
    use kdefs
    use gamtypes
    use helpertypes
    use strings, only : toUpper
    implicit none
    
    contains

    subroutine getPlanetParams(planet, xmlInp, pStrO)
        type(planet_T), intent(inout) :: planet
        type(XML_Input_T), intent(in) :: xmlInp
        character(len=*), intent(in), optional :: pStrO

        character(len=strLen) :: pID !Planet ID string
        logical :: doCorot
        real(rp) :: RIon
        

        if (present(pStrO)) then
            pID = pStrO
        else
            call xmlInp%Set_Val(pID,"/Kaiju/Gamera/prob/planet","Earth")
        endif

        select case (trim(toUpper(pID)))

        case("Earth","earth","EARTH")
            planet%rp_m = REarth
            planet%ri_m = RionE*1.e6  ! Defined in kdefs in 10000km
            planet%grav = 9.807
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", EarthM0g)
            planet%psiCorot = EarthPsi0
            planet%doGrav = .true.
        case("Saturn","saturn","SATURN")
            planet%rp_m = RSaturnXE*REarth
            planet%ri_m = 1.01*RSaturnXE*REarth
            planet%grav = 10.43
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", SaturnM0g)
            planet%psiCorot = 137.83*92 !kV
            planet%doGrav = .true.
        case("Jupiter", "jupiter", "JUPITER")
            planet%rp_m = RJupiterXE*REarth
            planet%ri_m = 1.01*RJupiterXE*REarth
            planet%grav = 24.79
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", JupiterM0g)
            planet%psiCorot = -2.5*1702.9*92.0 !kV
            planet%doGrav = .true.
        case("Mercury","mercury","MERCURY")
            planet%rp_m = RMercuryXE*REarth
            planet%ri_m = 1.05*RMercuryXE*REarth
            planet%grav = 0.0
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", MercuryM0g)
            planet%psiCorot = 0.0
            planet%doGrav = .false.
         case("Neptune","NEPTUNE")
            planet%rp_m = RNeptuneXE*REarth
            planet%ri_m = 1.01*RNeptuneXE*REarth
            planet%grav = 11.15
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", NeptuneM0g)
            planet%psiCorot = 10.024*92.0
            planet%doGrav = .false.
        case("Other","other","OTHER") ! Defaults to Earth values
            call xmlInp%Set_Val(planet%rp_m,"/Kaiju/Gamera/prob/x0",REarth)          ! [m]
            call xmlInp%Set_Val(Rion,"/Kaiju/Gamera/prob/Rion",RionE*1.e6/REarth)    ! [Rp]
            planet%ri_m =RIon*planet%rp_m                             ! [m]
            call xmlInp%Set_Val(planet%grav,"/Kaiju/Gamera/prob/G0",9.807)           ! [m/s2] 
            call xmlInp%Set_Val(planet%magMoment,"/Kaiju/Gamera/prob/M0",EarthM0g)   ! [gauss]
            call xmlInp%Set_Val(planet%psiCorot,"/Kaiju/Gamera/prob/Psi0",EarthPsi0) ! [kV]
            call xmlInp%Set_Val(planet%doGrav,"/Kaiju/Gamera/prob/doGrav",.true.)
        end select

        call xmlInp%Set_Val(doCorot,"/Kaiju/Gamera/prob/doCorot",.true.)
        if (.not. doCorot) then
            !Zero out corotation potential
            planet%psiCorot = 0.0
        endif

    end subroutine



end module planethelper
