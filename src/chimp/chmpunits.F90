!Unit information for CHIMP runs
!-------------------------------
!Input
!Grid data: L0 = Units of input corners [cm]
!Ie, for magnetospheric runs, L0 = Re = 6.38E+8 [cm]

!Input data: Velocity (V), B-Field (B), Time (T)
!Require scaling for input data to convert to cm/s, Gauss, and seconds
!Scaling for input: in2cms, in2G,in2s
!in2cms*V [cm/s], in2G*B [G], in2s*t [s]

!MHD Data unit conversion: Density (D), Pressure (P)
!inDScl*D = [#/cc]
!inPScl*P = [nPa]

!Output
!Output data: E-Field, B-Field, Time
!For now, output scaling is fixed to
!E-Field [mV/m], B-Field [nT], Time [s]

!Code units
!CGS w/ following normalization
!c = 1, Me (electron mass) = 1, L0 = 1, Qe (electron charge) = 1
module chmpunits
    use chmpdefs
    use cmidefs
    use xml_input
    use strings
    use planethelper

    implicit none

    !Globals for unit conversion
    real(rp) :: L0
    real(rp) :: inVScl,inBScl,inTScl !Converts input data to EB code units
    real(rp) :: inDScl=1.0,inPScl=1.0 !Units for MHD conversion
    real(rp) :: M0g = 0.0 !Magnetic moment [Gauss]
    real(rp) :: ebScl !CGS->Code EB units

    !Output units
    real(rp) :: oTScl !EB->s
    real(rp) :: oEScl !EB->mV/m
    real(rp) :: oBScl !EB->nT
    real(rp) :: oVScl !EB->km/s
    character(len=strLen) :: tStr="" !Time units string
    !Globals for background field data
    real(rp) :: MagM0=0 !Magnetic moment for dipoles

    !Input data to CGS conversions
    real(rp) :: in2cms,in2G,in2s

    !Radial value for field line endpoint to be considered closed
    real(rp) :: rClosed = 3.5

    !Constants
    real(rp), parameter, private :: V2mV = 1.0e+3 !V/m -> mV/m (duh)    

    contains

    !AMS 08052021
    !calculating scales in both functions identically
    !should maybe be moved to another subroutine so that any changes apply to both 

    !Set CHIMP units based on planet params
    subroutine setChimpUnitsVoltron(Model, planet, inpXML, uStrO)
        type(chmpModel_T), intent(inout) :: Model
        type(planet_T), intent(in) :: planet
        type(XML_Input_T), intent(in)    :: inpXML
        character(len=*), intent(in), optional :: uStrO
        
        real(rp) :: gv0, gB0, gP0

        !If not running in MAGE mode, assume planet doesn't exist and kick to setChimpUnits below
        if (.not. Model%isMAGE) then
            call setChimpUnits(Model,inpXML,uStrO)
        endif

        !Defaults
        call getGamPlanetNorms(planet, gv0=gv0, gB0=gB0, gP0=gP0)

        !Things we can get directly from planet params
        L0 = planet%rp_m*1.e2  ! m -> cm
        M0g = planet%magMoment  ! [Gauss]
        in2cms = 100*gv0 ! 100 km/s -> cm/s
        in2G   = gB0/G2nT
        in2s   = L0/in2cms
        inPScl = gP0

        !Special treatment for our solar system friends
        select case (trim(toUpper(planet%name)))
        case("EARTH")
            rClosed = 2.25  ! Correspond to EARTHCODE value below
        case("JUPITER")
            rClosed = 10.0
        case("SATURN")
            rClosed = 6.25
        case DEFAULT
            rClosed = 2.25
        end select

        !Now replace with user input if present
        call inpXML%Set_Val(rClosed,"domain/rClosed",rClosed)


        !Set main scaling values
        ebScl = (qe_cgs*L0/Me_cgs)/(vc_cgs**2.0)
        inTScl = in2s*vc_cgs/L0
        inVScl = in2cms/vc_cgs
        inBScl = in2G*ebScl

        MagM0 = -1.0*M0g*ebScl

        !Set output scaling values
        oTScl = (L0/vc_cgs)  !/(60*60.0) !ebtime->hrs
        tStr = "[Seconds]"

        oBScl = G2nT/ebScl !eb->nT
        oEScl = (G2T*vc_mks*V2mV)/ebScl !eb->mV/m
        oVScl = vc_cgs*1.0e-5 !eb (1/c) -> km/s
        write(*,*) '------------'
        write(*,*) 'CHIMP Units'
        write(*,*) 'inTScl = ', inTScl
        write(*,*) 'inBScl = ', inBScl
        write(*,*) 'inVScl = ', inVScl
        write(*,*) '------------'

    end subroutine setChimpUnitsVoltron

    !Use inpXML to set units unless uStrO is specified
    subroutine setChimpUnits(Model,inpXML,uStrO)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(in)    :: inpXML
        character(len=*), intent(in), optional :: uStrO

        real(rp) :: gv0, gB0, gP0

        if (present(uStrO)) then
            Model%uID = uStrO
        else
            call inpXML%Set_Val(Model%uID,"units/uid","Earth")
        endif

        !Defaults
        gv0 = defV0
        gB0 = defB0
        gP0 = defP0

        select case (trim(toUpper(Model%uID)))

        !Note, separating CODE units (during calculation) from output units (written to H5 files)
        case("EARTH")
            !Gamera units for Earth IO
            !Grid: Re
            !Velocity : 1 km/s
            !Field : 1 nT
            L0 = Re_cgs
            in2cms = gv0 ! 1 km/s -> 100 cm/s
            in2G   = 1.0/G2nT
            in2s   = 1.0
            M0g    = EarthM0g
            inPScl = 1.0 !Converted to nPa on output
            rClosed = 3.5
        case("EARTHCODE")
            !Gamera CODE units for Earth
            !Grid: Re
            !Velocity : 100 km/s
            !Field : 4.581 nT
            L0 = Re_cgs
            in2cms = 100*gv0 ! 100 km/s -> cm/s
            in2G   = gB0/G2nT
            in2s   = L0/in2cms
            M0g    = EarthM0g
            inPScl = gP0 !Gamera pressure -> nPa
            rClosed = 2.25
        case("JUPITER")
            !Gamera units for Jupiter
            L0 = RJupiterXE*Re_cgs
            in2cms = gv0 ! 1 km/s -> 100 cm/s
            in2G   = 1.0/G2nT
            in2s   = 1.0
            M0g    = JupiterM0g
            inPScl = 1.0 !Converted to nPa on output
            rClosed = 8.0
        case("JUPITERCODE")
            !Gamera units for Jupiter
            L0 = RJupiterXE*Re_cgs
            in2cms = 100*gv0 ! 100 km/s -> cm/s
            in2G   = gB0/G2nT
            in2s   = L0/in2cms
            M0g    = JupiterM0g
            inPScl = gP0 !Gamera pressure -> nPa
            rClosed = 10.0
        case("SATURN")
            L0 = RSaturnXE*Re_cgs
            in2cms = gv0 ! 1 km/s -> 100 cm/s
            in2G   = 1.0/G2nT
            in2s   = 1.0
            M0g = SaturnM0g
            inPScl = 1.0 !Converted to nPa on output
            rClosed = 6.25 !Inner boundary for Saturn
        case("SATURNCODE")
            L0 = RSaturnXE*Re_cgs
            in2cms = 100*gv0 ! 100 km/s -> cm/s
            in2G   = gB0/G2nT
            in2s   = L0/in2cms
            M0g = SaturnM0g
            inPScl = gP0 !Gamera pressure -> nPa
            rClosed = 6.25 !Inner boundary for Saturn
        case("HELIO")
            !Grid: Rs
            !Velocity : 150 km/s
            !Field: 100 nT
            !Gamera units for heliosphere runs
            L0 = 6.955e+10 !Rs in cm 
            in2cms = 1.0e-3/sqrt(4*PI*200*Mp_cgs) !150e+5 cm/s
            in2G   = 1.0e-3 !in [G]
            in2s   = L0/in2cms ! time in s 
            M0g = 0.0 
            inPScl = 1.0e-6*1.0e+8/4/pi  !Pressure  unit B[G]^2/4pi *1.e8 in [nPa]
            rClosed = 21.5 !Radius of inner boundary in units of grid length
        case("LFM")
            L0 = Re_cgs !Using scaled grid
            !Rest of units are already in cgs + Gauss
            in2cms = 1.0
            in2G = 1.0
            in2s = 1.0
            M0g = EarthM0g
            inPScl = 1.0 !Already converted from D/Cs -> #/cc & nPa
            rClosed = 3.5
        case("LFMJUPITER")
            L0 = 11.209*Re_cgs !Using scaled grid
            !Rest of units are already in cgs + Gauss
            in2cms = 1.0
            in2G = 1.0
            in2s = 1.0
            M0g = 4.28
            inPScl = 1.0 !Already converted from D/Cs -> #/cc & nPa
        case default
            write(*,*) '<Unknown units, using unscaled values!>'
            L0 = 1.0
            in2cms = 1.0
            in2G   = 1.0
            in2s   = 1.0
            M0g = 0.0

        end select

        !Now replace with user input if present
        call inpXML%Set_Val(rClosed,"domain/rClosed",rClosed)

        !Set main scaling values
        ebScl = (qe_cgs*L0/Me_cgs)/(vc_cgs**2.0)
        inTScl = in2s*vc_cgs/L0
        inVScl = in2cms/vc_cgs
        inBScl = in2G*ebScl

        MagM0 = -1.0*M0g*ebScl

        !Set output scaling values
        oTScl = (L0/vc_cgs)  !/(60*60.0) !ebtime->hrs
        tStr = "[Seconds]"

        select case (trim(toUpper(Model%uID)))
            case("HELIO")
                oTScl = (L0/vc_cgs)/in2s
                tStr = "[Dimensionless]"
                write(*,*) 'L0, in2cms, in2s', L0, in2cms, in2s
        end select


        oBScl = G2nT/ebScl !eb->nT
        oEScl = (G2T*vc_mks*V2mV)/ebScl !eb->mV/m
        oVScl = vc_cgs*1.0e-5 !eb (1/c) -> km/s
        write(*,*) '------------'
        write(*,*) 'CHIMP Units'
        write(*,*) 'inTScl = ', inTScl
        write(*,*) 'inBScl = ', inBScl
        write(*,*) 'inVScl = ', inVScl
        write(*,*) '------------'
        
    end subroutine setChimpUnits


    !Returns species data
    !Set charge (in |e|), mass (in Me)
    subroutine getSpecies(species,mass,charge)
        character(len=*), intent(in) :: species
        real(rp), intent(out) :: mass,charge

        real(rp) :: mScl

        mScl = Mp_cgs/Me_cgs

        !Set charge (in |e|), mass (in Me)
        select case(species)
            case ('H','H+','p','p+')
                !Hydrogen
                write(*,*) '<Using species H+>'
                charge = +1
                mass   = 1*mScl
            case ('O','Op','O+')
                !Oxygen
                write(*,*) '<Using species O+>'
                charge = +1
                mass   = 16*mScl
            case ('O6+','O6','o6','O6p')
                !Oxygen
                write(*,*) '<Using species O6+>'
                charge = +6
                mass   = 16*mScl
            case ('e','e-')
                !Electron
                write(*,*) '<Using species e->'
                charge  = -1
                mass    = 1
            case ('e+','ep')
                !Positron
                write(*,*) '<Using species e+>'
                charge  = +1
                mass    = 1
            case ('e++','epp')
                write(*,*) '<Using species e++>'
                !Doubly charged positron (nobody asked you)
                charge  = +2
                mass    = 1
            case ('epppp','e4p')
                write(*,*) '<Using species e++++>'
                !Quad charge positron
                charge  = +4
                mass    = 1
            case ('He+','Hep','hep')
                !He+
                write(*,*) '<Using species He+>'
                charge  = +1
                mass    = 4*mScl
            case ('He++','He','Hepp','hepp','He2')
                !He++
                write(*,*) '<Using species He++>'
                charge  = +2
                mass    = 4*mScl
            case default
                write(*,*) 'Unknown species/unspecified species, bailing ...'
                stop
        end select

    end subroutine getSpecies


end module chmpunits
