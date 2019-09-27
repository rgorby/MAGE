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

module chmpunits
    use chmpdefs
    use xml_input
    use strings

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
    !Gamera magnetosphere defaults
    !These are true as long as density is in #/cc and V is in 100km/s

    real(rp), parameter, private :: gamB0 = 4.58103037171500
    real(rp), parameter, private :: gamP0 = 1.670000007587940E-002
    contains

    !Use inpXML to set units unless uStrO is specified
    subroutine setUnits(Model,inpXML,uStrO)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(in)    :: inpXML
        character(len=*), intent(in), optional :: uStrO

        if (present(uStrO)) then
            Model%uID = uStrO
        else
            call inpXML%Set_Val(Model%uID,"units/uid","Earth")
        endif

        select case (trim(toUpper(Model%uID)))

        !Note, separating CODE units (during calculation) from output units (written to H5 files)
        case("EARTH")
            !Gamera units for Earth IO
            !Grid: Re
            !Velocity : 1 km/s
            !Field : 1 nT
            L0 = Re_cgs
            in2cms = 1.0e+5 ! km/s -> cm/s
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
            in2cms = 100*1.0e+5 ! 100 km/s -> cm/s
            in2G   = gamB0/G2nT
            in2s   = L0/in2cms
            M0g    = EarthM0g
            inPScl = gamP0 !Gamera pressure -> nPa
            rClosed = 3.5
        case("JUPITER")
            !Gamera units for Jupiter
            L0 = RJupiterXE*Re_cgs
            in2cms = 100*1.0e+5 ! 100 km/s -> cm/s
            in2G   = gamB0/G2nT
            in2s   = L0/in2cms
            M0g    = JupiterM0g
            inPScl = gamP0 !Gamera pressure -> nPa
            rClosed = 10.0
        case("SATURN")
            L0 = RSaturnXE*Re_cgs
            in2cms = 100*1.0e+5 ! 100 km/s -> cm/s
            in2G   = gamB0/G2nT
            in2s   = L0/in2cms
            M0g = SaturnM0g
            inPScl = gamP0 !Gamera pressure -> nPa
            rClosed = 5.0 !Inner boundary for Saturn
        case("HELIO")
            !E: Add code units here
            !Gamera units for heliosphere runs
            L0 = 1.0
            in2cms = 1.0
            in2G   = 1.0
            in2s   = 1.0
            M0g = 0.0
            inPScl = 1.0  !Gamera pressure -> nPa
            rClosed = 5.0 !Radius of inner boundary in units of grid length
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
        
    end subroutine SetUnits

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
                write(*,*) 'Unknown species, using H'
                charge  = +1
                mass    = 1*mScl
        end select

    end subroutine getSpecies
    
end module chmpunits