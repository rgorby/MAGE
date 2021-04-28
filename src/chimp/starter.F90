module starter
    use chmpdefs
    use chmpunits
    use tptypes
    use ebtypes
    use strings
    use xml_input
    use ebinit
    use gridloc
    use chmpdbz
    use gridinterp
    use usertpic
    use userebic
    use kaiomp
    
    implicit none

!----------------------------------------
!Initialization routines
    contains

    !Initialize Model, set IC pointers for particles/fields and call them
    subroutine goApe(Model,ebState,tpState,iXML,optFilename)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T),   intent(inout), optional :: ebState
        type(tpState_T),   intent(inout), optional :: tpState
        type(XML_Input_T), intent(out),   optional :: iXML
        character(len=*) , intent(in)   , optional   :: optFilename
 
        character(len=strLen) :: xmlStr
        type(XML_Input_T) :: inpXML
        write(*,*) "goApe: Create all data structures"
 
        CalcGC => DirectGC
       
    !Create input XML object
        if (present(optFilename)) then
            xmlStr = optFilename
        else
            call getIDeckStr(xmlStr)
        endif
        inpXML = New_XML_Input(trim(xmlStr),"Chimp",.true.)

        !Send XML object back if optional arg present
        if (present(iXML)) iXML = inpXML
        
    !----------------------------
    !Get information for model data structure
        write(*,*) '----------------------------'
        write(*,*) 'Initializing model ...'
        call setUnits (Model,inpXML)
        call initModel(Model,inpXML)
        call setBackground(Model,inpXML) !Background field info

    !----------------------------
    !Set and call EB field IC
        if (present(ebState)) then
            write(*,*) '----------------------------'
            write(*,*) 'Initializing fields ...'
        
            call initFields(Model,ebState,inpXML)
            
            !Also do localization init
            call InitLoc(Model,ebState%ebGr,inpXML)
        endif
    !----------------------------
    !Set and call TP ICs
        if (present(tpState)) then
            write(*,*) '----------------------------'
            write(*,*) 'Initializing particles ...'

            call initParticles(Model,ebState,tpState,inpXML)
        endif
    end subroutine goApe

    !Initialize model variable
    !Assume time values are in input units, ie code = inTScl*input time
    subroutine initModel(Model,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: inpXML

        integer :: Narg,nSeed,dR0
        integer, dimension(:), allocatable :: vSeed
        character(len=strLen) :: buffer
        real(rp) :: T0Out

        !Set all default outputs off
        Model%doEBOut = .false.
        Model%doTPOut = .false.
        Model%doFLOut = .false.

    !Time block
        call inpXML%Set_Val(Model%T0  ,'time/T0',0.0)
        call inpXML%Set_Val(Model%tFin,'time/tFin',60.0)
        call inpXML%Set_Val(Model%dt  ,'time/dt',1.0)
        T0Out = Model%T0 !Before scaling
        !Scale times
        Model%T0   = inTScl*Model%T0  
        Model%tFin = inTScl*Model%tFin
        Model%dt   = inTScl*Model%dt

    !Interpolation options
        call setInterpolation(Model,inpXML)

    !Numerical knobs for pusher
        call inpXML%Set_Val(Model%epsht,'pusher/epsht',5.0e-2)
        call inpXML%Set_Val(Model%epsgc,'pusher/epsgc',5.0e-2)
        call inpXML%Set_Val(Model%do2D ,'pusher/do2D',.false.)
        
        !Get integrator type
        call setIntegrator(Model,inpXML)

    !TP injection over time (streaming)
        call inpXML%Set_Val(Model%doStream,'stream/doStream',.false.)
        
    !Output block
        call inpXML%Set_Val(Model%tsOut,'output/tsOut',10)
        call inpXML%Set_Val(Model%dtOut,'output/dtOut',10.0)
        call inpXML%Set_Val(T0Out,'output/T0Out',T0Out)
        call inpXML%Set_Val(Model%doTimer,'output/timer',.false.)
        call inpXML%Set_Val(Model%doEQProj,'output/doEQProj',.false.)
        call inpXML%Set_Val(Model%doTrc,'output/doTrc',.false.)
        call inpXML%Set_Val(Model%doSlim,'output/doSlim',.false.)
        call inpXML%Set_Val(Model%doFat ,'output/doFat' ,.false.)
        call inpXML%Set_Val(Model%doLL  ,'output/doLL' ,.false.)

        Model%dtOut = inTScl*Model%dtOut
        T0Out = inTScl*T0Out

    !Run info
        call inpXML%Set_Val(Model%RunID,'sim/runid',"Sim")

    !Guiding center iteration
        call inpXML%Set_Val(Model%MaxIter,'gciter/maxIter',20)
        call inpXML%Set_Val(Model%TolGC  ,'gciter/TolGC',1.0e-3)

    !Field info
        call inpXML%Set_Val(Model%doNumB0 ,'fields/doNumB0',.false.)
        call inpXML%Set_Val(Model%doEBFix ,'fields/doEBFix',.true.)
        call inpXML%Set_Val(Model%doPureB0,'fields/doPureB0',.false.)
        call inpXML%Set_Val(Model%doEBOut,'fields/ebOut',.false.)
        call inpXML%Set_Val(Model%doMHD,'fields/doMHD',.false.)
        call inpXML%Set_Val(Model%doJ,  'fields/doJ'  ,.false.)

    !Tracer
        call inpXML%Set_Val(Model%epsds,'tracer/epsds',1.0e-2)

    !WPI
        call inpXML%Set_Val(Model%doWPI,'wpi/doWPI',.false.)
        call inpXML%Set_Val(Model%doEQScat,'wpi/doEQScat',.false.)
        !Outer boundary of scattering
        call inpXML%Set_Val(Model%reqScat, 'wpi/reqScat' ,HUGE)

    !Plasmapause
        call inpXML%Set_Val(Model%doPP,'pp/doPP',.false.)
        !Basic setup
        Model%t    = Model%T0
        Model%tOut = T0Out
        Model%nOut = 0

        !Set random seed information, use block #
        Narg = command_argument_count()
        if (Narg < 2) then
            write(*,*) '<No seed offset specified, defaulting to 1>'
            dR0 = 1
        else
            call get_command_argument(2,buffer)
            read(buffer,*) dR0
            write(*,'(a,I0,a)') '<Using random seed offset ', dR0,'>'
        endif
        call random_seed(size=nSeed)
        allocate(vSeed(nSeed))
        vSeed(:) = rseed0+dR0
        call random_seed(put=vSeed)
        Model%Nblk = dR0
        
    !Threading
        call SetOMP(inpXML)
    end subroutine initModel
    
    !Set integrator type
    subroutine setIntegrator(Model,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: inpXML

        character(len=strLen) :: iStr

        Model%isDynamic = .false.
        call inpXML%Set_Val(iStr,'pusher/imeth',"FO")
        select case(trim(toUpper(iStr)))

            case("FO","FULLORB","LORENTZ")
                Model%imeth = IFO
            case("GC","GUIDINGCENTER")
                Model%imeth = IGC
            case("DYNAMIC","DYN")
                Model%imeth = IDYN
                Model%isDynamic = .true.
        end select
    end subroutine setIntegrator

    !Set interpolation type, ie map weight pointers
    subroutine setInterpolation(Model,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: inpXML

        character(len=strLen) :: iStr

        !For now always using TSC for derivatives
        Wgt1Dp => tsc1Dp

        call inpXML%Set_Val(iStr,'interp/wgt',"TSC")
        select case(trim(toUpper(iStr)))
            case("TSC")
                Wgt1D => tsc1D
            case("LIN","BILIN","LINEAR","MLIN")
                Wgt1D => lin1D
            case("QUAD","BIQUAD","PARA")
                Wgt1D => quad1D
        end select
    end subroutine setInterpolation

end module starter
