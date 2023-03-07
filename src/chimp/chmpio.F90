module chmpio

    use clocks
    use xml_input
    use chmpdefs
    use chmpunits
    use tptypes
    use ebtypes
    use lineio
    use particleio
    use sliceio
    implicit none

    logical, private :: InitOut = .false. !Has IO been initialized

    contains
    
    subroutine cOutput(Model,ebState,tpState)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(tpState_T), intent(inout)   :: tpState
        
        real(rp) :: wTime,tpUps
        integer :: i,Nfo,Ngc,Ntot
        integer :: Npgc,Npfo
        type(prt_t) :: prt
        real(rp) :: Kevs(tpState%Np)

        associate( TPs=>tpState%TPs, NpT=>tpState%NpT,Np=>tpState%Np )
        !Calculate performance data
        wTime = readClock("Chimp")
        if (Model%ts >0) then
            !Count substeps
            Nfo = sum(TPs(:)%Nfo)
            Ngc = sum(TPs(:)%Ngc)
            Ntot = Nfo+Ngc
            !Reset substep counts
            TPs(:)%Nfo = 0
            TPs(:)%Ngc = 0

            !Particle updates per second
            tpUps = 1.0*Ntot/wTime
        else
            tpUps = 0.0
            Nfo = 0
            Ngc = 0
            Ntot = 0
        endif

        Npgc = count(TPs(:)%isGC .and. TPs(:)%isIn)
        Npfo = count( (.not. tpState%TPs(:)%isGC) .and. TPs(:)%isIn )

        !Calculate energies
        !$OMP PARALLEL DO
        do i=1,Np
            Kevs(i) = prt2kev(Model,TPs(i))
        enddo

        write(*,*) '----------------------------'
        write(*,'(a,f12.3,a,f12.3,a)')        'Sim Time   = ', Model%t*oTScl,' / ', Model%tFin*oTScl, ' ' // trim(tStr)
        write(*,'(a,I8)')                   '        ts = ', Model%ts

        if (Np>1) then
            
            write(*,'(a,I0,a,I0)')              '        TPs (Active/Total) = ', NpT,' / ',Np

            if (Model%doStream) then
                write(*,'(a,I0,a,I0)')              '        TPs (  Born/Total) = ', count(TPs(:)%isInit),' / ',Np
            endif
            if (Model%isDynamic) then
                write(*,'(a,I0,a,I0)')              '        TPs (    FO/GC   ) = ', Npfo,' / ',Npgc
            endif

            if (Ntot>0) then
                write(*,'(a,f12.3,a,f12.3,a)') '        Substeps (FO/GC) = ', Nfo*1.0e-6, 'M / ', Ngc*1.0e-6,'M'
            endif
            if (NpT>0) then
                write(*,'(a,f12.3)')     '     Avg K = ', sum(Kevs,mask=TPs(:)%isIn)/NpT
                write(*,'(a,f12.3)')     '     Max K = ', maxval(Kevs,mask=TPs(:)%isIn)!, ' @ Pid = ', maxloc(Kevs,mask=TPs(:))
                write(*,'(a,f12.3)')     '     Min K = ', minval(Kevs,mask=TPs(:)%isIn)!, ' @ Pid = ', minloc(Kevs,mask=TPs(:))
                write(*,'(a,f12.3)')     '     kPUps = ', tpUps*1.0e-3
            endif  
        else
            !Only 1 particle, do zoom output
            prt = TPs(1)
            write(*,'(a,i6)')         '     Particle ID              = ', prt%id
            write(*,'(a,3f8.3)')      '     Position                 = ', prt%Q(XPOS:ZPOS)
            write(*,'(a,f8.3)')       '     Particle Energy [keV]    = ', prt2kev(Model,prt)
            write(*,'(a,f8.3)')       '     Pitch Angle [o]          = ', prt%alpha*rad2deg
            write(*,'(a,2f8.3)')      '     Last Eq-XY               = ', prt%Qeq(EQX:EQY)
            write(*,'(a,f8.3)')       '     EQ-PA [o]                = ', prt%Qeq(EQALP)*rad2deg
            if (prt%isGC) then
                write(*,'(a,f8.3)')       '     p11/|P|                  = ', prt%Q(P11GC)/sqrt( prt%Q(GAMGC)**2.0 - 1 )
            endif
            write(*,'(a,f6.1,a,f6.1,a)') '        Substeps (FO/GC) = ', Nfo*1.0e-3, 'k / ', Ngc*1.0e-3,'k'
        endif
        write(*,*) '----------------------------'
        end associate
    end subroutine cOutput

    subroutine fOutput(Model,ebState,tpState,fLines)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in), optional   :: ebState
        type(tpState_T), intent(inout), optional   :: tpState
        type(fLine_T),   intent(in), optional   :: fLines(:)

        character(len=strLen) :: gStr
        integer :: Nfl

        if (.not. InitOut) then
            !if (present(ebState) .and. Model%doEBOut) call initEBio(Model,ebState)
            if (present(tpState)) call initPio (Model,tpState)
            if (present(fLines )) call initFLio(Model,fLines)

            !call InitOutput(Model,ebState,tpState)
            InitOut = .true.
            write(*,*) ''
        endif

        write(gStr,'(A,I0)') "Step#", Model%IO%nOut

        !write(*,'(a,f8.3,a)') '<Writing HDF5 output @ t = ', Model%t, ' >'
        
        if (Model%doEBOut .and. present(ebState)) then
            call Tic("ebOut")
            call writeEB(Model,ebState,gStr)
            call Toc("ebOut")
        endif
        
        if ( Model%doTPOut .and. present(tpState) .and. present(ebState) ) then
            call Tic("tpOut")
            call writeTP(Model,ebState,tpState,gStr)

            if (Model%doWPI) then
                tpState%TPs(:)%dAwpi = 0.0 !reset wpi diagnotsitcs every output
                tpState%TPs(:)%dKwpi = 0.0
            endif
            call Toc("tpOut")
        endif

        if (Model%doFLOut .and. present(fLines) .and. present(ebState) ) then
            call Tic("flOut")
            Nfl = size(fLines)
            call writeLines(Model,ebState%ebGr,fLines,Nfl)
            call Toc("flOut")
        endif

        !Setup for next output
        Model%tOut = Model%tOut + Model%dtOut
        Model%nOut = Model%nOut + 1

        Model%IO%tOut = Model%IO%tOut + Model%IO%dtOut
        Model%IO%nOut = Model%IO%nOut + 1

    end subroutine fOutput

    subroutine resOutput(Model,ebState,tpState)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(tpState_T), intent(inout)   :: tpState

        character(len=strLen) :: ResF,oStr,lnResF !Name of restart file
        logical :: fExist

        write (ResF, '(A,A,I0,A,I0.5,A)') trim(Model%RunID), "." , Model%Nblk, ".Res.", Model%IO%nRes, ".h5"

        call CheckAndKill(ResF)

        if (Model%t>1.0e-2) then
            write(oStr,'(f9.3,a,a)' ) Model%t*oTScl, ' ', trim(tStr)
        else
            write(oStr,'(es9.2,a,a)') Model%t*oTScl, ' ', trim(tStr)
        endif
        write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing HDF5 RESTART @ t = ', trim(oStr), ' >', ANSIRESET


        call writeTPres(Model,ebState,tpState,ResF)

        !Setup for next restart
        Model%IO%tRes = Model%IO%tRes + Model%IO%dtRes
        Model%IO%nRes = Model%IO%nRes + 1

        write (lnResF, '(A,A,I0,A,A,A)') trim(Model%RunID), "." , Model%Nblk, ".Res.", "XXXXX", ".h5"

        call MapSymLink(ResF,lnResF)

    end subroutine resOutput
end module chmpio
