module psdinit
    use chmpdefs
    use chmpunits
    use psdtypes
    use xml_input
    use particleio
    use psdutils
    use pdfuns
    
    use ioH5

    implicit none
    
    contains

    !Create phase space from XML info
    !Units: Radius (Re), phi (degrees), K (kev), alpha (degrees)

    subroutine genPS(Model,psGr,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(PSEq_T), intent(inout) :: psGr
        type(XML_Input_T), intent(inout) :: inpXML

        logical :: doLogs(NVARPS)
        character(len=strLen) :: sStr
        
        psGr%time = Model%t

        !Get phase space bounds
        !Set default cells/bounds (will convert degrees to radians later)
        psGr%Nr = 20
        psGr%Np = 8
        psGr%Nk = 15
        psGr%Na = 10
        !Default dimension bounds
        psGr%dimBds(PSRAD,:)   = [2.5,20.0]
        psGr%dimBds(PSPHI,:)   = [0.0,360.0]
        psGr%dimBds(PSKINE,:)  = [1.0,100.0]
        psGr%dimBds(PSALPHA,:) = [0.0,90.0]
        doLogs(PSRAD  ) = .true.
        doLogs(PSPHI  ) = .false.
        doLogs(PSKINE ) = .true.
        doLogs(PSALPHA) = .false.


        call getDim(Model,inpXML,"radius",psGr%Nr,psGr%rI,psGr%rC,psGr%rD,psGr%dimBds(PSRAD  ,1),psGr%dimBds(PSRAD  ,2),doLogs(PSRAD  ))
        call getDim(Model,inpXML,"phi"   ,psGr%Np,psGr%pI,psGr%pC,psGr%pD,psGr%dimBds(PSPHI  ,1),psGr%dimBds(PSPHI  ,2),doLogs(PSPHI  ))
        call getDim(Model,inpXML,"energy",psGr%Nk,psGr%kI,psGr%kC,psGr%kD,psGr%dimBds(PSKINE ,1),psGr%dimBds(PSKINE ,2),doLogs(PSKINE ))
        call getDim(Model,inpXML,"alpha" ,psGr%Na,psGr%aI,psGr%aC,psGr%aD,psGr%dimBds(PSALPHA,1),psGr%dimBds(PSALPHA,2),doLogs(PSALPHA))

        !Convert degrees to radians
        psGr%pI = (1/rad2deg)*psGr%pI
        psGr%aI = (1/rad2deg)*psGr%aI
        psGr%pC = (1/rad2deg)*psGr%pC
        psGr%aC = (1/rad2deg)*psGr%aC
        psGr%pD = (1/rad2deg)*psGr%pD
        psGr%aD = (1/rad2deg)*psGr%aD

        !Reset dimension bounds to grid values
        psGr%dimBds(PSRAD  ,:) = [minval(psGr%rI),maxval(psGr%rI)]
        psGr%dimBds(PSPHI  ,:) = [minval(psGr%pI),maxval(psGr%pI)]
        psGr%dimBds(PSKINE ,:) = [minval(psGr%kI),maxval(psGr%kI)]
        psGr%dimBds(PSALPHA,:) = [minval(psGr%aI),maxval(psGr%aI)]

        call inpXML%Set_Val(psGr%doShape,"sim/doShape",.true.)

        !Output grid data
        write(*,*) '--------------'
        write(*,'(a)') 'Phase space grid (N/Min/Max)'
        write(*,'(a,I5,a,f8.3,a,f8.3)') 'Radius: ', psGr%Nr,' / ',        psGr%dimBds(PSRAD  ,1),' / ',        psGr%dimBds(PSRAD  ,2)
        write(*,'(a,I5,a,f8.3,a,f8.3)') 'Phi   : ', psGr%Np,' / ',rad2deg*psGr%dimBds(PSPHI  ,1),' / ',rad2deg*psGr%dimBds(PSPHI  ,2)
        write(*,'(a,I5,a,f8.3,a,f8.3)') 'Energy: ', psGr%Nk,' / ',        psGr%dimBds(PSKINE ,1),' / ',        psGr%dimBds(PSKINE ,2)
        write(*,'(a,I5,a,f8.3,a,f8.3)') 'Alpha : ', psGr%Na,' / ',rad2deg*psGr%dimBds(PSALPHA,1),' / ',rad2deg*psGr%dimBds(PSALPHA,2)
        write(*,*) '--------------'
        
        !Create main arrays
        allocate(psGr%dVb(psGr%Nr,psGr%Np,psGr%Na))
        allocate(psGr%isClosed(psGr%Nr,psGr%Np))
        allocate(psGr%dG (psGr%Nr,psGr%Np,psGr%Nk,psGr%Na))
        psGr%dVb = 0.0
        psGr%dG = 0.0

        !Create flux tube holder
        allocate(psGr%bLns(psGr%Nr,psGr%Np))
        
        !Get species info
        call inpXML%Set_Val(sStr,"tps/species","X")
        call getSpecies(sStr,Model%m0,Model%q0)

        allocate(psGr%Qrp (psGr%Nr,psGr%Np,NVARMHD))
        allocate(psGr%kTeq(psGr%Nr,psGr%Np))
        allocate(psGr%Vreq(psGr%Nr,psGr%Np))

        if (.not. Model%doMHD) then
            !Not using MHD density/temperature, set defaults
            call inpXML%Set_Val(rho0,"fields/rho0",1.0_rp)
            call inpXML%Set_Val(kt0 ,"fields/kT0" ,1.0_rp)
        else
            rho0 = 0.0
            kT0 = 0.0
        endif

        psGr%Qrp(:,:,DEN) = rho0
        psGr%Qrp(:,:,VELX:VELZ) = 0.0
        psGr%QRP(:,:,PRESSURE) = 0.0
        psGr%kTeq = kT0
        psGr%Vreq = 0.0

    end subroutine genPS

    !Create population (collection of h5p files) for PSD calc.
    !Read from <population popid="XXX" ns="1" ne="4"/>

    subroutine genPop(Model,psPop,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(psdPop_T), intent(inout) :: psPop
        type(XML_Input_T), intent(inout) :: inpXML

        integer :: n,Np,Nt,NpTot
        real(rp) :: T0,m0,q0,dtStp
        logical :: fExist
        character(len=strLen) :: fIn

        !Get popid and range
        call inpXML%Set_Val(psPop%popid,"population/popid","chimp")
        call inpXML%Set_Val(psPop%ns,"population/ns",1)
        call inpXML%Set_Val(psPop%ne,"population/ne",4)

        !Get parameters for PSD IC
        call inpXML%Set_Val(psPop%kTScl,"population/kTScl",1.0_rp)
        call SetPSD0(Model,inpXML)

        NpTot = 0
        !Loop through files and check they exist
        do n=psPop%ns,psPop%ne
            write(fIn,'(a,a,I0.6,a)') trim(adjustl(psPop%popid)),'.',n,'.h5part'
            inquire(file=fIn,exist=fExist)
            if (.not. fExist) then
                write(*,*) 'Error opening file ', trim(fIn)
                stop
            endif
            call TPinfo(fIn,Np,Nt,dtStp,T0,m0,q0)
            NpTot = NpTot+Np
        enddo
        write(*,'(a,I9,a,I6,a)') 'Found ', NpTot, ' particles in ',psPop%ne-psPop%ns+1,' files.'

        psPop%NumTP = NpTot
        psPop%dtStp = dtStp
        psPop%T0 = T0
        allocate(psPop%TPs(NpTot,NVARPS))
        psPop%TPs = 0.0
        
        allocate(psPop%isWgt(NpTot))
        allocate(psPop%isIn (NpTot))
        psPop%isWgt = .false.
        psPop%isIn = .false.

        allocate(psPop%wgt(NpTot))
        allocate(psPop%tx (NpTot))
        psPop%wgt = 0.0 !Particle weights
        psPop%tx  = 0.0 !Particle marking time

        !For now just set Model%dt and Model%dtOut to dtStp
        write(*,*) "Setting dt/dtOut to dtStp ..."
        Model%dt    = dtStp
        Model%dtOut = dtStp

        if (Model%doStream) then
            !call inpXML%Set_Val(psPop%dTau,"stream/dTau",dtStp/inTScl)
            write(*,*) 'Setting tau to dtStp ...'
            psPop%dTau = dtStp
            call inpXML%Set_Val(psPop%dShell,"stream/dShell",1.0)
        endif
    end subroutine genPop

end module psdinit