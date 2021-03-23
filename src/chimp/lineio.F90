module lineio

    use chmpdefs
    use chmpunits
    use ebtypes
    use streamline
    use ioH5
    use files
    use ebtabutils
    implicit none

    character(len=strLen) :: flOutF
    integer, parameter :: MAXFLVS = 30

    contains

    !Initialize output
    subroutine initFLio(Model,fLines)
        type(chmpModel_T), intent(inout) :: Model
        type(fLine_T), intent(in), optional :: fLines(:)

        !Create filename
        write(flOutF,'(2a)') trim(adjustl(Model%RunID)),'.fl.h5'

        call CheckAndKill(flOutF)

    end subroutine initFLio
    
    !Write out block of streamlines
    subroutine writeLines(Model,ebGr,fLns,N)
        integer, intent(in) :: N
        type(chmpModel_T), intent(inout) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), dimension(N), intent(in) :: fLns

        type(IOVAR_T), dimension(MAXFLVS) :: IOVars

        character(len=strLen) :: gStr,lnStr
        integer :: i
        
        !Create group and write base data
        write(gStr,'(A,I0)') "Step#", Model%nOut
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        !call AddOutVar(IOVars,"MJD",MJDAt(ebState%ebTab,Model%t))
        
        !xxx, any others
        call WriteVars(IOVars,.true.,flOutF,gStr)
        call ClearIO(IOVars)

        !Now loop through and create subgroup for each line
        do i=1,N
            write(lnStr,'(A,I0)') "Line#", i-1
            call LineOut(Model,ebGr,fLns(i),flOutF,gStr,lnStr)
        enddo

    end subroutine writeLines
    
    !Write out individual line data
    subroutine LineOut(Model,ebGr,fL,fOut,gStr,lnStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: fL
        character(len=strLen), intent(in) :: fOut,gStr,lnStr

        type(IOVAR_T), dimension(MAXFLVS) :: IOVars
        integer :: i,Np,Nv
        real(rp) :: bS,bdV,OCb

        Np = fL%Nm + fL%Np + 1

        !Build output chain
        !write(*,*) 'Field line size ', Np
        !write(*,*) 'Bounds = ', fL%Nm,fL%Np

        !Test for degenerate chain (size = 1)
        if (Np<=1) then
            call ClearIO(IOVars)
            write(*,*) 'Degenerate field line, skipping ...'
            write(*,*) 'Seed = ', fL%xyz(0,:)
            return
        endif
        call AddOutVar(IOVars,"Np",Np)
        call AddOutVar(IOVars,"xyz",transpose(fL%xyz))
        call AddOutVar(IOVars,"n0",fL%Nm)

        !Add various information as attributes
        OCb = FLTop(Model,ebGr,fL)*1.0
        bS = FLEntropy(Model,ebGr,fL)
        bdV = FLVol(Model,ebGr,fL)

        call AddOutVar(IOVars,"FLTop",OCb)
        call AddOutVar(IOVars,"FLdV",bdV)
        call AddOutVar(IOVars,"FLEnt",bS)
        
        !Record seed point
        call AddOutVar(IOVars,"x0",fL%x0(XDIR))
        call AddOutVar(IOVars,"y0",fL%x0(YDIR))
        call AddOutVar(IOVars,"z0",fL%x0(ZDIR))

        if (Model%doMHD) then
            Nv = NumVFL
        else
            Nv = 0
        endif

        do i=0,Nv
            call AddOutVar(IOVars,fL%lnVars(i)%idStr,fL%lnVars(i)%V)
        enddo

        !Write output chain
        call WriteVars(IOVars,.true.,fOut,gStr,lnStr)
        call ClearIO(IOVars)

    end subroutine LineOut
end module lineio
