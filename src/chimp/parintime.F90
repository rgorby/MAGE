module parintime
	use chmpdefs
	use xml_input
	use chmpunits

	implicit none
	integer, private :: NumB = 0
	logical, private :: doParInT = .false. !// in time

	contains

	subroutine InitParInTime(Model,inpXML,OutID,OutF)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: inpXML
        character(len=*) , intent(in) :: OutID
        character(len=strLen), intent(out) :: OutF

        real(rp) :: dtB,T0
    !Possible // in time
        call inpXML%Set_Val(NumB,'parintime/NumB',NumB)
        if (NumB > 1) then
            doParInT = .true.
            if ( (Model%Nblk>NumB) .or. (Model%Nblk<1) ) then
                write(*,*) "This block outside of acceptable bounds"
                write(*,*) "Block = ",Model%Nblk
                write(*,*) "Bounds = ",1,NumB
                write(*,*) "Bailing ..."
                stop
            endif
            !Reset time bounds
            T0 = Model%T0
            dtB = (Model%tFin-Model%T0)/NumB
            write(*,*) '------'
            write(*,*) 'Resetting T0/TFin = ',Model%T0*oTScl,Model%tFin*oTScl
            write(*,*) 'Using block ', Model%Nblk
            Model%T0 = (Model%Nblk-1)*dtB + T0
            Model%tFin = Model%T0 + dtB
            write(*,*) 'To        T0/TFin = ',Model%T0*oTScl,Model%tFin*oTScl
            if (Model%Nblk < NumB) then
                !Cut off a bit from TFin to avoid overlap w/ start of next
                Model%tFin = Model%tFin-0.01*dtB
            endif

            !Don't bother offsetting, let the concatenating script handle it
            write(OutF,'(a,a,I0.4,a,a,a)') trim(adjustl(Model%RunID)),'.',Model%Nblk,'.',trim(adjustl(OutID)),'.h5'
            write(*,*) '------'
            
        else
            doParInT = .false.
            NumB = 0
            write(OutF,'(4a)') trim(adjustl(Model%RunID)),'.',trim(adjustl(OutID)),'.h5'
        endif
        write(*,*) "Writing output to ", trim(OutF)
	end subroutine InitParInTime


end module parintime