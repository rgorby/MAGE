module kaimagio
	use ioH5
	use kaimagdefs
	use kaimagtypes
	use clocks

	implicit none

	integer, parameter :: MAXIOVAR=30
	type(IOVAR_T), dimension(MAXIOVAR), private :: IOVars
	character(len=strLen), private :: h5File
	logical, private :: InitIO = .false.

	contains

	subroutine initKaimagIO(Model,Grid)
		type(kaimagModel_T), intent(in) :: Model 
		type(kaimagGrid_T), intent(in) :: Grid

		h5File = trim(Model%RunID) // "kaimag.h5"
		call CheckAndKill(h5File)

		call ClearIO(IOVars)
	    call AddOutVar(IOVARS,"XC",Grid%xc)
	    call AddOutVar(IOVARS,"YC",Grid%yc)
	    call AddOutVar(IOVARS,"ZC",Grid%zc)
	    call AddOutVar(IOVARS,"THETAC",Grid%thetac)
	    call AddOutVar(IOVARS,"PHIC",Grid%phic)
	    call WriteVars(IOVARS,.false.,h5File)

	end subroutine initkaimagIO

	subroutine writeKaimag(Model,Grid,State)
		type(kaimagModel_T), intent(in) :: Model 
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(in) :: State

		character(len=strLen) :: gStr 

		call Tic("IO")
		if (.not. initIO) then
			call initKaimagIO(Model,Grid)
			initIO = .true.
		endif


		write(gStr,'(A,I0)') "Step#",Model%IO%nOut
		write(*,*) 'Writing ', trim(gStr), State%time
		call ClearIO(IOVars)
		call AddOutVar(IOVARS,"time",State%time)
		call AddOutVar(IOVARS,"dt",Model%dt)
		call AddOutVar(IOVARS,"RHO",State%rho)
    	call AddOutVar(IOVARS,"UC",State%uc)
    	call AddOutVar(IOVARS,"VC",State%vc)
    	call WriteVars(IOVars,.true.,h5file,gStr)
    	call Toc("IO")

    end subroutine writeKaimag

    subroutine consoleOutputKaimag(Model,Grid,State)
		type(kaimagModel_T), intent(inout) :: Model 
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(in) :: State

		write(*,*) 'Console output step',Model%ts

		! Update the next timestep to output to console
		Model%IO%tsNext = Model%ts + Model%IO%tsOut

	end subroutine consoleOutputKaimag

	subroutine fOutputKaimag(Model,Grid,State)
		type(kaimagModel_T), intent(inout) :: Model 
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(in) :: State

		write(*,*) 'Inside fOutput',State%time

		call writeKaimag(Model,Grid,State)

		! Update the next time to output file
		Model%IO%tOut = Model%IO%tOut + Model%IO%dtOut
		Model%IO%nOut = Model%IO%nOut + 1

	end subroutine fOutputKaimag

	subroutine resOutputKaimag(Model,Grid,State)
		type(kaimagModel_T), intent(inout) :: Model 
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(in) :: State

		write(*,*) 'Inside resOutPut', State%time

		! Update the next time output a resatart file
		Model%IO%tRes = Model%IO%tRes + Model%IO%dtRes
		Model%IO%nRes = Model%IO%nRes + 1
	end subroutine resOutputKaimag


end module kaimagio