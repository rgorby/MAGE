module kaimagmain
	use kaimagtypes
	use kaimaginit
	use kaimagsolver  
	use xml_input
	use kaimagprob
	use kaimagclaw

	implicit none

	type kaimagApp_T
		type(KaimagGrid_T) :: Grid
		type(kaimagModel_T) :: Model
		type(kaimagState_T) :: State, oState
	end type kaimagApp_T

	contains

	subroutine init_kaimag(kaimagApp,optFilename,updateV)
		type(kaimagApp_T), intent(inout) :: kaimagApp
		character(len=*), optional, intent(in) :: optFilename
		procedure(kaimagupdateV_T), pointer, intent(inout) :: updateV 

		character(len=strLen) :: inpXML
		character(len=strLen) :: vStr
		type(XML_Input_T) :: xmlInp 

		inpXML = optFilename
		write(*,*) 'Reading input deck from ', trim(inpXML)
		xmlInp = New_XML_Input(trim(inpXML),'Kaiju/kaimag',.true.)

		call initKaimagModel(kaimagApp%Model,xmlInp)
		call createUniformGrid(kaimagApp%Model,kaimagApp%Grid,xmlInp)

		call xmlInp%Set_Val(vStr,"sim/vType","BELL")
		call kaimagSetUpdateV_T(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State,&
								xmlInp,updateV,vStr)

		call setInitialConditions(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State,&
								  xmlInp,updateV)


	end subroutine init_kaimag

	subroutine stepKaimag(kaimagApp,updateV)
		type(kaimagApp_T), intent(inout) :: kaimagApp
		procedure(kaimagupdateV_T), pointer, intent(inout) :: updateV

		kaimagApp%Model%dt = CalcDT(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State)

		call UpdateStateData(kaimagApp,updateV)

		call EnforceBCs(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State)

	end subroutine stepKaimag

	subroutine UpdateStateData(kaimagApp,updateV)
		type(kaimagApp_T), intent(inout) :: kaimagApp
		procedure(kaimagupdateV_T), pointer, intent(inout) :: updateV
		kaimagApp%oState = kaimagApp%State

		call Advance2DSph(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State,kaimagApp%oState,&
						  updateV)

		!call AdvanceClaw(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State,kaimagApp%oState,&
		!				  updateV)		

		kaimagApp%Model%ts = kaimagApp%Model%ts + 1
		kaimagApp%Model%t = kaimagApp%Model%t + kaimagApp%Model%dt
		kaimagApp%State%time = kaimagApp%Model%t 

	end subroutine UpdateStateData


end module kaimagmain