

module sifapp

	use xml_input
	use volttypes
	use sifstarter

	implicit none

	type, extends(innerMagBase_T) :: sifIMAG_T

		type(sifModel_T) :: Model
		type(sifGrid_T)  :: Grid
		type(sifState_T) :: State


        contains

        ! over-ride the base functions with SIF versions
        procedure :: doInit => initSIF
        !procedure :: doAdvance => advanceSIF
        !procedure :: doEval => EvalSIF
        !procedure :: doIO => doSIFIO
        !procedure :: doConIO => doSIFConIO
        !procedure :: doRestart => doSIFRestart

    end type

    contains


    subroutine initSIF(imag,iXML,isRestart,vApp)
    	class(sifIMAG_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        type(voltApp_T), intent(inout) :: vApp

        character(len=strLen) :: tmpStr
        !character(len=strLen) :: RunID
        !real(rp) :: t0

        ! Start init
        ! Model
        call sifInitModel(imag%Model, vApp%planet, iXML)
        ! Grid and config info
        call iXML%Set_Val(tmpStr, "grid/gType","UNISPH")

    end subroutine initSIF

end module sifapp