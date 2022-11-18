

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

        ! over-ride the base functions with RCM versions
        procedure :: doInit => initSIF
        procedure :: doAdvance => advanceSIF
        procedure :: doEval => EvalSIF
        procedure :: doIO => doSIFIO
        procedure :: doConIO => doSIFConIO
        procedure :: doRestart => doSIFRestart

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

        call iXML%Set_Val(tmpStr, "grid/gType","UNISPH")

        ! Start init
        if (tmpStr .eq. "SHGRID") then
        	call sifInitModel(imag%Model, imag%Grid, vApp%planet, iXML, vApp%shGrid)
        else
        	call sifInitModel(imag%Model, imag%Grid, vApp%planet, iXML)
        endif

    end subroutine initSIF

end module sifapp