module mixparams
  use mixdefs
  use mixtypes
  use xml_input

  implicit none

  contains
    subroutine initMIXParams(Params, optFilename)
      type(mixParams_T), intent(out) :: Params
      character(len=*), optional, intent(in) :: optFilename

      character(len=strLen) :: inpXML,tmpStr
      integer :: Narg
      logical :: fExist
      type(XML_Input_T) :: xmlInp
     
      if(present(optFilename)) then
         ! read from the prescribed file
         inpXML = optFilename
      else
         !Find input deck
         Narg = command_argument_count()
         if (Narg .eq. 0) then
            write(*,*) 'No input deck specified, defaulting to mixparams.xml'
            inpXML = "mixparams.xml"
         else
            call get_command_argument(1,inpXML)
         endif
      endif

      !First set input deck reader
      xmlInp = New_XML_Input(trim(inpXML),'REMIX',.true.)

        ! =========== CONDUCTANCE MODEL PARAMTERS =================== !
        ! EUV_MODEL_TYPE
        call xmlInp%Set_Val(tmpStr,"conductance/euv_model_type","AMIE")
        select case (tmpSTR)
           case ("AMIE")
              Params%euv_model_type = AMIE
           case ("MOEN_BREKKE")
              Params%euv_model_type = MOEN_BREKKE
           case default 
              stop "The EUV model type entered is not supported."              
        end select

        ! ET_MODEL_TYPE
        call xmlInp%Set_Val(tmpStr,"conductance/et_model_type","NONE")
        select case (tmpSTR)
           case ("NONE")
              Params%et_model_type = NONE
           case ("ADHOC")
              Params%et_model_type = ADHOC
           case ("CMIT")
              Params%et_model_type = CMIT
           case default 
              stop "The ET model type entered is not supported."              
        end select

        ! Numerical constants
        call xmlInp%Set_Val(Params%alpha,"conductance/alpha",1.0332467_rp)
        call xmlInp%Set_Val(Params%beta,"conductance/beta",0.4362323_rp)
        call xmlInp%Set_Val(Params%R,"conductance/R",0.083567956_rp)
        call xmlInp%Set_Val(Params%F107,"conductance/F107",120._rp)
        call xmlInp%Set_Val(Params%pedmin,"conductance/pedmin",2.0_rp)
        call xmlInp%Set_Val(Params%hallmin,"conductance/hallmin",1.0_rp)
        call xmlInp%Set_Val(Params%sigma_ratio,"conductance/sigma_ratio",3.0_rp)
        call xmlInp%Set_Val(Params%ped0,"conductance/ped0",10.0_rp)
        call xmlInp%Set_Val(Params%const_sigma,"conductance/const_sigma",.false.)
        call xmlInp%Set_Val(Params%doRamp,"conductance/doRamp",.true.)
        call xmlInp%Set_Val(Params%doChill,"conductance/doChill",.false.)
        call xmlInp%Set_Val(Params%doStarlight,"conductance/doStarlight",.false.)        
        call xmlInp%Set_Val(Params%doMR,"conductance/doMR",.false.)        
        call xmlInp%Set_Val(Params%apply_cap,"conductance/apply_cap",.true.)
        call xmlInp%Set_Val(Params%doGCM,"conductance/doGCM",.false.)
        call xmlInp%Set_Val(Params%doGCM2way,"conductance/doGCM2way",.false.)
        ! =========== CONDUCTANCE MODEL PARAMTERS =================== !

        ! =========== SOLVER PARAMTERS =================== !
        call xmlInp%Set_Val(Params%maxitr,"solver/maxitr",400)
        call xmlInp%Set_Val(Params%mr,"solver/mr",30)
        call xmlInp%Set_Val(Params%tol_abs,"solver/tol_abs",0.000001_rp)
        call xmlInp%Set_Val(Params%tol_rel,"solver/tol_rel",0.000001_rp)
        call xmlInp%Set_Val(Params%llbc_value,"solver/llbc_value",0.0_rp)
        ! =========== SOLVER PARAMTERS =================== !

        ! =========== GRID PARAMTERS =================== !
        call xmlInp%Set_Val(Params%Np,"grid/Np",256)
        call xmlInp%Set_Val(Params%Nt,"grid/Nt",128)
        call xmlInp%Set_Val(Params%LowLatBoundary,"grid/LowLatBoundary",45.0_rp)
        ! =========== GRID PARAMTERS =================== !

        ! =========== INIT PARAMTERS =================== !
        call xmlInp%Set_Val(Params%init_from_file,"grid/init_from_file",.false.)
        ! =========== INIT PARAMTERS =================== !

        ! =========== IO PARAMTERS =================== !
        call xmlInp%Set_Val(Params%dtOut,"output/dtOut",1.0_rp)
        ! =========== IO PARAMTERS =================== !

    end subroutine initMIXParams
end module mixparams
