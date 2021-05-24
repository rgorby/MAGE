module mixparams
  use mixdefs
  use mixtypes
  use xml_input

  implicit none

  !Grid default parameter type
  type MixGrid0_T
    integer  :: Np=256,Nt=128
    real(rp) :: LowLatBC = 45.0
  end type MixGrid0_T

  type(MixGrid0_T), private :: MixGrid0

  contains

    !Set default mix grid parameters based on inner MHD boundary and number of k cells
    subroutine SetMixGrid0(Rin,Nk)
      real(rp), intent(in) :: Rin
      integer , intent(in) :: Nk

      real(rp) :: dDeg

      MixGrid0%LowLatBC = asin(sqrt(1.0/Rin))*180.0/PI
      !Set default grid resolution (in degrees)
      if (Nk<=128) then
        !QUAD
        dDeg = 1.0
      else if (Nk<=256) then
        !OCT
        dDeg = 0.5
      else
        !HEX or above
        dDeg = 0.25
      endif

      MixGrid0%Np = nint(360.0/dDeg)
      MixGrid0%Nt = nint(MixGrid0%LowLatBC/dDeg)

    end subroutine SetMixgrid0

    subroutine initMIXParams(Params, optFilename)
      type(mixParams_T), intent(out) :: Params
      character(len=*), optional, intent(in) :: optFilename

      character(len=strLen) :: inpXML,tmpStr
      integer :: Narg
      logical :: fExist,doSuicide
      type(XML_Input_T) :: xmlInp
      
      doSuicide = .false.

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

        ! AURORA_MODEL_TYPE
        call xmlInp%Set_Val(tmpStr,"precipitation/aurora_model_type","FEDDER")
        select case (tmpSTR)
           case ("FEDDER")
              Params%aurora_model_type = FEDDER
           case ("ZHANG")
              Params%aurora_model_type = ZHANG
           case ("RCMONO")
              Params%aurora_model_type = RCMONO
           case default 
              stop "The aurora model type entered is not supported (Available options: FEDDER, ZHANG, RCMONO)."
        end select

        ! Numerical constants
        call xmlInp%Set_Val(Params%alpha,"precipitation/alpha",1.0332467_rp)
        call xmlInp%Set_Val(Params%beta,"precipitation/beta",0.4362323_rp)
        call xmlInp%Set_Val(Params%R,"precipitation/R",0.083567956_rp)
        call xmlInp%Set_Val(Params%doAuroralSmooth,"precipitation/doAuroralSmooth",.false.)        
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
        ! =========== CONDUCTANCE MODEL PARAMTERS =================== !

        ! =========== SOLVER PARAMTERS =================== !
        call xmlInp%Set_Val(Params%maxitr,"solver/maxitr",400)
        call xmlInp%Set_Val(Params%mr,"solver/mr",30)
        call xmlInp%Set_Val(Params%tol_abs,"solver/tol_abs",0.000001_rp)
        call xmlInp%Set_Val(Params%tol_rel,"solver/tol_rel",0.000001_rp)
        call xmlInp%Set_Val(Params%llbc_value,"solver/llbc_value",0.0_rp)
        ! =========== SOLVER PARAMTERS =================== !

        ! =========== GRID PARAMTERS =================== !
        call xmlInp%Set_Val(Params%Np,"grid/Np",MixGrid0%Np)
        call xmlInp%Set_Val(Params%Nt,"grid/Nt",MixGrid0%Nt)
        call xmlInp%Set_Val(Params%LowLatBoundary,"grid/LowLatBoundary",MixGrid0%LowLatBC)

        ! =========== GRID PARAMTERS =================== !

        ! =========== INIT PARAMTERS =================== !
        call xmlInp%Set_Val(Params%init_from_file,"grid/init_from_file",.false.)
        ! =========== INIT PARAMTERS =================== !

        ! =========== IO PARAMTERS =================== !
        call xmlInp%Set_Val(Params%dtOut,"output/dtOut",1.0_rp)
        call xmlInp%Set_Val(Params%nRes,"/gamera/restart/nRes",-1)
        ! =========== IO PARAMTERS =================== !

        ! =========== DEBUG PARAMETERS ================ !
        call xmlInp%Set_Val(Params%mklmsglvl,"debug/mklmsglvl", 0)
        ! =========== DEBUG PARAMETERS   ============== !

        ! Check for old-style precipitation quantities
        write(*,*) ANSIRED
        if (xmlInp%Exists("conductance/alpha")) then
            write(*,*) "conductance/alpha is gone. You're entering a world of pain. Use precipitation/alpha"
            doSuicide=.true.
        else if (xmlInp%Exists("conductance/beta")) then
            write(*,*) "conductance/beta is gone. You're entering a world of pain. Use precipitation/beta"
            doSuicide=.true.
        else if (xmlInp%Exists("conductance/R")) then
            write(*,*) "conductance/R is gone. You're entering a world of pain. Use precipitation/R"
            doSuicide=.true.
        else if (xmlInp%Exists("conductance/aurora_model_type")) then
            write(*,*) "conductance/aurora_model_type is gone. You're entering a world of pain. Use precipitation/aurora_model_type"
            doSuicide=.true.
        end if
        write(*,*) ANSIRESET, ''
        if (doSuicide) stop

    end subroutine initMIXParams
end module mixparams
