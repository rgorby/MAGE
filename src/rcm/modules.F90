!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CONSTANTS
  USE rcm_precision
  USE kdefs, ONLY : EarthPsi0,Re_cgs,Me_cgs,Mp_cgs,Mu0,Kbltz,eCharge
  use rcmdefs, ONLY : isize,jsize,jwrap
  REAL(rprec),PARAMETER :: radius_earth_m = Re_cgs*1.0e-2 ! Earth's radius in meters
  REAL(rprec),PARAMETER :: boltz = Kbltz*1.0e-7
  REAl(rprec),PARAMETER :: mass_proton =Mp_cgs*1.0e-3
  REAL(rprec),PARAMETER :: mass_electron=Me_cgs*1.0e-3
  REAL(rprec),PARAMETER :: ev=eCharge
  REAL(rprec),PARAMETER :: gamma=5.0/3.0
  REAL(rprec),PARAMETER :: big_vm = -1.0e5
  REAL(rprec),PARAMETER :: nt = 1.0e-9
  !REAL(rprec),PARAMETER :: tiote = 7.8
  REAL(rprec),PARAMETER :: tiote = 4.0 !Changed by K 5/9/20
  REAL(rprec),PARAMETER :: RCMCorot = EarthPsi0*1.0e+3 ! Convert corotation to V

END MODULE CONSTANTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE rice_housekeeping_module
  USE kdefs, ONLY : strLen
  USE rcm_precision, only : iprec,rprec
  use xml_input
  use kronos
  use strings
  use earthhelper, ONLY : SetKp0
  use rcmdefs, ONLY : DenPP0, ELOSS_FDG, ELOSS_SS, ELOSS_C05, ELOSS_C19, ELOSS_DW
  
  IMPLICIT NONE
  
  LOGICAL :: L_write_rcmu          = .false., &
             L_write_rcmu_torcm    = .false., &
             L_write_tracing_debug = .false., &
             L_write_vars_debug    = .false.
  
  INTEGER(iprec) :: nSubstep = 4
  INTEGER(iprec) :: rcm_record
  REAL(rprec) :: HighLatBD,LowLatBD
  LOGICAL :: doLatStretch = .false.
  LOGICAL :: doFLCLoss = .true. !Use FLC losses
  LOGICAL :: doNewCX = .true. !Use newer CX loss estimate
  LOGICAL :: doSmoothDDV = .true. !Whether to smooth ij deriv of residual FTV
  LOGICAL :: doSmoothBNDLOC = .true. !Whether to do bndloc smoothing
  LOGICAL :: doPPSmooth = .true. !Try to smooth plasmapause

! set this to true to tilt the dipole, must turn off corotation also
  LOGICAL :: rcm_tilted = .false.
! set this to false to turn off the dynamic plasmasphere  07242020  sbao
  LOGICAL :: dp_on = .true.
  LOGICAL, PARAMETER :: use_plasmasphere = .true.
  LOGICAL :: doAvg2MHD = .true.
  LOGICAL :: doPPRefill = .false.!Whether to refill plasmasphere
  LOGICAL :: doRelax    = .true. !Whether to relax energy distribution
  LOGICAL :: doQ0 = .true. !Whether to include implicit cold ions in tomhd moments

  INTEGER(iprec) :: ELOSSMETHOD        
  INTEGER(iprec) :: InitKp = 1, NowKp
  LOGICAL :: doFLOut = .false. !Whether to output field lines (slow)
  INTEGER(iprec) :: nSkipFL = 8 !Stride for outputting field lines

  LOGICAL :: doKapDef = .true. !Whether to do kappa by default
  LOGICAL :: doRescaleDef = .true. !Whether to rescale D,P=>eta by default
  REAL(rprec) :: staticR = 0.0
  REAL(rprec) :: LowLatMHD = 0.0
  REAL(rprec) :: rcm_pFloor = 0.0 !nPa
  REAL(rprec) :: rcm_dFloor = 0.0 !#/cc

  REAL(rprec) :: epsPk = 1.0e-3

  type RCMEllipse_T
      !Ellipse parameters
      real(rprec) :: xSun=12.5,xTail=-15.0,yDD=15.0
      logical  :: isDynamic=.true. !Whether to update parameters  
  end type RCMEllipse_T

  type EWMTauIn_T !electron lifetime wave model input
      integer(iprec) :: Nm=24, Nl=20, Nk=7 ,Ne=100
      real(rprec), ALLOCATABLE :: MLTi(:), Li(:), Kpi(:), Eki(:)
      real(rprec), ALLOCATABLE :: tau1i(:,:,:,:), tau2i(:,:,:,:)
  end type EWMTauIn_T 

  type(EWMTauIn_T) :: EWMTauInput
  type(RCMEllipse_T) :: ellBdry
  type(TimeSeries_T), private :: KpTS

  CONTAINS

      !Get RCM params from Kaiju-style XML file
      subroutine RCM_MHD_Params_XML(iXML)
        type(XML_Input_T), intent(in), optional :: iXML
        character(len=strLen) :: inpXML,tmpStr
        type(XML_Input_T) :: xmlInp

        if(present(iXML)) then
          call iXML%GetFileStr(inpXML)
        else
          !Find input deck filename
          call getIDeckStr(inpXML)
        endif

        !Create new XML reader w/ RCM as root
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju/RCM',.true.)

        !Read various parameters
        call xmlInp%Set_Val(L_write_rcmu_torcm,"output/toRCM",L_write_rcmu_torcm)
        call xmlInp%Set_Val(L_write_rcmu,"output/toMHD",L_write_rcmu)
        call xmlInp%Set_Val(L_write_vars_debug,"output/debug",L_write_vars_debug)
        call xmlInp%Set_Val(nSkipFL,"output/nSkipFL",nSkipFL)
        call xmlInp%Set_Val(doFLOut,"output/doFLOut",doFLOut)
        call xmlInp%Set_Val(rcm_tilted,"tilt/isTilt",rcm_tilted)

        !Grid bounds
        call xmlInp%Set_Val(HighLatBD,"grid/HiLat" ,75.0_rprec)
        call xmlInp%Set_Val(LowLatBD ,"grid/LowLat",30.0_rprec)
        call xmlInp%Set_Val(doLatStretch ,"grid/doLatStretch",.false.)

        !Ellipse parameters
        call xmlInp%Set_Val(ellBdry%xSun ,"ellipse/xSun" ,ellBdry%xSun )
        call xmlInp%Set_Val(ellBdry%xTail,"ellipse/xTail",ellBdry%xTail)
        call xmlInp%Set_Val(ellBdry%yDD  ,"ellipse/yDD"  ,ellBdry%yDD  )
        call xmlInp%Set_Val(ellBdry%isDynamic,"ellipse/isDynamic"  ,.true.)
        
        !Dynamic plasmaspehre parameters
        call xmlInp%Set_Val(dp_on      ,"plasmasphere/isDynamic",dp_on)
        call xmlInp%Set_Val(InitKp     ,"plasmasphere/initKp",InitKp) 
        call xmlInp%Set_Val(staticR    ,'plasmasphere/staticR',staticR)
        call xmlInp%Set_Val(doPPRefill ,'plasmasphere/doRefill',doPPRefill)
        call xmlInp%Set_Val(DenPP0     ,'plasmasphere/DenPP0',DenPP0)
        call xmlInp%Set_Val(doPPSmooth ,'plasmasphere/doPPSmooth',doPPSmooth)

        call SetKp0(InitKp)
        NowKp = InitKp

        !Loss options
        call xmlInp%Set_Val(doFLCLoss,"loss/doFLCLoss",doFLCLoss) 
        call xmlInp%Set_Val(tmpStr,"loss/eLossMethod","FDG")
        select case (tmpSTR)
           case ("FDG")
              ELOSSMETHOD = ELOSS_FDG
           case ("SS")
              ELOSSMETHOD = ELOSS_SS
           case ("C05")
              ELOSSMETHOD = ELOSS_C05
           case ("C19")
              ELOSSMETHOD = ELOSS_C19
           case ("DW")
              ELOSSMETHOD = ELOSS_DW
           case default
              stop "The electron loss type entered is not supported (Available options: FDG, SS, C05, C19)."
        end select

        call xmlInp%Set_Val(doNewCX  ,"loss/doNewCX"  ,doNewCX  )
        call xmlInp%Set_Val(doRelax  ,"loss/doRelax"  ,doRelax  )

        !Tomhd parameters
        call xmlInp%Set_Val(doAvg2MHD ,"tomhd/doAvg2MHD" ,doAvg2MHD )
        call xmlInp%Set_Val(doRelax   ,"tomhd/doRelax"   ,doRelax   )
        call xmlInp%Set_Val(doQ0      ,"tomhd/doQ0"      ,doQ0      )
        

        !Torcm parameters
        call xmlInp%Set_Val(doKapDef      ,"torcm/doKappa"        ,doKapDef      )
        call xmlInp%Set_Val(doRescaleDef  ,"torcm/doRescale" ,doRescaleDef)
        call xmlInp%Set_Val(doSmoothBNDLOC,"torcm/doSmoothBNDLOC" ,doSmoothBNDLOC)

        !Advance parameters
        call xmlInp%Set_Val(nSubstep,"sim/nSubstep", nSubstep)

        !Advection
        call xmlInp%Set_Val(doSmoothDDV,"advect/doSmoothDDV",doSmoothDDV)
        call xmlInp%Set_Val(epsPk      ,"advect/epsPk      ",epsPk      )

      !Initialize Kp (and maybe other indices) time series
        call xmlInp%Set_Val(KpTS%wID,"/Kaiju/Gamera/wind/tsfile","NONE")
        call KpTS%initTS("Kp",doLoudO=.false.)

      !Get floors from gamera part of XML
        ! call xmlInp%Set_Val(rcm_pFloor,"/Kaiju/gamera/floors/pFloor",rcm_pFloor)
        ! call xmlInp%Set_Val(rcm_dFloor,"/Kaiju/gamera/floors/dFloor",rcm_dFloor)


      end subroutine RCM_MHD_Params_XML

      !Update any indices in RCM that may be necessary
      subroutine UpdateRCMIndices(time)
        REAL(rprec), intent(in) :: time
        INTEGER(iprec) :: n
        REAL(rprec)    :: t0,t,KpMax

        NowKp = InitKp
        if (time<=0) return

        t0 = time
        KpMax = 0.0
        !Loop over +/- 15min, find max Kp
        do n=-1,+1
          t = t0 + 15.0*60.0*n
          KpMax = max(KpMax,KpTS%evalAt(t))
        enddo
        NowKp = nint(KpMax) !Cast to integer
      
      end subroutine UpdateRCMIndices

END MODULE rice_housekeeping_module
