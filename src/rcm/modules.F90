!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CONSTANTS
  USE rcm_precision
  USE kdefs, ONLY : EarthPsi0,Me_cgs,Mp_cgs,Mu0
  use rcmdefs, ONLY : isize,jsize,jwrap
  REAL(rprec),PARAMETER :: radius_earth_m = 6380.e3 ! Earth's radius in meters
  REAL(rprec),PARAMETER :: radius_iono_m  = 6380.e3 + 100.e3 + 20.e3 ! ionosphere radius in meters
  REAL(rprec),PARAMETER :: boltz = 1.38E-23
  REAl(rprec),PARAMETER :: mass_proton =Mp_cgs*1.0e-3
  REAL(rprec),PARAMETER :: mass_electron=Me_cgs*1.0e-3
  !REAl(rprec),PARAMETER :: mass_proton=1.6726e-27
  !REAL(rprec),PARAMETER :: mass_electron=9.1094e-31
  REAL(rprec),PARAMETER :: ev=1.6022e-19
  REAL(rprec),PARAMETER :: gamma=1.6667
  REAL(rprec),PARAMETER :: one_over_gamma=0.6
  !REAL(rprec),PARAMETER :: mu0 = 4.0e-7*3.14159
  REAL(rprec),PARAMETER :: big_vm = -1.0e5
  REAL(rprec),PARAMETER :: nt = 1.0e-9
  !REAL(rprec),PARAMETER :: tiote = 7.8
  REAL(rprec),PARAMETER :: tiote = 4.0 !Changed by K 5/9/20
  REAL(rprec),PARAMETER :: pressure_factor = 2./3.*ev/radius_earth_m*nt
  REAL(rprec),PARAMETER :: density_factor = nt/radius_earth_m
  REAL(rprec),PARAMETER :: RCMCorot = EarthPsi0*1.0e+3 ! Convert corotation to V
END MODULE CONSTANTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE rice_housekeeping_module
  USE rcm_precision, only : iprec,rprec,strLen
  use xml_input
  use strings

  IMPLICIT NONE
  
  LOGICAL :: L_write_rcmu          = .false., &
             L_write_rcmu_torcm    = .false., &
             L_write_tracing_debug = .false., &
             L_write_vars_debug    = .false., &
             L_write_int_grid_debug= .true.
  
  !INTEGER(iprec) :: Idt_overwrite         = 1
  INTEGER(iprec) :: Idt_overwrite         = 5 !K: Setting this to avoid as much unnecessary subcycling
  INTEGER(iprec) :: rcm_record
  REAL(rprec) :: HighLatBD,LowLatBD
  LOGICAL :: doLatStretch = .true.

! set this to true to tilt the dipole, must turn off corotation also
  LOGICAL :: rcm_tilted = .false.
! set this to false to turn off the dynamic plasmasphere  07242020  sbao
  LOGICAL :: dp_on = .true.
  LOGICAL, PARAMETER :: use_plasmasphere = .true.
  INTEGER(iprec) :: InitKp = 1
  REAL(rprec) :: staticR = 0.0
  REAL(rprec) :: LowLatMHD = 0.0
  type RCMEllipse_T
      !Ellipse parameters
      real(rprec) :: xSun=10.0,xTail=-15.0,yDD=10.0
      logical  :: isDynamic=.true. !Whether to update parameters
          
  end type RCMEllipse_T
  type(RCMEllipse_T) :: ellBdry

  CONTAINS
  
      SUBROUTINE Read_rcm_mhd_params

      IMPLICIT NONE

      INTEGER, PARAMETER :: Lun = 10
      LOGICAL :: L_flag

      INQUIRE (FILE='rcm_mhd.params',EXIST=L_flag)

      IF (.NOT.L_flag) THEN
         WRITE (*,*) ' RCM_MHD: no rcm_mhd.params file found, default values will be used'

      ELSE
         OPEN (LUN, FILE='rcm_mhd.params', STATUS='OLD')
         READ (LUN,*) L_write_rcmu_torcm
         READ (LUN,*) L_write_rcmu
         READ (LUN,*) L_write_vars_debug
         READ (LUN,*) Idt_overwrite
         READ (LUN,*) rcm_tilted
      END IF

      WRITE (*,*)
      WRITE (*,'(A,L7)') ' RCM_MHD:  rcmu_torcm.dat file      (in TORCM) will be written?_____', L_write_rcmu_torcm
      WRITE (*,'(A,L7)') ' RCM_MHD:  rcmu.dat  file      (in TOMHD) will be written?_____', L_write_rcmu    
      WRITE (*,'(A,L7)') ' RCM_MHD:  debug vars print to stdout (in TORCM,TOMHD)   ?_____', L_write_vars_debug
      WRITE (*,'(A,I5)') ' RCM_MHD:  Internal RCM time step (in s)  will be set to ?_____', Idt_overwrite
      WRITE (*,'(A,L7)') ' RCM_MHD:  Running RCM in tilted mode   ?______________________',rcm_tilted
      WRITE (*,*)

      RETURN
      END SUBROUTINE Read_rcm_mhd_params

      !Get RCM params from Kaiju-style XML file
      subroutine RCM_MHD_Params_XML(iXML)
        type(XML_Input_T), intent(in), optional :: iXML

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp

        if(present(iXML)) then
          call iXML%GetFileStr(inpXML)
        else
          !Find input deck filename
          call getIDeckStr(inpXML)
        endif

        !Create new XML reader w/ RCM as root
        xmlInp = New_XML_Input(trim(inpXML),'RCM',.true.)

        !Read various parameters
        call xmlInp%Set_Val(L_write_rcmu_torcm,"output/toRCM",L_write_rcmu_torcm)
        call xmlInp%Set_Val(L_write_rcmu,"output/toMHD",L_write_rcmu)
        call xmlInp%Set_Val(L_write_vars_debug,"output/debug",L_write_vars_debug)
        call xmlInp%Set_Val(rcm_tilted,"tilt/isTilt",rcm_tilted)

        !Grid bounds
        call xmlInp%Set_Val(HighLatBD,"grid/HiLat" ,75.0_rprec)
        call xmlInp%Set_Val(LowLatBD ,"grid/LowLat",15.0_rprec)
        call xmlInp%Set_Val(doLatStretch ,"grid/doLatStretch",.true.)

        !Ellipse parameters
        call xmlInp%Set_Val(ellBdry%xSun ,"ellipse/xSun" ,ellBdry%xSun )
        call xmlInp%Set_Val(ellBdry%xTail,"ellipse/xTail",ellBdry%xTail)
        call xmlInp%Set_Val(ellBdry%yDD  ,"ellipse/yDD"  ,ellBdry%yDD  )
        call xmlInp%Set_Val(ellBdry%isDynamic,"ellipse/isDynamic"  ,.true.)
        
        !Dynamic plasmaspehre parameters
        call xmlInp%Set_Val(dp_on,"plasmasphere/isDynamic",dp_on)
        call xmlInp%Set_Val(InitKp ,"plasmasphere/initKp",InitKp) 
        call xmlInp%Set_Val(staticR ,'plasmasphere/staticR',staticR)
        !For now just using default Idt_overwrite

      end subroutine RCM_MHD_Params_XML

END MODULE rice_housekeeping_module
