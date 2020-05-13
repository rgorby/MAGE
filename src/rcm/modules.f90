!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CONSTANTS
  USE rcm_precision
  USE kdefs, ONLY : EarthPsi0
  USE Rcm_mod_subs, ONLY : isize,jsize,jwrap
  REAL(rprec),PARAMETER :: radius_earth_m = 6380.e3 ! Earth's radius in meters
  REAL(rprec),PARAMETER :: radius_iono_m  = 6380.e3 + 100.e3 + 20.e3 ! ionosphere radius in meters
  REAL(rprec),PARAMETER :: boltz = 1.38E-23
  REAl(rprec),PARAMETER :: mass_proton=1.6726e-27
  REAL(rprec),PARAMETER :: mass_electron=9.1094e-31
  REAL(rprec),PARAMETER :: ev=1.6022e-19
  REAL(rprec),PARAMETER :: gamma=1.6667
  REAL(rprec),PARAMETER :: one_over_gamma=0.6
  REAL(rprec),PARAMETER :: mu0 = 4.0e-7*3.14159
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
!  USE Rcm_mod_subs, ONLY : rprec,iprec
  USE rcm_precision, only : iprec,rprec,strLen
  use xml_input
  use strings

  IMPLICIT NONE
  
  LOGICAL :: L_write_rcmu          = .false., &
             L_write_rcmu_torcm    = .false., &
             L_write_tracing_debug = .false., &
             L_write_vars_debug    = .false., &
             L_write_int_grid_debug= .true.
  
  INTEGER(iprec) :: Idt_overwrite         = 1
  INTEGER(iprec) :: rcm_record
  REAL(rprec) :: HighLatBD,LowLatBD
  
! set this to true to tilt the dipole, must turn off corotation also
  LOGICAL :: rcm_tilted = .false.

    type RCMEllipse_T
        !Ellipse parameters
        real(rprec) :: xSun=7.5,xTail=-15.0,yDD=10.0
        logical  :: isDynamic=.false. !Whether to update parameters
        real(rprec) :: dRadMHD = 2.0
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

        !Ellipse parameters
        call xmlInp%Set_Val(ellBdry%xSun ,"ellipse/xSun" ,7.5_rprec)
        call xmlInp%Set_Val(ellBdry%xTail,"ellipse/xTail",-15.0_rprec)
        call xmlInp%Set_Val(ellBdry%yDD  ,"ellipse/yDD"  ,10.0_rprec)
        call xmlInp%Set_Val(ellBdry%isDynamic,"ellipse/isDynamic"  ,.true.)

        call xmlInp%Set_Val(ellBdry%dRadMHD ,"ellipse/dRadMHD" ,ellBdry%dRadMHD)

        !For now just using default Idt_overwrite

      end subroutine RCM_MHD_Params_XML

END MODULE rice_housekeeping_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE ionosphere_exchange
  use rcm_mhd_interfaces
  use rcm_mod_subs, ONLY: isize, jsize, jwrap, colat, aloct
  
  contains 
    
    !> Allocate Ionosphere Grid variables and read Ion grid from "RCM-ion.dat".
    !! Don't forget to deallocate!
    SUBROUTINE setupIon(RM)
      IMPLICIT NONE
      type(rcm_mhd_T),intent(inout) ::RM
      integer(iprec) :: lat,lon

      rm%nLat_ion = isize
      rm%nLon_ion = jsize-jwrap+1

      ALLOCATE( rm%gcolat(rm%nLat_ion) )
      ALLOCATE( rm%glong(rm%nLon_ion) )

      ALLOCATE( rm%pot(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%eng_avg(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%flux(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%fac(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Pave(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Nave(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Vol(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Bmin(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%X_bmin(rm%nLat_ion, rm%nLon_ion, 3) )
      ALLOCATE( rm%iopen(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Prcm(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Npsph(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Nrcm(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%beta_average(rm%nLat_ion, rm%nLon_ion))
      ALLOCATE( rm%sigmap(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%sigmah(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%latc(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%lonc(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Lb  (rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Tb  (rm%nLat_ion, rm%nLon_ion) )

      ALLOCATE( rm%toMHD(rm%nLat_ion, rm%nLon_ion) )
      
      rm%gcolat (:) = colat (:,1)
      rm%glong  (:) = aloct (1,jwrap:jsize)
      if (rm%glong(rm%nLon_ion) < pi) rm%glong(rm%nLon_ion) = rm%glong(rm%nLon_ion) + 2*pi
    END SUBROUTINE setupIon

    !> Deallocate any variables allocated by setupIon.
    SUBROUTINE tearDownIon(rm)
      type(rcm_mhd_T),intent(inout) ::RM

      if (ALLOCATED(rm%pot)) DEALLOCATE(rm%pot)
      if (ALLOCATED(rm%sigmap)) DEALLOCATE(rm%sigmap)
      if (ALLOCATED(rm%sigmah)) DEALLOCATE(rm%sigmah)
      if (ALLOCATED(rm%gcolat)) DEALLOCATE(rm%gcolat)
      if (ALLOCATED(rm%glong)) DEALLOCATE(rm%glong)
      if (ALLOCATED(rm%flux)) DEALLOCATE(rm%flux)
      if (ALLOCATED(rm%fac)) DEALLOCATE(rm%fac)
      if (ALLOCATED(rm%Pave)) DEALLOCATE(rm%Pave)
      if (ALLOCATED(rm%Nave)) DEALLOCATE(rm%Nave)
      if (ALLOCATED(rm%Vol)) DEALLOCATE(rm%Vol)
      if (ALLOCATED(rm%Bmin)) DEALLOCATE(rm%Bmin)
      if (ALLOCATED(rm%X_bmin)) DEALLOCATE(rm%X_bmin)
      if (ALLOCATED(rm%iopen)) DEALLOCATE(rm%iopen)
      if (ALLOCATED(rm%Prcm)) DEALLOCATE(rm%Prcm)
      if (ALLOCATED(rm%Npsph)) DEALLOCATE(rm%Npsph)
      if (ALLOCATED(rm%Nrcm)) DEALLOCATE(rm%Nrcm)
      if (ALLOCATED(rm%beta_average)) DEALLOCATE(rm%beta_average)
      if (ALLOCATED(rm%latc)) DEALLOCATE(rm%latc)
      if (ALLOCATED(rm%lonc)) DEALLOCATE(rm%lonc)
      if (ALLOCATED(rm%toMHD)) DEALLOCATE(rm%toMHD)

    END SUBROUTINE tearDownIon

  END MODULE ionosphere_exchange


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE conversion_module
!USE Rcm_mod_subs, ONLY : rprec,iprec
  USE rcm_precision
  IMPLICIT NONE
  INTEGER(iprec) :: idim,jdim,kdim
  REAL(rprec), ALLOCATABLE :: bndloc_old(:),almmin(:),almmax(:),almdel(:),&
       eta_midnight(:)
  REAL(rprec), ALLOCATABLE :: x0(:,:),y0(:,:),z0(:,:)
  REAL(rprec), ALLOCATABLE :: x0_sm(:,:),y0_sm(:,:),z0_sm(:,:)
  REAL(rprec), ALLOCATABLE :: te(:,:),ti(:,:),to(:,:),&
       eetabnd(:,:),&
       den(:,:),press(:,:),&
       deno(:,:),presso(:,:),&
       beta_average(:,:)

  REAL(rprec), ALLOCATABLE :: eeta_new(:,:,:)
  INTEGER(iprec), ALLOCATABLE :: iopen(:,:),imin_j_old(:),inner_bndy(:)
END MODULE conversion_module