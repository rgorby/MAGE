module usergamic
    use kdefs
    use gamtypes
    use gambctypes
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use ioH5
    use cmetypes
    use glutils
    use glsolution
    use helioutils
    use clocks

    implicit none

    enum, bind(C)
       ! variables passed via innerbc file
       ! Br, Vr, Rho, Temperature, Br @ kface, Vr @ kface
       enumerator :: BRIN=1,VRIN,RHOIN,TIN,BRKFIN,VRKFIN
    endenum 

    enum, bind(C)
        ! variables passed via innerbc file
        ! Br, Bt, Bp, Vr, Inside_Mask
        enumerator :: CMEBR = 1, CMEBT, CMEBP, CMEVR, CMEMASK
    endenum 

    integer, private, parameter :: NVARSIN=6 ! SHOULD be the same as the number of vars in the above enumerator
    real(rp), dimension(:,:,:,:), allocatable :: ibcVars
    ! Br, Bt, Bp, Vr, Inside_Mask
    integer, private, parameter :: NCMEVARS=5
    real(rp) :: Den_CME, T_CME, Pres_CME

    !Various global would go here
    real (rp) :: Rho0, P0, Vslow,Vfast, wScl, Cs0, B0, MJD_c

    ! things we keep reusing
    real(rp), dimension(NDIM) :: xyz,xyz0,rHat,phiHat
    real(rp) :: Rfactor

    ! global grid
    integer :: gNkp

    ! FIXME
    ! setting it here temporarily. eventually need to read from HDF or something
    ! note also that this is slightly incorrect, since Rbc below is used as the radius
    ! of the center of the first ghost cell.

    !for inner heliosphere
    real(rp) :: Rbc = 21.5

    real(rp) :: Tsolar ! Solar rotation period, defined in apps/helioutils.F90

    !DeltaT sets a rotation of the wsa map relative to +X at simulation time t=0 
    !DeltaT = 0 then the map joint is at +X at t=0
    !DeltaT = Tsolar/2 then the map joint is at -X at t=0 (facing Earth)
    real(rp) :: DeltaT 

    character(len=strLen) :: wsaFile

    ! CME Models
    class(baseCMEModel_T), allocatable :: cccmeModel, kfcmeModel, jfcmeModel
    class(baseCMEState_t), allocatable :: cccmeState, kfcmeState, jfcmeState
    class(baseCMESolution_T), allocatable :: cccmeSolution, kfCMESolution, jfCMESolution
    real(rp) :: emerge_lastP, t_smooth
    real(rp) :: emergence_times(5)
    type(XML_Input_T) :: inpXMLCME
    character(len=strLen) :: inpXMLStr

    ! use this to fix the Efield at the inner boundary
    !    real(rp), allocatable :: inEijk(:,:,:,:)

    ! type for solar wind BC
    type, extends(innerIBC_T) :: SWInnerBC_T

        !Main electric field structures
        real(rp), allocatable, dimension(:,:,:,:) :: inEijk,inExyz
        ! CME Model
        contains

        !        procedure :: doInit => InitIonInner
        ! TODO: this shoudl be made generic (wsa, mas, etc.) How
        procedure :: doBC => wsaBC
    end type SWInnerBC_T

    contains

    !>
    !>
    !>
    !>
    subroutine initUser(Model,Grid,State,inpXML)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML
        procedure(GasIC_T), pointer :: Wxyz
        procedure(HackStep_T), pointer :: tsHack
        procedure(HackE_T), pointer :: eHack

        integer :: i,j,k,nvar,nr,d
        integer :: n1, n2
        integer :: kg, ke, kb, jg
        integer, dimension(6) :: bounds
        real(rp) :: a
        real(rp) :: Tsolar_synodic

        real(rp) :: ibcVarsStatic !For br only
        logical :: isSpSymSW
        ! if (.not.allocated(inEijk)) allocate(inEijk(1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM))

        ! set units and other thins, like Tsolar
        call setHeliosphere(Model,inpXML,Tsolar)

        ! grab inner 
        call inpXML%Set_Val(wsaFile,"prob/wsaFile","innerbc.h5" )
        call inpXML%Set_Val(isSpSymSW, "prob/isSpSymSW", .false.)
        call inpXML%Set_Val(DeltaT, "prob/DeltaT", 0.0)
        ! compute global Nkp
        gNkp = Grid%Nkp*Grid%NumRk

        ! initial conditions
        ! TODO: change using norm. units in Model, set in helioutils
        Cs0   = 0.267  ! 40 km/s
        Vslow = 1.33   ! 200 km/s
        Vfast = 5.33   ! 800 km/s
        !for inner helio
        B0    = 2.0    ! 200 nT
        Rho0  = 1.0    ! 200/cc
        P0    = 1.0e-4*Rho0*Cs0**2.0/Model%gamma

        ![OHelio] for 1-10 au helio
        !Cs0 = 0.78    !27 km/s
        !B0 = 1.       ! 5nT
        !Rho0 = 1.     ! 10/cc
        !Vslow = 8.5   !300 km/s 

        ! deallocate default BCs
        ! required because defaults are triply-periodic
        ! need to wipe
        call WipeBCs(Model,Grid)

        !Set BCs
        allocate(SWInnerBC_T        :: Grid%externalBCs(1)%p)
        allocate(helioOuterIBC_T    :: Grid%externalBCs(2)%p)
        allocate(helioInnerJBC_T    :: Grid%externalBCs(3)%p)
        allocate(helioOuterJBC_T    :: Grid%externalBCs(4)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(5)%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(6)%p)

        ! Use CME Model
        if(Model%doCME) then
            ! Grid%isg,Grid%is-1, Grid%js:Grid%je, Grid%ks:Grid%ke
            bounds(1) = 1
            bounds(2) = Model%Ng
            bounds(3) = Grid%js
            bounds(4) = Grid%je
            bounds(5) = Grid%ks
            bounds(6) = Grid%ke
            call getIDeckStr(inpXMLStr)
            inpXMLCME = New_XML_Input(trim(inpXMLStr),'Kaiju/CME',.true.)
            call inpXMLCME%Set_Val(Model%cmeModel,"sim/model","GL") 
            select case (Model%cmeModel) 
                case ("GL")
                    allocate(glModel_T :: cccmeModel)
                    allocate(glState_T :: cccmeState)
                    allocate(glSolution_T :: ccCMESolution)
                    allocate(glModel_T :: jfcmeModel)
                    allocate(glState_T :: jfcmeState)
                    allocate(glSolution_T :: jfCMESolution)
                    allocate(glModel_T :: kfcmeModel)
                    allocate(glState_T :: kfcmeState)
                    allocate(glSolution_T :: kfCMESolution)
                    select type (ccCMESolution)
                        type is (glSolution_T)
                            select type (cccmeState)
                                type is(glState_T)
                                    select type (cccmeModel)
                                        type is(glModel_T)
                                            cccmeState%Ranki = Grid%Ri
                                            cccmeState%Rankj = Grid%Rj
                                            cccmeState%Rankk = Grid%Rk
                                            call initGLInterface(Grid%xyzcc(Grid%isg:Grid%is-1, Grid%js:Grid%je, Grid%ks:Grid%ke, :),cccmeModel, cccmeState, ccCMESolution, inpXMLCME, bounds)
                                            emergence_times = calcEmergenceTimes(cccmeModel,cccmeState, Model%Units%gT0)
                                            emerge_lastP = maxval(emergence_times)
                                    end select
                            end select
                    end select
                    select type (jfCMESolution)
                        type is (glSolution_T)
                            select type (jfcmeState)
                                type is(glState_T)
                                    select type (jfcmeModel)
                                        type is(glModel_T)
                                            jfcmeState%Ranki = Grid%Ri
                                            jfcmeState%Rankj = Grid%Rj
                                            jfcmeState%Rankk = Grid%Rk
                                            call initGLInterface(Grid%xfc(Grid%isg:Grid%is-1, Grid%js:Grid%je, Grid%ks:Grid%ke, :, JDIR), jfcmeModel, jfcmeState, jfCMESolution, inpXMLCME, bounds)                              
                                    end select
                            end select
                    end select
                    select type (kfCMESolution)
                        type is (glSolution_T)
                            select type (kfcmeState)
                                type is(glState_T)
                                    select type (kfcmeModel)
                                        type is(glModel_T)
                                            kfcmeState%Ranki = Grid%Ri
                                            kfcmeState%Rankj = Grid%Rj
                                            kfcmeState%Rankk = Grid%Rk
                                            call initGLInterface(Grid%xfc(Grid%isg:Grid%is-1, Grid%js:Grid%je, Grid%ks:Grid%ke, :, KDIR), kfcmeModel, kfcmeState, kfCMESolution, inpXMLCME, bounds)                                     
                                    end select
                            end select
                    end select            
            end select
            call inpXMLCME%Set_Val(Den_CME,"prob/Den_CME",0. )
            call inpXMLCME%Set_Val(T_CME,"prob/T_CME",0. )
            t_smooth = emerge_lastP + 3600./Model%Units%gT0
            Pres_CME = (Den_CME/Model%Units%gD0)*Kbltz*T_CME/Model%Units%gP0
        end if

        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

        ! Add gravity
        !        eHack  => EFix
        !        Model%HackE => eHack
        !        tsHack => PerStep
        !        Model%HackStep => tsHack
   
         !Write MJD_center_of_WSA_map to the root of H5 output 
         Model%HackIO_0 => writeMJDcH5Root

         if (isSpSymSW) then
            ! spherically symmetric solar wind
            call sphereSymSW(Model, Grid, inpXML)
        else
            ! everybody reads WSA data
            call readIBC(wsaFile)
        end if
        !MJD0 is MJD_center_of_WSA_map - Tsolar_synodic/2; Tsolar_synodic = 27.28
        ![EP] TO DO: add synodic Tsolar to constants
        Tsolar_synodic = 27.28
        Model%MJD0 = MJD_c - Tsolar_synodic/2.

        !if not restart set State%time according the tSpin
        State%time  = Model%t

        !Map IC to grid
        Wxyz => GasIC
        call GasIC2State(Model,Grid,State,Wxyz)

        ! NOTE, not filling ghost cells here
        ! relying on J, K boundary conditions
        ! first zero everything out
        State%magFlux = 0.0

        do k=Grid%ks,Grid%ke
           kg = k+Grid%ijkShift(KDIR)
           ! map rotating to static grid
           call mapKinit(kg,ke,kb,a)

           do j=Grid%js,Grid%je
              jg = j+Grid%ijkShift(JDIR)
              do i=Grid%isg,Grid%ieg+1

                 ibcVarsStatic = 1._rp

                 if ( (jg>=Grid%js).and.(jg<=size(ibcVars,2)) ) then
                     ! interpolate BRIN linearly from rotating to inertial frame
                     ibcVarsStatic = a*ibcVars(Model%Ng,jg,kb,BRIN)+(1-a)*ibcVars(Model%Ng,jg,ke,BRIN)
                 endif

                 ! FIXME: calc Rbc appropriately above!
                 Rfactor = Rbc/norm2(Grid%xfc(Grid%is,j,k,:,IDIR))
                 !write(*,*) shape(Grid%xfc)
                 ! note scaling br to the first active face
                 ! this is opposite to what we did on the python side
                 ! not elegant, but otherwise we'd have to store both br and br_iface in teh innerbc.h5 file
                 !
                 ! note also that ibcVars(Model%Ng,:,:,:) corresponds to the center of the first ghost cell (just below the boundary)
                 !State%magFlux(i,j,k,IDIR) = ibcVars(Model%Ng,j+Grid%ijkShift(JDIR),k+Grid%ijkShift(KDIR),BRIN)*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
                  State%magFlux(i,j,k,IDIR) = ibcVarsStatic*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
              enddo
           enddo
        enddo

        !Local functions
        !NOTE: Don't put BCs here as they won't be visible after the initialization call

        contains
            subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
                real (rp), intent(in) :: x,y,z
                real (rp), intent(out) :: D,Vx,Vy,Vz,P
                real (rp) :: r, theta, phi
                real (rp) :: r_unit(NDIM)

                D = Rho0
                P = P0

                !Calculate radial hat vector
                r = sqrt(x**2.0 + y**2.0 + z**2.0)
                theta = acos(z/R)
                phi = atan2(y,x)
                r_unit = [x,y,z]/r

                !Set primitives (already have D/P)
                Vx = r_unit(XDIR)*Vslow
                Vy = r_unit(YDIR)*Vslow
                Vz = r_unit(ZDIR)*Vslow
            end subroutine GasIC

            subroutine mapKinit(k,ke,kb,a)
                ! find the lower and upper k-index of the rotating cell on the static grid
                integer, intent(in) :: k
                integer, intent(out) :: ke,kb  ! upper (ke) and lower (kb) indices
                real(rp), intent(out):: a      ! interp coefficient

                real(rp) :: kprime

                kprime = modulo(k-gNkp*(State%time+DeltaT)/Tsolar,real(gNkp,kind=rp))

                if ((kprime.ge.0).and.(kprime.lt.1)) then
                    kb=gNkp
                    ke=1
                else
                    kb=floor(kprime)
                    ke=ceiling(kprime)
                endif
                a = ke-kprime

            end subroutine mapKinit

    end subroutine initUser

    !> Inner-I BC for WSA-Gamera
    !>
    !>
    subroutine wsaBC(bc,Model,Grid,State)
      class(SWInnerBC_T), intent(inout) :: bc
      type(Model_T), intent(in) :: Model
      type(Grid_T), intent(in) :: Grid
      type(State_T), intent(inout) :: State

      ! local variables
      integer :: i,j,k, kb, ke, knew
      integer :: kg, jg, ig ! global indices (in preparation for MPI)
      integer :: var ! ibcVar variable number
      real(rp) :: a, dphi
      real(rp) :: ibcVarsStatic(NVARSIN)

      !including ghost cells in k
      real(rp) :: ccCMEVars(Model%Ng, Grid%js:Grid%je, Grid%ks:Grid%ke, NCMEVARS)
      real(rp) :: jfCMEVars(Model%Ng, Grid%js:Grid%je, Grid%ks:Grid%ke, NCMEVARS)  ! j-face
      real(rp) :: kfCMEVars(Model%Ng, Grid%js:Grid%je, Grid%ks:Grid%ke, NCMEVARS)  ! k-face

      real(rp) :: ccCMEVarsStatic(NCMEVARS)
      real(rp) :: jfCMEVarsStatic(NCMEVARS)
      real(rp) :: kfCMEVarsStatic(NCMEVARS)
      real(rp) :: R, Theta, Phi
      real(rp) :: Theta_kf, R_kf ! kface
      real(rp), dimension(NVAR) :: conVar, pVar

      ! Set cme vars to zero initially
      ccCMEVars = 0.
      jfCMEVars = 0.
      kfCMEVars = 0.
      ccCMEVarsStatic = 0.
      jfCMEVarsStatic = 0.
      kfCMEVarsStatic = 0.
      ! After Initial transient passage, generate CME solution
      if (Model%doCME .and. (Model%t >= cccmeModel%Tstart_transient/Model%Units%gT0)) then
            call Tic("CME Update")
            ! cell-center
            !call cccmeModel%updateModelTime(Model%t, Model%Units%gT0)
            if(Model%rotateCME) then 
                ! Set new longitude in this rank for model based on solar rotation
                dphi = 2.*PI*(Model%t+DeltaT)/Tsolar
                cccmeState%currentLongitude = cccmeModel%longitude + dphi
                jfcmeState%currentLongitude = jfcmeModel%longitude + dphi
                kfcmeState%currentLongitude = kfcmeModel%longitude + dphi
                ! Update model grid based on the above longitude update
                call cccmeState%updateGrid(Grid%xyzcc(Grid%isg:Grid%is-1, Grid%js:Grid%je, Grid%ks:Grid%ke, :), cccmeModel)
                call jfcmeState%updateGrid(Grid%xfc(Grid%isg:Grid%is-1, Grid%js:Grid%je, Grid%ks:Grid%ke, :, JDIR), jfcmeModel)
                call kfcmeState%updateGrid(Grid%xfc(Grid%isg:Grid%is-1, Grid%js:Grid%je, Grid%ks:Grid%ke, :, KDIR), kfcmeModel)                               
            end if
            cccmeModel%time = (Model%t - cccmeModel%Tstart_transient/Model%Units%gT0)*Model%Units%gT0
            if (cccmeModel%isDebug) then
                write(*,"(1X,A14,2X,F)") "Sim time: ", Model%t 
                write(*,"(1X,A14,2X,F)") "Updated CME time: ", cccmeModel%time 
            end if
            call ccCMESolution%generateSolution(ccCMESolution, cccmeModel, cccmeState)
            ! J-Faces
            !call jfcmeModel%updateModelTime(Model%t, Model%Units%gT0)
            jfcmeModel%time = (Model%t - jfcmeModel%Tstart_transient/Model%Units%gT0)*Model%Units%gT0
            call jfCMESolution%generateSolution(jfCMESolution, jfcmeModel, jfcmeState)
            ! K-Faces
            !call kfcmeModel%updateModelTime(Model%t, Model%Units%gT0)
            kfcmeModel%time = (Model%t - kfcmeModel%Tstart_transient/Model%Units%gT0)*Model%Units%gT0
            call kfCMESolution%generateSolution(kfCMESolution, kfcmeModel, kfcmeState)
            if (cccmeModel%isDebug) then
                write(*,"(1X,A14,2X,3F)") "Max Vr: ", maxval(ccCMESolution%v(:,:,:,XDIR)), maxval(jfCMESolution%v(:,:,:,XDIR)), maxval(kfCMESolution%v(:,:,:,XDIR))
                write(*,"(1X,A14,2X,3F)") "Max Br: ", maxval(ccCMESolution%b(:,:,:,XDIR)), maxval(jfCMESolution%b(:,:,:,XDIR)), maxval(kfCMESolution%b(:,:,:,XDIR))
            end if
            ccCMEVars(:, :, :, CMEBR:CMEBP) = ccCMESolution%b
            ccCMEVars(:, :, :, CMEVR) = ccCMESolution%v(:, :, :, XDIR)
            ccCMEVars(:, :, :, CMEMASK) = ccCMESolution%inside_mask
            jfCMEVars(:, :, :, CMEBR:CMEBP) = jfCMESolution%b
            jfCMEVars(:, :, :, CMEVR) = jfCMESolution%v(:, :, :, XDIR)
            jfCMEVars(:, :, :, CMEMASK) = jfCMESolution%inside_mask
            kfCMEVars(:, :, :, CMEBR:CMEBP) = kfCMESolution%b
            kfCMEVars(:, :, :, CMEVR) = kfCMESolution%v(:, :, :, XDIR)
            kfCMEVars(:, :, :, CMEMASK) = kfCMESolution%inside_mask
            call Toc("CME Update")
      end if

      !i-boundaries (IN)
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,k,jg,kg,ke,kb,knew,a,var,xyz,R,Theta,Phi,rHat,phiHat) &
      !$OMP private(ibcVarsStatic,ccCMEVarsStatic,jfCMEVarsStatic,kfCMEVarsStatic) &
      !$OMP private(pVar,conVar,xyz0,R_kf,Theta_kf)
      do k=Grid%ksg,Grid%keg+1  ! note, going all the way to last face for mag fluxes
         kg = k+Grid%ijkShift(KDIR)
         ! map rotating to static grid
         call mapK(kg,ke,kb,a)

         do j=Grid%js,Grid%je+1 ! note, going all the way to last face for mag fluxes
            jg = j+Grid%ijkShift(JDIR)

            do ig=1,Model%Ng
                i=ig-Model%Ng
                ! note, ibcVars are not defined in the global j-corners or k-corners
                ! access to k-corners is fine because mapK function will always map into active domain (FIXME: check for k=nk+1 face!)
                ! access to j-corners gets into unallocated space for ibcVars. Use this trick instead:
                ! set everything arbitrarily to 1. in the global corners (on low and high j boundaries)
                ! we then apply the j-boundary after i boundary anyway, so the corners will be overwritten
                ibcVarsStatic = 1._rp  
                ! otherwise
                ! interpolate linearly from rotating to inertial frame 
                if ((jg>=Grid%js) .and. (jg<=size(ibcVars,2))) then
                    do var=1,NVARSIN
                        ibcVarsStatic(var) = a*ibcVars(ig,jg,kb,var) + (1-a)*ibcVars(ig,jg,ke,var)
                    end do
                end if

                ! Set CME local values
                if(Model%doCME .and. (Model%t >= cccmeModel%Tstart_transient/Model%Units%gT0)) then
                    !filling ghost cells in k with for CME Solution
                    !if (k .gt. Grid%ke) then
                    !    knew = modulo(k, Grid%ke)
                    !    ccCMEVars(ig, j, k, CMEBR:CMEBP) = ccCMESolution%b(ig, j, knew, :)
                    !    ccCMEVars(ig, j, k, CMEVR) = ccCMESolution%v(ig, j, knew, XDIR)
                    !    ccCMEVars(ig, j, k, CMEMASK) = ccCMESolution%inside_mask(ig, j, knew)
                    !    jfCMEVars(ig, j, k, CMEBR:CMEBP) = jfCMESolution%b(ig, j, knew, :)
                    !    jfCMEVars(ig, j, k, CMEVR) = jfCMESolution%v(ig, j, knew, XDIR)
                    !    jfCMEVars(ig, j, k, CMEMASK) = jfCMESolution%inside_mask(ig, j, knew)
                    !    kfCMEVars(ig, j, k, CMEBR:CMEBP) = kfCMESolution%b(ig, j, knew, :)
                    !    kfCMEVars(ig, j, k, CMEVR) = kfCMESolution%v(ig, j, knew, 1)
                    !    kfCMEVars(ig, j, k, CMEMASK) = kfCMESolution%inside_mask(ig, j, knew)
                    !end if
                    !if (k .lt. Grid%ks) then
                    !    knew = Grid%ke + k
                    !    ccCMEVars(ig, j, k, CMEBR:CMEBP) = ccCMESolution%b(ig, j, knew, :)
                    !    ccCMEVars(ig, j, k, CMEVR) = ccCMESolution%v(ig, j, knew, XDIR)
                    !    ccCMEVars(ig, j, k, CMEMASK) = ccCMESolution%inside_mask(ig, j, knew)
                    !    jfCMEVars(ig, j, k, CMEBR:CMEBP) = jfCMESolution%b(ig, j, knew, :)
                    !    jfCMEVars(ig, j, k, CMEVR) = jfCMESolution%v(ig, j, knew, XDIR)
                    !    jfCMEVars(ig, j, k, CMEMASK) = jfCMESolution%inside_mask(ig, j, knew)
                    !    kfCMEVars(ig, j, k, CMEBR:CMEBP) = kfCMESolution%b(ig, j, knew, :)
                    !    kfCMEVars(ig, j, k, CMEVR) = kfCMESolution%v(ig, j, knew, XDIR)
                    !    kfCMEVars(ig, j, k, CMEMASK) = kfCMESolution%inside_mask(ig, j, knew)
                    !end if
                    !if ((j<=Grid%je) .and. (k>=Grid%ks) .and. (k<=Grid%ke)) then
                    cccmeVarsStatic = ccCMEVars(ig, j, k, :)
                    jfcmeVarsStatic = jfCMEVars(ig, j, k, :)
                    kfcmeVarsStatic = kfCMEVars(ig, j, k, :)
                    !end if
                end if
                
                ! do cell centered things for cell-centers only
                if ( (j/=Grid%je+1).and.(k/=Grid%keg+1) ) then
                    ! various geometrical quantities for the cell center
                    xyz  = Grid%xyzcc(i,j,k,:)
                    R    = norm2(xyz)
                    Theta = acos(xyz(3)/R)
                    Phi = atan2(xyz(2),xyz(1))
                    rHat = xyz/R
                    phiHat = [-sin(phi),cos(phi),0._rp]

                    ! NOTE, WSA data were already scaled appropriately in the python code
                    ! TODO: save them in the innerbc hdf file and set in helioutils appropriately
                    if (Model%doCME .and. (ccCMEVarsStatic(CMEMASK) .gt. 0.1) .and. (Model%t >= cccmeModel%Tstart_transient/Model%Units%gT0)) then !inside the bubble mask is 1
                        if(cccmeModel%isDebug) write(*,"(1X,A14,2X,1F,2X,A14,2X,3I)") "Inside bubble at time: ", Model%t, "i,j,k = ", i, j, k
			            if (Model%t <= emerge_lastP) then
                            if(ccCMEModel%isDebug) write(*,"(1X,A14,2X,4F)") "CME VR Before emerge_lastP: time, unscaled VR, scaled VR, WSA VRIN", Model%t, ccCMEVarsStatic(CMEVR), ccCMEVarsStatic(CMEVR)/(1.d-5*Model%Units%gv0), ibcVarsStatic(VRIN)
                            pVar(DEN) = Den_CME/Model%Units%gD0
                            pVar(PRESSURE) = Pres_CME
                            !pVar(VELX:VELZ) = rHat*max(ccCMEVarsStatic(CMEVR)/(1.d-5*Model%Units%gv0), ibcVarsStatic(VRIN))
                            pVar(VELX:VELZ) = rHat*ccCMEVarsStatic(CMEVR)/(1.d-5*Model%Units%gv0)
                        else if ((Model%t > emerge_lastP) .and. (Model%t <= t_smooth)) then !last point passed -- smoothly transition to WSA via linear changing from CME values to WSA over 1 hour
                            if(ccCMEModel%isDebug) write(*,"(1X,A14,2X,4F)") "CME VR After emerge_lastP: time, unscaled VR, scaled VR, WSA VRIN:", Model%t, ccCMEVarsStatic(CMEVR), ccCMEVarsStatic(CMEVR)/(1.d-5*Model%Units%gv0), ibcVarsStatic(VRIN)
                            pVar(DEN) = Den_CME/Model%Units%gD0 + &
                                            (Model%t - emerge_lastP)/(3600./Model%Units%gt0)*(ibcVarsStatic(RHOIN) - Den_CME/Model%Units%gD0)
                            pVar(PRESSURE)  = Pres_CME + &
                                                (Model%t - emerge_lastP)/(3600./Model%Units%gt0)*(ibcVarsStatic(RHOIN)*Model%Units%gD0*Kbltz*ibcVarsStatic(TIN)/(Model%Units%gP0) - Pres_CME)
                            !pVar(VELX:VELZ) = rHat*max(ccCMEVarsStatic(CMEVR)/(1.d-5*Model%Units%gv0), ibcVarsStatic(VRIN))
                            !pVar(VELX:VELZ) = rHat*ccCMEVarsStatic(CMEVR)/(1.d-5*Model%Units%gv0)
                            pVar(VELX:VELZ) = rHat*ibcVarsStatic(VRIN)
                        else ! when Model%t is larger than t_smooth use WSA
                            if(cccmeModel%isDebug) write(*,"(1X,A14,2X,2F)") "Time is greater then t_smooth: ", Model%t, t_smooth
			                pVar(DEN) = ibcVarsStatic(RHOIN)
                            pVar(PRESSURE) = ibcVarsStatic(RHOIN)*Model%Units%gD0*Kbltz*ibcVarsStatic(TIN)/(Model%Units%gP0)
                            pVar(VELX:VELZ) = rHat*ibcVarsStatic(VRIN)
                            !pVar(VELX:VELZ) = rHat*ccCMEVarsStatic(CMEVR)/(1.d-5*Model%Units%gv0)
                        end if
                        if(cccmeModel%isDebug) write(*,"(1X,A34,2X,5F,3I)") "Inside pVar(Dens, Pres, VELX:VELZ): ", pVar(DEN), pVar(PRESSURE), pVar(VELX:VELZ), i, j, k
                    else if (Model%doCME .and. (Model%t >= cccmeModel%Tstart_transient/Model%Units%gT0)) then !outside the bubble mask is 0
                        pVar(DEN) = ibcVarsStatic(RHOIN)
                        pVar(PRESSURE) = ibcVarsStatic(RHOIN)*Model%Units%gD0*Kbltz*ibcVarsStatic(TIN)/(Model%Units%gP0)
                        pVar(VELX:VELZ) = rhat*ibcVarsStatic(VRIN)
                        !pVar(VELX:VELZ) = rHat*ccCMEVarsStatic(CMEVR)/(1.d-5*Model%Units%gv0)
                    else !Spin-up wind
                        !Set primitives
                        pVar(DEN)       = ibcVarsStatic(RHOIN)
                        ! note conversion to my units with B0^2/4pi in the denominator
                        pVar(PRESSURE)  = ibcVarsStatic(RHOIN)*Model%Units%gD0*Kbltz*ibcVarsStatic(TIN)/(Model%Units%gP0)
                        pVar(VELX:VELZ) = rHat*ibcVarsStatic(VRIN) 
                        !if(cccmeModel%isDebug) write(*,"(1X,A36,2X,5F,3I)") "Outside pVar(Dens, Pres, VELX:VELZ): ", pVar(DEN), pVar(PRESSURE), pVar(VELX:VELZ), i, j, k
                    end if
                    !if(cccmeModel%isDebug) write(*,"(1X,A14,2X,1F)") "Current Time: ", Model%t
                    !if(cccmeModel%isDebug) write(*,"(1X,A36,2X,5F)") "pVar(Dens, Pres, VELX:VELZ): ", pVar(DEN), pVar(PRESSURE), pVar(VELX:VELZ)
                    !Swap prim->con in ghost variables
                    call CellP2C(Model,pVar,conVar)
                    State%Gas(i,j,k,:,BLK) = conVar
                    ! note, don't need cc Bxyz because we run flux2field through ghosts
                end if

                ! also need theta at k-face for k-flux
                ! although we're assuming theta and R don't change from cell to k-face center,
                ! xyzcc used under if statement above is only defined for cell centers so need to define it here
                xyz0 = Grid%xfc(i,j,k,:,KDIR) !just reusing a temp var here.
                R_kf = norm2(xyz0)
                Theta_kf = acos(xyz0(ZDIR)/R_kf)

                ! note scaling for Bi. See note above in InitUser
                !if(norm2(Grid%xfc(Grid%is,j,k,:,IDIR)) == 0.0) write(*,"(1X,A20,2X,3I,3F)") "Grid xfc zero!: ", Grid%is, j, k, Grid%xfc(Grid%is,j,k,:,IDIR)
                Rfactor = Rbc/norm2(Grid%xfc(Grid%is,j,k,:,IDIR))
                if (Model%doCME .and. (ccCMEVarsStatic(CMEMASK) .gt. 0.1) .and. (Model%t >= cccmeModel%Tstart_transient/Model%Units%gT0)) then
                    State%magFlux(i,j,k,IDIR) = ibcVarsStatic(BRIN)*Rfactor**2*Grid%face(Grid%is,j,k,IDIR) + &
                                                ccCMEVarsStatic(CMEBR)*Grid%face(i,j,k,IDIR)/Model%Units%gB0
                    State%magFlux(i,j,k,JDIR) = jfCMEVarsStatic(CMEBT)*Grid%face(i,j,k,JDIR)/Model%Units%gB0
                    State%magFlux(i,j,k,KDIR) = - 2*PI/Tsolar*R_kf*sin(Theta_kf)/ibcVarsStatic(VRKFIN)*ibcVarsStatic(BRKFIN)*Grid%face(i,j,k,KDIR) + &
                                                kfCMEVarsStatic(CMEBP)*Grid%face(i,j,k,KDIR)/Model%Units%gB0
                    if(cccmeModel%isDebug) write(*,"(1X,A20,2X,3F,3I)") "Inside Mag Flux: ", State%magFlux(i,j,k,XDIR),State%magFlux(i,j,k,YDIR),State%magFlux(i,j,k,ZDIR), i, j, k
                else
                    State%magFlux(i,j,k,IDIR) = ibcVarsStatic(BRIN)*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
                    State%magFlux(i,j,k,JDIR) = 0.0
                    State%magFlux(i,j,k,KDIR) = - 2*PI/Tsolar*R_kf*sin(Theta_kf)/ibcVarsStatic(VRKFIN)*ibcVarsStatic(BRKFIN)*Grid%face(i,j,k,KDIR)
                    !if(cccmeModel%isDebug) write(*,"(1X,A20,2X,3F,3I)") "Outside Mag Flux: ", State%magFlux(i,j,k,XDIR),State%magFlux(i,j,k,YDIR),State%magFlux(i,j,k,ZDIR), i, j, k
               end if
            end do
         end do
      end do

    contains
      subroutine mapK(k,ke,kb,a)
        ! find the lower and upper k-index of the rotating cell on the static grid
        integer, intent(in) :: k
        integer, intent(out) :: ke,kb  ! upper (ke) and lower (kb) indices
        real(rp), intent(out):: a      ! interp coefficient

        real(rp) :: kprime

        kprime = modulo(k-gNkp*(State%time+DeltaT)/Tsolar,real(gNkp,kind=rp)) 

        if ((kprime.ge.0).and.(kprime.lt.1)) then
           kb=gNkp 
           ke=1
        else
           kb=floor(kprime)
           ke=ceiling(kprime)
        endif
        a = ke-kprime

      end subroutine mapK

    end subroutine wsaBC

    subroutine eFix(Model,Gr,State)
      type(Model_T), intent(in) :: Model
      type(Grid_T), intent(in) :: Gr
      type(State_T), intent(inout) :: State

      ! see example of how to do this in voltron/ICs/earthcmi.F90
      ! Grid%externalBCs(1)%p
      ! State%Efld(Gr%is  ,:,:,JDIR:KDIR) = inEijk(1,:,:,JDIR:KDIR)*Gr%edge(Gr%is  ,:,:,JDIR:KDIR)


    end subroutine eFix

    !> Read in innerbc file for helio background
    !>
    !>
    subroutine readIBC(ibcH5)
        character(len=*), intent(in) :: ibcH5
        logical :: fExist
        integer :: i,nvar,dims(3)
        integer, parameter :: MAXIOVAR = 50
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars


        !Reset IO chain
        call ClearIO(IOVars)

        inquire(file=ibcH5,exist=fExist)
        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open innerbc file, exiting'
            stop
        endif

        !Setup input chain
        call AddInVar(IOVars,"vr")
        call AddInVar(IOVars,"vr_kface")
        call AddInVar(IOVars,"rho")
        call AddInVar(IOVars,"temp")
        call AddInVar(IOVars,"br")
        call AddInVar(IOVars,"br_kface")
        call AddInVar(IOVars,"MJD")

        call ReadVars(IOVars,.false.,ibcH5) !Don't use io precision

        ! NOTE, assuming they all have the same dimesnions here (NO2,NJ,NK)
        ! see wsa2gamera
        dims=IOVars(1)%dims(1:3) ! i,j,k
        write(*,"(1x,A16,2x,3I)") "innerbc dims: ", dims
        if (.not.allocated(ibcVars)) allocate(ibcVars(dims(1),dims(2),dims(3),NVARSIN))

        do i=1,NVARSIN
            select case (i)
            case (BRIN)
            nvar= FindIO(IOVars,"br")
            case (VRIN)
            nvar= FindIO(IOVars,"vr")
            case (RHOIN)
            nvar= FindIO(IOVars,"rho")
            case (TIN)
            nvar= FindIO(IOVars,"temp")
            case (BRKFIN)
            nvar= FindIO(IOVars,"br_kface")
            case (VRKFIN)
            nvar= FindIO(IOVars,"vr_kface")
            end select

            ibcVars(:,:,:,i) = reshape(IOVars(nvar)%data,dims)
        end do
        !reading modified julian date from innerbc
        MJD_c = GetIOReal(IOVars,"MJD")
   
    end subroutine readIBC

    !> Set IBC vars to spherically symmetric solar wind
    !>
    !>
    subroutine sphereSymSW(Model, Grid, inpXML)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(XML_Input_T), intent(in) :: inpXML
        real(rp), dimension(:), allocatable :: R_gc
        real(rp), dimension(NVARSIN) :: helioVarsIn
        real(rp) :: R, Theta, ratio_sq
        integer :: i, j, k, ig, jg, kg
        character(len=strLen) :: swmodel

        write(*,"(1x,A16,2x,3I)") "innerbc dims: ", Model%Ng, Grid%Njp*Grid%NumRj, Grid%Nkp*Grid%NumRk
        !write(*,"(1x,A16,2x,6I)") "innerbc dims: ", Grid%isg, Grid%is-1, Grid%js, Grid%je + 1, Grid%ks, Grid%ke + 1
        if (.not. allocated(ibcVars)) allocate (ibcVars(Model%Ng, Grid%Njp*Grid%NumRj, Grid%Nkp*Grid%NumRk, NVARSIN))
        allocate (R_gc(Grid%isg:Grid%is - 1))
        ibcVars = 0.0
        call inpXML%Set_Val(helioVarsIn(TIN), "helio/tin", 3.d5 )
        call inpXML%Set_Val(helioVarsIn(BRIN), "helio/brin", 0.001)
        call inpXML%Set_Val(helioVarsIn(BRKFIN), "helio/brkfin", 0.001)
        call inpXML%Set_Val(helioVarsIn(RHOIN), "helio/rhoin", 800.)
        call inpXML%Set_Val(helioVarsIn(VRIN), "helio/vrin", 400. )
        call inpXML%Set_Val(helioVarsIn(VRKFIN), "helio/vrkin", 400. )
        call inpXML%Set_Val(swmodel, "prob/swmodel", "monopole" )
        write(*,"(1x,A14,6F)") "Helio Vars: ", helioVarsIn
        ! R_gc(-3:0); Grid%is = 1; Grid%isg = -3
        ! calculate radii of ghost cells
        !do i = -3, 0
        do i = Grid%isg, Grid%is - 1
            R_gc(i) = norm2(Grid%xyzcc(i, 1, 1, :))
        end do
        write(*,"(1x,A14,4F)") "Rgc: ", R_gc

        ! ibcVars(1:4, ....)
        do kg = 1, Grid%Nkp*Grid%NumRk
            k = modulo(kg, Grid%Nkp)
            do jg = 1, Grid%Njp*Grid%NumRj
                j = modulo(jg, Grid%Njp)
                do ig = 1, Model%Ng
                    i = ig - Model%Ng

                    xyz = Grid%xyzcc(i, j, k, :)
                    R = norm2(xyz)
                    Theta = acos(xyz(ZDIR)/R)
                    !calculate ratio to set 1/r^2 fall of density and br
                    ratio_sq = R_gc(Grid%is-1)*R_gc(Grid%is-1)/R_gc(i)/R_gc(i)

                    if( swmodel == "dipole") then
                        ibcVars(ig, jg, kg, BRIN) = helioVarsIn(BRIN)/Model%Units%gB0*ratio_sq*cos(Theta)/abs(cos(Theta))
                        ibcVars(ig, jg, kg, BRKFIN) = helioVarsIn(BRKFIN)/Model%Units%gB0*ratio_sq*cos(Theta)/abs(cos(Theta))
                        ibcVars(ig, jg, kg, RHOIN) = helioVarsIn(RHOIN)/Model%Units%gD0*ratio_sq
                    elseif ( swmodel == "monopole") then                    
                        ibcVars(ig, jg, kg, BRIN) = helioVarsIn(BRIN)/Model%Units%gB0*ratio_sq
                        ibcVars(ig, jg, kg, BRKFIN) = helioVarsIn(BRKFIN)/Model%Units%gB0*ratio_sq
                        ibcVars(ig, jg, kg, RHOIN) = helioVarsIn(RHOIN)/Model%Units%gD0*ratio_sq
                    else 
                        write(*,*) "Please specify a valid solar wind model, stoppin..."
                        stop
                    end if
                    ibcVars(ig, jg, kg, VRIN) = helioVarsIn(VRIN)/(1.d-5*Model%Units%gv0)
                    ibcVars(ig, jg, kg, VRKFIN) = helioVarsIn(VRKFIN)/(1.d-5*Model%Units%gv0)
                    ibcVars(ig, jg, kg, TIN) = helioVarsIn(TIN)
                end do
            end do
        end do

        MJD_c = 0.0

    end subroutine sphereSymSW

    !>
    !>
    !>
    subroutine writeMJDcH5Root(Model,Grid,IOVars)
        type(Model_T), intent(in)    :: Model
        type(Grid_T) , intent(in)    :: Grid
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars

        call AddOutVar(IOVars,"MJDc", MJD_c)

    end subroutine writeMJDcH5Root

end module usergamic
