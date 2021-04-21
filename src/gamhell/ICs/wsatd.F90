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
    use helioutils
    ![EP] for TD 
    use ebinit
    use ebtypes
    use volttypes
    use chmpfields

    implicit none

    !enum, bind(C)
    !   ! variables passed via innerbc file
    !   ! Br, Vr, Rho, Temperature, Br @ kface, Vr @ kface
    !   enumerator :: BRIN=1,VRIN,RHOIN,TIN,VRKFIN,VRKFIN
    !endenum 

    !For TD
    enum, bind(C)
       ! variables passed via innerbc file
       ! Br, Bp_kface, Bt_jface, Vr, Vt, Vp, Rho, Cs (TIN)
       enumerator :: BRIN=1,VRIN,VTIN,VPIN,RHOIN,TIN,BTJFIN,BPKFIN,BRKFIN,VRKFIN
    endenum 

    integer, private, parameter :: NVARSIN=6 ! SHOULD be the same as the number of vars in the above enumerator
    !for TD
    integer, parameter :: NVARSINTD = 8
    integer, parameter :: NVARSCC = 6 
    !real(rp), dimension(:,:,:,:), allocatable :: ibcVars

    !Various global would go here
    real (rp) :: Rho0, P0, Vslow, Vfast, wScl, Cs0, B0, MJD_c

    !Scaling from innerbc to CGS
    real (rp) :: km2cm = 1.e5
    real (rp) :: nT2Gs = 1.e-5

    ! things we keep reusing
    real(rp), dimension(NDIM) :: xyz,xyz0,rHat,phiHat
    real(rp) :: Rfactor

    ! global grid
    integer :: gNkp

    ! FIXME
    ! setting it here temporarily. eventually need to read from HDF or something
    ! note also that this is slightly incorrect, since Rbc below is used as the radius
    ! of the center of the first ghost cell.
    real(rp) :: Rbc = 21.5   

    real(rp) :: Tsolar ! Solar rotation period, defined in apps/helioutils.F90
    
    character(len=strLen) :: wsaFile

    ! use this to fix the Efield at the inner boundary
    !real(rp), allocatable :: inEijk(:,:,:,:)
    type :: ibcMap_T
        real(rp), dimension(:,:,:,:), allocatable :: ccVars !6 cell-center vars
        real(rp), dimension(:,:,:), allocatable :: jfVars !1 j-face
        real(rp), dimension(:,:,:,:), allocatable :: kfVars !3 k-face Bp, Br, Vr
       !real(rp), dimension(:,:,:), allocatable :: Rho, Temp, Vr, Vt, Vp, Br, Bt_jf, Bp_kf, Br_kf, Vr_kf
    !   !6 cell-centered variables Rho, Temp, Vr, Vt, Vp, Br have dimensions (Ni, Nt, Nk)
    !   !Bt at j-faces (Ni, Nt+1, Nk)
    !   !Bp, Br, Vr at k-faces (Ni, Nt, Nk+1)
    end type ibcMap_T

    ![EP] data structure for TD
    type :: wsaData_T 

        type(ebTab_T)   :: ebTab
        logical :: doStatic = .true.
        integer :: Nr,Nt,Np !dimensions
        !real(rp), dimension(:,:), allocatable :: X,Y
        real(rp) :: wsaT1,wsaT2 !Times of two data slices
        integer  :: wsaN1,wsaN2 !Indices of two data slices
        !vars from time slice in innerbc.h5
        type(ibcMap_T) :: ibcMap1
        type(ibcMap_T) :: ibcMap2

    end type wsaData_T

    ! type for solar wind BC
    type, extends(innerIBC_T) ::WSAInnerBC_T

        type(wsaData_T) :: wsaData
        !Main electric field structures
        real(rp), allocatable, dimension(:,:,:,:) :: inEijk,inExyz

        contains

        procedure :: doInit => InitwsaBC
          ! TODO: this shoudl be made generic (wsa, mas, etc.) How
        procedure :: doBC => wsaBC
    end type WSAInnerBC_T

    contains

    subroutine initUser(Model,Grid,State,inpXML)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML
        procedure(GasIC_T), pointer :: Wxyz
        procedure(HackStep_T), pointer :: tsHack
        procedure(HackE_T), pointer :: eHack

        real(rp), dimension(:,:,:), allocatable :: ibcBr


        integer :: i,j,k,nvar,nr,d


!        if (.not.allocated(inEijk)) allocate(inEijk(1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM))

        ! set units and other things, like Tsolar
        call setHeliosphere(Model,inpXML,Tsolar)

        ! grab inner 
        call inpXML%Set_Val(wsaFile,"prob/wsaFile","innerbc.h5" )

        ![EP] get the WSA map for current State%time
        ibcBr = 0.0
        call initTSlice(wsaFile,inpXML,Model,State,ibcBr)
        !now we have WSA "map" for a current code time
        ![EP] TODO: For restart add boundary cases in time

        ! compute global Nkp
        gNkp = Grid%Nkp*Grid%NumRk

        ! initial conditions
        ! TODO: change using norm. units in Model, set in helioutils
        Cs0   = 0.267  ! 40 km/s
        Vslow = 1.33   ! 200 km/s
        Vfast = 5.33   ! 800 km/s
        B0    = 2.0    ! 200 nT
        Rho0  = 1.0    ! 200/cc
        P0    = 1.0e-4*Rho0*Cs0**2.0/Model%gamma

        ! deallocate default BCs
        ! required because defaults are triply-periodic
        ! need to wipe
        call WipeBCs(Model,Grid)

        !Set BCs
        allocate(WSAInnerBC_T        :: Grid%externalBCs(1)%p)
        allocate(helioOuterIBC_T    :: Grid%externalBCs(2)%p)
        allocate(helioInnerJBC_T    :: Grid%externalBCs(3)%p)
        allocate(helioOuterJBC_T    :: Grid%externalBCs(4)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(5)%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(6)%p)

        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

        ! Add gravity
!        tsHack => PerStep
!        Model%HackStep => tsHack
!        eHack  => EFix
!        Model%HackE => eHack

        ! everybody reads WSA data
        !call readIBC(wsaFile)

        Model%MJD0 = MJD_c

        !Map IC to grid
        Wxyz => GasIC
        call GasIC2State(Model,Grid,State,Wxyz)

        ! NOTE, not filling ghost cells here
        ! relying on J, K boundary conditions
        ! first zero everything out
        State%magFlux = 0.0

        do k=Grid%ks,Grid%ke
           do j=Grid%js,Grid%je
              do i=Grid%isg,Grid%ieg+1
                 ! FIXME: calc Rbc appropriately above!
                 Rfactor = Rbc/norm2(Grid%xfc(Grid%is,j,k,:,IDIR))
                 ! note scaling br to the first active face
                 ! this is opposite to what we did on the python side
                 ! not elegant, but otherwise we'd have to store both br and br_iface in teh innerbc.h5 file
                 !
                 ! note also that ibcVars(Model%Ng,:,:,:) corresponds to the center of the first ghost cell (just below the boundary)
                 State%magFlux(i,j,k,IDIR) = ibcBr(Model%Ng,j+Grid%ijkShift(JDIR),k+Grid%ijkShift(KDIR))*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
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

    end subroutine initUser

    !Initialization for WSA BCs
    ! Here we get times from innerbc, ghost grid dimensions (rdTab func)
    ! Determine time bracket, reading WSA maps from those times
    ! Fill out ibcVars
    subroutine InitwsaBC(bc,Model,Grid,State,xmlInp)
      class(WSAInnerBC_T), intent(inout) :: bc
      type(Model_T), intent(inout) :: Model
      type(Grid_T), intent(in) :: Grid
      type(State_T), intent(in) :: State
      type(XML_Input_T), intent(in) :: xmlInp

      real(rp) :: w1, w2
      integer :: n1, n2
      integer :: Nr, Nt, Np
      integer :: dims(3)

      !call inpXML%Set_Val(wsaFile,"prob/wsaFile","innerbc.h5" )

      bc%wsaData%ebTab%bStr = wsaFile

      ![EP] read times, convert to times to code units, read grid dimensions
      call rdTab(bc%wsaData%ebTab,xmlInp,wsaFile,doTSclO=.false.)
      ![EP] ebTab%N is a number of time steps
      if (bc%wsaData%ebTab%N>1) then
            bc%wsaData%doStatic = .false.
      endif
      write(*,*) '[EP in InitwsaBC] doStatic = ', bc%wsaData%doStatic
      ![EP] check
      write(*,*) '[EP in InitwsaBC] rdTab check ', bc%wsaData%ebTab%N, bc%wsaData%ebTab%dNi, bc%wsaData%ebTab%dNj, bc%wsaData%ebTab%dNk
      !!!add a check for steady state and time-dep
      write(*,*) '[EP in InitwsaBC] times in innerbc ', bc%wsaData%ebTab%times

      ![EP] dimensions of i-ghost grid
      bc%wsaData%Nr = bc%wsaData%ebTab%dNi
      bc%wsaData%Nt = bc%wsaData%ebTab%dNj
      bc%wsaData%Np = bc%wsaData%ebTab%dNk
      write(*,*) '[EP in InitwsaBC] Dimensions ', bc%wsaData%Nr, bc%wsaData%Nt, bc%wsaData%Np
      write(*,*) '[EP in InitwsaBC] State%time ', State%time

      !initialization
      ![EP] TD: find bounding time slices from ebTab file
      call findSlc(bc%wsaData%ebTab,State%time*Model%Units%gT0,n1,n2)
      write(*,*) '[EP in InitwsaBC] Bounding slices ', n1, n2

      !read map from Step#n1
      call rdWSAMap(bc%wsaData,Model,n1,bc%wsaData%ibcMap1)
      bc%wsaData%wsaN1 = n1
      bc%wsaData%wsaT1 = bc%wsaData%ebTab%times(n1) ! in sec

      !read map from Step#2
      call rdWSAMap(bc%wsaData,Model,n2,bc%wsaData%ibcMap2)
      bc%wsaData%wsaN2 = n2
      bc%wsaData%wsaT2 = bc%wsaData%ebTab%times(n2) !in sec
      write(*,*)'[EP in InitwsaBC] Bounding times ', bc%wsaData%wsaT1, bc%wsaData%wsaT2

      !![EP]interpolation (a) calculate weights (b) interpolate in time 
      !call tCalcWeights(bc%wsaData,State%time*Model%Units%gT0,w1,w2)
      !write(*,*) '[EP in InitwsaBC] Bounding weights ', w1, w2

      !dims = [bc%wsaData%Nr,bc%wsaData%Nt+1,bc%wsaData%Np+1] !i,j,k
      !write(*,*)'[EP in InitwsaBC] dims in initSlice for ibcVars', dims
      !if (.not.allocated(ibcVars)) allocate(ibcVars(dims(1),dims(2),dims(3),NVARSINTD))

      !!INTERPOLATE ALL VARS that we get from innerbc.h5
      !ibcVars(:,:,:,1:VPIN) = w1*bc%wsaData%ibcVarsW1(:,:,:,:) + w2*bc%wsaData%ibcVarsW2(:,:,:,:)
      !write(*,*)'[EP in InitwsaBC] filled out ibcVars'  

    end subroutine InitwsaBC


    !Inner-I BC for WSA-Gamera
    subroutine wsaBC(bc,Model,Grid,State)
      class(WSAInnerBC_T), intent(inout) :: bc
      type(Model_T), intent(in) :: Model
      type(Grid_T), intent(in) :: Grid
      type(State_T), intent(inout) :: State

      type(ibcMap_T) :: ibcVars

      ! local variables
      integer :: i,j,k, kb, ke
      integer :: kg, jg, ig ! global indices (in preparation for MPI)
      integer :: var, vark ! ibcVar variable number
      real(rp) :: a
      real(rp) :: ibcVarsStatic(10)
      real(rp) :: R, Theta, Phi
      real(rp) :: Theta_kf, R_kf ! kface
      real(rp), dimension(NVAR) :: conVar, pVar
      real(rp) :: w1, w2
      integer :: n1, n2
      integer :: Nr, Nt, Np
      integer :: dims_cc(3), dims_jf(3), dims_kf(3)

      !write(*,*)'[EP in wsaBC] Begin. Time: ', State%time*Model%Units%gT0, bc%wsaData%wsaT1, bc%wsaData%wsaT2
  
      if (.not.((State%time*Model%Units%gT0 >= bc%wsaData%wsaT1) .and. (State%time*Model%Units%gT0 <= bc%wsaData%wsaT2))) then
         !find bounding slices
         call findSlc(bc%wsaData%ebTab,State%time*Model%Units%gT0,n1,n2)
         write(*,*) '[EP in wsaBC] getting in a new time backet', n1, n2

        !n1=1
        !n2=2
        !read a map from Step#n1
        call rdWSAMap(bc%wsaData,Model,n1,bc%wsaData%ibcMap1)
        bc%wsaData%wsaN1 = n1
        !time from Step#n1
        bc%wsaData%wsaT1 = bc%wsaData%ebTab%times(n1)
        

        !read a map from Step#2
        call rdWSAMap(bc%wsaData,Model,n2,bc%wsaData%ibcMap2)
        bc%wsaData%wsaN2 = n2
        !time from Step#n2
        bc%wsaData%wsaT2 = bc%wsaData%ebTab%times(n2)
        
      endif

      ![EP]interpolation (a) calculate weights (b) interpolate in time 
      call tCalcWeights(bc%wsaData,State%time*Model%Units%gT0,w1,w2)
      write(*,*) '[EP in wsaBC] Bounding weights ', w1, w2
      
      !Always first map
      !w1 = 1.
      !w2 = 0.
      dims_cc = [bc%wsaData%Nr,bc%wsaData%Nt,bc%wsaData%Np] !i,j,k
      dims_jf = [bc%wsaData%Nr,bc%wsaData%Nt+1,bc%wsaData%Np]
      dims_kf = [bc%wsaData%Nr,bc%wsaData%Nt,bc%wsaData%Np+1]

      ! allocate ibcVars 
      if (.not.allocated(ibcVars%ccVars)) allocate(ibcVars%ccVars(dims_cc(1),dims_cc(2),dims_cc(3),NVARSCC))
      if (.not.allocated(ibcVars%jfVars)) allocate(ibcVars%jfVars(dims_jf(1),dims_jf(2),dims_jf(3)))
      if (.not.allocated(ibcVars%kfVars)) allocate(ibcVars%kfVars(dims_kf(1),dims_kf(2),dims_kf(3),3))

      ! interpolate in time
      ibcVars%ccVars(:,:,:,:) = w1*bc%wsaData%ibcMap1%ccVars(:,:,:,:) + w2*bc%wsaData%ibcMap2%ccVars(:,:,:,:)
      ibcVars%jfVars(:,:,:) = w1*bc%wsaData%ibcMap1%jfVars(:,:,:) + w2*bc%wsaData%ibcMap2%jfVars(:,:,:)
      ibcVars%kfVars(:,:,:,1) = w1*bc%wsaData%ibcMap1%kfVars(:,:,:,1) + w2*bc%wsaData%ibcMap2%kfVars(:,:,:,1)
      !write(*,*)'[EP in wsaBC] filled ibcVars '
      
      !To obtain Vr and Br values at k-faces do simple interpolation in a uniform grid
      !Br_kface for kfaces 1:Nk
      !TODO: Make enumerator for k-Face Vars 1 - BpK, 2 - BrK, 3 - VrK
      ibcVars%kfVars(:,:,2:bc%wsaData%Np,2) = 0.5*(ibcVars%ccVars(:,:,1:bc%wsaData%Np-1,BRIN)+ ibcVars%ccVars(:,:,2:bc%wsaData%Np,BRIN))
      !boundary faces k=1 and Np+1
      !periodic conditions
      ibcVars%kfVars(:,:,1,2) = 0.5*(ibcVars%ccVars(:,:,1,BRIN) + ibcVars%ccVars(:,:,bc%wsaData%Np,BRIN))
      ibcVars%kfVars(:,:,bc%wsaData%Np+1,2) = ibcVars%kfVars(:,:,1,2)

      !Vr_kface
      ibcVars%kfVars(:,:,2:bc%wsaData%Np,3) = 0.5*(ibcVars%ccVars(:,:,1:bc%wsaData%Np-1,VRIN)+ ibcVars%ccVars(:,:,2:bc%wsaData%Np,VRIN))
      ibcVars%kfVars(:,:,1,3) = 0.5*(ibcVars%ccVars(:,:,1,VRIN) + ibcVars%ccVars(:,:,bc%wsaData%Np,VRIN))
      ibcVars%kfVars(:,:,bc%wsaData%Np+1,3) = ibcVars%kfVars(:,:,1,3)

      !write(*,*)'[EP in wsaBC] Added VRKFIN and BRKFIN to ibcVars. Beginning a big loop. '
      !i-boundaries (IN)
      !!$OMP PARALLEL DO default(shared) &
      !!$OMP private(i,j,k,jg,kg,ke,kb,a,vark,var,xyz,R,Theta,Phi,rHat,phiHat) &
      !!$OMP private(ibcVarsStatic,pVar,conVar,xyz0,R_kf,Theta_kf) &
      !!$OMP private(w1, w2, n1, n2, Nr, Nt, Np) 
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
               !if ( (jg>=Grid%js).and.(jg<=size(ibcVars%ccVars,2) )) then
               if (.not.((Grid%hasUpperBC(JDIR)==.True.).and.(j==Grid%je+1))) then 
                  ibcVarsStatic(1:NVARSCC) = a*ibcVars%ccVars(ig,jg,kb,:)+(1-a)*ibcVars%ccVars(ig,jg,ke,:)
                  ibcVarsStatic(BPKFIN:VRKFIN) = a*ibcVars%kfVars(ig,jg,kb,:)+(1-a)*ibcVars%kfVars(ig,jg,ke,:)
                endif

                ibcVarsStatic(BTJFIN) = a*ibcVars%jfVars(ig,jg,kb)+(1-a)*ibcVars%jfVars(ig,jg,ke)

               ! do cell centered things for cell-centers only
               if ( (j/=Grid%jeg+1).and.(k/=Grid%keg+1) ) then
                  ! various geometrical quantities for the cell center
                  xyz  = Grid%xyzcc(i,j,k,:)
                  R    = norm2(xyz)
                  Theta = acos(xyz(3)/R)
                  Phi = atan2(xyz(2),xyz(1))
                  rHat = xyz/R
                  phiHat = [-sin(phi),cos(phi),0._rp]

                  ! [VGM] NOTE, WSA data were already scaled appropriately in the python code
                  ! [VGM] TODO: save them in the innerbc hdf file and set in helioutils appropriately

                  !Set primitives
                  pVar(VELX:VELZ) = rHat*ibcVarsStatic(VRIN)
                   !pVar(VELX:VELZ) = rHat*3.0
                  ! note conversion to my units with B0^2/4pi in the denominator
                  pVar(PRESSURE)  = ibcVarsStatic(RHOIN)*Model%Units%gD0*Kbltz*ibcVarsStatic(TIN)/(Model%Units%gP0)
                  !pVar(PRESSURE)  = 5.* Model%Units%gD0*Kbltz*1.e6/(Model%Units%gP0)
                  pVar(DEN)      = ibcVarsStatic(RHOIN)
                   !pVar(DEN)      = 5.
                  !Swap prim->con in ghost variables
                  call CellP2C(Model,pVar,conVar)
                  State%Gas(i,j,k,:,BLK) = conVar

                  ! note, don't need cc Bxyz because we run flux2field through ghosts
               end if

               if (i==0.and.j==64.and.k==64) then
                  write(*,*) '[EP in wsaBC] ibcVarsStatic', ibcVarsStatic(:)
               end if

               ! also need theta at k-face for k-flux
               ! although we're assuming theta and R don't change from cell to k-face center,
               ! xyzcc used under if statement above is only defined for cell centers so need to define it here
               xyz0 = Grid%xfc(i,j,k,:,KDIR) !just reusing a temp var here.
               R_kf = norm2(xyz0)
               Theta_kf = acos(xyz0(ZDIR)/R_kf)

               ! note scaling for Bi. See note above in InitUser
               !no scaling for Bt_jface
               !BRIN=1,VRIN,RHOIN,TIN,BPKFIN,BTJFIN,VTIN,VPIN
               Rfactor = Rbc/norm2(Grid%xfc(Grid%is,j,k,:,IDIR))
               ![EP] at all i-faces flux should be the same. We re-normalize from cc to face
               !FIX that
               State%magFlux(i,j,k,IDIR) = ibcVarsStatic(BRIN)*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
               !State%magFlux(i,j,k,IDIR) = 1.*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
               State%magFlux(i,j,k,JDIR) = ibcVarsStatic(BTJFIN)*Grid%face(i,j,k,JDIR)
               !State%magFlux(i,j,k,JDIR) = 0.0
               !multiplying by Rbc is not correct. Need to multiply to r_cc
               State%magFlux(i,j,k,KDIR) = - 2*PI/Tsolar*R_kf*sin(Theta_kf)/ibcVarsStatic(VRKFIN)*ibcVarsStatic(BRKFIN)*Grid%face(i,j,k,KDIR) &
                                           + ibcVarsStatic(BPKFIN)*Rbc/norm2(Grid%xfc(i,j,k,:,KDIR))*Grid%face(i,j,k,KDIR)
               !State%magFlux(i,j,k,KDIR) = - 2*PI/Tsolar*R_kf*sin(Theta_kf)/ibcVarsStatic(VRKFIN)*ibcVarsStatic(BRKFIN)*Grid%face(i,j,k,KDIR)
                !State%magFlux(i,j,k,KDIR) = 0.0
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

        kprime = modulo(k-gNkp*State%time/Tsolar,real(gNkp,kind=rp)) 

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


    ! !Do per timestep, includes lazy gravitational force term
    ! subroutine PerStep(Model,Gr,State)
    !     type(Model_T), intent(in) :: Model
    !     type(Grid_T), intent(inout) :: Gr
    !     type(State_T), intent(inout) :: State

    !     integer :: i,j,k

    !     real(rp), dimension(NDIM) :: xyz, Vxyz, rHat
    !     real(rp), dimension(NVAR) :: pW,pCon
    !     real(rp) :: D,IntE,r
    !     real(rp) :: GM0

    !     !Scaling for gravitational force
    !     GM0 = UN/(UL**3*UB**2)*6.67408*1.99*4*pi*1.67/6.955/10  ! 2.74e4cm/s^2

    !     !Add grav force
    !     !$OMP PARALLEL DO default(shared) &
    !     !$OMP private(i,j,k,xyz,rHat,Vxyz,pW,pCon,r,D,IntE)
    !     do k=Gr%ksg,Gr%keg
    !         do j=Gr%jsg,Gr%jeg
    !             do i=Gr%isg,Gr%ieg
    !                 xyz = Gr%xyzcc(i,j,k,:)
    !                 r = norm2(xyz)
    !                 rHat = xyz/r

    !                 pCon = State%Gas(i,j,k,:,BLK)
    !                 call CellC2P(Model,pCon,pW)

    !                 D = pW(DEN)
    !                 IntE = pW(PRESSURE)/(Model%gamma-1)
    !                 Vxyz = pW(VELX:VELZ)
    !                 Vxyz = Vxyz - Model%dt*GM0*rHat/(r*r)

    !                 !Reset conserved State
    !                 pCon(DEN) = D
    !                 pCon(MOMX:MOMZ) = D*Vxyz
    !                 pCon(ENERGY) = IntE + 0.5*D*dot_product(Vxyz,Vxyz)
    !                 State%Gas(i,j,k,:,BLK) = pCon
    !             enddo
    !         enddo
    !     enddo
    ! end subroutine PerStep

    subroutine eFix(Model,Gr,State)
      type(Model_T), intent(in) :: Model
      type(Grid_T), intent(in) :: Gr
      type(State_T), intent(inout) :: State

      ! see example of how to do this in voltron/ICs/earthcmi.F90
       !Grid%externalBCs(1)%p
!      State%Efld(Gr%is  ,:,:,JDIR:KDIR) = inEijk(1,:,:,JDIR:KDIR)*Gr%edge(Gr%is  ,:,:,JDIR:KDIR)


    end subroutine eFix

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
      !if (.not.allocated(ibcVars)) allocate(ibcVars(dims(1),dims(2),dims(3),NVARSIN))

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

         !ibcVars(:,:,:,i) = reshape(IOVars(nvar)%data,dims)
      end do
         !reading modified julian date from innerbc
         MJD_c = GetIOReal(IOVars,"MJD")
   
    end subroutine readIBC

    subroutine initTSlice(wsaFile,inpXML,Model,State,ibcBr)
       type(XML_Input_T), intent(in) :: inpXML
       type(Model_T), intent(in) :: Model
       type(State_T), intent(in) :: State
       character(len=strLen), intent(in) :: wsaFile
       real(rp), dimension(:,:,:), intent(out), allocatable :: ibcBr

       type(wsaData_T) :: wsaData
       integer :: n1, n2
       real(rp) :: w1, w2
       integer :: dims(3)


       wsaData%ebTab%bStr = wsaFile

       ![EP] read times, convert to times to code units, read grid dimensions
       call rdTab(wsaData%ebTab,inpXML,wsaFile,doTSclO=.false.)
       ![EP] ebTab%N is a number of time steps
        if (wsaData%ebTab%N>1) then
            wsaData%doStatic = .false.
        endif
        ![EP] debug
        write(*,*) '[EP in initTSlice] doStatic = ', wsaData%doStatic
        write(*,*) '[EP in initTSlice] rdTab check ', wsaData%ebTab%N, wsaData%ebTab%dNi, wsaData%ebTab%dNj, wsaData%ebTab%dNk
        !!!add a check for steady state and time-dep
        write(*,*) '[EP in initTSlice] times in innerbc ', wsaData%ebTab%times

        ![EP] dimensions of i-ghost grid
        wsaData%Nr = wsaData%ebTab%dNi
        wsaData%Nt = wsaData%ebTab%dNj
        wsaData%Np = wsaData%ebTab%dNk
        write(*,*) '[EP in initTSlice] Dimensions ', wsaData%Nr, wsaData%Nt, wsaData%Np
        write(*,*) '[EP in initTSlice] State%time ', State%time
        !TODO allocate ibcVars

        !initialization
        ![EP] TD: find bounding time slices from ebTab file
        call findSlc(wsaData%ebTab,State%time*Model%Units%gT0,n1,n2)
        write(*,*) '[EP in initTSlice] Bounding slices ', n1, n2

        !read map from Step#n1
        call rdWSAMap(wsaData,Model,n1,wsaData%ibcMap1)
        wsaData%wsaN1 = n1
        wsaData%wsaT1 = wsaData%ebTab%times(n1)

        !read map from Step#2
        call rdWSAMap(wsaData,Model,n2,wsaData%ibcMap2)
        wsaData%wsaN2 = n2
        wsaData%wsaT2 = wsaData%ebTab%times(n2)
        write(*,*)'[EP in initTSlice] Bounding times ', wsaData%wsaT1, wsaData%wsaT2 

        ![EP]interpolation (a) calculate weights (b) interpolate in time 
        call tCalcWeights(wsaData,State%time*Model%Units%gT0,w1,w2)
         write(*,*) '[EP in initTSlice] Bounding weights ', w1, w2
   

        dims = [wsaData%Nr,wsaData%Nt,wsaData%Np] !i,j,k cell-centered
        write(*,*)'[EP in initTSlice] dims in initSlice for ibcBr', dims
        if (.not.allocated(ibcBr)) allocate(ibcBr(dims(1),dims(2),dims(3)))

        !INTERPOLATE BRIN that we get from innerbc.h5
        ibcBr(:,:,:) = w1*wsaData%ibcMap1%ccVars(:,:,:,BRIN) + w2*wsaData%ibcMap2%ccVars(:,:,:,BRIN)
        !write(*,*)'[EP in initTSlice] Filled out ibcBr, checking ibcBr(4,:,1), ibc(4,64,:)', ibcBr(4,:,1), ibcBr(4,64,:)

        write(*,*)'[EP in initTSlice] filled out ibcBr'

    end subroutine initTSlice


    ![EP] reads WSA map for a given time step Step#n in innerbc.h5
    subroutine rdWSAMap(wsaData,Model,n,ibcMap)
        type(wsaData_T), intent(inout) :: wsaData
        type(Model_T), intent(in) :: Model
        integer, intent(in) :: n !timeslice
        !real(rp), dimension(:,:,:,:), intent(out), allocatable :: W
        
        integer, parameter :: MAXIOVAR = 50
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        type(ibcMap_T),intent(out) :: ibcMap


        !integer, parameter :: NIOVAR = 2
        !type(IOVAR_T), dimension(NIOVAR) :: IOVars
        !integer :: dims(2)
        integer :: i,nvar, dims_cc(3), dims_jf(3), dims_kf(3)

        !reading group for this time step n
        write(*,*) 'Reading file/group = ', &
             trim(wsaData%ebtab%bStr),'/',trim(wsaData%ebTab%gStrs(n))

        !Reset IO chain
        call ClearIO(IOVars)

        !Setup input chain
        call AddInVar(IOVars,"br",vSclO=nT2Gs/Model%Units%gB0)
        call AddInVar(IOVars,"vr",vSclO=km2cm/Model%Units%gv0)
        call AddInVar(IOVars,"vt",vSclO=km2cm/Model%Units%gv0)
        call AddInVar(IOVars,"vp",vSclO=km2cm/Model%Units%gv0)
        call AddInVar(IOVars,"rho",vSclO=1/Model%Units%gD0)
        call AddInVar(IOVars,"T") !Temp in K 
        call AddInVar(IOVars,"bt_jface",vSclO=nT2Gs/Model%Units%gB0)
        call AddInVar(IOVars,"bp_kface",vSclO=nT2Gs/Model%Units%gB0)
        call AddInVar(IOVars,"MJD")
        call AddInVar(IOVars,"time",vSclO=1/Model%Units%gT0)

        call ReadVars(IOVars,.false.,wsaData%ebTab%bStr,wsaData%ebTab%gStrs(n)) 

        !do i=1,10
        !   write(*,*)'[EP in rdWSAMap]  IOVars%dims', IOVars(i)%dims(1:3)
        !end do 

        !add +1 in j and k because of field components at faces
        !Bp_kface has dim 257, 128, 4
        !Bt_jface has dim 256, 129, 4
        dims_cc = [wsaData%Nr,wsaData%Nt,wsaData%Np] !i,j,k
        dims_jf = [wsaData%Nr,wsaData%Nt+1,wsaData%Np]
        dims_kf = [wsaData%Nr,wsaData%Nt,wsaData%Np+1]
        !dims=IOVars(1)%dims(1:3) ! i,j,k
        write(*,*) '[EP in rdWSAMap] Dimensions of W in rdWSAMap', dims_cc, dims_jf, dims_kf
        !Allocate W
        if (.not.allocated(ibcMap%ccVars)) allocate(ibcMap%ccVars(dims_cc(1),dims_cc(2),dims_cc(3),NVARSCC))
        if (.not.allocated(ibcMap%jfVars)) allocate(ibcMap%jfVars(dims_jf(1),dims_jf(2),dims_jf(3)))
        if (.not.allocated(ibcMap%kfVars)) allocate(ibcMap%kfVars(dims_kf(1),dims_kf(2),dims_kf(3),3))

        !zero out
        ibcMap%ccVars = 0.
        ibcMap%jfVars = 0.
        ibcMap%kfVars = 0.

        !BRIN=1,VRIN,VTIN,VPIN,RHOIN,TIN,BTJFIN,BPKFIN,BRKFIN,VRKFIN  
        !NVARSINTD=8

        do i=1,NVARSINTD
         select case (i)
         case (BRIN)
            nvar= FindIO(IOVars,"br")
         case (VRIN)
            nvar= FindIO(IOVars,"vr")
         case (VTIN)
            nvar= FindIO(IOVars,"vt")
         case (VPIN)
            nvar= FindIO(IOVars,"vp")
         case (RHOIN)
            nvar= FindIO(IOVars,"rho")
         case (TIN)
            nvar= FindIO(IOVars,"T")
         case (BTJFIN)
            nvar= FindIO(IOVars,"bt_jface")
         case (BPKFIN)
            nvar= FindIO(IOVars,"bp_kface")
         end select

         !Scale from innerbc units (km/s, cm-3, nT, K) to code dim-less units
         if (i<=NVARSCC) then
            ibcMap%ccVars(:,:,:,i) = reshape(IOVars(nvar)%data*IOVars(nvar)%scale,dims_cc)
         else if (i == BTJFIN) then
            ibcMap%jfVars(:,:,:) = reshape(IOVars(nvar)%data*IOVars(nvar)%scale,dims_jf)
         else if (i == BPKFIN) then
            ibcMap%kfVars(:,:,:,1) = reshape(IOVars(nvar)%data*IOVars(nvar)%scale,dims_kf)
         end if

         !W(:,:,:,i) = reshape(IOVars(nvar)%data*IOVars(nvar)%scale,dims)
       end do

         !reading MJD
         MJD_c = GetIOReal(IOVars,"MJD")
        
    end subroutine rdWSAMap

    subroutine tCalcWeights(wsaData,t,w1,w2)
        type(wsaData_T), intent(inout) :: wsaData
        real(rp), intent(in)  :: t
        real(rp), intent(out) :: w1,w2
        real(rp) :: dt
        if (wsaData%doStatic) then
            w1 = 1.0
            w2 = 0.0
            return
        endif

        if (t > wsaData%wsaT2) then
            w1 = 0.0
            w2 = 1.0
        else if (t < wsaData%wsaT1) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = wsaData%wsaT2-wsaData%wsaT1
            if (dt>TINY) then
                w1 = (wsaData%wsaT2-t)/dt
                w2 = (t - wsaData%wsaT1)/dt
            else
                w1 = 0.5
                w2 = 0.5
            endif
        endif !Weights
    end subroutine tCalcWeights

end module usergamic
