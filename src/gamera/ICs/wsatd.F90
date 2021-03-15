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

    implicit none

    !enum, bind(C)
    !   ! variables passed via innerbc file
    !   ! Br, Vr, Rho, Temperature, Br @ kface, Vr @ kface
    !   enumerator :: BRIN=1,VRIN,RHOIN,TIN,BRKFIN,VRKFIN
    !endenum 

    !For TD
    enum, bind(C)
       ! variables passed via innerbc file
       ! Br, Bp_kface, Bt_jface, Vr, Vt, Vp, Rho, Cs (TIN)
       enumerator :: BRIN=1,VRIN,RHOIN,TIN,BPKFIN,BTJFIN,VTIN,VPIN
    endenum 

    ![EP] data structure for TD
    ![EP] should that be 'type, extends(baseBC_T)' ??
    type :: wsaData_T 

        type(ebTab_T)   :: ebTab
        logical :: doStatic = .true.
        integer :: Nr,Nt,Np !dimensions
        !real(rp), dimension(:,:), allocatable :: X,Y
        real(rp) :: wsaT1,wsaT2 !Times of two data slices
        integer  :: wsaN1,wsaN2 !Indices of two data slices
        real(rp), dimension(:,:,:,:), allocatable :: ibcVarsW1,ibcVarsW2

    end type wsaData_T


    integer, private, parameter :: NVARSIN=6 ! SHOULD be the same as the number of vars in the above enumerator
    !for TD
    integer, parameter :: NVARSINTD = 8
    real(rp), dimension(:,:,:,:), allocatable :: ibcVars

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
    real(rp) :: Rbc = 21.5   

    real(rp) :: Tsolar ! Solar rotation period, defined in apps/helioutils.F90
    
    character(len=strLen) :: wsaFile

    ! use this to fix the Efield at the inner boundary
!    real(rp), allocatable :: inEijk(:,:,:,:)

    ! type for solar wind BC
    type, extends(innerIBC_T) :: SWInnerBC_T

        !Main electric field structures
        real(rp), allocatable, dimension(:,:,:,:) :: inEijk,inExyz

        contains

!        procedure :: doInit => InitIonInner
          ! TODO: this shoudl be made generic (wsa, mas, etc.) How
        procedure :: doBC => wsaBC
    end type SWInnerBC_T

    contains

    subroutine initUser(wsaData,Model,Grid,State,inpXML)
        type(wsaData_T), intent(inout) :: wsaData
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML
        procedure(GasIC_T), pointer :: Wxyz
        procedure(HackStep_T), pointer :: tsHack
        procedure(HackE_T), pointer :: eHack

        integer :: i,j,k,nvar,nr,d
        integer :: Nrr, Nt, Np
        integer :: n1, n2
        real(rp) :: w1, w2

!        if (.not.allocated(inEijk)) allocate(inEijk(1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM))

        ! set units and other thins, like Tsolar
        call setHeliosphere(Model,inpXML,Tsolar)

        ! grab inner 
        call inpXML%Set_Val(wsaFile,"prob/wsaFile","innerbc.h5" )

        ![EP] we need to covert time in seconds to code units so doTSclO=true
        wsaData%ebTab%bStr = wsaFile
        ![EP] read times, convert to times to code units, read grid dimensions
        call rdTab(wsaData%ebTab,inpXML,wsaFile,doTSclO=.true.)
        ![EP] ebTab%N is a number of time steps
        if (wsaData%ebTab%N>1) then
            wsaData%doStatic = .false.
        endif
        write(*,*) wsaData%doStatic
        ![EP] check
        write(*,*) wsaData%ebTab%N, wsaData%ebTab%dNi, wsaData%ebTab%dNj, wsaData%ebTab%dNk
        !!!add a check for steady state and time-dep

        ![EP] diemsions of i-ghost grid
        Nr = wsaData%ebTab%dNi
        Nt = wsaData%ebTab%dNj
        Np = wsaData%ebTab%dNk
        write(*,*) Nrr, Nt, Np

        !allocate ibcVars

        !initialization
        ![EP] TD: find bounding time slices from ebTab file
        call findSlc(wsaData%ebTab,State%time,n1,n2)
        write(*,*) n1, n2

        !read map from Step#n1
        call rdWSAMap(wsaData,n1,wsaData%ibcVarsW1)
        wsaData%wsaN1 = n1
        wsaData%wsaT1 = wsaData%ebTab%times(n1)

        !read map from Step#2
        call rdWSAMap(wsaData,n2,wsaData%ibcVarsW2)
        wsaData%wsaN2 = n2
        wsaData%wsaT2 = wsaData%ebTab%times(n2)

        ![EP]interpolation (a) calculate weights (b) interpolate in time 
        call tCalcWeights(wsaData,State%time,w1,w2)

        !(interp BR) INTERPOLATE ALL VARS
        ibcVars(:,:,:,:) = w1*wsaData%ibcVarsW1(:,:,:,:) + w2*wsaData%ibcVarsW2(:,:,:,:)
        !now we have WSA "map" for a current code time

        ![EP] For restart add boundary cases in time


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
        allocate(SWInnerBC_T        :: Grid%externalBCs(1)%p)
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
                 State%magFlux(i,j,k,IDIR) = ibcVars(Model%Ng,j+Grid%ijkShift(JDIR),k+Grid%ijkShift(KDIR),BRIN)*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
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

    !Inner-I BC for WSA-Gamera
    subroutine wsaBC(wsaData,bc,Model,Grid,State)
      type(wsaData_T), intent(inout) :: wsaData
      class(SWInnerBC_T), intent(inout) :: bc
      type(Model_T), intent(in) :: Model
      type(Grid_T), intent(in) :: Grid
      type(State_T), intent(inout) :: State

      ! local variables
      integer :: i,j,k, kb, ke
      integer :: kg, jg, ig ! global indices (in preparation for MPI)
      integer :: var ! ibcVar variable number
      real(rp) :: a
      real(rp) :: ibcVarsStatic(NVARSIN)
      real(rp) :: R, Theta, Phi
      real(rp) :: Theta_kf, R_kf ! kface
      real(rp), dimension(NVAR) :: conVar, pVar
      real(rp) :: w1, w2, n1, n2
  
      !i-boundaries (IN)
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,k,jg,kg,ke,kb,a,var,xyz,R,Theta,Phi,rHat,phiHat) &
      !$OMP private(ibcVarsStatic,pVar,conVar,xyz0,R_kf,Theta_kf)
      
      if (.not.((State%time >= wsaData%wsaT1) .and. (State%time <= wsaData%wsaT2))) then
         !find bounding slices
         call findSlc(wsaData%ebTab,State%time,n1,n2)
         write(*,*) n1, n2

        !read a map from Step#n1
        call rdWSAMap(wsaData,n1,wsaData%ibcVarsW1)
        wsaData%wsaN1 = n1
        !time from Step#n1
        wsaData%wsaT1 = wsaData%ebTab%times(n1)

        !read a map from Step#2
        call rdWSAMap(wsaData,n2,wsaData%ibcVarsW2)
        wsaData%wsaN2 = n2
        !time from Step#n2
        wsaData%wsaT2 = wsaData%ebTab%times(n2)

        ![EP]interpolation (a) calculate weights (b) interpolate in time 
        call tCalcWeights(wsaData,State%time,w1,w2)

        !interpolate in time between two WSA maps (all vars)
        ibcVars(:,:,:,:) = w1*wsaData%ibcVarsW1(:,:,:,:) + w2*wsaData%ibcVarsW2(:,:,:,:)

      endif



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
               if ( (jg>=Grid%js).and.(jg<=size(ibcVars,2)) ) then
                  do var=1,NVARSINTD
                     ibcVarsStatic(var) = a*ibcVars(ig,jg,kb,var)+(1-a)*ibcVars(ig,jg,ke,var)
                  end do
               end if

               ! do cell centered things for cell-centers only
               if ( (j/=Grid%jeg+1).and.(k/=Grid%keg+1) ) then
                  ! various geometrical quantities for the cell center
                  xyz  = Grid%xyzcc(i,j,k,:)
                  R    = norm2(xyz)
                  Theta = acos(xyz(3)/R)
                  Phi = atan2(xyz(2),xyz(1))
                  rHat = xyz/R
                  phiHat = [-sin(phi),cos(phi),0._rp]

                  ! NOTE, WSA data were already scaled appropriately in the python code
                  ! TODO: save them in the innerbc hdf file and set in helioutils appropriately

                  !Set primitives
                  pVar(VELX:VELZ) = rHat*ibcVarsStatic(VRIN)
                  ! note conversion to my units with B0^2/4pi in the denominator
                  pVar(PRESSURE)  = ibcVarsStatic(RHOIN)*Model%Units%gD0*Kbltz*ibcVarsStatic(TIN)/(Model%Units%gP0)
                  pVar(DEN)       = ibcVarsStatic(RHOIN)

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
               !no scaling for Bt_jface
               !BRIN=1,VRIN,RHOIN,TIN,BPKFIN,BTJFIN,VTIN,VPIN
               Rfactor = Rbc/norm2(Grid%xfc(Grid%is,j,k,:,IDIR))
               State%magFlux(i,j,k,IDIR) = ibcVarsStatic(BRIN)*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
               State%magFlux(i,j,k,JDIR) = ibcVarsStatic(BTJFIN)*Grid%face(i,j,k,JDIR)
               State%magFlux(i,j,k,KDIR) = - 2*PI/Tsolar*R_kf*sin(Theta_kf)/ibcVarsStatic(VRKFIN)*ibcVarsStatic(BRKFIN)*Grid%face(i,j,k,KDIR) &
                                           + ibcVarsStatic(BPKFIN)*Rbc/norm2(Grid%xfc(i,j,k,:,KDIR))*Grid%face(i,j,k,KDIR)
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

    ![EP] reads WSA map for a given time step Step#n in innerbc.h5
    subroutine rdWSAMap(wsaData,n,W)
        type(wsaData_T), intent(inout) :: wsaData
        integer, intent(in) :: n !timeslice
        real(rp), dimension(:,:,:,:), intent(out) :: W
        !real(rp), dimension(:,:,:,:), allocatable :: ibcVars

        integer, parameter :: MAXIOVAR = 50
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars


        !integer, parameter :: NIOVAR = 2
        !type(IOVAR_T), dimension(NIOVAR) :: IOVars
        !integer :: dims(2)
        integer :: i,nvar,dims(3)

        !reading group for this time step n
        write(*,*) 'Reading file/group = ', &
             trim(wsaData%ebtab%bStr),'/',trim(wsaData%ebTab%gStrs(n))

        !Reset IO chain
        call ClearIO(IOVars)

        !Setup input chain
        call AddInVar(IOVars,"vr")
        call AddInVar(IOVars,"vt")
        call AddInVar(IOVars,"vp")
        call AddInVar(IOVars,"rho")
        !SOUND SPEED
        call AddInVar(IOVars,"cs")
        call AddInVar(IOVars,"br")
        call AddInVar(IOVars,"bp_kface")
        call AddInVar(IOVars,"bt_jface")
        call AddInVar(IOVars,"MJD")
        call AddInVar(IOVars,"time")

        call ReadVars(IOVars,.false.,wsaData%ebTab%bStr,wsaData%ebTab%gStrs(n)) 

        dims = [wsaData%Nr,wsaData%Nt,wsaData%Np] !i,j,k
        write(*,*) dims
        !Allocate W
        !Bp_kface has dim 257, 128, 4
        !Bt_jface has dim 256, 129, 4
        if (.not.allocated(W)) allocate(W(dims(1),dims(2),dims(3),NVARSINTD))

        !BRIN=1,VRIN,RHOIN,TIN,BPKFIN,BTJFIN,VTIN,VPIN

        do i=1,NVARSINTD
         select case (i)
         case (BRIN)
            nvar= FindIO(IOVars,"br")
         case (VRIN)
            nvar= FindIO(IOVars,"vr")
         case (RHOIN)
            nvar= FindIO(IOVars,"rho")
         case (TIN)
            nvar= FindIO(IOVars,"cs")
         case (BPKFIN)
            nvar= FindIO(IOVars,"bp_kface")
         case (BTJFIN)
            nvar= FindIO(IOVars,"bt_jface")
         case (VTIN)
            nvar= FindIO(IOVars,"vt")
         case (VPIN)
            nvar= FindIO(IOVars,"vp")
         end select

         W(:,:,:,i) = reshape(IOVars(nvar)%data,dims)
        end do
         !reading MJD
         MJD_c = GetIOReal(IOVars,"MJD")
        
    end subroutine rdEmpMap

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
                w1 = (wsaData%empT2-t)/dt
                w2 = (t - wsaData%empT1)/dt
            else
                w1 = 0.5
                w2 = 0.5
            endif
        endif !Weights
    end subroutine tCalcWeights

end module usergamic
