! Various utilities for heliospheric runs

module helioutils
    use gamtypes
    use gamutils
    use math
    use gridutils
    use output

    implicit none

    ! normalization
    real(rp), private :: gD0, gB0, gx0, gT0, gv0, gP0

    type, extends(baseBC_T) :: helioInnerJBC_T
        contains
        procedure :: doBC => helio_ibcJ
    end type helioInnerJBC_T

    type, extends(baseBC_T) :: helioOuterJBC_T
        contains
        procedure :: doBC => helio_obcJ
    end type helioOuterJBC_T

    type, extends(baseBC_T) :: helioOuterIBC_T
        contains
        procedure :: doBC => helio_obcI
    end type helioOuterIBC_T

    contains

      subroutine setHeliosphere(Model,inpXML,Tsolar)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: inpXML
        real(rp),intent(out) :: Tsolar ! Solar rotation period

        ! normalization
        gD0=200.         ! [/cc]
        gB0=1.e-3        ! [Gs], 100 nT
        gx0=Rsolar*1.e5  ! [cm], solar radius

        ! get the necessary units 
        gv0 = gB0/sqrt(4*pi*gD0*mp_cgs) ! [cm/s] ~ 154km/s for gD0=200. and gB0 = 1.e-3
        gT0 = gx0/gv0                   ! [s] ~ 1.25 hour for above values
        gP0 = gB0**2/(4*pi)               ! [erg/cm3]   

        ! Use gamma=1.5 for SW calculations (set in xml, but defaults to 1.5 here)
        call inpXML%Set_Val(Model%gamma,"physics/gamma",1.5_rp)
        call inpXML%Set_Val(Tsolar,"prob/Tsolar",25.38_rp)    ! in days

        ! convert Tsolar to code units
        Tsolar = Tsolar*24.*3600./gt0
      
        !Add gravity if required
        ! TODO: turn gravity on later
        Model%doGrav = .false.
        if (Model%doGrav) then
            !Force spherical gravity (zap non-radial components)
!            Model%doSphGrav = .true.
!            Model%Phi => PhiGrav
        endif

        !Change console output pointer
        ! don't use for now
!        timeString => helioTime 
        
        if (Model%isLoud) then
            write(*,*) '---------------'
            write(*,*) 'Heliosphere normalization'
            write(*,*) 'T0   [hr]       = ', gT0/3600.
            write(*,*) 'x0   [Rsolar]   = ', gx0
            write(*,*) 'v0   [km/s]  = '   , gv0*1.e-5
            write(*,*) 'P0   [erg/cm3]  = ', gP0
            write(*,*) 'B0   [nT]   = '    , gB0*1.e5
            write(*,*) '---------------'
        endif
        
        !Save scaling to gUnits_T structure in Model
        Model%Units%gT0 = gT0
        Model%Units%gx0 = gx0
        Model%Units%gv0 = gv0
        Model%Units%gD0 = gD0
        Model%Units%gP0 = gP0
        Model%Units%gB0 = gB0
!         Model%Units%gG0 = gG0   ! unused for helio, defaults to 0.

        ! without setting the scaling below, it defaults to 1. 
        !Add normalization/labels to output slicing
        ! Model%gamOut%tScl = gT0   !/3600. 
        ! Model%gamOut%dScl = gD0
        ! Model%gamOut%vScl = gv0   !*1.0e-5 !km/s
        ! Model%gamOut%pScl = gP0
        ! Model%gamOut%bScl = gB0   !*1.e5

        ! Model%gamOut%tID = 'Helio'
        ! Model%gamOut%tID = 's'  !'hr'
        ! Model%gamOut%dID = '#/cc'
        ! Model%gamOut%vID = 'km/s'
        ! Model%gamOut%pID = 'erg/cm3'
        ! Model%gamOut%bID = 'nT'

      end subroutine setHeliosphere

      subroutine helioTime(T,tStr)
        real(rp), intent(in) :: T
        character(len=strLen), intent(out) :: tStr
        
        write(tStr,'(f9.3,a)' ) T*gT0/3600.0, ' [hr]'
      end subroutine helioTime


      subroutine helio_ibcJ(bc,Model,Grid,State)
        ! improved versions of Kareems zeroGrad_(i,o)bcJ
        ! improvements:
        ! 1. Do true zero gradient (js+n-1 -> js-n, etc)
        ! 2. Apply to Br field, rather than Bi flux
        ! 3. Apply to Vr rather than V-vector
        
        class(helioInnerJBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        
        real(rp), dimension(NDIM) :: xyz,rhat
        real(rp) :: Vr,Vt,Vp,Vx,Vy,Vz
        integer :: n,i,j,k
        
        !j-boundaries (IN)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,i,k,Vx,Vy,Vz,Vr,xyz,rhat)
        do k=Grid%ksg,Grid%keg
           do n=1,Model%Ng  
              do i=Grid%isg,Grid%ieg
                 State%Gas(i,Grid%js-n,k,:,:) = State%Gas(i,Grid%js+n-1,k,:,:)
                   
                 ! for velocity do radial only
                 Vx = State%Gas(i,Grid%js+n-1,k,MOMX,BLK)
                 Vy = State%Gas(i,Grid%js+n-1,k,MOMY,BLK)
                 Vz = State%Gas(i,Grid%js+n-1,k,MOMZ,BLK)

                 !Calculate radial velocity
                 xyz  = Grid%xyzcc(i,Grid%js+n-1,k,:)
                 rhat = xyz/norm2(xyz)
                 Vr   = dot_product([Vx,Vy,Vz],rhat)
 
                 !Calculate radial hat vector
                 xyz  = Grid%xyzcc(i,Grid%js-n,k,:)
                 rhat = xyz/norm2(xyz)
                 State%Gas(i,Grid%js-n,k,MOMX,:) = Vr*rhat(XDIR)
                 State%Gas(i,Grid%js-n,k,MOMY,:) = Vr*rhat(YDIR)
                 State%Gas(i,Grid%js-n,k,MOMZ,:) = Vr*rhat(ZDIR)

                 State%magFlux(i,Grid%js-n,k,IDIR) = State%magFlux(i,Grid%js+n-1,k,IDIR)/Grid%Face(i,Grid%js+n-1,k,IDIR)*Grid%Face(i,Grid%js-n,k,IDIR)
                 State%magFlux(i,Grid%js-n,k,KDIR) = State%magFlux(i,Grid%js+n-1,k,KDIR)
                 State%magFlux(i,Grid%js-n,k,JDIR) = State%magFlux(i,Grid%js+n,k,JDIR)
              enddo
           enddo
        enddo
      end subroutine helio_ibcJ

      subroutine helio_obcJ(bc,Model,Grid,State)
        ! see comment above in helio_ibcJ
        class(helioOuterJBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        
        real(rp), dimension(NDIM) :: xyz,rhat
        real(rp) :: Vr,Vt,Vp,Vx,Vy,Vz
        integer :: n,i,j,k
        
        !j-boundaries (OUT)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,i,k,Vx,Vy,Vz,Vr,xyz,rhat)
        do k=Grid%ksg,Grid%keg
           do n=1,Model%Ng  
              do i=Grid%isg,Grid%ieg
                 State%Gas(i,Grid%je+n,k,:,:) = State%Gas(i,Grid%je-n+1,k,:,:)

                 ! for velocity do radial only
                 Vx = State%Gas(i,Grid%je-n+1,k,MOMX,BLK)
                 Vy = State%Gas(i,Grid%je-n+1,k,MOMY,BLK)
                 Vz = State%Gas(i,Grid%je-n+1,k,MOMZ,BLK)
                   
                 !Calculate radial velocity
                 xyz  = Grid%xyzcc(i,Grid%je-n+1,k,:)
                 rhat = xyz/norm2(xyz)
                 Vr   = dot_product([Vx,Vy,Vz],rhat)
 
                 !Calculate radial hat vector
                 xyz  = Grid%xyzcc(i,Grid%je+n,k,:)
                 rhat = xyz/norm2(xyz)
                 State%Gas(i,Grid%je+n,k,MOMX,:) = Vr*rhat(XDIR)
                 State%Gas(i,Grid%je+n,k,MOMY,:) = Vr*rhat(YDIR)
                 State%Gas(i,Grid%je+n,k,MOMZ,:) = Vr*rhat(ZDIR)

                 State%magFlux(i,Grid%je+n,k  ,IDIR) = State%magFlux(i,Grid%je-n+1,k,IDIR)/Grid%Face(i,Grid%je-n+1,k,IDIR)*Grid%Face(i,Grid%je+n,k,IDIR)
                 State%magFlux(i,Grid%je+n,k  ,KDIR) = State%magFlux(i,Grid%je-n+1,k,KDIR)
                 State%magFlux(i,Grid%je+n+1,k,JDIR) = State%magFlux(i,Grid%je-n+1,k ,JDIR)
              enddo
           enddo
        enddo
      end subroutine helio_obcJ

      subroutine helio_obcI(bc,Model,Grid,State)
        ! A version of Kareem's zerGrad_obcI (see ../bcs.F90)
        ! works better for superfast helio
        class(helioOuterIBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,ig,j,k
        real(rp), dimension(NDIM) :: Bxyz
        real(rp), dimension(NDIM) :: xyz
        real(rp) :: R, R0

        !i-boundaries (OUT)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,j,k,ig,xyz,R0,R)
        do k=Grid%ksg,Grid%keg
           do j=Grid%jsg,Grid%jeg
              !Get Cartesian field in last physical cell
              do n=1,Model%Ng
                 ig = Grid%ie+n

                 xyz = Grid%xyzcc(Grid%ie,j,k,:)
                 R0 = norm2(xyz)
                 xyz = Grid%xyzcc(ig,j,k,:)
                 R = norm2(xyz)

                 State%Gas(ig,j,k,:,:)    = State%Gas(Grid%ie,j,k,:,:)
                 ! is this necessary? I seem to remember it worked better. Test again?
                 State%Gas(ig,j,k,DEN,:)  = (R0/R)**2*State%Gas(Grid%ie,j,k,DEN,:) ! drop density to suck stuff out
                 ! note, we're lazier here with Bj,Bk fluxes than in the j-BCs above. 
                 State%magFlux(ig+1,j,k,IDIR) = State%magFlux(Grid%ie,j,k,IDIR) 
                 State%magFlux(ig  ,j,k,JDIR) = State%magFlux(Grid%ie,j,k,JDIR) 
                 State%magFlux(ig  ,j,k,KDIR) = State%magFlux(Grid%ie,j,k,KDIR)
              enddo
           enddo
        enddo
      end subroutine helio_obcI

end module helioutils
