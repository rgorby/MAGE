!Routines for boundary conditions
module bcs
    use types
    use gamutils
    use math
    use gridutils
    use ringutils
    use basebc

    implicit none

    type, extends(baseBC_T) :: lazyBC_T
        contains
        procedure :: doBC => lazyBC
    end type lazyBC_T

    type, extends(baseBC_T) :: periodicInnerIBC_T
        contains
        procedure :: doBC => periodic_ibcI
    end type periodicInnerIBC_T

    type, extends(baseBC_T) :: periodicOuterIBC_T
        contains
        procedure :: doBC => periodic_obcI
    end type periodicOuterIBC_T

    type, extends(baseBC_T) :: periodicInnerJBC_T
        contains
        procedure :: doBC => periodic_ibcJ
    end type periodicInnerJBC_T

    type, extends(baseBC_T) :: periodicOuterJBC_T
        contains
        procedure :: doBC => periodic_obcJ
    end type periodicOuterJBC_T

    type, extends(baseBC_T) :: periodicInnerKBC_T
        contains
        procedure :: doBC => periodic_ibcK
    end type periodicInnerKBC_T

    type, extends(baseBC_T) :: periodicOuterKBC_T
        contains
        procedure :: doBC => periodic_obcK
    end type periodicOuterKBC_T

    type, extends(baseBC_T) :: zeroExtensionInnerIBC_T
        contains
        procedure :: doBC => zeroExten_ibcI
    end type zeroExtensionInnerIBC_T

    type, extends(baseBC_T) :: zeroExtensionOuterIBC_T
        contains
        procedure :: doBC => zeroExten_obcI
    end type zeroExtensionOuterIBC_T

    type, extends(baseBC_T) :: zeroExtensionInnerJBC_T
        contains
        procedure :: doBC => zeroExten_ibcJ
    end type zeroExtensionInnerJBC_T

    type, extends(baseBC_T) :: zeroExtensionOuterJBC_T
        contains
        procedure :: doBC => zeroExten_obcJ
    end type zeroExtensionOuterJBC_T

    type, extends(baseBC_T) :: zeroExtensionInnerKBC_T
        contains
        procedure :: doBC => zeroExten_ibcK
    end type zeroExtensionInnerKBC_T

    type, extends(baseBC_T) :: zeroExtensionOuterKBC_T
        contains
        procedure :: doBC => zeroExten_obcK
    end type zeroExtensionOuterKBC_T

    type, extends(baseBC_T) :: cartesianReflectingInnerIBC_T
        contains
        procedure :: doBC => cartReflect_ibcI
    end type cartesianReflectingInnerIBC_T

    type, extends(baseBC_T) :: cartesianReflectingOuterIBC_T
        contains
        procedure :: doBC => cartReflect_obcI
    end type cartesianReflectingOuterIBC_T

    type, extends(baseBC_T) :: cartesianReflectingInnerJBC_T
        contains
        procedure :: doBC => cartReflect_ibcJ
    end type cartesianReflectingInnerJBC_T

    type, extends(baseBC_T) :: cartesianReflectingOuterJBC_T
        contains
        procedure :: doBC => cartReflect_obcJ
    end type cartesianReflectingOuterJBC_T

    type, extends(baseBC_T) :: cartesianReflectingInnerKBC_T
        contains
        procedure :: doBC => cartReflect_ibcK
    end type cartesianReflectingInnerKBC_T

    type, extends(baseBC_T) :: cartesianReflectingOuterKBC_T
        contains
        procedure :: doBC => cartReflect_obcK
    end type cartesianReflectingOuterKBC_T

    type, extends(baseBC_T) :: zeroGradientInnerIBC_T
        contains
        procedure :: doBC => zeroGrad_ibcI
    end type zeroGradientInnerIBC_T

    type, extends(baseBC_T) :: zeroGradientOuterIBC_T
        contains
        procedure :: doBC => zeroGrad_obcI
    end type zeroGradientOuterIBC_T

    type, extends(baseBC_T) :: zeroGradientInnerJBC_T
        contains
        procedure :: doBC => zeroGrad_ibcJ
    end type zeroGradientInnerJBC_T

    type, extends(baseBC_T) :: zeroGradientOuterJBC_T
        contains
        procedure :: doBC => zeroGrad_obcJ
    end type zeroGradientOuterJBC_T

    type, extends(baseBC_T) :: zeroGradientInnerKBC_T
        contains
        procedure :: doBC => zeroGrad_ibcK
    end type zeroGradientInnerKBC_T

    type, extends(baseBC_T) :: zeroGradientOuterKBC_T
        contains
        procedure :: doBC => zeroGrad_obcK
    end type zeroGradientOuterKBC_T

    type, extends(baseBC_T) :: sphereInBC_T
        contains
        procedure :: doBC => SphIn
    end type sphereInBC_T

    type, extends(baseBC_T) :: sphereOutBC_T
        contains
        procedure :: doBC => SphOut
    end type sphereOutBC_T

    type, extends(baseBC_T) :: cylindricalPoleBC_T
        contains
        procedure :: doBC => cylpole
    end type cylindricalPoleBC_T

    type, extends(baseBC_T) :: lfmInBC_T
        contains
        procedure :: doBC => lfmIn
    end type lfmInBC_T

    type, extends(baseBC_T) :: lfmOutBC_T
        contains
        procedure :: doBC => lfmOut
    end type lfmOutBC_T

    contains

    !Below here are various predefined BCs
!--------------------------------------------
    !Lazy, do-nothing boundary
    subroutine lazyBC(bc,Model,Grid,State)
        class(lazyBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        return
    end subroutine lazyBC

!--------------------------------------------    
    !Standard periodic BCs
    subroutine periodic_ibcI(bc,Model,Grid,State)
        class(periodicInnerIBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        integer :: n,i,j,k

        !i-boundaries (IN)
            do k=Grid%ksg,Grid%keg
                do j=Grid%jsg,Grid%jeg
                    do n=1,Model%Ng
                        State%Gas(Grid%is -n,j,k,:,:)  = State%Gas(Grid%ie-n+1,j,k,:,:)
                        State%Bxyz(Grid%is-n,j,k,:)  = State%Bxyz(Grid%ie-n+1,j,k,:)
                        State%magFlux(Grid%is-n,j,k,YDIR:ZDIR) = State%magFlux(Grid%ie -n+1,j,k,YDIR:ZDIR)
                        State%magFlux(Grid%is-n,j,k,XDIR) = State%magFlux(Grid%ie -n+1,j,k,XDIR)
                    enddo
                enddo
            enddo

    end subroutine periodic_ibcI

    subroutine periodic_obcI(bc,Model,Grid,State)
        class(periodicOuterIBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        integer :: n,i,j,k

            do k=Grid%ksg,Grid%keg
                do j=Grid%jsg,Grid%jeg
                    do n=1,Model%Ng
                        State%Gas(Grid%ie+n,j,k,:,:) = State%Gas(Grid%is +n-1,j,k,:,:)
                        State%Bxyz(Grid%ie+n,j,k,:) = State%Bxyz(Grid%is +n-1,j,k,:)
                        State%magFlux(Grid%ie+n,j,k,YDIR:ZDIR) = State%magFlux(Grid%is +n-1,j,k,YDIR:ZDIR)
                        State%magFlux(Grid%ie+n+1,j,k,XDIR) = State%magFlux(Grid%is+n,j,k,XDIR)
                    enddo
                enddo
            enddo    
    end subroutine periodic_obcI

    subroutine periodic_ibcJ(bc,Model,Grid,State)
        class(periodicInnerJBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        integer :: n,i,j,k

            do k=Grid%ksg,Grid%keg
                do n=1,Model%Ng  
                    do i=Grid%isg,Grid%ieg
                        State%Gas(i,Grid%js-n,k,:,:) = State%Gas(i,Grid%je-n+1,k,:,:)
                        State%Bxyz(i,Grid%js-n,k,:) = State%Bxyz(i,Grid%je-n+1,k,:)
                        State%magFlux(i,Grid%js-n,k,:) = State%magFlux(i,Grid%je-n+1,k,:)         
                    enddo
                enddo
            enddo    
    end subroutine periodic_ibcJ

    subroutine periodic_obcJ(bc,Model,Grid,State)
        class(periodicOuterJBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        integer :: n,i,j,k

            do k=Grid%ksg,Grid%keg
                do n=1,Model%Ng  
                    do i=Grid%isg,Grid%ieg
                        State%Gas(i,Grid%je+n,k,:,:) = State%Gas(i,Grid%js+n-1,k,:,:)
                        State%Bxyz(i,Grid%je+n,k,:) = State%Bxyz(i,Grid%js+n-1,k,:)
                        State%magFlux(i,Grid%je+n,k,XDIR) = State%magFlux(i,Grid%js+n-1,k,XDIR)
                        State%magFlux(i,Grid%je+n,k,ZDIR) = State%magFlux(i,Grid%js+n-1,k,ZDIR)
                        State%magFlux(i,Grid%je+n+1,k,YDIR) = State%magFlux(i,Grid%js+n  ,k,YDIR)
                    enddo
                enddo
            enddo

    end subroutine periodic_obcJ

    subroutine periodic_ibcK(bc,Model,Grid,State)
        class(periodicInnerKBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        integer :: n,i,j,k

        if (Model%do25D) then
            !Just replicate ks slice
            do n=1,Model%Ng
                
                State%Gas(:,:,Grid%ks-n,:,:) = State%Gas(:,:,Grid%ks,:,:)
                State%Bxyz(:,:,Grid%ks-n,:) = State%Bxyz(:,:,Grid%ks,:)
                State%magFlux(:,:,Grid%ks-n,:) = State%magFlux(:,:,Grid%ks,:)   
            enddo
        else
                !$OMP PARALLEL DO default(shared) &
                !$OMP private(i,j,n)
                do j=Grid%jsg,Grid%jeg
                    do i=Grid%isg,Grid%ieg
                        do n=1,Model%Ng          
                            State%Gas (i,j,Grid%ks-n,:,:)  = State%Gas (i,j,Grid%ke-n+1,:,:)
                            State%Bxyz(i,j,Grid%ks-n,:)    = State%Bxyz(i,j,Grid%ke-n+1,:)
                            State%magFlux(i,j,Grid%ks-n,:) = State%magFlux(i,j,Grid%ke-n+1,:)                
                        enddo
                    enddo
                enddo
        endif

    end subroutine periodic_ibcK

    subroutine periodic_obcK(bc,Model,Grid,State)
        class(periodicOuterKBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        integer :: n,i,j,k

        if (Model%do25D) then
            !Just replicate ks slice
            do n=1,Model%Ng
                State%Gas(:,:,Grid%ks+n,:,:) = State%Gas(:,:,Grid%ks,:,:)
                State%Bxyz(:,:,Grid%ks+n,:) = State%Bxyz(:,:,Grid%ks,:)
                State%magFlux(:,:,Grid%ks+n,:) = State%magFlux(:,:,Grid%ks,:)  
            enddo
        else
                !$OMP PARALLEL DO default(shared) &
                !$OMP private(i,j,n)        
                do j=Grid%jsg,Grid%jeg
                    do i=Grid%isg,Grid%ieg  
                        do n=1,Model%Ng         
                            State%Gas(i,j,Grid%ke+n,:,:)  = State%Gas(i,j,Grid%ks+n-1,:,:)
                            State%Bxyz(i,j,Grid%ke+n,:)  = State%Bxyz(i,j,Grid%ks+n-1,:)
                            State%magFlux(i,j,Grid%ke+n,XDIR:YDIR)  = State%magFlux(i,j,Grid%ks+n-1,XDIR:YDIR)
                            State%magFlux(i,j,Grid%ke+n+1,ZDIR)  = State%magFlux(i,j,Grid%ks+n,ZDIR)
            
                        enddo
                    enddo
                enddo
        endif

    end subroutine periodic_obcK

!--------------------------------------------
    !Zeroth extension BCs
    subroutine zeroExten_ibcI(bc,Model,Grid,State)
        class(zeroExtensionInnerIBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !i-boundaries (IN)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do n=1,Model%Ng
                    ! gas variables equal to i=is
                    State%Gas(Grid%is -n,j,k,:,:)  = State%Gas(Grid%is,j,k,:,:)
                    ! magnetic fields equal to i=is
                    State%Bxyz(Grid%is -n,j,k,:)  = State%Bxyz(Grid%is,j,k,:)
                    ! i,j,k-magflux equals to is
                    State%magFlux(Grid%is-n,j,k,:) = State%magFlux(Grid%is,j,k,:)
                enddo
            enddo
        enddo

    end subroutine zeroExten_ibcI

    subroutine zeroExten_obcI(bc,Model,Grid,State)
        class(zeroExtensionOuterIBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !i-boundaries (OUT)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do n=1,Model%Ng
                    ! gas variables equal to i=ie
                    State%Gas(Grid%ie+n,j,k,:,:)  = State%Gas(Grid%ie,j,k,:,:)
                    ! magnetic fields equal to i=ie
                    State%Bxyz(Grid%ie+n,j,k,:)  = State%Bxyz(Grid%ie,j,k,:)
                    ! i-magflux equals to ie+1
                    State%magFlux(Grid%ie+n+1,j,k,IDIR) = State%magFlux(Grid%ie+1,j,k,IDIR)
                    ! j,k-magFlux equals to ie
                    State%magFlux(Grid%ie+n,j,k,JDIR:KDIR)  = State%Bxyz(Grid%ie,j,k,JDIR:KDIR)
                enddo
            enddo
        enddo    
    end subroutine zeroExten_obcI

    subroutine zeroExten_ibcJ(bc,Model,Grid,State)
        class(zeroExtensionInnerJBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k
        
        !j-boundaries (IN)
        do k=Grid%ksg,Grid%keg
            do n=1,Model%Ng
                do i=Grid%isg,Grid%ieg
                    ! gas variables equal to j=js
                    State%Gas(i,Grid%js-n,k,:,:) = State%Gas(i,Grid%js,k,:,:)
                    ! magnetic fields equal to j=js
                    State%Bxyz(i,Grid%js-n,k,:) = State%Bxyz(i,Grid%js,k,:)
                    ! i,j,k-magflux equal to j = js
                    State%magFlux(i,Grid%js-n,k,:) = State%magFlux(i,Grid%js,k,:)
                enddo
            enddo
        enddo

    end subroutine zeroExten_ibcJ

    !Zero grad BCs
    subroutine zeroExten_obcJ(bc,Model,Grid,State)
        class(zeroExtensionOuterJBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !j-boundaries (OUT)
        do k=Grid%ksg,Grid%keg
            do n=1,Model%Ng  
                do i=Grid%isg,Grid%ieg
                    ! gas variables equal to j=je
                    State%Gas(i,Grid%je+n,k,:,:) = State%Gas(i,Grid%je,k,:,:)
                    ! magnetic fields equal to j = je
                    State%Bxyz(i,Grid%je+n,k,:) = State%Bxyz(i,Grid%je,k,:)
                    ! i, k-magflux equal to j=je
                    State%magFlux(i,Grid%je+n,k  ,IDIR) = State%magFlux(i,Grid%je,k  ,IDIR)
                    State%magFlux(i,Grid%je+n,k  ,KDIR) = State%magFlux(i,Grid%je,k  ,KDIR)
                    ! j-magflux equals to je+1
                    State%magFlux(i,Grid%je+n+1,k,JDIR) = State%magFlux(i,Grid%je+1,k,JDIR)
                enddo
            enddo
        enddo

    end subroutine zeroExten_obcJ

    subroutine zeroExten_ibcK(bc,Model,Grid,State)
        class(zeroExtensionInnerKBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !k-boundaries (IN)
        do n=1,Model%Ng
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg    
                    ! gas variables equal to k = ks       
                    State%Gas(i,j,Grid%ks-n,:,:)  = State%Gas(i,j,Grid%ks,:,:)
                    ! magnetic fields equal to k = ks
                    State%Bxyz(i,j,Grid%ks-n,:)  = State%Bxyz(i,j,Grid%ks,:)
                    ! i,j,k-magflux equal to k = ks
                    State%magFlux(i,j,Grid%ks-n,:)  = State%magFlux(i,j,Grid%ks,:)                
                enddo
            enddo
        enddo

    end subroutine zeroExten_ibcK

    subroutine zeroExten_obcK(bc,Model,Grid,State)
        class(ZeroExtensionOuterKBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !k-boundaries (OUT)
        do n=1,Model%Ng
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg  
                    ! gas variables equal to k = ke         
                    State%Gas(i,j,Grid%ke+n,:,:)  = State%Gas(i,j,Grid%ke,:,:)
                    ! magnetic fields equal to k = ke
                    State%Bxyz(i,j,Grid%ke+n,:)  = State%Bxyz(i,j,Grid%ke,:)
                    ! i,j-magflux equal to k = ks
                    State%magFlux(i,j,Grid%ke+n,IDIR:JDIR)  = State%magFlux(i,j,Grid%ke,IDIR:JDIR)
                    ! k-magFlux equals to k = ks+1
                    State%magFlux(i,j,Grid%ke+n+1,KDIR)  = State%magFlux(i,j,Grid%ke+1,KDIR)
                enddo
            enddo
        enddo

    end subroutine zeroExten_obcK

!---------------------------------------------------------
! Reflecting boundary only works for Cartesian grids
    subroutine cartReflect_ibcI(bc,Model,Grid,State)
        class(cartesianReflectingInnerIBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !i-boundaries (IN)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do n=1,Model%Ng
                    ! density, y,z-momentum and energy has zero gradient at i=is-1/2
                    State%Gas(Grid%is -n,j,k,DEN,:)  = State%Gas(Grid%is+n-1,j,k,DEN,:)
                    State%Gas(Grid%is -n,j,k,MOMY:MOMZ,:)  = State%Gas(Grid%is+n-1,j,k,MOMY:MOMZ,:)
                    State%Gas(Grid%is -n,j,k,ENERGY,:)  = State%Gas(Grid%is+n-1,j,k,ENERGY,:)
                    ! x-momentum has zero value at i=is-1/2
                    State%Gas(Grid%is -n,j,k,MOMX,:)  = -1.0*State%Gas(Grid%is+n-1,j,k,MOMX,:)
                    ! magnetic fields have zero gradient at i=is-1/2
                    State%Bxyz(Grid%is -n,j,k,:)  = State%Bxyz(Grid%is+n-1,j,k,:)
                    ! i-magflux has zero gradient at i=is
                    State%magFlux(Grid%is-n,j,k,IDIR) = State%magFlux(Grid%is+n,j,k,IDIR)
                    ! j,k-magflux has zero gradient at i=is-1/2
                    State%magFlux(Grid%is-n,j,k,JDIR:KDIR) = State%magFlux(Grid%is+n-1,j,k,JDIR:KDIR)
                enddo
            enddo
        enddo

    end subroutine cartReflect_ibcI

    subroutine cartReflect_obcI(bc,Model,Grid,State)
        class(cartesianReflectingOuterIBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,ig,j,k

        !i-boundaries (OUT)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do n=1,Model%Ng
                    ig = Grid%ie+n
                    ! density, y,z-momentum and energy has zero gradient at i=ie+1/2
                    State%Gas(Grid%ie+n,j,k,DEN,:)  = State%Gas(Grid%ie-n+1,j,k,DEN,:)
                    State%Gas(Grid%ie+n,j,k,MOMY:MOMZ,:)  = State%Gas(Grid%ie-n+1,j,k,MOMY:MOMZ,:)
                    State%Gas(Grid%ie+n,j,k,ENERGY,:)  = State%Gas(Grid%ie-n+1,j,k,ENERGY,:)
                    ! x-momentum has zero value at i=ie+1/2
                    State%Gas(Grid%ie+n,j,k,MOMX,:)  = -1*State%Gas(Grid%ie-n+1,j,k,MOMX,:)
                    ! i-magflux has zero gradient at i=ie+1
                    State%magFlux(Grid%ie+n+1,j,k,IDIR)  = State%magFlux(Grid%ie-n+1,j,k,IDIR)
                    ! j,k-magflux have zero gradient at i=ie+1/2
                    State%magFlux(Grid%ie+n,j,k,JDIR:KDIR)  = State%magFlux(Grid%ie-n+1,j,k,JDIR:KDIR)
                enddo
            enddo
        enddo    
    end subroutine cartReflect_obcI

    subroutine cartReflect_ibcJ(bc,Model,Grid,State)
        class(cartesianReflectingInnerJBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k
        
        !j-boundaries (IN)
        do k=Grid%ksg,Grid%keg
            do n=1,Model%Ng  
                do i=Grid%isg,Grid%ieg
                    ! fluid variables are cell-centered, only flip the normal component: j-dir
                    State%Gas(i,Grid%js-n,k,DEN,:) = State%Gas(i,Grid%js+n-1,k,DEN,:)
                    State%Gas(i,Grid%js-n,k,MOMX,:) = State%Gas(i,Grid%js+n-1,k,MOMX,:)
                    State%Gas(i,Grid%js-n,k,MOMY,:) = -1*State%Gas(i,Grid%js+n-1,k,MOMY,:)
                    State%Gas(i,Grid%js-n,k,MOMZ,:) = State%Gas(i,Grid%js+n-1,k,MOMZ,:)
                    State%Gas(i,Grid%js-n,k,ENERGY,:) = State%Gas(i,Grid%js+n-1,k,ENERGY,:)
                    ! mag-field and mag-flux are set to be zero gradient
                    ! mag-field components are cell-centered so treated the same way as fluid variables
                    State%Bxyz(i,Grid%js-n,k,:) = State%Bxyz(i,Grid%js+n-1,k,:)
                    ! mag-flux components are slightly more complicated
                    ! in the j-dir, i-flux and j-flux are cell-centered so treated the same way as fluid variables
                    !               j-flux is face-centered along j-dir so the index is shifted "right" by 1
                    State%magFlux(i,Grid%js-n,k,IDIR) = State%magFlux(i,Grid%js+n-1,k,IDIR)
                    State%magFlux(i,Grid%js-n,k,KDIR) = State%magFlux(i,Grid%js+n-1,k,KDIR) 
                    State%magFlux(i,Grid%js-n,k,JDIR) = State%magFlux(i,Grid%js+n,k,JDIR)        
                enddo
            enddo
        enddo

    end subroutine cartReflect_ibcJ

    subroutine cartReflect_obcJ(bc,Model,Grid,State)
        class(cartesianReflectingOuterJBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !j-boundaries (OUT)
        do k=Grid%ksg,Grid%keg
            do n=1,Model%Ng  
                do i=Grid%isg,Grid%ieg
                    ! fluid variables are cell-centered, only flip the normal component: j-dir
                    State%Gas(i,Grid%je+n,k,DEN,:) = State%Gas(i,Grid%je-n+1,k,DEN,:)
                    State%Gas(i,Grid%je+n,k,MOMX,:) = State%Gas(i,Grid%je-n+1,k,MOMX,:)
                    State%Gas(i,Grid%je+n,k,MOMY,:) = -1*State%Gas(i,Grid%je-n+1,k,MOMY,:)
                    State%Gas(i,Grid%je+n,k,MOMZ,:) = State%Gas(i,Grid%je-n+1,k,MOMZ,:)
                    State%Gas(i,Grid%je+n,k,ENERGY,:) = State%Gas(i,Grid%je-n+1,k,ENERGY,:)
                    ! mag-field and mag-flux are set to be zero gradient
                    ! mag-field components are cell-centered so treated the same way as fluid variables
                    State%Bxyz(i,Grid%je+n,k,:) = State%Bxyz(i,Grid%je-n+1,k,:)
                    ! mag-flux components are slightly more complicated
                    ! in the j-dir, i-flux and j-flux are cell-centered so treated the same way as fluid variables
                    !               j-flux is face-centered along j-dir so the index is shifted "left" by 1
                    State%magFlux(i,Grid%je+n,  k,IDIR) = State%magFlux(i,Grid%je-n+1,k,IDIR)
                    State%magFlux(i,Grid%je+n,  k,KDIR) = State%magFlux(i,Grid%je-n+1,k,KDIR)
                    State%magFlux(i,Grid%je+n+1,k,JDIR) = State%magFlux(i,Grid%je-n+1,k,JDIR)
                enddo
            enddo
        enddo

    end subroutine cartReflect_obcJ

    subroutine cartReflect_ibcK(bc,Model,Grid,State)
        class(cartesianReflectingInnerKBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !k-boundaries (IN)
        do n=1,Model%Ng
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg
                    ! fluid variables are cell-centered, only flip the normal component: k-dir
                    State%Gas(i,j,Grid%ks-n,DEN,:)  = State%Gas(i,j,Grid%ks+n-1,DEN,:)
                    State%Gas(i,j,Grid%ks-n,MOMX,:)  = State%Gas(i,j,Grid%ks+n-1,MOMX,:)
                    State%Gas(i,j,Grid%ks-n,MOMY,:)  = State%Gas(i,j,Grid%ks+n-1,MOMY,:)
                    State%Gas(i,j,Grid%ks-n,MOMZ,:)  = -1*State%Gas(i,j,Grid%ks+n-1,MOMZ,:)
                    State%Gas(i,j,Grid%ks-n,ENERGY,:)  = State%Gas(i,j,Grid%ks+n-1,ENERGY,:)
                    ! mag-field and mag-flux are set to be zero gradient
                    ! mag-field components are cell-centered so treated the same way as fluid variables
                    State%Bxyz(i,j,Grid%ks-n,:)  = State%Bxyz(i,j,Grid%ks+n-1,:)
                    ! mag-flux components are slightly more complicated
                    ! in the k-dir, i-flux and j-flux are cell-centered so treated the same way as fluid variables
                    !               k-flux is face-centered along k-dir so the index is shifted "right" by 1
                    State%magFlux(i,j,Grid%ks-n,IDIR) = State%magFlux(i,j,Grid%ks+n-1,IDIR) ! the reflect point is ks-1/2
                    State%magFlux(i,j,Grid%ks-n,JDIR) = State%magFlux(i,j,Grid%ks+n-1,JDIR) ! the reflect point is ks-1/2
                    State%magFlux(i,j,Grid%ks-n,KDIR) = State%magFlux(i,j,Grid%ks+n  ,KDIR) ! the reflect point is ks
                enddo
            enddo
        enddo

    end subroutine cartReflect_ibcK
   
    subroutine cartReflect_obcK(bc,Model,Grid,State)
        class(cartesianReflectingOuterKBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !k-boundaries (OUT)
        do n=1,Model%Ng
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg
                    ! fluid variables are cell-centered, only flip the normal component: k-dir
                    State%Gas(i,j,Grid%ke+n,DEN,:)  = State%Gas(i,j,Grid%ke-n+1,DEN,:)
                    State%Gas(i,j,Grid%ke+n,MOMX,:)  = State%Gas(i,j,Grid%ke-n+1,MOMX,:)
                    State%Gas(i,j,Grid%ke+n,MOMY,:)  = State%Gas(i,j,Grid%ke-n+1,MOMY,:)
                    State%Gas(i,j,Grid%ke+n,MOMZ,:)  = -1.0*State%Gas(i,j,Grid%ke-n+1,MOMZ,:)
                    State%Gas(i,j,Grid%ke+n,ENERGY,:)  = State%Gas(i,j,Grid%ke-n+1,ENERGY,:)
                    ! mag-field and mag-flux are set to be zero gradient
                    ! mag-field components are cell-centered so treated the same way as fluid variables
                    State%Bxyz(i,j,Grid%ke+n,:)  = State%Bxyz(i,j,Grid%ke-n+1,:)
                    ! in the k-dir, i-flux and j-flux are cell-centered so treated the same way as fluid variables
                    !               k-flux is face-centered along k-dir so the index is shifted "left" by 1
                    State%magFlux(i,j,Grid%ke+n,IDIR)  = State%magFlux(i,j,Grid%ke-n+1,IDIR)   ! the reflect point is ke+1/2
                    State%magFlux(i,j,Grid%ke+n,JDIR)  = State%magFlux(i,j,Grid%ke-n+1,JDIR)   ! the reflect point is ke+1/2
                    State%magFlux(i,j,Grid%ke+n+1,KDIR)  = State%magFlux(i,j,Grid%ke-n+1,KDIR) ! the reflect point is ke+1

                enddo
            enddo
        enddo

    end subroutine cartReflect_obcK

!--------------------------------------------
    !Zero grad BCs
    subroutine zeroGrad_ibcI(bc,Model,Grid,State)
        class(zeroGradientInnerIBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !i-boundaries (IN)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do n=1,Model%Ng
                    State%Gas(Grid%is -n,j,k,:,:)  = State%Gas(Grid%is,j,k,:,:)
                    State%magFlux(Grid%is-n,j,k,:) = State%magFlux(Grid%is,j,k,:)
                enddo
            enddo
        enddo


    end subroutine zeroGrad_ibcI

    subroutine zeroGrad_obcI(bc,Model,Grid,State)
        class(zeroGradientOuterIBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,ig,j,k
        real(rp), dimension(NDIM) :: Bxyz

        !i-boundaries (OUT)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                !Get Cartesian field in last physical cell
                Bxyz = CellBxyz(Model,Grid,State%magFlux,Grid%ie,j,k)
                do n=1,Model%Ng
                    ig = Grid%ie+n
                    State%Gas(ig,j,k,:,:)  = State%Gas(Grid%ie,j,k,:,:)
                    State%magFlux(ig+1,j,k,IDIR) = Grid%face(ig+1,j,k,IDIR)*dot_product(Bxyz,Grid%Tf(ig+1,j,k,NORMX:NORMZ,IDIR))
                    State%magFlux(ig  ,j,k,JDIR) = Grid%face(ig  ,j,k,JDIR)*dot_product(Bxyz,Grid%Tf(ig  ,j,k,NORMX:NORMZ,JDIR))
                    State%magFlux(ig  ,j,k,KDIR) = Grid%face(ig  ,j,k,KDIR)*dot_product(Bxyz,Grid%Tf(ig  ,j,k,NORMX:NORMZ,KDIR))
                enddo
            enddo
        enddo    
    end subroutine zeroGrad_obcI

    subroutine zeroGrad_ibcJ(bc,Model,Grid,State)
        class(zeroGradientInnerJBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k
        
        !j-boundaries (IN)
        do k=Grid%ksg,Grid%keg
            do n=1,Model%Ng
                do i=Grid%isg,Grid%ieg
                    State%Gas(i,Grid%js-n,k,:,:) = State%Gas(i,Grid%js,k,:,:)
                    State%magFlux(i,Grid%js-n,k,IDIR) = State%magFlux(i,Grid%js,k,IDIR)
                    State%magFlux(i,Grid%js-n,k,JDIR) = State%magFlux(i,Grid%js,k,JDIR)
                    State%magFlux(i,Grid%js-n,k,KDIR) = State%magFlux(i,Grid%js,k,KDIR)
                enddo
            enddo
        enddo

    end subroutine zeroGrad_ibcJ

    !Zero grad BCs
    subroutine zeroGrad_obcJ(bc,Model,Grid,State)
        class(zeroGradientOuterJBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !j-boundaries (OUT)
        do k=Grid%ksg,Grid%keg
            do n=1,Model%Ng  
                do i=Grid%isg,Grid%ieg
                    State%Gas(i,Grid%je+n,k,:,:) = State%Gas(i,Grid%je,k,:,:)
                    State%magFlux(i,Grid%je+n,k  ,IDIR) = State%magFlux(i,Grid%je,k   ,IDIR)
                    State%magFlux(i,Grid%je+n,k  ,KDIR) = State%magFlux(i,Grid%je,k   ,KDIR)
                    State%magFlux(i,Grid%je+n+1,k,JDIR) = State%magFlux(i,Grid%je+1,k ,JDIR)
                enddo
            enddo
        enddo

    end subroutine zeroGrad_obcJ

    subroutine zeroGrad_ibcK(bc,Model,Grid,State)
        class(zeroGradientInnerKBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !k-boundaries (IN)
        do n=1,Model%Ng
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg           
                    State%Gas(i,j,Grid%ks-n,:,:)  = State%Gas(i,j,Grid%ks,:,:)
                    State%magFlux(i,j,Grid%ks-n,:)  = State%magFlux(i,j,Grid%ks,:)                
                enddo
            enddo
        enddo

    end subroutine zeroGrad_ibcK

    subroutine zeroGrad_obcK(bc,Model,Grid,State)
        class(zeroGradientOuterKBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,j,k

        !k-boundaries (OUT)
        do n=1,Model%Ng
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg           
                    State%Gas(i,j,Grid%ke+n,:,:)  = State%Gas(i,j,Grid%ke,:,:)
                    State%magFlux(i,j,Grid%ke+n,XDIR:YDIR)  = State%magFlux(i,j,Grid%ke,IDIR:JDIR)
                    State%magFlux(i,j,Grid%ke+n+1,ZDIR)  = State%magFlux(i,j,Grid%ke+1,JDIR)
    
                enddo
            enddo
        enddo

    end subroutine zeroGrad_obcK

!--------------------------------------------
    !Various geometry-specific BCs
    subroutine SphIn(bc,Model,Grid,State)
        class(sphereInBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,ig,j,k

        !i-boundaries (IN)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                !Get Cartesian field in last physical cell
                do n=1,Model%Ng
                    ig = Grid%is-n
                    State%Gas(ig,j,k,:,:)  = State%Gas(Grid%is,j,k,:,:)
                    State%magFlux(ig  ,j,k,IDIR) = State%magFlux(Grid%is,j,k,IDIR) 
                    State%magFlux(ig  ,j,k,JDIR) = State%magFlux(Grid%is,j,k,JDIR) 
                    State%magFlux(ig  ,j,k,KDIR) = State%magFlux(Grid%is,j,k,KDIR)
                    State%Bxyz(ig,j,k,:) = State%Bxyz(Grid%is,j,k,:)
                enddo
            enddo
        enddo

    end subroutine SphIn
    
    subroutine SphOut(bc,Model,Grid,State)   
        class(sphereOutBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,ig,j,k

        !i-boundaries (OUT)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                !Get Cartesian field in last physical cell
                do n=1,Model%Ng
                    ig = Grid%ie+n
                    State%Gas(ig,j,k,:,:)  = State%Gas(Grid%ie,j,k,:,:)
                    State%magFlux(ig+1,j,k,IDIR) = State%magFlux(Grid%ie,j,k,IDIR) 
                    State%magFlux(ig  ,j,k,JDIR) = State%magFlux(Grid%ie,j,k,JDIR) 
                    State%magFlux(ig  ,j,k,KDIR) = State%magFlux(Grid%ie,j,k,KDIR)
                    State%Bxyz(ig,j,k,:) = State%Bxyz(Grid%ie,j,k,:)
                enddo
            enddo
        enddo

    end subroutine SphOut

!--------------------------------------------
    !Various pole BCs
    !Inner-I BC for cylindrical geometry w/ pole
    subroutine cylpole(bc,Model,Grid,State)
        class(cylindricalPoleBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k,n
        integer :: Np,Np2,jp

        Np = Model%Ring%Np
        Np2 = Np/2 !Halfway about pole

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(k,j,n,jp)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do n=1,Model%Ng
                    jp = j+Np2
                    if (jp > Np) jp = jp-Np !Wrap back
                    !Set through-pole BCs
                    State%Gas(Grid%is -n,j,k,:,:)  = State%Gas(Grid%is+n-1,jp,k,:,:)
                    State%magFlux(Grid%is-n,j,k,IDIR) = -1.0*State%magFlux(Grid%is+n  ,jp,k,IDIR)
                    State%magFlux(Grid%is-n,j,k,JDIR) = -1.0*State%magFlux(Grid%is+n-1,jp,k,JDIR)
                    State%magFlux(Grid%is-n,j,k,KDIR) =      State%magFlux(Grid%is+n-1,jp,k,KDIR)                
                    State%Bxyz(Grid%is -n,j,k,:)    = State%Bxyz(Grid%is+n-1,jp,k,:)
                 enddo
             enddo
         enddo


    end subroutine cylpole

    !Routines for X-axis pole on LFM grid
    subroutine lfmIn(bc,Model,Grid,State)
        class(lfmInBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,k
        integer :: ig,jg,kg,ip,jp,kp
        integer :: Np,Np2

        !i,jg,k = ghost cell
        !i,jp,kp = opposite cell
        Np = Model%Ring%Np
        Np2 = Np/2 !Halfway about pole
        
        !j-boundaries (IN)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,i,k,ig,jg,kg,ip,jp,kp)
        do k=Grid%ksg,Grid%keg
            do i=Grid%isg,Grid%ieg
                do n=1,Model%Ng
                    !Set ghost cells
                    ig = i
                    jg = Grid%js-n
                    kg = k
                    call lfmIJK(Model,Grid,ig,jg,kg,ip,jp,kp)

                    State%Gas    (ig,jg,kg,:,:)  =  State%Gas    (ip,jp  ,kp,:,:)
                    State%Bxyz   (ig,jg,kg,:)    =  State%Bxyz   (ip,jp  ,kp,:)
                    State%magFlux(ig,jg,kg,IDIR) =  State%magFlux(ip,jp  ,kp,IDIR)
                    State%magFlux(ig,jg,kg,JDIR) = -State%magFlux(ip,jp+1,kp,JDIR)
                    State%magFlux(ig,jg,kg,KDIR) = -State%magFlux(ip,jp  ,kp,KDIR)

                enddo
            enddo
        enddo
        
    end subroutine lfmIn

    subroutine lfmOut(bc,Model,Grid,State)
        class(lfmOutBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,k
        integer :: ig,jg,kg,ip,jp,kp
        integer :: Np,Np2

        !i,jg,k = ghost cell
        !i,jp,kp = opposite cell

        Np = Model%Ring%Np
        Np2 = Np/2 !Halfway about pole

        !j-boundaries (OUT)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,i,k,ig,jg,kg,ip,jp,kp)
        do k=Grid%ksg,Grid%keg
            do i=Grid%isg,Grid%ieg
                do n=1,Model%Ng
                    !Set ghost cells
                    ig = i
                    jg = Grid%je+n
                    kg = k
                    call lfmIJK(Model,Grid,ig,jg,kg,ip,jp,kp)

                    State%Gas    (ig,jg  ,kg,:,:)  = State%Gas     (ip,jp,kp,:,:)
                    State%Bxyz   (ig,jg  ,kg,:)    = State%Bxyz    (ip,jp,kp,:)
                    State%magFlux(ig,jg  ,kg,IDIR) =  State%magFlux(ip,jp,kp,IDIR)
                    State%magFlux(ig,jg+1,kg,JDIR) = -State%magFlux(ip,jp,kp,JDIR)
                    State%magFlux(ig,jg  ,kg,KDIR) = -State%magFlux(ip,jp,kp,KDIR)

                enddo
            enddo
        enddo

    end subroutine lfmOut
    
!--------------------------------------------
    !Various exotic BCs

end module bcs
