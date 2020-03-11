module gamdebug
	use kdefs
	use gamtypes

    implicit none

    logical, private :: doVerboseDB = .true.
    logical, private :: didFix = .false.
    logical, private :: firstFix = .false.
    contains


    !Checks magfluxes on K-periodic boundary
    subroutine ChkEFieldLFM(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State

        integer :: i,j,k
        real(rp) :: dEi,dEj

        if (.not. firstFix) return
        
        if (doVerboseDB) write(*,*) 'Checking E field ...'
        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dEi = abs(State%Efld(i,j,Grid%ks,IDIR) - State%Efld(i,j,Grid%ke+1,IDIR))
                dEj = abs(State%Efld(i,j,Grid%ks,JDIR) - State%Efld(i,j,Grid%ke+1,JDIR))

                if ( (dEi+dEj)> TINY ) then
                    if (doVerboseDB) write(*,*) 'E: i/j, dEi/dEj = ',i,j,dEi,dEj
                endif

            enddo
        enddo

    end subroutine ChkEFieldLFM

    !Checks magfluxes on K-periodic boundary
    subroutine ChkFluxLFM(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State

        integer :: i,j,k
        real(rp) :: dB

        if (.not. firstFix) return
        if (doVerboseDB) write(*,*) 'Checking Bk ...'
        do i=Grid%is,Grid%ie
            do j=Grid%js,Grid%je
                dB = abs(State%magFlux(i,j,Grid%ks,KDIR)-State%magFlux(i,j,Grid%ke+1,KDIR))
                if ( dB > TINY ) then
                    if (doVerboseDB) write(*,*) 'Bk: i/j, Bk = ',i,j,dB
                endif
            enddo
        enddo

    end subroutine ChkFluxLFM

    !Checks magfluxes on K-periodic boundary
    subroutine ChkGasFluxLFM(Model,Grid,gFlx,mFlx)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), intent(inout)           :: gFlx(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(inout), optional :: mFlx(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM,1:NDIM)
        

        integer :: i,j,k
        real(rp) :: dFg,dFm,dB

        if (.not. firstFix) return

        if (doVerboseDB) write(*,*) 'Checking stresses ...'
        do i=Grid%is,Grid%ie
            do j=Grid%js,Grid%je
                dFg = norm2( gFlx(i,j,Grid%ks,:,KDIR,:) - gFlx(i,j,Grid%ke+1,:,KDIR,:) )
                dFm = norm2( mFlx(i,j,Grid%ks,:,KDIR  ) - mFlx(i,j,Grid%ke+1,:,KDIR  ) )

                dB = dFg+dFm
                if ( dB > TINY ) then
                    if (doVerboseDB) write(*,*) 'Fluxes: i/j, dFg,dFm = ',i,j,dFg,dFm

                endif
            enddo
        enddo

    end subroutine ChkGasFluxLFM

    !Fixes gas fluxes on K-periodic boundary
    subroutine FixGasFluxLFM(Model,Grid,gFlx,mFlx)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), intent(inout)           :: gFlx(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(inout), optional :: mFlx(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM,1:NDIM)

        integer :: i,j,k
        real(rp), dimension(NDIM) :: mF
        real(rp), dimension(NVAR,BLK:Model%nSpc) :: gF
        logical :: doMaxwell

        if (present(mFlx)) doMaxwell = .true.

        if ( Grid%hasLowerBC(KDIR) .and. Grid%hasUpperBC(KDIR) ) then

            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,mF,gF)            
            do i=Grid%is,Grid%ie
                do j=Grid%js,Grid%je
                    !Gas fluxes
                    gF = 0.5*( gFlx(i,j,Grid%ks,:,KDIR,:) + gFlx(i,j,Grid%ke+1,:,KDIR,:) )
                    gFlx(i,j,Grid%ks  ,:,KDIR,:) = gF
                    gFlx(i,j,Grid%ke+1,:,KDIR,:) = gF
                    
                    if (doMaxwell) then
                        mF = 0.5*( mFlx(i,j,Grid%ks,:,KDIR) + mFlx(i,j,Grid%ke+1,:,KDIR) )
                        mFlx(i,j,Grid%ks  ,:,KDIR) = mF
                        mFlx(i,j,Grid%ke+1,:,KDIR) = mF
                    endif

                enddo
            enddo
        endif !Both low-K and hi-K

    end subroutine FixGasFluxLFM

    !Fixes magfluxes on K-periodic boundary
    subroutine FixFluxLFM(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k
        real(rp) :: dB

        if (firstFix) return

        do i=Grid%is,Grid%ie
            do j=Grid%js,Grid%je
                dB = 0.5*(State%magFlux(i,j,Grid%ks,KDIR) + State%magFlux(i,j,Grid%ke+1,KDIR))
                State%magFlux(i,j,Grid%ks  ,KDIR) = dB
                State%magFlux(i,j,Grid%ke+1,KDIR) = dB

             enddo
        enddo
        firstFix = .true.

    end subroutine FixFluxLFM

    !Force agreement of E fields across K-periodic boundary
    subroutine FixEFieldLFM(Model,Gr,E)
        type(Model_T), intent(in)    :: Model
        type(Grid_T) , intent(in)    :: Gr
        real(rp)     , intent(inout) :: E(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM)

        integer :: i,j,k
        real(rp) :: Ei,Ej

        if ( Gr%hasLowerBC(KDIR) .and. Gr%hasUpperBC(KDIR) ) then

            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,Ei,Ej)            
            do i=Gr%is,Gr%ie+1
                do j=Gr%js,Gr%je+1
                    Ei = 0.5*(E(i,j,Gr%ks,IDIR) + E(i,j,Gr%ke+1,IDIR))
                    Ej = 0.5*(E(i,j,Gr%ks,JDIR) + E(i,j,Gr%ke+1,JDIR))
                    
                    E(i,j,Gr%ks  ,IDIR) = Ei
                    E(i,j,Gr%ke+1,IDIR) = Ei
                    E(i,j,Gr%ks  ,JDIR) = Ej
                    E(i,j,Gr%ke+1,JDIR) = Ej
                enddo
            enddo

        endif !Both low-K and hi-K

    end subroutine FixEFieldLFM

	!Checks periodicity for metric terms on LFM-style grid
    subroutine ChkMetricLFM(Model,Grid)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid

        integer :: i,j,k

        real(rp) :: dEi,dEj
        real(rp) :: dFk

        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dEi = sum(abs(Grid%Te(i,j,Grid%ks,:,IDIR) - Grid%Te(i,j,Grid%ke+1,:,IDIR)))
                dEj = sum(abs(Grid%Te(i,j,Grid%ks,:,JDIR) - Grid%Te(i,j,Grid%ke+1,:,JDIR)))

                if ( (dEi+dEj)> TINY ) then
                    if (doVerboseDB) write(*,*) 'Te: i/j, dEi/dEj = ',i,j,dEi,dEj
                endif
            enddo
        enddo

        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dEi = sum(abs(Grid%Teb(i,j,Grid%ks,:,IDIR) - Grid%Teb(i,j,Grid%ke+1,:,IDIR)))
                dEj = sum(abs(Grid%Teb(i,j,Grid%ks,:,JDIR) - Grid%Teb(i,j,Grid%ke+1,:,JDIR)))

                if ( (dEi+dEj)> TINY ) then
                    if (doVerboseDB) write(*,*) 'Teb: i/j, dEi/dEj = ',i,j,dEi,dEj
                endif
            enddo
        enddo

        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dEi = abs(Grid%edge(i,j,Grid%ks,IDIR) - Grid%edge(i,j,Grid%ke+1,IDIR))
                dEj = abs(Grid%edge(i,j,Grid%ks,JDIR) - Grid%edge(i,j,Grid%ke+1,JDIR))

                if ( (dEi+dEj)> TINY ) then
                    if (doVerboseDB) write(*,*) 'edge: i/j, dEi/dEj = ',i,j,dEi,dEj
                endif
            enddo
        enddo

        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dFk = abs(Grid%face(i,j,Grid%ks,KDIR) - Grid%face(i,j,Grid%ke+1,KDIR))
                if (dFk > TINY) then
                    if (doVerboseDB) write(*,*) 'fA: i/j, dFk = ',i,j,dFk
                endif
            enddo
        enddo

        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dFk = sum(abs(Grid%Tf(i,j,Grid%ks,:,KDIR) - Grid%Tf(i,j,Grid%ke+1,:,KDIR)))
                if (dFk > TINY) then
                    if (doVerboseDB) write(*,*) 'Tf: i/j, dFk = ',i,j,dFk
                endif
            enddo
        enddo

        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dFk = sum(abs(Grid%xfc(i,j,Grid%ks,:,KDIR) - Grid%xfc(i,j,Grid%ke+1,:,KDIR)))
                
                if (dFk > TINY) then
                    if (doVerboseDB) write(*,*) 'xfc: i/j, dFk = ',i,j,dFk
                endif
            enddo
        enddo

        if (Model%doBackground) then
            do i=Grid%is,Grid%ie+1
                do j=Grid%js,Grid%je+1
                    dFk = sum(abs(Grid%fcB0(i,j,Grid%ks,:,KDIR) - Grid%fcB0(i,j,Grid%ke+1,:,KDIR)))
                    if (dFk > TINY) then
                        if (doVerboseDB) write(*,*) 'fcB0: i/j, dFk = ',i,j,dFk
                    endif
                enddo
            enddo

            do i=Grid%is,Grid%ie+1
                do j=Grid%js,Grid%je+1
                    dFk = abs(Grid%bFlux0(i,j,Grid%ks,KDIR) - Grid%bFlux0(i,j,Grid%ke+1,KDIR))
                    if (dFk > TINY) then
                        if (doVerboseDB) write(*,*) 'bFlx0: i/j, dFk = ',i,j,dFk
                    endif
                enddo
            enddo

        endif      
    end subroutine ChkMetricLFM

end module gamdebug