module chmpfields
    use chmpdefs
    use chmpunits
    use ebtypes
    use ioH5
    use math

    implicit none

    logical, parameter :: doRering = .true.
    integer, parameter :: NumRR = 1

    integer, parameter :: EBIOVARS = 10
    type(IOVAR_T), dimension(EBIOVARS) :: ebIOs

    logical :: chkStatic = .false. !Whether to check for static->time interval switch

    !Generic field initialization routine
    !Update ebState to time t
    abstract interface
        subroutine updateField_T(Model,ebState,t)
            Import :: rp,chmpModel_T,ebState_T
            type(chmpModel_T), intent(in)    :: Model
            type(ebState_T), intent(inout)   :: ebState 
            real(rp), intent(in) :: t    

        end subroutine updateField_T

    end interface

    procedure(updateField_T), pointer :: updateFields !=> ebUpdate

    contains

    !Standard eb update
    subroutine ebUpdate(Model,ebState,t)
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T), intent(inout)   :: ebState 
        real(rp), intent(in) :: t

        logical :: skipUpdate
        integer :: i1,i2
        !Check if doStatic should be flipped
        if (chkStatic .and. (t >= ebState%eb1%time) ) then            
            ebState%doStatic = .false.
            chkStatic = .false.
        endif
        !write(*,*) 'Updating to ', t
        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        skipUpdate = (t >= ebState%eb1%time .and. t <= ebState%eb2%time) .or. ebState%doStatic
        if (t < ebState%eb1%time) then
            !Set to static until time is within interval
            ebState%doStatic = .true.
            skipUpdate = .true.
            chkStatic = .true.
        endif

        if (skipUpdate) return

        !Test for out of time
        if (t>=maxval(ebTab%times(:))) then
            write(*,*) 'Out of time data, switching to static fields ...'
            ebState%doStatic = .true.

            !Copy eb2->eb1
            eb1%time = eb2%time
            eb1%dB   = eb2%dB
            eb1%E    = eb2%E
            if (Model%doMHD) eb1%W = eb2%W

            !Bail out
            return        
        endif

        !Do work if still here
        !Go ahead and reread both (lazy way of avoiding corner cases)
        call findSlc(ebState%ebTab,t,i1,i2)

        !Read eb1
        call readEB(Model,ebGr,ebTab,eb1,ebTab%gStrs(i1))

        !Read eb2
        call readEB(Model,ebGr,ebTab,eb2,ebTab%gStrs(i2))

        end associate

    end subroutine ebUpdate

    !Finds bounding slices from ebTab file
    !NOTE: findloc isn't supported by most gfortran versions, so this is a lazy workaround
    !When gfortran gets its act together, just comment out last few lines using findloc
    subroutine findSlc(ebTab,t,i1,i2)
        type(ebTab_T), intent(in) :: ebTab
        real(rp), intent(in) :: t
        integer, intent(out) :: i1,i2

        integer :: n
    
        !Work-around code        
        do n=1,ebTab%N
            if (ebTab%times(n)>t) exit
        enddo
        i1 = n-1
        i2 = i1+1
        i1 = max(1,i1)
        i2 = min(ebTab%N,i2)

        !Old code
        ! i1 = findloc(ebTab%times .le. t,.true.,dim=1,back=.true.)
        ! i2 = findloc(ebTab%times .gt. t,.true.,dim=1)
        ! i1 = max(1,i1)
        ! i2 = min(ebTab%N,i2)
        

        if (i2 == i1) i2=i1+1 !Possible if none of the tab slices are in range
        !write(*,*) 'i1 / i2 = ', i1,i2
        !write(*,*) 'T(i1) / T / T(i2) = ', oTScl*ebTab%times(i1),oTScl*t,oTScl*ebTab%times(i2)

    end subroutine findSlc

    !Reads specific slice given by group string
    !File: bStr, Group: gStr
    !NOTE: Calculating convective electric field to avoid diffusive terms
    subroutine readEB(Model,ebGr,ebTab,ebF,gStr)
        type(chmpModel_T), intent(in)     :: Model
        type(   ebGrid_T), intent(in)     :: ebGr
        type(ebTab_T), intent(in) :: ebTab
        type(  ebField_T), intent(inout)  :: ebF
        character(len=strLen), intent(in) :: gStr

        character(len=strLen) :: ebFile

        real(rp), dimension(:,:,:), allocatable :: Bx,By,Bz,Vx,Vy,Vz,D,P
        integer :: Nip,Njp,Nkp
        integer :: i,j,k,is,ie,js,je,ks,ke
        integer :: dN(NDIM)
        real(rp), dimension(NDIM) :: B0xyz, Vxyz,Bxyz,Exyz,xcc

        write(*,'(5a)') '<Reading eb from ', trim(ebTab%bStr), '/', trim(gStr), '>'

        !------------
        !Init holders (global # of active cells)
        Nip = ebGr%Nip
        Njp = ebGr%Njp
        Nkp = ebGr%Nkp

        allocate(Bx(Nip,Njp,Nkp))
        allocate(By(Nip,Njp,Nkp))
        allocate(Bz(Nip,Njp,Nkp))

        allocate(Vx(Nip,Njp,Nkp))
        allocate(Vy(Nip,Njp,Nkp))
        allocate(Vz(Nip,Njp,Nkp))
        if (Model%doMHD) then
            allocate(D(Nip,Njp,Nkp))
            allocate(P(Nip,Njp,Nkp))
        endif

        !------------
        !Get data from files
        dN = [ebTab%dNi,ebTab%dNj,ebTab%dNk] !Individual blocks
        do k=1,ebTab%Rk
            do j=1,ebTab%Rj
                do i=1,ebTab%Ri
                    ebFile = genName(ebTab,i,j,k)

                    !Get piece from file
                    call ClearIO(ebIOs)
                    call AddInVar(ebIOs,"Bx")
                    call AddInVar(ebIOs,"By")
                    call AddInVar(ebIOs,"Bz")
                    call AddInVar(ebIOs,"Vx")
                    call AddInVar(ebIOs,"Vy")
                    call AddInVar(ebIOs,"Vz")
                    
                    if (Model%doMHD) then
                        call AddInVar(ebIOs,"D")
                        call AddInVar(ebIOs,"P")
                    endif

                    call AddInVar(ebIOs,"time",vTypeO=IOREAL)
                    call ReadVars(ebIOs,.true.,ebFile,gStr)

                    !Push piece to grid
                    is = (i-1)*ebTab%dNi + 1
                    js = (j-1)*ebTab%dNj + 1
                    ks = (k-1)*ebTab%dNk + 1
                    ie = is + ebTab%dNi - 1
                    je = js + ebTab%dNj - 1
                    ke = ks + ebTab%dNk - 1

                    Bx(is:ie,js:je,ks:ke) = inBScl*reshape(ebIOs(1)%data,dN)
                    By(is:ie,js:je,ks:ke) = inBScl*reshape(ebIOs(2)%data,dN)
                    Bz(is:ie,js:je,ks:ke) = inBScl*reshape(ebIOs(3)%data,dN)
                    Vx(is:ie,js:je,ks:ke) = inVScl*reshape(ebIOs(4)%data,dN)
                    Vy(is:ie,js:je,ks:ke) = inVScl*reshape(ebIOs(5)%data,dN)
                    Vz(is:ie,js:je,ks:ke) = inVScl*reshape(ebIOs(6)%data,dN)
                    if (Model%doMHD) then
                        !Convert incoming data to [#/cc] and [nPa] using defined scaling params
                        D(is:ie,js:je,ks:ke) = inDScl*reshape(ebIOs(7)%data,dN)
                        P(is:ie,js:je,ks:ke) = inPScl*reshape(ebIOs(8)%data,dN)
                    endif
                enddo
            enddo
        enddo
        !Get time data from last file
        i = FindIO(ebIOs,"time")
        ebF%time = inTScl*ebIOs(i)%data(1)

        !------------
        !Pull data into ebField
        ebF%dB = 0.0
        ebF%E  = 0.0
        if (Model%doMHD) then
            ebF%W = 0.0
        endif

        !FIXME: Lazily assuming is=1,ie=Nxp
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(B0xyz,Bxyz,Vxyz,Exyz)    
        do k=1,Nkp
            do j=1,Njp
                do i=1,Nip
                    !Get residual field
                    B0xyz = ebGr%B0cc(i,j,k,:) !Background field

                    Bxyz = [Bx(i,j,k),By(i,j,k),Bz(i,j,k)]
                    Vxyz = [Vx(i,j,k),Vy(i,j,k),Vz(i,j,k)]

                    ebF%dB(i,j,k,:) = Bxyz-B0xyz

                    !Get electric field from V/B
                    Exyz = -cross(Vxyz,Bxyz)
                    ebF%E(i,j,k,:) = Exyz
                    if (Model%doMHD) then
                        ebF%W(i,j,k,DEN) = D(i,j,k)
                        ebF%W(i,j,k,VELX:VELZ) = Vxyz
                        ebF%W(i,j,k,PRESSURE) = P(i,j,k)
                    endif
                enddo
            enddo
        enddo
        !------------
        !Fill in ghosts
        call ebGhosts(Model,ebGr,ebF)
        
        !------------
        !Clean up
        call ClearIO(ebIOs)
        deallocate(Bx,By,Bz,Vx,Vy,Vz)
        if (Model%doMHD) deallocate(D,P)

    end subroutine readEB

    !Fill in ghost cells for EB data
    subroutine ebGhosts(Model,ebGr,ebF)
        type(chmpModel_T), intent(in)    :: Model
        type(   ebGrid_T), intent(in)    :: ebGr
        type(  ebField_T), intent(inout) :: ebF

        integer :: n,d,i,j,k,ip,jp,kp
        integer :: Nk,Nk2
        logical :: iGh,jGh,kGh

        select case(ebGr%GrID)
        case(LFMGRID,EGGGRID)
            if (doRering) then
                !Average ring
                !NumRR = # of rings to average over
                Nk = ebGr%Nkp
                do n=0,NumRR-1
                    !$OMP PARALLEL DO default(shared) &
                    !$OMP private(i,d)
                    do i=ebGr%is,ebGr%ie
                        !Fields
                        do d=1,NDIM
                            ebF%dB(i,ebGr%js+n,ebGr%ks:ebGr%ke,d) = sum(ebF%dB(i,ebGr%js+n,ebGr%ks:ebGr%ke,d))/Nk
                            ebF%E (i,ebGr%js+n,ebGr%ks:ebGr%ke,d) = sum(ebF%E (i,ebGr%js+n,ebGr%ks:ebGr%ke,d))/Nk

                            ebF%dB(i,ebGr%je-n,ebGr%ks:ebGr%ke,d) = sum(ebF%dB(i,ebGr%je-n,ebGr%ks:ebGr%ke,d))/Nk
                            ebF%E (i,ebGr%je-n,ebGr%ks:ebGr%ke,d) = sum(ebF%E (i,ebGr%je-n,ebGr%ks:ebGr%ke,d))/Nk
                        enddo
                        !MHD variables
                        if (Model%doMHD) then
                            do d=1,NVARMHD
                                ebF%W (i,ebGr%js+n,ebGr%ks:ebGr%ke,d) = sum(ebF%W (i,ebGr%js+n,ebGr%ks:ebGr%ke,d))/Nk
                                ebF%W (i,ebGr%je-n,ebGr%ks:ebGr%ke,d) = sum(ebF%W (i,ebGr%je-n,ebGr%ks:ebGr%ke,d))/Nk
                            enddo

                        endif
                    enddo !i-loop
                enddo
            endif

            !Now do ghosts, lazily loop through *EVERYTHING*
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,ip,jp,kp,iGh,jGh,kGh)
            do k=ebGr%ksg,ebGr%keg
                do j=ebGr%jsg,ebGr%jeg
                    do i=ebGr%isg,ebGr%ieg
                        iGh = (i<ebGr%is) .or. (i>ebGr%ie)
                        jGh = (j<ebGr%js) .or. (j>ebGr%je)
                        kGh = (k<ebGr%ks) .or. (k>ebGr%ke)
                        if (iGh .or. jGh .or. kGh) then
                            call ijk2Active(Model,ebGr,i,j,k,ip,jp,kp)
                            ebF%dB(i,j,k,:) = ebF%dB(ip,jp,kp,:)
                            ebF%E (i,j,k,:) = ebF%E (ip,jp,kp,:)
                            if (Model%doMHD) then
                                ebF%W (i,j,k,:) = ebF%W (ip,jp,kp,:)
                            endif
                        endif
                    enddo
                enddo
            enddo
                    
        case(CARTGRID)
            write(*,*) 'Cartesian not implemented'
            stop
        end select

    end subroutine ebGhosts


    !Takes i,j,k cell index and returns active cell ip,jp,kp of mirror
    !Map in i,k,j order
    subroutine ijk2Active(Model,Grid,i,j,k,ip,jp,kp)
        type(chmpModel_T), intent(in)    :: Model
        type(   ebGrid_T), intent(in)    :: Grid
        integer, intent(in) :: i,j,k
        integer, intent(out) :: ip,jp,kp

        integer :: Np,Np2

        Np  = Grid%Nkp
        Np2 = Grid%Nkp/2

        !Start w/ i index, do mirror back into active
        if (i < Grid%is) then
            ip = Grid%is + (Grid%is-i) - 1
        elseif (i > Grid%ie) then
            ip = Grid%ie - (i-Grid%ie) + 1
        else
            ip = i
        endif

        !Next do k, map via periodicity
        if (k < Grid%ks) then
            kp = Grid%ke - (Grid%ks-k) + 1
        elseif (k > Grid%ke) then
            kp = Grid%ks + (k-Grid%ke) - 1
        else
            kp = k
        endif

        !Finally do j
        if (j < Grid%js) then
            jp = Grid%js + (Grid%js-j) - 1
            kp = k+Np2
            if (kp>Np) kp = kp-Np
        elseif (j > Grid%je) then
            jp = Grid%je - (j-Grid%je) + 1
            kp = k+Np2
            if (kp>Np) kp = kp-Np
        else
            jp = j
        endif
    end subroutine ijk2Active

    !Generate name of output file based on tiling
    function genName(ebTab,i,j,k) result(fName)
        type(ebTab_T), intent(in) :: ebTab
        integer, intent(in) :: i,j,k
        character(len=strLen) :: fName

        character(len=strLen) :: fHd,fRn,fijk,fTl
        integer :: n

        if (ebTab%isMPI) then
            n = (j-1) + (i-1)*ebTab%Rj + (k-1)*ebTab%Ri*ebTab%Rj
            write(fHd ,'(a,a)') trim(ebTab%bStr), '_'
            write(fRn ,'(I0.4,a,I0.4,a,I0.4,a)') ebTab%Ri,'_',ebTab%Rj,'_',ebTab%Rk,'_'
            write(fijk,'(I0.4,a,I0.4,a,I0.4,a)') i-1,'_',j-1,'_',k-1,'_'
            write(fTl ,'(I0.12,a)') n,'.h5'
            
            fName = trim(fHd) // trim(fRn) // trim(fijk) // trim(fTl)
        else
            fName = ebTab%bStr
        endif
        !write(*,*) 'ijk / file = ',i,j,k,trim(fName)
    end function genName
end module chmpfields
