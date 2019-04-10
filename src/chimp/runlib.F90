!Driver for stand-alone field line tracing
module runlib_chimp
    use clocks
    use chmpdefs
    use chmpio
    use ebtypes
    use starter
    use streamline
    use chmpfields
    use math

#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    !Main data structures
    type(chmpModel_T) :: Model
    type(ebState_T)   :: ebState
    type(XML_Input_T) :: inpXML

    integer :: n,NumP

    real(rp), allocatable :: Xs(:,:)
    ! store this so we don't have to back engineer the unrolled index vs i,j,k
    integer,  allocatable :: mapind(:,:)    ! unrolled index,DIR (i,j,k)
    real(rp),  allocatable :: eqmap(:,:,:,:)     ! equatorial mapping of cell ijk (i,j,k,DIR)
    type(fLine_T), allocatable :: fLs(:)
    real(rp) :: wT

    contains 


    subroutine setModel()
      character(len=strLen) :: xmlStr

      !Create input XML object
      call getIDeckStr(xmlStr)
      inpXML = New_XML_Input(trim(xmlStr),"Chimp",.true.)
      
      !----------------------------
      !Get information for model data structure
      write(*,*) '----------------------------'
      write(*,*) 'CHIMP: Initializing model ...'
      call setUnits (Model,inpXML)
      call initModel(Model,inpXML)
      call setBackground(Model,inpXML) !Background field info
    end subroutine setModel

    subroutine setEBState(gridCorners)
      real(rp), dimension(:,:,:,:), intent(in) :: gridCorners

      write(*,*) '----------------------------'
      write(*,*) 'CHIMP: Initializing fields ...'
        
      call ebInit_fromMHDGrid(Model,ebState,inpXML,gridCorners)
      
      !Also do localization init
      call InitLoc(Model,ebState%ebGr,inpXML)
    end subroutine setEBState

    subroutine Btot2ebState(Btot)
      real(rp), dimension(:,:,:,:), intent(in) :: Btot
      
      integer :: i,j,k
      real(rp) :: B0xyz(NDIM)

      associate( ebGr=>ebState%ebGr,eb1=>ebState%eb1,eb2=>ebState%eb2 )

      !$OMP PARALLEL DO default(shared) collapse(2) &
      !$OMP private(B0xyz) 
      do k=ebGr%ks,ebGr%ke
         do j=ebGr%js,ebGr%je
            do i=ebGr%is,ebGr%ie
               !Get residual field
               B0xyz = ebGr%B0cc(i,j,k,:) !Background field
               eb1%dB(i,j,k,:) = Btot(i,j,k,:) - B0xyz
               eb2%dB(i,j,k,:) = Btot(i,j,k,:) - B0xyz  ! NOTE, assuming static
            enddo
         enddo
      enddo
        
      !NOTE: Need to fill in ghosts here
      call ebGhosts(Model,ebGr,eb1)
      call ebGhosts(Model,ebGr,eb2)

      end associate
    end subroutine Btot2ebState

    ! intended to replace getChimp
    ! assumes data passed from MHD for now
    ! needs modification to adapt for standalone CHIMP
    subroutine ChimpInit(gridCorners)
      real(rp), dimension(:,:,:,:), intent(in) :: gridCorners

      call setModel()
      call setEBState(gridCorners) ! note, this points to a specific implementation of initFields -- modify for standalone Chimp

      print *,'---- DONE CHIMP INIT ----'

      !Setup timers
      call initClocks()

      Model%doFLOut = .false.

      call genPtsFromGrid()

      ! debug stuff
      ! print *,'------- CHIMP INIT ------'
      ! print *,'xmin,xmax: ', minval(gridCorners(:,:,:,1)),maxval(gridCorners(:,:,:,1))
      ! print *,'ymin,ymax: ', minval(gridCorners(:,:,:,2)),maxval(gridCorners(:,:,:,2))
      ! print *,'zmin,zmax: ', minval(gridCorners(:,:,:,3)),maxval(gridCorners(:,:,:,3))

      ! print *,'xmin,xmax: ', minval(ebState%ebGr%xyz(:,:,:,1)),maxval(ebState%ebGr%xyz(:,:,:,1))
      ! print *,'ymin,ymax: ', minval(ebState%ebGr%xyz(:,:,:,2)),maxval(ebState%ebGr%xyz(:,:,:,2))
      ! print *,'zmin,zmax: ', minval(ebState%ebGr%xyz(:,:,:,3)),maxval(ebState%ebGr%xyz(:,:,:,3))

      Model%t = Model%T0
    end subroutine ChimpInit

    subroutine ChimpProject()
      real(rp) :: xout(NDIM)
      integer :: i,j,k
      integer :: is,ie,js,je,ks,ke  ! aliases for brevity

      is = ebState%ebGr%is
      ie = ebState%ebGr%ie
      js = ebState%ebGr%js
      je = ebState%ebGr%je
      ks = ebState%ebGr%ks
      ke = ebState%ebGr%ke

      ! clean map first
      !$OMP PARALLEL DO default(shared) collapse(2) 
      do k=ks,ke
         do j=js,je
            do i=is,ie
               eqmap(i,j,k,:) = [-999.,-999.,-999.]
            enddo
         enddo
      enddo

      !$OMP PARALLEL DO default(shared) &
      !$OMP& private(xout) &
      !$OMP& schedule(dynamic)
      do n=1,NumP
         ! trace only below certain colatitude
         if ( abs(acos(Xs(n,ZDIR)/norm2(Xs(n,:)))*180./pi) >= 30.  ) then
            call getProjection(Model,ebState,Xs(n,:),Model%t,xout)
            eqmap(mapind(n,IDIR),mapind(n,JDIR),mapind(n,KDIR),:) = xout
         end if
      enddo

      ! after the above is done, eqmap has either been project to equator or set to [-999.,-999.,-999.] 
      ! the latter can happen for 3 reasons: 
      ! 1. it's outside of seed domain (set by genPtsFromGrid() function). 
      ! 2. it was not traced (e.g., by the adhoc seed point limiter above tracing only poleward of 30 deg)
      ! 3. it did not project to equator (terminated at domain boundary, i.e., the FL is open)
    end subroutine ChimpProject


    subroutine ChimpGetStream()

    !Do main loop
    do while (Model%t<=Model%tFin)
        call Tic("Omega")

        call Tic("Step")
        !Update fields to current time
        call updateFields(Model,ebState,Model%t)
        call Toc("Step")

        
        !Trace field lines
        call Tic("Tracer")
        !$OMP PARALLEL DO
        do n=1,NumP
            call genStream(Model,ebState,Xs(n,:),Model%t,fLs(n))
        enddo
        call Toc("Tracer")

        call Tic("Output")

        !Write lines
        call fOutput(Model,fLines=fLs)

        if (modulo(Model%ts,Model%tsOut) ==0) then
            wT = readClock("Tracer")
            write(*,'(a,f8.3,a)') 'T = ', Model%t*oTScl, ' ' // trim(tStr)
            write(*,'(a,f8.3)') '   kFLps = ', 1.0e-3*NumP*Model%tsOut/wT
        endif

        call Toc("Output")


        call Toc("Omega")
        !Timing book keeping
        if (modulo(Model%ts,Model%tsOut) ==0) then
            if (Model%doTimer) call printClocks()
            call cleanClocks()
        endif

        !Update time
        Model%t = Model%t + Model%dt
        Model%ts = Model%ts+1

    enddo

    end subroutine ChimpGetStream

    !Generate points to use for tracing
    !points/grType = RP/TP (polar/spherical shell)
    subroutine genPts()
      
      real(rp) :: rIn
      integer :: i
      real(rp), allocatable, dimension(:) :: x1,x2
      character(len=strLen) :: grType
      !Guess at minimum radius
      rIn = norm2(ebState%ebGr%xyz(2,1,1,:))
      
      !Start w/ number of points (use global NumP)
      call inpXML%Set_Val(NumP,'points/Np',100)
      
      if (allocated(Xs)) deallocate(Xs)
      allocate(Xs(NumP,NDIM))
      
      !Generate based on points/grType
      call inpXML%Set_Val(grType,'points/grType',"RP")
      write(*,*) 'Grid type = ', trim(toUpper(grType))

      select case(trim(toUpper(grType)))
      case("RP","POLAR")
         write(*,*) 'Generating polar grid ...'
         call getSample(Model,inpXML,"radius",NumP,x1,4.0_rp, 10.0_rp)
         call getSample(Model,inpXML,"phi"   ,NumP,x2,0.0_rp,360.0_rp)
         x2 = x2/rad2deg
         do i=1,NumP
            Xs(i,XDIR) = x1(i)*cos(x2(i))
            Xs(i,YDIR) = x1(i)*sin(x2(i))
            Xs(i,ZDIR) = 0.0
         enddo
      case("TP","SPHERICAL")
         call getSample(Model,inpXML,"theta" ,NumP,x1,0.0_rp,180.0_rp)
         call getSample(Model,inpXML,"phi"   ,NumP,x2,0.0_rp,360.0_rp)
         x1 = x1/rad2deg
         x2 = x2/rad2deg
         do i=1,NumP
            Xs(i,XDIR) = rIn*cos(x2(i))*sin(x1(i))
            Xs(i,YDIR) = rIn*sin(x2(i))*sin(x1(i))
            Xs(i,ZDIR) = rIn*cos(x1(i))
         enddo
      end select
      
    end subroutine genPts

    ! a version of the above but without randomizing the seed points
    ! only works for a wedge of points in equatorial plane
    subroutine genOrderedPts()
      real(rp) :: rmin,rmax,pmin,pmax
      integer :: i, j, np, nr 
      real(rp), allocatable, dimension(:) :: x1,x2
      character(len=strLen) :: grType

      !Start w/ number of points (use global NumP)
      call inpXML%Set_Val(np,'points/Np',100)  ! #points in azimuth
      call inpXML%Set_Val(nr,'points/Nr',100)  ! #points in azimuth
      NumP = nr*np

      if (allocated(Xs)) deallocate(Xs)
      allocate(Xs(NumP,NDIM))

      if(.not.(allocated(x1))) allocate(x1(nr))
      if(.not.(allocated(x2))) allocate(x2(np))

      !Generate based on points/grType
      write(*,*) 'Generating polar grid ...'

      call inpXML%Set_Val(rmin,'radius/min',5.0)
      call inpXML%Set_Val(rmax,'radius/max',25.0)
      call inpXML%Set_Val(pmin,'phi/min',0.0)
      call inpXML%Set_Val(pmax,'phi/max',360.0)

      x1 = rmin + (rmax-rmin)/nr*[(i, i = 0, nr)]
      x2 = pmin + (pmax-pmin)/np*[(i, i = 0, np)]
      x2 = x2/rad2deg

      do i=1,np
         Xs((i-1)*nr+1:i*nr,XDIR) = x1*cos(x2(i))
         Xs((i-1)*nr+1:i*nr,YDIR) = x1*sin(x2(i))
         Xs((i-1)*nr+1:i*nr,ZDIR) = 0.0
      enddo

    end subroutine genOrderedPts

    subroutine genPtsFromGrid()
      real(rp) :: rmin,rmax
      integer :: i,j,k,n
      integer :: ni, nj, nk, ind
      integer :: is,ie,js,je,ks,ke,isg,ieg,jsg,jeg,ksg,keg  ! aliases for brevity

      ni = ebState%ebGr%Nip
      nj = ebState%ebGr%Njp
      nk = ebState%ebGr%Nkp

      is = ebState%ebGr%is
      ie = ebState%ebGr%ie
      js = ebState%ebGr%js
      je = ebState%ebGr%je
      ks = ebState%ebGr%ks
      ke = ebState%ebGr%ke

      isg = ebState%ebGr%isg
      ieg = ebState%ebGr%ieg
      jsg = ebState%ebGr%jsg
      jeg = ebState%ebGr%jeg
      ksg = ebState%ebGr%ksg
      keg = ebState%ebGr%keg

      ! allocate storage for output and passage to Gamera
      if (.not.allocated(eqmap)) allocate(eqmap(is:ie,js:je,ks:ke,NDIM))

      call inpXML%Set_Val(rmin,'radius/min',5.0)
      call inpXML%Set_Val(rmax,'radius/max',25.0)

      Model%rmin = rmin
      Model%rmax = rmax

      ! find nightside i-location just inside the domain
      ! note minloc assumes indices start from 1, thus the -1+isg at the end
      ind=minloc(rmax+ebState%ebGr%xyz(:,je+1,ks,XDIR),dim=1,mask=rmax+ebState%ebGr%xyz(:,je+1,ks,XDIR)>=0.)-1+isg  
      ni = ind-is+1
      NumP = ni*nj*nk

      if (.not.allocated(mapind)) allocate(mapind(NumP,NDIM))
      if (.not.allocated(Xs)) allocate(Xs(NumP,NDIM))

      do i=is,ni
         do j=js,je
            do k = ks,ke
               n = (i-is)*nj*nk+(j-js)*nk+k
               Xs(n, :) = ebState%ebGr%xyzcc(i,j,k,:)
               ! store this so we don't have to backengineer later
               mapind(n,:) = [i,j,k]
            end do
         end do
      end do
    end subroutine genPtsFromGrid
end module runlib_chimp
