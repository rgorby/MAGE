!Driver for stand-alone field line tracing
program projectx
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

    real(rp), allocatable :: Xs(:,:), Xe(:,:)
    real(rp) :: xout(NDIM)
    type(fLine_T), allocatable :: fLs(:)
    real(rp) :: wT

!    write(*,*) "Num threads=",omp_get_max_threads()
    !------------
    !Setup timers
    call initClocks()

    !----------------------------
    !Initialize model and fields
    call goApe(Model,ebState,iXML=inpXML)
    Model%doFLOut = .true.
    
    !Get tracing points and field line holders
!    call genPts()
!    call genOrderedPts()
    call genPtsFromGrid()
    allocate(fLs(NumP))

    !Loop from T0 -> tFin
    Model%t = Model%T0

    write(*,*) ''
    write(*,'(a,I0,a)')  'Tracing ', NumP, ' points'
    write(*,'(a,f8.3,a,f8.3)') 'Time Interval = ', Model%T0*oTScl,' / ', Model%tFin*oTScl
    write(*,'(a,f8.3)') 'dt = ', Model%dt*oTScl

    open (unit = 111, file = "output")

    !Do main loop
    do while (Model%t<=Model%tFin)
        call Tic("Omega")

        call Tic("Step")
        !Update fields to current time
        call updateFields(Model,ebState,Model%t)
        call Toc("Step")

        
        !Trace field lines
        call Tic("Tracer")
!        !$OMP PARALLEL DO
        !$OMP PARALLEL DO default(shared) &
        !$OMP& private(xout) &
        !$OMP& schedule(dynamic)
        do n=1,NumP
           ! trace only below certain colatitude
           if ( abs(acos(Xs(n,ZDIR)/norm2(Xs(n,:)))*180./pi) >= 30.  ) then
              call getEquatorProjection(Model,ebState,Xs(n,:),Model%t,xout)
              Xe(n,:) = xout
           else
              Xe(n,:) = [-999.,-999.,-999.]
           end if
        enddo

        ! do n=1,NumP
        !     write(111,"(I5,6F10.2)") n,Xs(n,1),Xs(n,2),Xs(n,3),Xe(n,1),Xe(n,2),Xe(n,3)
        ! enddo

        call Toc("Tracer")


        call Tic("Output")

        !Write lines
!        call fOutput(Model,fLines=fLs)

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

    close(111)
    contains 

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
          integer :: i,j,k
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
          
          call inpXML%Set_Val(rmin,'radius/min',5.0)
          call inpXML%Set_Val(rmax,'radius/max',25.0)

          ! find nightside i-location just inside the domain
          ! note minloc assumes indices start from 1, thus the -1+isg at the end
          ind=minloc(rmax+ebState%ebGr%xyz(:,je+1,ks,XDIR),dim=1,mask=rmax+ebState%ebGr%xyz(:,je+1,ks,XDIR)>=0.)-1+isg  
          ni = ind-is+1
          NumP = ni*nj*nk
          if (allocated(Xs)) deallocate(Xs)
          allocate(Xs(NumP,NDIM))
          if (allocated(Xe)) deallocate(Xe)
          allocate(Xe(NumP,NDIM))

          do i=is,ni
             do j=js,je
                do k = ks,ke
                   Xs( (i-is)*nj*nk+(j-js)*nk+k, :) = ebState%ebGr%xyzcc(i,j,k,:)
                end do
             end do
          end do
        end subroutine genPtsFromGrid

end program projectx
