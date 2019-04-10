!Driver for stand-alone field line tracing
program tracex
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
    type(fLine_T), allocatable :: fLs(:)
    real(rp) :: wT

    !write(*,*) "Num threads=",omp_get_max_threads()
    !------------
    !Setup timers
    call initClocks()

    !----------------------------
    !Initialize model and fields
    call getChimp(Model,ebState,iXML=inpXML)
    Model%doFLOut = .true.
    
    !Get tracing points and field line holders
!    call genPts()
    call genOrderedPts()
    !call genPtsFromGrid()
    allocate(fLs(NumP))

    !Loop from T0 -> tFin
    Model%t = Model%T0

    write(*,*) ''
    write(*,'(a,I0,a)')  'Tracing ', NumP, ' points'
    write(*,'(a,f8.3,a,f8.3)') 'Time Interval = ', Model%T0*oTScl,' / ', Model%tFin*oTScl
    write(*,'(a,f8.3)') 'dt = ', Model%dt*oTScl

    !Do main loop
    do while (Model%t<=Model%tFin)
        call Tic("Omega")

        call Tic("Step")
        !Update fields to current time
        call updateFields(Model,ebState,Model%t)
        call Toc("Step")

        
        !Trace field lines
        call Tic("Tracer")
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic)
        do n=1,NumP
            call genStream(Model,ebState,Xs(n,:),Model%t,fLs(n))
        enddo
        call Toc("Tracer")

        call Tic("Output")

        !Write lines
        call fOutput(Model,ebState=ebState,fLines=fLs)

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
            real(rp) :: rmin,rmax,pmin,pmax,tmin,tmax
            real(rp) :: x1min,x1max,x2min,x2max
            integer :: i, j, k,Nx1,Nx2,Nx3
            real(rp), allocatable, dimension(:) :: x1,x2,x3
            character(len=strLen) :: grType

            !Generate based on points/grType
            call inpXML%Set_Val(grType,'points/grType',"RP")
            write(*,*) 'Grid type = ', trim(toUpper(grType))

            call inpXML%Set_Val(Nx1,'points/Nx1',100)
            call inpXML%Set_Val(Nx2,'points/Nx2',100)
            call inpXML%Set_Val(Nx3,'points/Nx3',1)

            NumP = Nx1*Nx2*Nx3
            write(*,*) "Number of seed points = ", NumP
            
            if (allocated(Xs)) deallocate(Xs)
            allocate(Xs(NumP,NDIM))
           
            if(.not.(allocated(x1))) allocate(x1(Nx1))
            if(.not.(allocated(x2))) allocate(x2(Nx2))
            if(.not.(allocated(x3))) allocate(x3(Nx3))

            call inpXML%Set_Val(x1min,'points/x1min',-10.0)
            call inpXML%Set_Val(x1max,'points/x1max', 10.0)
            call inpXML%Set_Val(x2min,'points/x2min',-10.0)
            call inpXML%Set_Val(x2max,'points/x2max', 10.0)

            x1 = x1min + (x1max-x1min)/Nx1*[(i, i = 0, Nx1)]
            x2 = x2min + (x2max-x2min)/Nx2*[(i, i = 0, Nx2)]

            !Generate based on points/grType
            select case(trim(toUpper(grType)))
            !---------
            case("RP")
              write(*,*) 'Generating polar grid ...'
              call inpXML%Set_Val(rmin,'radius/min',5.0)
              call inpXML%Set_Val(rmax,'radius/max',25.0)
              call inpXML%Set_Val(pmin,'phi/min',0.0)
              call inpXML%Set_Val(pmax,'phi/max',360.0)

              x1 = rmin + (rmax-rmin)/Nx1*[(i, i = 0, Nx1)]
              x2 = pmin + (pmax-pmin)/Nx2*[(i, i = 0, Nx2)]
              x2 = x2/rad2deg
              n = 1
              do i=1,Nx1
                do j=1,Nx2
                  Xs(n,XDIR) = x1(i)*cos(x2(j))
                  Xs(n,YDIR) = x1(i)*sin(x2(j))
                  Xs(n,ZDIR) = 0.0
                  n = n+1
                enddo
                
              enddo
              !---------
              case("XY")
                write(*,*) 'Generating XY grid ...'
                n = 1
                do i=1,Nx1
                  do j=1,Nx2
                    Xs(n,XDIR) = x1(i)
                    Xs(n,YDIR) = x2(j)
                    Xs(n,ZDIR) = 0.0
                    n = n+1
                  enddo
                enddo
              !---------
              case("RTP","SPHERICAL")
                write(*,*) 'Generating spherical grid ...'
                call inpXML%Set_Val(rmin,'radius/min',5.0)
                call inpXML%Set_Val(rmax,'radius/max',25.0)
                call inpXML%Set_Val(pmin,'phi/min',0.0)
                call inpXML%Set_Val(pmax,'phi/max',360.0)
                call inpXML%Set_Val(tmin,'theta/min',0.0)
                call inpXML%Set_Val(tmax,'theta/max',90.0)
                x1 = rmin + (rmax-rmin)/Nx1*[(i, i = 0, Nx1)]
                x2 = pmin + (pmax-pmin)/Nx2*[(i, i = 0, Nx2)]
                x3 = tmin + (tmax-tmin)/Nx3*[(i, i = 0, Nx3)]
                x2 = x2/rad2deg
                x3 = x3/rad2deg

                n = 1
                do i=1,Nx1
                  do j=1,Nx2
                    do k=1,Nx3
                      Xs(n,XDIR) = x1(i)*cos(x2(j))*sin(x3(k))
                      Xs(n,YDIR) = x1(i)*sin(x2(j))*sin(x3(k))
                      Xs(n,ZDIR) = x1(i)*cos(x3(k))
                      n = n+1
                    enddo
                  enddo
                enddo

            end select
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

          Model%rmin = rmin
          Model%rmax = rmax

          ! find nightside i-location just inside the domain
          ! note minloc assumes indices start from 1, thus the -1+isg at the end
          ind=minloc(rmax+ebState%ebGr%xyz(:,je+1,ks,XDIR),dim=1,mask=rmax+ebState%ebGr%xyz(:,je+1,ks,XDIR)>=0.)-1+isg  
          ni = ind-is+1
          NumP = ni*nj*nk
          if (allocated(Xs)) deallocate(Xs)
          allocate(Xs(NumP,NDIM))

          do i=is,ni
             do j=js,je
                do k = ks,ke
                   Xs( (i-is)*nj*nk+(j-js)*nk+k, :) = ebState%ebGr%xyzcc(i,j,k,:)
                end do
             end do
          end do
        end subroutine genPtsFromGrid

end program tracex
