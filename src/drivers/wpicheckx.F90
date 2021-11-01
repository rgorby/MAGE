!Driver for running quick tests validating wpi based subroutines 
program wpicheck
    use kdefs
    use chmpdefs
    use chmpunits
    use ioH5
    use tptypes
    use tputils
    use starter
    use wpicalc
    use wpitypes
    use xml_input

    implicit none

    !Main data structures
    type(chmpModel_T) :: Model
    type(wave_T) :: wave
    type(wModel_T)   :: wModel_1prt
    type(XML_Input_T) :: inpXML

    integer, parameter :: MAXIOVAR = 50
    type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

    character(len=strLen) :: sStr

    !Initialize model
    call goApe(Model,iXML=inpXML)

    call inpXML%Set_Val(sStr,"tps/species","e-")
    call getSpecies(sStr,Model%m0,Model%q0)

    !use xml to set up daa test
    call initWPI(wModel_1prt,wave,inpXML)

    ! call singlePrtDiff(Model,wave,wModel_1prt,inpXML)

    call diffCruveBounce(Model,wave,wModel_1prt,inpXML)

    call DiffCoef_summers05(Model,wave,inpXML)

    call diffCurve_summers98(Model,wave)

    call resTest(Model,wave,wModel_1prt)

    contains

        ! all particles are GC
        subroutine createPrts(Model,prt,K0,alpha,psi,B0)
            type(chmpModel_T), intent(in) :: Model
            type(prt_T), intent(inout) :: prt
            real(rp), intent(in) :: K0,alpha,psi,B0

            real(rp) :: gamma,p11,pMag,Mu
            integer :: pSgn

            if (.not.prt%isGC) then
                write(*,*) 'Partciles must be guiding center!'
                stop
            endif

            !ignoring spatial position 
            prt%Q(XPOS:ZPOS) = 0
            !Note always setting initial EQ values regardless of initial z
            prt%Qeq(EQX:EQY) = 0

            !Set particle energy
            gamma = (K0*1.0e-3)/(Model%m0*mec2) + 1.0
            prt%Q(GAMGC) = gamma

            !For now set equatorial energies (both) to K0
            prt%Qeq(EQKEV) = K0
            prt%Qeq(EQKEB) = K0 

            prt%Q(PSITP) = psi

            !Set particle pitch angle (not a dynamic quantity)
            prt%alpha = alpha
            prt%Qeq(EQALP) = alpha

            !set direction for p11
            if (alpha <= PI/2) then
                pSgn = +1
            else
                pSgn = -1
            endif

            pMag = Model%m0*sqrt(gamma**2.0 - 1.0)
            p11 = pSgn*pMag*sqrt( 1 - sin(alpha)**2.0 )

            prt%Q(P11GC) = p11
            Mu = ( (Model%m0**2)*(gamma**2.0 - 1.0) - p11*p11) / (2*Model%m0*B0)
            prt%Q(MUGC) = Mu

            !Initialize various quantities
            prt%ijk0 = 0
            prt%Ngc = 0
            prt%Nfo = 0
            prt%ddt = 0.0
        end subroutine createPrts

        subroutine DiffCoef_summers05(Model,wave,inpXML)
            type(chmpModel_T), intent(inout)    :: Model
            type(wave_T),      intent(in)       :: wave
            type(XML_Input_T), intent(inout)    :: inpXML
            type(wModel_T)                      :: wModel
            type(prt_T)                         :: prt
            type(tpState_T),dimension(:),allocatable :: tpState
              
            real(rp) :: B0 ![nT] uniform background field
            integer  :: Nk,Nas,Np,n,i,j
            real(rp) :: alpha0,psi0,Daa
            real(rp) :: xj,yj
            real(rp), dimension(:), allocatable   :: K0,astar !particle energies
            real(rp), dimension(:,:), allocatable :: daa_f1,pa_f1,xjs_f1,daa_f2,pa_f2
            real(rp) :: dpa,palpha, Kprt

            character(len=strLen) :: outH5 = "diffCoef.h5"

            !setting values
            B0     = 344.0
            alpha0 = 0.0
            psi0   = 0.0

            !normalizing variables
            B0 = B0*inBScl
            alpha0 = alpha0*deg2rad
            psi0 = psi0*deg2rad

            !manually set wave spectrum to match values used in Summers 
            wModel%xm = 0.35
            wModel%Dx = 0.3
            wModel%B1 = 0.1*inBScl
            
            !Initialize the particles
            K0 = [100.0,300.0,1000.0,3000.0] ![keV] normalized when converted to gamma in createPrts
  
            Nk = size(K0)
            allocate(tpState(Nk))
            !Get number
            call inpXML%Set_Val(Np,"tps/Np",500)
            tpState%Np = Np
            do i=1,Nk
                allocate(tpState(i)%TPs(Np))
                tpState(i)%TPs(:)%isIn = .true.
                tpState(i)%TPs(:)%isGC = .true.
            end do

            !get random sample of pitch-angles between 0 and 90 deg.
            ! call getSample(Model,inpXML,"alpha",Np,alpha,0.0_rp,90.0_rp)
            ! alpha = alpha*deg2rad
            dpa = (90.0/Np)*deg2rad

            do i=1,Nk
                do n=1,Np
                    palpha = alpha0+n*dpa
                    !Create particles 
                    tpState(i)%TPs(n)%id = n !+(i-1)*Np
                    call createPrts(Model,tpState(i)%TPs(n),K0(i),palpha,psi0,B0)
                    !call createPrts(Model,tpState(i)%TPs(n),K0(i),alpha(n),psi0,B0)
                end do 
            end do

            astar = [0.4444,0.16,0.04,0.0178,0.01] !plasma parameter
            Nas = size(astar)

            !allocate arrays to hold output
            allocate(daa_f1(Nk,Np))
            allocate(pa_f1(Nk,Np))
            allocate(xjs_f1(Nk,Np))
            allocate(daa_f2(Nas,Np))
            allocate(pa_f2(Nas,Np))

            !Figure 1 of Summers 2005  Daa vs pa vs K
            do i=1,Nk
                do n=1,Np
                    prt = tpState(i)%TPs(n)
                    Kprt = pGC2K(Model,prt)
                    call Resonance(wave,wModel,Model%m0,Kprt,prt%alpha,astar(2),xj,yj)
                    if (xj == 999) then
                        xjs_f1(i,n) = xj
                        daa_f1(i,n) = xj
                        pa_f1(i,n)  = xj
                    else
                        xjs_f1(i,n) = xj
                        Daa = DiffCoef(wave,wModel,Model%m0,Kprt,prt%alpha,astar(2),B0,xj,yj) 
                        daa_f1(i,n) = Daa/oTScl
                        pa_f1(i,n)  = prt%alpha*rad2deg
                    end if
                end do 
            end do

            !Figure 2 of Summers 2005 Daa vs pa vs astar
            do j=1,Nas
                do n=1,Np
                    prt = tpState(3)%TPs(n) !take 1 MeV particles
                    Kprt = pGC2K(Model,prt)
                    call Resonance(wave,wModel,Model%m0,Kprt,prt%alpha,astar(j),xj,yj)
                    if (xj == 999) then
                        daa_f2(j,n) = xj
                        pa_f2(j,n)  = xj
                    else
                        Daa = DiffCoef(wave,wModel,Model%m0,Kprt,prt%alpha,astar(j),B0,xj,yj) 
                        daa_f2(j,n) = Daa/oTScl
                        pa_f2(j,n)  = prt%alpha*rad2deg
                    end if
                end do 
            end do

            !write data to a file
            call CheckAndKill(outH5)
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"astar",astar)
            call AddOutVar(IOVars,"K",K0)
            call AddOutVar(IOVars,"Daa1",daa_f1)
            call AddOutVar(IOVars,"alpha1",pa_f1)
            call AddOutVar(IOVars,"xjs",xjs_f1)
            call AddOutVar(IOVars,"Daa2",daa_f2)
            call AddOutVar(IOVars,"alpha2",pa_f2)
            call WriteVars(IOVars,.true.,outH5)
            call ClearIO(IOVars)

        end subroutine DiffCoef_summers05

        subroutine diffCurve_summers98(Model,wave)
            type(chmpModel_T), intent(in) :: Model
            type(wave_T),      intent(in) :: wave
            type(tpState_T)               :: tpState
            type(wModel_T)                :: wModel
            type(prt_T)                   :: prt

            real(rp) :: B0  ![nT] uniform background field (arbitrary)
            real(rp) :: dAlim![rad] Limiting change in pith-angle to be below this value each wpi 
            real(rp) :: alpha0,psi0,astar !Summers '98 used astar=[1,3,10] in Fig. 3
            real(rp), dimension(5) :: K0
            integer :: Na,Np,n,i,j
            real(rp), dimension(:,:), allocatable :: p11,pperp,Ks, xjs
            real(rp) :: xj,yj,da,dp,Daa,ddt,dGDda !Daa, ddt, dGDda are not used, using const da
            real(rp) :: gamNew,aNew,pNew,Mu,p11Mag,Kprt
            integer :: pSgn=1

            character(len=strLen) :: outH5 = "diffCurve.h5"

            !setting values
            B0 = 344
            dAlim = -0.05
            alpha0 = 180.0
            psi0=0
            astar=1.0 ![1,3,10]
            Daa=0.0001 ! arbitrary, not used
            dGDda = 0.5 ! arbitrary, not used
            ddt=0.05

            !manually set diff curve test to allow resonance with all waves
            wModel%xm = 0.5
            wModel%Dx = 0.5
            wModel%B1 = 0.1*inBScl

            !setting up diffusion curve test
            !starting energies for each particle
            K0 = [50.0,100.0,180.0,350.0,700.0]

            Na = abs((PI/2)/dAlim)
            Np = size(K0)

            !allocate arrays to hold output
            allocate(p11(Np,Na))
            allocate(pperp(Np,Na))
            allocate(Ks(Np,Na))
            allocate(xjs(Np,Na))

            !normalizing variables
            B0 = B0*inBScl
            alpha0 = alpha0*deg2rad

            !generating test particles for diffusion curve comparison
            tpState%Np = Np
            allocate(tpState%TPs(Np))

            tpState%TPs(:)%isIn = .true.
            tpState%TPs(:)%isGC = .true.

            ! Create particles for diffusion curve test
            do n=1,Np
                !Create particles 
                tpState%TPs(n)%id = n
                call createPrts(Model,tpState%TPs(n),K0(n),alpha0,psi0,B0)
            enddo

            ! loop over particles and pitch-angles 
            do i=1,Np 
                prt = tpState%TPs(i)
                do j=1,Na
                    Kprt = pGC2K(Model,prt)
                    call Resonance(wave,wModel,Model%m0,Kprt,prt%alpha,astar,xj,yj)
                    prt%xj = xj
                    prt%yj = yj

                    !save particle state before Updating
                    p11(i,j) = prt%Q(P11GC)
                    pperp(i,j) = sqrt(abs(prt%Q(MUGC))*2*B0*Model%m0)
                    xjs(i,j) = xj
                    Ks(i,j) = prt2kev(Model,prt)/1000.0 ![MeV]

                    ! Calculate the resulting change in pitch angle and energy of the particle
                    call LangevinEq(wave,wModel,Model,prt,dGDda,Daa,ddt,astar,aNew,pNew,dAlim)

                    !Update the pitch angle
                    prt%alpha = aNew

                    if (aNew <= PI/2) then
                        pSgn = +1
                    else
                        pSgn = -1
                    endif
                    !Updating the particle momentum and energy
                    p11Mag = pSgn*pNew*sqrt( 1 - sin(aNew)**2.0 )
                    gamNew = sqrt(1+(pNew/Model%m0)**2.0)
                    Mu = ( (Model%m0**2)*(gamNew**2.0 - 1.0) - p11Mag*p11Mag) / (2*Model%m0*B0)
                    prt%Q(P11GC) = p11Mag
                    prt%Q(MUGC ) = Mu 
                    prt%Q(GAMGC) = gamNew

                end do
            end do

            !write data to a file
            call CheckAndKill(outH5)
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"astar",astar)
            call AddOutVar(IOVars,"p11",p11)
            call AddOutVar(IOVars,"pperp",pperp)
            call AddOutVar(IOVars,"xjs",xjs)
            call AddOutVar(IOVars,"K",Ks)
            call WriteVars(IOVars,.true.,outH5)
            call ClearIO(IOVars)

            !performing diffusion curve test
        end subroutine diffCurve_summers98

        subroutine diffCruveBounce(Model,wave,wModel,inpXML)
            type(chmpModel_T), intent(inout)    :: Model
            type(wave_T),      intent(in)       :: wave
            type(wModel_T),    intent(in)       :: wModel
            type(XML_Input_T), intent(inout)    :: inpXML
            type(tpState_T)                     :: tpState
            type(prt_T)   :: prtDC
            real(rp), dimension(:), allocatable :: time,B,astar,lat
            real(rp), dimension(:,:), allocatable :: kDC,aDC,xjDC,yjDC,p11DC,gamDC,muDC
            real(rp) :: B0,n0
            real(rp) :: K0,alpha0,psi0
            real(rp) :: currlat,Bm,Bcurr,Acurr
            integer :: Na,Np,n,i,j
            real(rp) :: wpe,inWScl,xj,yj,Daa=0.0001,dGDda=0.5 !Daa,dGDda are not used, using const da
            real(rp) :: dAlim = -0.00087
            real(rp) :: gamNew,aNew,pNew,Mu,p11Mag,Kprt
            integer  :: pSgn=1

            character(len=strLen) :: outH5 = "diffCurveBounce.h5"

            logical :: doAllWaves

            doAllWaves = .true.

            !pull K, alpha, psi for prt from xml file 
            call inpXML%Set_Val(K0,'energy/min',500.0)
            call inpXML%Set_Val(alpha0,'alpha/min',45.0)
            call inpXML%Set_Val(psi0,'psi/min',0.0)

            B0     = 156.0_rp
            n0     = 9.0_rp

            !normalizing variables
            B0 = B0*inBScl
            alpha0 = alpha0*deg2rad
            psi0 = psi0*deg2rad                   
            
            !Bmirror
            Bm = B0/(sin(alpha0)**2.0)

            ! Na = abs((PI/2)/dAlim)
            Na = int(abs(alpha0)/abs(dAlim))
            ! Na = int(abs(alpha0-PI/2)/dAlim)
            Np = 11 !1 degree resolution in lat from -10 to 0

            !generating test particles for diffusion curve comparison
            tpState%Np = Np
            allocate(tpState%TPs(Np))

            tpState%TPs(:)%isIn = .true.
            tpState%TPs(:)%isGC = .true.

            !allocate arrays to hold output
            allocate(B(Np))
            allocate(time(Np))
            allocate(lat(Np))
            allocate(astar(Np))
            allocate(kDC(Np,Na))
            allocate(aDC(Np,Na))
            allocate(xjDC(Np,Na))
            allocate(yjDC(Np,Na))
            allocate(p11DC(Np,Na))
            allocate(muDC(Np,Na))
            allocate(gamDC(Np,Na))

            ! Create particles for diffusion curve test
            currlat = -10.0_rp*deg2rad
            do n=1,Np
                ! loop through bounce motion to create particles and astar at each latitude 
                Bcurr = B0*sqrt(1.0+3.0*(sin(currlat)**2.0))/(cos(currlat)**6.0)
                Acurr = asin(sqrt(Bcurr/Bm))
                tpState%TPs(n)%id = n
                call createPrts(Model,tpState%TPs(n),K0,Acurr,psi0,Bcurr)
                Bcurr = B0
                write(*,*)"Bcurr,Bm,Acurr,currlat: ",Bcurr*oBScl,Bm*oBScl,Acurr*rad2deg,currlat*rad2deg
                !Calulating the ratio of nonrelativistic gyrofrequency to plasma frequency at prt's location
                inWScl = 114704.0207604 !(qe_cgs*L0/(vc_cgs*sqrt(Me_cgs)))^2 
                wpe = sqrt(4*PI*n0) 
                astar(n) = Bcurr**2/((wpe**2)*inWScl)
                B(n) = Bcurr
                time(n) = n
                lat(n) = currlat*rad2deg
                currlat = currlat+1.0_rp*deg2rad
            enddo

            ! loop over particles and pitch-angles 
            do i=1,Np 
                prtDC = tpState%TPs(i)
                do j=1, Na
                    write(*,*) "iteration: ",i,j
                    Kprt = pGC2K(Model,prtDC)
                    call Resonance(wave,wModel,Model%m0,Kprt,prtDC%alpha,astar(i),xj,yj,doAllWaves)
                    prtDC%xj = xj
                    prtDC%yj = yj
                    if(xj == 999) cycle

                    !save particle state before Updating
                    kDC(i,j) = prt2kev(Model,prtDC)
                    aDC(i,j) = prtDC%alpha*rad2deg
                    xjDC(i,j)= prtDC%xj
                    yjDC(i,j)= prtDC%yj
                    p11DC(i,j) = prtDC%Q(P11GC)
                    muDC(i,j)  = prtDC%Q(MUGC )
                    gamDC(i,j) = prtDC%Q(GAMGC)
                    write(*,*) "B,lat, astar: ",B(i)*oBScl,lat(i),astar(i)
                    write(*,*) "Energy, pa,xj,yj: ",kDC(i,j),aDC(i,j),xjDC(i,j),yjDC(i,j)
                    ! Calculate the resulting change in pitch angle and energy of the particle
                    call LangevinEq(wave,wModel,Model,prtDC,dGDda,Daa,Model%dt,astar(i),aNew,pNew,dAlim)

                    if (aNew >= PI/2) exit
                    !Update the pitch angle
                    prtDC%alpha = aNew

                    if (aNew <= PI/2) then
                        pSgn = +1
                    else
                        pSgn = -1
                    endif
                    !Updating the particle momentum and energy
                    p11Mag = pSgn*pNew*sqrt( 1 - sin(aNew)**2.0 )
                    gamNew = sqrt(1+(pNew/Model%m0)**2.0)
                    Mu = ( (Model%m0**2)*(gamNew**2.0 - 1.0) - p11Mag*p11Mag) / (2*Model%m0*B(i))
                    prtDC%Q(P11GC) = p11Mag
                    prtDC%Q(MUGC ) = Mu 
                    prtDC%Q(GAMGC) = gamNew
                enddo
            end do

            !write data to a file
            call CheckAndKill(outH5)
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"time",time)
            call AddOutVar(IOVars,"lat",lat)
            call AddOutVar(IOVars,"B",B*oBScl)
            call AddOutVar(IOVars,"kDC",kDC)
            call AddOutVar(IOVars,"aDC",aDC)
            call AddOutVar(IOVars,"xjDC",xjDC)
            call AddOutVar(IOVars,"yjDC",yjDC)
            call AddOutVar(IOVars,"p11DC",p11DC)
            call AddOutVar(IOVars,"muDC",muDC)
            call AddOutVar(IOVars,"gamDC",gamDC)
            call WriteVars(IOVars,.true.,outH5)
            call ClearIO(IOVars)

        end subroutine diffCruveBounce

        subroutine singlePrtDiff(Model,wave,wModel,inpXML)
            type(chmpModel_T), intent(inout)    :: Model
            type(wave_T),      intent(in)       :: wave
            type(wModel_T),    intent(in)       :: wModel
            type(XML_Input_T), intent(inout)    :: inpXML
            type(ebState_T) :: ebState
            type(tpState_T),dimension(:),allocatable :: tpState
            type(prt_T)   :: prt,prtDC

            character(len=strLen) :: outH5 = "singlePrtTest.h5"

            real(rp), dimension(:), allocatable :: kwpi,awpi,xjwpi,yjwpi,p11wpi,gamwpi,muwpi,Nwpi,dKwpi,dAwpi
            real(rp), dimension(:), allocatable :: Keq,Aeq,time,lat,B,z
            real(rp), dimension(:,:), allocatable :: kDC,aDC,xjDC,yjDC,p11DC,gamDC,muDC
            real(rp), dimension(:), allocatable :: astar
            real(rp) :: B0,n0,B1
            real(rp) :: K0,alpha0,psi0
            real(rp) :: t0,t1,dt
            real(rp) :: x,currlat,L,r,Bm,Bcurr,Acurr,gamma,pMag
            real(rp), dimension(NDIM) :: req
            integer  :: Nt,s,Na,n,i
            real(rp) :: wpe,inWScl,xj,yj,Daa=0.0001,dGDda=0.5 !Daa,dGDda are not used, using const da
            real(rp) :: dAlim = 0.00087
            real(rp) :: gamNew,aNew,pNew,Mu,p11Mag,Kprt
            integer  :: pSgn=1

            logical :: doAllWaves

            doAllWaves = .true.

            !pull K, alpha, psi for prt from xml file 
            call inpXML%Set_Val(K0,'energy/min',500.0)
            call inpXML%Set_Val(alpha0,'alpha/min',45.0)
            call inpXML%Set_Val(psi0,'psi/min',0.0)

            ! pulling timing info, already normalized
            t0 = Model%T0
            t1 = Model%tFin
            dt = Model%dt
            Nt = int((t1-t0)/dt)

            !setting ebState wave information (nothing else is needed in ebState)
            ebState%ebWave = wave
            ebState%ebWmodel = wModel

            !setting values at equator
            ! B0     = 341.3
            ! n0     = 7.0
            B0     = 156.0_rp
            n0     = 9.0_rp
            ! B0     = 312.0

            !normalizing variables
            B0 = B0*inBScl
            alpha0 = alpha0*deg2rad
            psi0 = psi0*deg2rad 

            prt%isIn = .true.
            prt%isGC = .true.
            call createPrts(Model,prt,K0,alpha0,psi0,B0)                    

            ! setting position of particle, putting at midnight so y=0
            prt%Q(XPOS:ZPOS) = [-6.0_rp, 0.0_rp, 0.0_rp] 
            prt%Qeq(EQX:EQY) = [-6.0_rp, 0.0_rp] 
            req = [prt%Qeq(EQX),prt%Qeq(EQY),0.0_rp]
            L = norm2(req)
            
            !Bmirror
            Bm = B0/(sin(alpha0)**2.0)

            ! particle initially at equator
            prt%Qeq(EQKEV) = prt2kev(Model,prt)
            prt%Qeq(EQALP) = prt%alpha

            prt%ddt = Model%dt

            !allocate arrays to hold output
            allocate(kwpi(Nt))
            allocate(awpi(Nt))
            allocate(xjwpi(Nt))
            allocate(yjwpi(Nt))
            allocate(p11wpi(Nt))
            allocate(gamwpi(Nt))
            allocate(muwpi(Nt))
            allocate(Nwpi(Nt))
            allocate(time(Nt))
            allocate(dKwpi(Nt))
            allocate(dAwpi(Nt))

            allocate(Keq(Nt))
            allocate(Aeq(Nt))
            allocate(lat(Nt))
            allocate(z(Nt))
            allocate(B(Nt))

            do s=1,Nt
                !moving particle along field line
                currlat =  sin(s/(20.0*PI))*(PI/7.79) !bounces between +/- 23˚ 
                Bcurr = B0*sqrt(1.0+3.0*(sin(currlat)**2.0))/(cos(currlat)**6.0)
                r = L*(cos(currlat)**2.0)
                z(s) = r*sin(currlat)
                x = r*cos(currlat)
                prt%Q(XPOS:ZPOS) = [-x, 0.0_rp, z(s)]
                ! Acurr = asin(sqrt(Bcurr/Bm))

                ! prt%alpha = Acurr

                ! if (Acurr <= PI/2) then
                !     pSgn = +1
                ! else
                !     pSgn = -1
                ! endif
                ! !Energy and mu unchange, assume no vExB 
                ! gamma = prt%Q(GAMGC)
                ! pMag = Model%m0*sqrt(gamma**2.0 - 1.0)
                ! p11Mag = pSgn*pMag*sqrt( 1 - sin(Acurr)**2.0 )
                ! prt%Q(P11GC) = p11Mag
                ! Mu = ( (Model%m0**2)*(ebGam**2.0 - 1.0) - p11Mag*p11Mag) / (2*Model%m0*MagB)
                ! prt%Q(MUGC ) = Mu 

                ! WPI calculations
                call PerformWPI(prt,Model%t,Model,ebState,Bcurr,n0)

                !check for equatoria crossing
                if ((s > 1) .and. (z(s-1)*z(s) <=0)) then
                    prt%Qeq(EQKEV) = prt2kev(Model,prt)
                    prt%Qeq(EQALP) = prt%alpha
                end if

                Model%t = Model%t+Model%dt

                kwpi(s) = prt2kev(Model,prt)
                awpi(s) = prt%alpha*rad2deg
                xjwpi(s)= prt%xj
                yjwpi(s)= prt%yj
                p11wpi(s) = prt%Q(P11GC)
                muwpi(s)  = prt%Q(MUGC )
                gamwpi(s) = prt%Q(GAMGC)
                Nwpi(s) = prt%Nwpi
                time(s) = oTScl*Model%t
                lat(s)  = currlat*rad2deg
                B(s)    = Bcurr*oBScl
                dKwpi(s)= prt%dKwpi
                dAwpi(s)= prt%dAwpi
                Keq(s) = prt%Qeq(EQKEV)
                Aeq(s) = prt%Qeq(EQALP)*rad2deg
            enddo

            ! diffusion curve calculation
            prtDC%isIn = .true.
            prtDC%isGC = .true.
            ! call createPrts(Model,prtDC,K0,alpha0,psi0,B0)

            Na = int(abs(alpha0-PI/2)/dAlim)

            !allocate arrays to hold output
            allocate(kDC(2,Na))
            allocate(aDC(2,Na))
            allocate(xjDC(2,Na))
            allocate(yjDC(2,Na))
            allocate(p11DC(2,Na))
            allocate(muDC(2,Na))
            allocate(gamDC(2,Na))
            allocate(astar(2))

            !Calulating the ratio of nonrelativistic gyrofrequency to plasma frequency at prt's location
            inWScl = 114704.0207604 !(qe_cgs*L0/(vc_cgs*sqrt(Me_cgs)))^2 
            wpe = sqrt(4*PI*n0) 
            astar(1) = B0**2/((wpe**2)*inWScl)
            B1 = B0*sqrt(1.0+3.0*(sin(PI/18)**2.0))/(cos(PI/18)**6.0)
            astar(2) = B1**2/((wpe**2)*inWScl)

            write(*,*) "SinglePrtTest:: B0, n0, astar: ",B0*oBScl,n0,astar

            do i=1,2
                call createPrts(Model,prtDC,K0,alpha0,psi0,B0)
                do n=1, Na
                    Kprt = pGC2K(Model,prtDC)
                    call Resonance(wave,wModel,Model%m0,Kprt,prtDC%alpha,astar(i),xj,yj,doAllWaves)
                    prtDC%xj = xj
                    prtDC%yj = yj

                    !save particle state before Updating
                    kDC(i,n) = prt2kev(Model,prtDC)
                    aDC(i,n) = prtDC%alpha*rad2deg
                    xjDC(i,n)= prtDC%xj
                    yjDC(i,n)= prtDC%yj
                    p11DC(i,n) = prtDC%Q(P11GC)
                    muDC(i,n)  = prtDC%Q(MUGC )
                    gamDC(i,n) = prtDC%Q(GAMGC)

                    ! Calculate the resulting change in pitch angle and energy of the particle
                    call LangevinEq(wave,wModel,Model,prtDC,dGDda,Daa,Model%dt,astar(i),aNew,pNew,dAlim)

                    !Update the pitch angle
                    prtDC%alpha = aNew

                    if (aNew <= PI/2) then
                        pSgn = +1
                    else
                        pSgn = -1
                    endif
                    !Updating the particle momentum and energy
                    p11Mag = pSgn*pNew*sqrt( 1 - sin(aNew)**2.0 )
                    gamNew = sqrt(1+(pNew/Model%m0)**2.0)
                    Mu = ( (Model%m0**2)*(gamNew**2.0 - 1.0) - p11Mag*p11Mag) / (2*Model%m0*B0)
                    prtDC%Q(P11GC) = p11Mag
                    prtDC%Q(MUGC ) = Mu 
                    prtDC%Q(GAMGC) = gamNew
                enddo
            enddo

            !write data to a file
            call CheckAndKill(outH5)
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"time",time)
            call AddOutVar(IOVars,"lat",lat)
            call AddOutVar(IOVars,"B",B)
            call AddOutVar(IOVars,"Z",z)
            call AddOutVar(IOVars,"Keq",Keq)
            call AddOutVar(IOVars,"Aeq",Aeq)
            call AddOutVar(IOVars,"kwpi",kwpi)
            call AddOutVar(IOVars,"awpi",awpi)
            call AddOutVar(IOVars,"Nwpi",Nwpi)
            call AddOutVar(IOVars,"dKwpi",dKwpi)
            call AddOutVar(IOVars,"dAwpi",dAwpi)
            call AddOutVar(IOVars,"xjwpi",xjwpi)
            call AddOutVar(IOVars,"yjwpi",yjwpi)
            call AddOutVar(IOVars,"p11wpi",p11wpi)
            call AddOutVar(IOVars,"muwpi",muwpi)
            call AddOutVar(IOVars,"gamwpi",gamwpi)
            call AddOutVar(IOVars,"kDC",kDC)
            call AddOutVar(IOVars,"aDC",aDC)
            call AddOutVar(IOVars,"xjDC",xjDC)
            call AddOutVar(IOVars,"yjDC",yjDC)
            call AddOutVar(IOVars,"p11DC",p11DC)
            call AddOutVar(IOVars,"muDC",muDC)
            call AddOutVar(IOVars,"gamDC",gamDC)
            call WriteVars(IOVars,.true.,outH5)
            call ClearIO(IOVars)

        end subroutine singlePrtDiff

        subroutine resTest(Model,wave,wModel)
            type(chmpModel_T), intent(in) :: Model
            type(wave_T),      intent(in) :: wave
            type(wModel_T),    intent(in) :: wModel
            type(prt_T)                   :: prt

            real(rp) :: B0 ![nT] uniform background field (arbitrary)
            real(rp) :: K0, alpha0,psi0,astar
            real(rp) :: xj,yj
            real(rp) :: p11,pperp,alpha,K,Daa,Kprt
            real(rp), dimension(:), allocatable :: xjs, yjs,x90,y90
            integer :: n,Na
            real(rp) :: dAlim
            real(rp) :: pa
            real(rp), dimension(:), allocatable :: xroots,yroots,resid,alphas,astars,Ks
            complex(rp), dimension(:,:), allocatable :: allRoots
            complex(rp), dimension(NROOTS+1) :: coef
            logical :: doAllWaves

            character(len=strLen) :: outH5 = "rootTest.h5"

            doAllWaves = .true.
            
            !setting values
            ! B0     = 344.0
            K0     = 100.0_rp
            alpha0 = 89.0_rp
            psi0   = 0.0_rp
            ! astar  = 0.16_rp
            B0     = 156.0_rp
            ! K0 = 136.21044791583898
            ! alpha0 = 89.999648419152763
            ! alpha0 = 89.999
            astar = 0.026!1!0.16!0.026

            !normalizing variables
            B0 = B0*inBScl
            alpha0 = alpha0*deg2rad
            psi0 = psi0*deg2rad

            prt%isIn = .true.
            prt%isGC = .true.

            call createPrts(Model,prt,K0,alpha0,psi0,B0)

            Kprt = pGC2K(Model,prt)

            ! do n=1, 100000
            call Resonance(wave,wModel,Model%m0,Kprt,prt%alpha,astar,xj,yj)

            ! call res90deg(wave,Model%m0,Kprt,astar,x90,y90)
            ! write(*,*) "ResTest::alpha=90:: xjs,yjs: ",x90,y90
            ! call ResRoots(wave,Model%m0,Kprt,prt%alpha,astar,xjs,yjs)
            ! write(*,*) "ResTest::All roots:: xjs,yjs: ",xjs,yjs

            Daa = DiffCoef(wave,wModel,Model%m0,Kprt,prt%alpha,astar,B0,xj,yj)

            p11 = prt%Q(P11GC)
            alpha = prt%alpha/deg2rad
            pperp = sqrt(abs(prt%Q(MUGC))*2*B0*Model%m0)
            K = prt2kev(Model,prt)!/1000.0 ![MeV]

            write(*,*) "ResTest:: K,p11,pperp,alpha: ",K,p11,pperp,alpha
            write(*,*) "ResTest:: xj,yj: ",xj,yj
            write(*,*) "ResTest:: astar, Daa: ",astar,Daa, Daa/oTScl
            ! enddo

            ! testing quartic root solver
            dAlim = 0.0000087_rp
            Na = int(abs(alpha0-PI/2)/dAlim)
            allocate(allRoots(Na,NROOTS))
            allocate(xroots(Na))
            allocate(yroots(Na))
            allocate(resid(Na))
            allocate(alphas(Na))

            do n=1, Na
                pa = alpha0+(n-1)*dAlim
                call ResRoots(wave,Model%m0,Kprt,pa,astar,xjs,yjs,coef)
                call fastQuarticSolver(coef,allRoots(n,:))
                call selectRoot(wModel,pa,xjs,yjs,xj,yj,doAllWaves)
                alphas(n) = pa
                xroots(n) = xj
                yroots(n) = yj
                resid(n) = quarticResidual(xj,coef)
                ! write(*,*) "Resonant Energy and pitch angle: ",K,pa*rad2deg
                ! write(*,*) allRoots(n,:)
                ! write(*,*) xj,yj
                ! write(*,*) " "
            enddo

            allocate(Ks(2))
            allocate(astars(2))
            Ks(:) = K
            astars(:) = astar

           !write data to a file
            call CheckAndKill(outH5)
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"xroots",xroots)
            call AddOutVar(IOVars,"yroots",yroots)
            call AddOutVar(IOVars,"resid",resid)
            call AddOutVar(IOVars,"alphas",alphas)
            call AddOutVar(IOVars,"Ks",Ks)
            call AddOutVar(IOVars,"astars",astars)
            ! call AddOutVar(IOVars,"allRoots",allRoots)
            call WriteVars(IOVars,.true.,outH5)
            call ClearIO(IOVars)

        end subroutine resTest

        function quarticResidual(x,coef) result(resid)
            real(rp), intent(in) :: x
            complex(rp), dimension(NROOTS+1), intent(in) :: coef

            real(rp) :: resid

            !calculating value of quartic at a given x for array of coef
            resid = coef(5)+x*(coef(4)+x*(coef(3)+x*(coef(2)+x*coef(1))))

        end function quarticResidual


end program wpicheck