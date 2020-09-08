!Driver for running quick tests validating wpi based subroutines 
program wpicheck
    use kdefs
    use chmpdefs
    use chmpunits
    use ioH5
    use tptypes
    use starter
    use wpicalc
    use wpitypes
    use xml_input

    implicit none

    !Main data structures
    type(chmpModel_T) :: Model
    type(tpState_T)   :: tpState_daa
    type(wave_T) :: wave
    type(wModel_T)   :: wModel_daa
    type(XML_Input_T) :: inpXML

    integer, parameter :: MAXIOVAR = 50
    type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

    character(len=strLen) :: sStr

    !Initialize model
    call goApe(Model,iXML=inpXML)

    call inpXML%Set_Val(sStr,"tps/species","e-")
    call getSpecies(sStr,Model%m0,Model%q0)

    !use xml to set up daa test
    call initWPI(Model,wModel_daa,wave,inpXML)

    call DiffCoef_summers05(Model,wave,wModel_daa,inpXML)

    call diffCurve_summers98(Model,wave)

    call resTest(Model,wave,wModel_daa)

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

        subroutine DiffCoef_summers05(Model,wave,wModel,inpXML)
            type(chmpModel_T), intent(inout)    :: Model
            type(wave_T),      intent(in)       :: wave
            type(wModel_T),    intent(in)       :: wModel
            type(XML_Input_T), intent(inout)    :: inpXML
            type(prt_T)                         :: prt
            type(tpState_T),dimension(:),allocatable :: tpState
            ! type(tpState_T),allocatable :: tpState1,tpState2,tpState3,tpState4,tpState5
            
            real(rp) :: B0 ![nT] uniform background field
            integer  :: Nk,Nas,Np,n,i,j
            real(rp) :: alpha0,psi0,Daa
            real(rp) :: xj,yj
            real(rp), dimension(:), allocatable   :: K0,astar !particle energies
            real(rp), dimension(:), allocatable   :: alpha
            real(rp), dimension(:,:), allocatable :: daa_f1,pa_f1,xjs_f1,daa_f2,pa_f2

            character(len=strLen) :: outH5 = "diffCoef.h5"

            !setting values
            B0     = 344.0
            alpha0 = 180.0
            psi0   = 0.0
            
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
            call getSample(Model,inpXML,"alpha",Np,alpha,0.0_rp,90.0_rp)
            alpha = alpha/rad2deg

            do i=1,Nk
                do n=1,Np
                    !Create particles 
                    tpState(i)%TPs(n)%id = n !+(i-1)*Np
                    call createPrts(Model,tpState(i)%TPs(n),K0(i),alpha(n),psi0,B0)
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
                   call Resonance(Model,wave,wModel,prt,astar(2),xj,yj)
                   xjs_f1(i,n) = xj
                   Daa = DiffCoef(Model,wave,wModel,prt,astar(2),B0,xj,yj) 
                   daa_f1(i,n) = Daa/oTScl
                   pa_f1(i,n)  = prt%alpha*rad2deg
                end do 
            end do

            !Figure 2 of Summers 2005 Daa vs pa vs astar
            do j=1,Nas
                do n=1,Np
                   prt = tpState(3)%TPs(n) !take 1 MeV particles
                   call Resonance(Model,wave,wModel,prt,astar(j),xj,yj)
                   Daa = DiffCoef(Model,wave,wModel,prt,astar(j),B0,xj,yj) 
                   daa_f2(j,n) = Daa/oTScl
                   pa_f2(j,n)  = prt%alpha*rad2deg
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
            real(rp) :: xj,yj,da,dp,Daa,ddt !Daa and ddt not used, using const da
            real(rp) :: gamOld,gamNew,aNew,pMag,Mu,p11Mag
            integer :: pSgn=1

            character(len=strLen) :: outH5 = "diffCurve.h5"

            !setting values
            B0 = 344
            dAlim = -0.05
            alpha0 = 180.0
            psi0=0
            astar=1.0 ![1,3,10]
            Daa=0.0001 ! arbitrary, not used
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
                    call Resonance(Model,wave,wModel,prt,astar,xj,yj)
                    !save particle state before Updating
                    p11(i,j) = prt%Q(P11GC)
                    pperp(i,j) = sqrt(abs(prt%Q(MUGC))*2*B0*Model%m0)
                    xjs(i,j) = xj
                    Ks(i,j) = prt2kev(Model,prt)/1000.0 ![MeV]

                    ! Calculate the resulting change in pitch angle and energy of the particle
                    call DiffCurve(Model,wave,prt,Daa,ddt,xj,yj,da,dp,dAlim)

                    !Update the pitch angle
                    aNew = prt%alpha + da
                    prt%alpha = aNew

                    if (aNew <= PI/2) then
                        pSgn = +1
                    else
                        pSgn = -1
                    endif
                    !Updating the particle momentum and energy
                    gamOld = prt2Gam(prt,Model%m0)
                    pMag = Model%m0*sqrt(gamOld**2.0 - 1.0)+dp !dp is scalar and change in total momentum
                    p11Mag = pSgn*pMag*sqrt( 1 - sin(aNew)**2.0 )
                    gamNew = sqrt(1+(pMag/Model%m0)**2.0)
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

        subroutine resTest(Model,wave,wModel)
            type(chmpModel_T), intent(in) :: Model
            type(wave_T),      intent(in) :: wave
            type(wModel_T),    intent(in) :: wModel
            type(prt_T)                   :: prt

            real(rp) :: B0 ![nT] uniform background field (arbitrary)
            real(rp) :: K0, alpha0,psi0,astar
            real(rp) :: xj,yj !Daa and ddt not used, using const da
            real(rp) :: p11,pperp,alpha,K,Daa

            !setting values
            B0     = 344.0
            K0     = 100.0
            alpha0 = 180.0
            psi0   = 0.0
            astar  = 1.0

            !normalizing variables
            B0 = B0*inBScl
            alpha0 = alpha0*deg2rad
            psi0 = psi0*deg2rad

            prt%isIn = .true.
            prt%isGC = .true.

            call createPrts(Model,prt,K0,alpha0,psi0,B0)

            call Resonance(Model,wave,wModel,prt,astar,xj,yj)

            Daa = DiffCoef(Model,wave,wModel,prt,astar,B0,xj,yj)

            p11 = prt%Q(P11GC)
            alpha = prt%alpha/deg2rad
            pperp = sqrt(abs(prt%Q(MUGC))*2*B0*Model%m0)
            K = prt2kev(Model,prt)!/1000.0 ![MeV]

            write(*,*) "ResTest:: K,p11,pperp,alpha: ",K,p11,pperp,alpha
            write(*,*) "ResTest:: xj,yj: ",xj,yj
            write(*,*) "ResTest:: astar, Daa: ",astar,Daa


        end subroutine resTest


end program wpicheck