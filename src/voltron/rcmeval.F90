!Routines for interpolating (TSC) on RCM grid
module rcmeval
    use volttypes
    use rcm_mhd_interfaces
    use rcmdefs, only : DenPP0,PSPHKT
    use earthhelper
    use gdefs, only : dFloor,pFloor

    implicit none
    logical  :: doWolfLim   = .true. !Whether to do wolf-limiting
    logical  :: doWolfNLim  = .false.  !If wolf-limiting whether to do wolf-limiting on density as well
    real(rp) :: nBounce     = 2.0 !Scaling factor for Alfven transit
    real(rp) :: maxBetaLim  = 6.0/5.0 !Largest beta to use in weighting formula
	logical  :: doBounceDT = .true. !Whether to use Alfven bounce in dt-ingest

    !Points for interpolation
    integer, parameter, private :: Np = 9

    contains

    !Enforce Wolf-limiting on an MHD/RCM thermodynamic state
    !Density [#/cc], pressure [nPa]
    subroutine WolfLimit(nrc,prc,npsph,nmhd,pmhd,beta,nlim,plim)
        real(rp), intent(in)  :: nrc,prc,npsph,nmhd,pmhd,beta
        real(rp), intent(out) :: nlim,plim

        real(rp) :: nrcm,prcm,ppsph
        real(rp) :: alpha,blim,dVoV,wRCM,wMHD
        logical :: doRC,doPP

        nlim = 0.0
        plim = 0.0
        nrcm = 0.0
        prcm = 0.0

        !Get a low but non-zero pressure for plasmasphere
        ppsph = DkT2P(npsph,PSPHKT) !Using ~eV default plasmasphere temperature
    !Incorporate RC/PP contributions
        !Test RC/PP contribution
        doRC = (prc   >= TINY  )
        doPP = (npsph >= DenPP0)

        if (doRC) then
            !Incorporate RC contribution
            nrcm = nrcm + nrc
            prcm = prcm + prc
        endif
        if (doPP) then
            !Incorporate plasmasphere contribution
            nrcm = nrcm + npsph
            prcm = prcm + ppsph
        endif
        !Now have total density/pressure contributions from RC+PP
    !Think about bailing
        !Return raw values if not limiting, and don't limit if there's no RC
        if ( (.not. doWolfLim) .or. (.not. doRC) ) then
            nlim = nrcm
            plim = prcm
            return
        endif

    !If still here we've gotta wolf limit
        !Experiment w/ limiting max value of beta allowed
        blim = min(beta,maxBetaLim)
        !Get scaling term
        alpha = 1.0 + blim*5.0/6.0

        wRCM = 1.0/alpha
        wMHD = (alpha-1.0)/alpha ! = 1 - wRCM
        plim = wRCM*prcm + wMHD*pmhd

        !Check whether to limit density
        if (doWolfNLim .and. (nrcm>TINY)) then
            !n_R V = (n_M + dn)(V + dV)
            !nlim = n_M + dn, Drop dn*dV =>
            dVoV = 0.5*(blim/alpha)*(prcm-pmhd)/pmhd
            nlim = nrcm - nmhd*dVoV
            if (nlim <= dFloor) then
                !Something went bad, nuke everything
                nlim = 0.0
                plim = 0.0
            endif !nlim
        else
            nlim = nrcm !Raw density
        endif !doWolfNLim
        
    end subroutine WolfLimit

    !Interpolate state at lat/lon
    subroutine InterpRCM(RCMApp,lat,lon,t,imW,isEdible)
    	type(rcm_mhd_T), intent(in)  :: RCMApp
        real(rp)       , intent(in)  :: lat,lon,t
        real(rp)       , intent(out) :: imW(NVARIMAG)
        logical        , intent(out) :: isEdible

        real(rp) :: colat,nrcm,prcm,npp,pScl,beta,pmhd,nmhd,wIM
        real(rp) :: plim,nlim,Tb
        integer, dimension(2) :: ij0
        integer :: Ni,Nj
        integer , dimension(Np,2) :: IJs
        real(rp), dimension(Np) :: Ws
        logical , dimension(Np) :: isGs

        Ni = RCMApp%nLat_ion
        Nj = RCMApp%nLon_ion

        !Set defaults
        imW(:) = 0.0
        imW(IMDEN ) = 0.0
        imW(IMPR  ) = 0.0
        imW(IMTSCL) = 0.0
        isEdible = .false.

        colat = PI/2 - lat

        !Do 1st short cut tests
        isEdible =  (colat >= RCMApp%gcolat(1)) .and. (colat <= RCMApp%gcolat(RCMApp%nLat_ion)) &
                    .and. (lat > TINY)

        if (.not. isEdible) return

        !If still here, find mapping (i,j) on RCM grid of point
        call GetRCMLoc(lat,lon,ij0)

        !Do second short cut tests
        isEdible = RCMApp%toMHD(ij0(1),ij0(2))
        if (.not. isEdible) return

        call GetInterp(lat,lon,ij0,IJs,Ws,isGs)
        !Do last short cut
        if (.not. all(isGs)) return

        prcm = rcmPScl*AvgQ(RCMApp%Prcm ,IJs,Ws,Ni,Nj)
        nrcm = rcmNScl*AvgQ(RCMApp%Nrcm ,IJs,Ws,Ni,Nj)
        npp  = rcmNScl*AvgQ(RCMApp%Npsph,IJs,Ws,Ni,Nj)
        pmhd = rcmPScl*AvgQ(RCMApp%Pave ,IJs,Ws,Ni,Nj)
        nmhd = rcmNScl*AvgQ(RCMApp%Nave ,IJs,Ws,Ni,Nj)
        beta = AvgQ(RCMApp%beta_average ,IJs,Ws,Ni,Nj)
        wIM  = AvgQ(RCMApp%wImag        ,IJs,Ws,Ni,Nj)
        Tb   = AvgQ(RCMApp%Tb           ,IJs,Ws,Ni,Nj)

        nlim = 0.0
        plim = 0.0

        if (doWolfLim) then
            call WolfLimit(nrcm,prcm,npp,nmhd,pmhd,beta,nlim,plim)
        else
            !Just lazyily use same function w/ beta=0
            call WolfLimit(nrcm,prcm,npp,nmhd,pmhd,0.0_rp,nlim,plim)
        endif

        !Store values
        imW(IMDEN) = nlim
        imW(IMPR)  = plim

        if (doBounceDT) then
            !Use Alfven bounce timescale
            imW(IMTSCL) = nBounce*RCMApp%Tb(ij0(1),ij0(2))
        endif

        imW(IMX1)   = rad2deg*lat
        imW(IMX2)   = rad2deg*lon

        contains

        !Get ij's of stencil points and weights
        subroutine GetInterp(lat,lon,ij0,IJs,Ws,isGs)
            real(rp), intent(in)  :: lat,lon
            integer , intent(in)  :: ij0(2)
            integer , intent(out) :: IJs(Np,2)
            real(rp), intent(out) :: Ws(Np)
            logical , intent(out) :: isGs(Np)

            integer :: i0,j0,n,di,dj,ip,jp
            real(rp) :: colat,dcolat,dlon,eta,zeta
            real(rp), dimension(-1:+1) :: wE,wZ

            !Single point
            isGs = .true.
            IJs(:,:) = 1
            Ws = 0.0
            IJs(1,:) = [ij0]
            Ws (1  ) = 1.0
            associate(gcolat=>RCMApp%gcolat,glong=>RCMApp%glong, &
            	      nLat=>RCMApp%nLat_ion,nLon=>RCMApp%nLon_ion, &
            	      toMHD=>RCMApp%toMHD)

            i0 = ij0(1)
            j0 = ij0(2)

            if ( (i0==1) .or. (i0==nLat) ) return !Don't bother if you're next to lat boundary
            
            !Get index space mapping: eta,zeta in [-0.5,0.5]
            colat = PI/2 - lat
            dcolat = ( gcolat(i0+1)-gcolat(i0-1) )/2
            dlon  = glong(2)-glong(1) !Assuming constant spacing

            eta  = ( colat - gcolat(i0) )/ dcolat
            zeta = ( lon - glong(j0) )/dlon

            !Clamp mappings
            call ClampMap(eta)
            call ClampMap(zeta)
            !Calculate weights
            call weight1D(eta,wE)
            call weight1D(zeta,wZ)

            n = 1
            do dj=-1,+1
                do di=-1,+1
                    ip = i0+di
                    jp = j0+dj
                    !Wrap around boundary, repeated point at 1/isize
                    if (jp<1)    jp = nLon-1
                    if (jp>nLon) jp = 2
                    IJs(n,:) = [ip,jp]
                    Ws(n) = wE(di)*wE(dj)
                    isGs(n) = toMHD(ip,jp)
                    n = n + 1
                enddo
            enddo !dj

            end associate            
        end subroutine GetInterp

        !1D triangular shaped cloud weights
        !1D weights for triangular shaped cloud interpolation
        !Assuming on -1,1 reference element, dx=1
        !Check for degenerate cases ( |eta| > 0.5 )
        subroutine weight1D(eta,wE)
            real(rp), intent(in)  :: eta
            real(rp), intent(out) :: wE(-1:1)

            wE(-1) = 0.5*(0.5-eta)**2.0
            wE( 1) = 0.5*(0.5+eta)**2.0
            wE( 0) = 0.75 - eta**2.0

        end subroutine weight1D

        !Clamps mapping in [-0.5,0.5]
        subroutine ClampMap(ez)
          REAL(rprec), intent(inout) :: ez
          if (ez<-0.5) ez = -0.5
          if (ez>+0.5) ez = +0.5
        end subroutine ClampMap

        function AvgQ(Q,IJs,Ws,Ni,Nj) 
            integer , intent(in) :: Ni,Nj
            integer , intent(in) :: IJs(Np,2)
            real(rp), intent(in) :: Ws(Np)
            real(rp), intent(in) :: Q(Ni,Nj)

            real(rp) :: AvgQ
            integer :: n,i0,j0
            real(rp) :: Qs(Np)
            AvgQ = 0.0

            do n=1,Np
                i0 = IJs(n,1)
                j0 = IJs(n,2)
                Qs(n) = Q(i0,j0)
            enddo
            AvgQ = dot_product(Qs,Ws)
            
        end function AvgQ

        subroutine GetRCMLoc(lat,lon,ij0)
            real(rp), intent(in) :: lat,lon
            integer, intent(out) :: ij0(2)

            integer :: iX,jX,iC,n
            real(rp) :: colat,dp,dcol,dI,dJ

            associate(gcolat=>RCMApp%gcolat,glong=>RCMApp%glong, &
                      nLat=>RCMApp%nLat_ion,nLon=>RCMApp%nLon_ion)

            !Assuming constant lon spacing
            dp = glong(2) - glong(1)

            !Get colat point
            colat = PI/2 - lat
!Use findloc w/ intel for speed
#ifdef __INTEL_COMPILER
            iC = findloc(gcolat >= colat,.true.,dim=1) - 1
#else 
!Bypass as findloc does not work for gfortran<9    
           !Work-around code        
            do n=1,nLat
                if (gcolat(n) >= colat) exit
            enddo
            iC = n-1
#endif
            dcol = gcolat(iC+1)-gcolat(iC)
            dI = (colat-gcolat(iC))/dcol
            if (dI <= 0.5) then
                iX = iC
            else
                iX = iC+1
            endif

            !Get lon point
            dJ = lon/dp
            if ( (dJ-floor(dJ)) <= 0.5 ) then
                jX = floor(dJ)+1
            else
                jX = floor(dJ)+2
            endif
            
            !Impose bounds just in case
            iX = max(iX,1)
            iX = min(iX,nLat)
            jX = max(jX,1)
            jX = min(jX,nLon)

            ij0 = [iX,jX]

            end associate
        end subroutine GetRCMLoc

    end subroutine InterpRCM    
end module rcmeval