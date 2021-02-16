!Routines to map native source grids to Bios-Savart (ground system) grids
module calcdbremap
	use kdefs
	use chmpdefs
	use ebtypes
    use calcdbtypes
    use ebinterp
    use chmpfields

	implicit none

	integer, parameter, private :: i0 = 1

	contains

	!Get magBS,ionBS,facBS using current MJD and source grids
	subroutine remapBS(Model,t,ebState,ionGrid,facGrid,magBS,ionBS,facBS)
        type(chmpModel_T), intent(in) :: Model
        real(rp)         , intent(in) :: t
        type(ebState_T)  , intent(in) :: ebState
        type(ionGrid_T)  , intent(in) :: ionGrid
        type(facGrid_T)  , intent(in) :: facGrid
        type(BSGrid_T), intent(inout) :: magBS,ionBS,facBS

        real(rp) :: mjd,w1,w2,dV
        integer :: i,j,k,l,n
        real(rp), dimension(:,:,:,:), allocatable :: Jxyz
        real(rp), dimension(NDIM) :: smX,smJ,geoX,geoJ

        mjd = MJDAt(ebState%ebTab,t) !Current MJD

        !CALCDB-TODO: Add any initialization for geopack at a given time here

        !CALCDB-TODO: Make sm2geo actually do something

  	!----
  	!Magnetospheric part
  		
  	!Start by getting Jxyz at current time
  		!Jxyz are SM currents at time t
  		associate(ebGr=>ebState%ebGr)
  		allocate(Jxyz(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        call GetTWgts(Model,ebState,t,w1,w2)
        !$OMP PARALLEL WORKSHARE
    	Jxyz = w1*ebState%eb1%Jxyz + w2*ebState%eb2%Jxyz
        !$OMP END PARALLEL WORKSHARE
        
    !Now do remap and store into BSGrid
    	do k=ebGr%ks,ebGr%ke
			do j=ebGr%js,ebGr%je
				do i=ebGr%is,ebGr%ie
					!Get n-d => 1D index
					n = ijk2n(i,j,k,ebGr%is,ebGr%ie,ebGr%js,ebGr%je,ebGr%ks,ebGr%ke)

					smX = ebGr%xyz(i,j,k,:) !Cell-centered MHD grid coordinates (SM)
					smJ =     Jxyz(i,j,k,:) !SM current at cell-center
					dV  = ebGr%dV(i,j,k) !Volume element, Re^3

					if (i <= ebGr%is+i0-1) then
						dV = 0.0
						!Zap this contribution
					endif

					call sm2geo(smX,smJ,mjd,geoX,geoJ)
					magBS%XYZcc(n,XDIR:ZDIR) = geoX
					magBS%Jxyz (n,XDIR:ZDIR) = geoJ
					magBS%dV(n) = dV

				enddo
			enddo
		enddo

		end associate

  	!----
  	!Ionospheric part
  		do k=1,2 !Hemisphere
  			do j=1,ionGrid%Nth !Theta
  				do i=1,ionGrid%Np !phi
                    !Get n-d => 1D index
                    n = ijk2n(i,j,k,1,ionGrid%Np,1,ionGrid%Nth,1,2)

  					smX = ionGrid%XYZcc(i,j,k,:)
  					smJ = ionGrid%Jxyz (i,j,k,:)
  					dV  = ionGrid%dS(i,j,k)
  					call sm2geo(smX,smJ,mjd,geoX,geoJ)
					ionBS%XYZcc(n,XDIR:ZDIR) = geoX
					ionBS%Jxyz (n,XDIR:ZDIR) = geoJ
					ionBS%dV(n) = dV

  				enddo
  			enddo
  		enddo

  	!----
  	!FAC part

  		do l=1,2 !Hemisphere
  			do k=1,facGrid%rSegs
  				do j=1,facGrid%Nth
  					do i=1,facGrid%Np
                        !Get n-d => 1D index
                        n = ijkl2n(i,j,k,l,1,ionGrid%Np,1,ionGrid%Nth,1,facGrid%rSegs,1,2)

  						smX = facGrid%XYZcc(i,j,k,l,:)
  						smJ = facGrid%Jxyz (i,j,k,l,:)
  						dV  = facGrid%dV   (i,j,k,l)
  						call sm2geo(smX,smJ,mjd,geoX,geoJ)
						facBS%XYZcc(n,XDIR:ZDIR) = geoX
						facBS%Jxyz (n,XDIR:ZDIR) = geoJ
						facBS%dV(n) = dV
					enddo !i
				enddo
			enddo
		enddo !l

	end subroutine remapBS


	!Convert coordinate and vector from SM to GEO
	!CALCDB-TODO: Make sm2geo actually do something, right now it's just doing identity
	subroutine sm2geo(smX,smJ,mjd,geoX,geoJ)
		real(rp), intent(in ) :: smX(NDIM),smJ(NDIM),mjd
		real(rp), intent(out) :: geoX(NDIM),geoJ(NDIM)
		geoX = smX
		geoJ = smJ
	end subroutine sm2geo

    !TODO: Rewrite these routines to be smarter
    function ijk2n(i,j,k,is,ie,js,je,ks,ke) result(n)
        integer, intent(in) :: i,j,k,is,ie,js,je,ks,ke
        integer :: n

        integer :: ip,jp,kp
        integer :: Ni,Nj

        ip = i-is+1
        jp = j-js+1
        kp = k-ks+1
        Ni = ie-is+1
        Nj = je-js+1

        !Convert to n
        n = ip + (jp-1)*Ni + (kp-1)*Ni*Nj
    end function ijk2n

    function ijkl2n(i,j,k,l,is,ie,js,je,ks,ke,ls,le) result(n)
        integer, intent(in) :: i,j,k,l,is,ie,js,je,ks,ke,ls,le
        integer :: n

        integer :: ip,jp,kp,lp
        integer :: Ni,Nj,Nk

        ip = i-is+1
        jp = j-js+1
        kp = k-ks+1
        lp = l-ls+1
        Ni = ie-is+1
        Nj = je-js+1
        Nk = ke-ks+1

        !Convert to n
        n = ip + (jp-1)*Ni + (kp-1)*Ni*Nj + (lp-1)*Ni*Nj*Nk
    end function ijkl2n

end module calcdbremap