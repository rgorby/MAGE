MODULE ionosphere_exchange
  use rcm_mhd_interfaces
  use rcm_mod_subs, ONLY: isize, jsize, jwrap, colat, aloct
  use kdefs, ONLY : PI
  contains 
    
    !> Allocate Ionosphere Grid variables and read Ion grid from "RCM-ion.dat".
    !! Don't forget to deallocate!
    SUBROUTINE setupIon(RM)
      IMPLICIT NONE
      type(rcm_mhd_T),intent(inout) ::RM
      integer(iprec) :: lat,lon

      rm%nLat_ion = isize
      rm%nLon_ion = jsize-jwrap+1

      ALLOCATE( rm%gcolat(rm%nLat_ion) )
      ALLOCATE( rm%glong(rm%nLon_ion) )

      ALLOCATE( rm%pot(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%nflx(rm%nLat_ion, rm%nLon_ion, 2) )
      ALLOCATE( rm%eng_avg(rm%nLat_ion, rm%nLon_ion, 2) )
      ALLOCATE( rm%flux(rm%nLat_ion, rm%nLon_ion, 2) )
      ALLOCATE( rm%fac(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Pave(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Nave(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Vol(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Bmin(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%X_bmin(rm%nLat_ion, rm%nLon_ion, 3) )
      ALLOCATE( rm%iopen(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Prcm(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Npsph(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Nrcm(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%beta_average(rm%nLat_ion, rm%nLon_ion))
      ALLOCATE( rm%sigmap(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%sigmah(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%latc(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%lonc(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Lb  (rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Tb  (rm%nLat_ion, rm%nLon_ion) )

      ALLOCATE( rm%toMHD(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%losscone(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%radcurv(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%wIMAG(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%oxyfrac(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%Percm(rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%errD (rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%errP (rm%nLat_ion, rm%nLon_ion) )
      ALLOCATE( rm%nTrc (rm%nLat_ion, rm%nLon_ion) )

      rm%gcolat (:) = colat (:,1)
      rm%glong  (:) = aloct (1,jwrap:jsize)
      if (rm%glong(rm%nLon_ion) < pi) rm%glong(rm%nLon_ion) = rm%glong(rm%nLon_ion) + 2*pi

      ! initialize all other variables to 0
      rm%pot = 0.0
      rm%nflx = 0.0
      rm%eng_avg = 0.0
      rm%flux = 0.0
      rm%fac = 0.0
      rm%Pave = 0.0
      rm%Nave = 0.0
      rm%Vol = 0.0
      rm%Bmin = 0.0
      rm%X_bmin = 0.0
      rm%iopen = 0
      rm%Prcm = 0.0
      rm%Npsph = 0.0
      rm%Nrcm = 0.0
      rm%beta_average = 0.0
      rm%sigmap = 0.0
      rm%sigmah = 0.0
      rm%latc = 0.0
      rm%lonc = 0.0
      rm%Lb = 0.0
      rm%Tb = 0.0
      rm%toMHD = .False.
      rm%losscone = 0.0
      rm%radcurv = 0.0
      rm%wIMAG = 0.0
      rm%oxyfrac = 0.0
      rm%Percm = 0.0
      rm%errD = 0.0
      rm%errP = 0.0
      rm%nTrc = 0

    END SUBROUTINE setupIon

    !> Deallocate any variables allocated by setupIon.
    SUBROUTINE tearDownIon(rm)
      type(rcm_mhd_T),intent(inout) ::RM

      if (ALLOCATED(rm%pot))          DEALLOCATE(rm%pot)
      if (ALLOCATED(rm%sigmap))       DEALLOCATE(rm%sigmap)
      if (ALLOCATED(rm%sigmah))       DEALLOCATE(rm%sigmah)
      if (ALLOCATED(rm%gcolat))       DEALLOCATE(rm%gcolat)
      if (ALLOCATED(rm%glong))        DEALLOCATE(rm%glong)
      if (ALLOCATED(rm%flux))         DEALLOCATE(rm%flux)
      if (ALLOCATED(rm%fac))          DEALLOCATE(rm%fac)
      if (ALLOCATED(rm%Pave))         DEALLOCATE(rm%Pave)
      if (ALLOCATED(rm%Nave))         DEALLOCATE(rm%Nave)
      if (ALLOCATED(rm%Vol))          DEALLOCATE(rm%Vol)
      if (ALLOCATED(rm%Bmin))         DEALLOCATE(rm%Bmin)
      if (ALLOCATED(rm%X_bmin))       DEALLOCATE(rm%X_bmin)
      if (ALLOCATED(rm%iopen))        DEALLOCATE(rm%iopen)
      if (ALLOCATED(rm%Prcm))         DEALLOCATE(rm%Prcm)
      if (ALLOCATED(rm%Npsph))        DEALLOCATE(rm%Npsph)
      if (ALLOCATED(rm%Nrcm))         DEALLOCATE(rm%Nrcm)
      if (ALLOCATED(rm%beta_average)) DEALLOCATE(rm%beta_average)
      if (ALLOCATED(rm%latc))         DEALLOCATE(rm%latc)
      if (ALLOCATED(rm%lonc))         DEALLOCATE(rm%lonc)
      if (ALLOCATED(rm%toMHD))        DEALLOCATE(rm%toMHD)
      if (ALLOCATED(rm%losscone))     DEALLOCATE(rm%losscone)
      if (ALLOCATED(rm%oxyfrac))      DEALLOCATE(rm%oxyfrac)
      if (ALLOCATED(rm%Percm))        DEALLOCATE(rm%Percm)
      if (ALLOCATED(rm%radcurv))      DEALLOCATE(rm%radcurv)
      if (ALLOCATED(rm%wIMAG))        DEALLOCATE(rm%wIMAG)
      if (ALLOCATED(rm%errD))         DEALLOCATE(rm%errD)
      if (ALLOCATED(rm%errP))         DEALLOCATE(rm%errP)
      if (ALLOCATED(rm%nTrc))         DEALLOCATE(rm%nTrc)

    END SUBROUTINE tearDownIon

  END MODULE ionosphere_exchange

