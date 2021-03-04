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

      rm%gcolat (:) = colat (:,1)
      rm%glong  (:) = aloct (1,jwrap:jsize)
      if (rm%glong(rm%nLon_ion) < pi) rm%glong(rm%nLon_ion) = rm%glong(rm%nLon_ion) + 2*pi
    END SUBROUTINE setupIon

    !> Deallocate any variables allocated by setupIon.
    SUBROUTINE tearDownIon(rm)
      type(rcm_mhd_T),intent(inout) ::RM

      if (ALLOCATED(rm%pot)) DEALLOCATE(rm%pot)
      if (ALLOCATED(rm%sigmap)) DEALLOCATE(rm%sigmap)
      if (ALLOCATED(rm%sigmah)) DEALLOCATE(rm%sigmah)
      if (ALLOCATED(rm%gcolat)) DEALLOCATE(rm%gcolat)
      if (ALLOCATED(rm%glong)) DEALLOCATE(rm%glong)
      if (ALLOCATED(rm%flux)) DEALLOCATE(rm%flux)
      if (ALLOCATED(rm%fac)) DEALLOCATE(rm%fac)
      if (ALLOCATED(rm%Pave)) DEALLOCATE(rm%Pave)
      if (ALLOCATED(rm%Nave)) DEALLOCATE(rm%Nave)
      if (ALLOCATED(rm%Vol)) DEALLOCATE(rm%Vol)
      if (ALLOCATED(rm%Bmin)) DEALLOCATE(rm%Bmin)
      if (ALLOCATED(rm%X_bmin)) DEALLOCATE(rm%X_bmin)
      if (ALLOCATED(rm%iopen)) DEALLOCATE(rm%iopen)
      if (ALLOCATED(rm%Prcm)) DEALLOCATE(rm%Prcm)
      if (ALLOCATED(rm%Npsph)) DEALLOCATE(rm%Npsph)
      if (ALLOCATED(rm%Nrcm)) DEALLOCATE(rm%Nrcm)
      if (ALLOCATED(rm%beta_average)) DEALLOCATE(rm%beta_average)
      if (ALLOCATED(rm%latc)) DEALLOCATE(rm%latc)
      if (ALLOCATED(rm%lonc)) DEALLOCATE(rm%lonc)
      if (ALLOCATED(rm%toMHD)) DEALLOCATE(rm%toMHD)
      if (ALLOCATED(rm%losscone)) DEALLOCATE(rm%losscone)
      if (ALLOCATED(rm%oxyfrac)) DEALLOCATE(rm%oxyfrac)
      if (ALLOCATED(rm%Percm)) DEALLOCATE(rm%Percm)
      if (ALLOCATED(rm%radcurv)) DEALLOCATE(rm%radcurv)
      if (ALLOCATED(rm%wIMAG)) DEALLOCATE(rm%wIMAG)

    END SUBROUTINE tearDownIon

  END MODULE ionosphere_exchange

