module mixstate
  use mixdefs
  use mixtypes

  implicit none
  
  contains
    subroutine init_state_fromfile(G,St,pot_in,j_in,sigmap_in,sigmah_in,rho_in,cs_in)
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      real(rp), dimension(:,:) :: pot_in,j_in,sigmap_in,sigmah_in,rho_in,cs_in

      ! set up variable state
      if (.not.allocated(St%Vars)) allocate(St%Vars(G%Np,G%Nt,nVars))
      St%Vars(:,:,POT) = pot_in
      St%Vars(:,:,FAC) = j_in
      St%Vars(:,:,SIGMAP) = sigmap_in
      St%Vars(:,:,SIGMAH) = sigmah_in
      St%Vars(:,:,DENSITY) = rho_in
      St%Vars(:,:,SOUND_SPEED) = cs_in
    end subroutine init_state_fromfile

    subroutine init_state(G,St)
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      ! set up variable state
      if (.not.allocated(St%Vars)) allocate(St%Vars(G%Np,G%Nt,nVars))
      St%Vars(:,:,:) = 0.0D0
    end subroutine init_state
end module mixstate
