module gamdebug
	use kdefs
	use gamtypes
    contains

	!Checks periodicity for metric terms on LFM-style grid
    subroutine ChkMetricLFM(Model,Grid)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid

        integer :: i,j,k

        real(rp) :: dEi,dEj
        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dEi = sum(abs(Grid%Te(i,j,Grid%ks,:,IDIR) - Grid%Te(i,j,Grid%ke+1,:,IDIR)))
                dEj = sum(abs(Grid%Te(i,j,Grid%ks,:,JDIR) - Grid%Te(i,j,Grid%ke+1,:,JDIR)))

                if ( (dEi+dEj)> TINY ) then
                    write(*,*) 'Te: i/j, dEi/dEj = ',i,j,dEi,dEj
                endif
            enddo
        enddo

        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dEi = sum(abs(Grid%Teb(i,j,Grid%ks,:,IDIR) - Grid%Teb(i,j,Grid%ke+1,:,IDIR)))
                dEj = sum(abs(Grid%Teb(i,j,Grid%ks,:,JDIR) - Grid%Teb(i,j,Grid%ke+1,:,JDIR)))

                if ( (dEi+dEj)> TINY ) then
                    write(*,*) 'Teb: i/j, dEi/dEj = ',i,j,dEi,dEj
                endif
            enddo
        enddo

        do i=Grid%is,Grid%ie+1
            do j=Grid%js,Grid%je+1
                dEi = abs(Grid%edge(i,j,Grid%ks,IDIR) - Grid%edge(i,j,Grid%ke+1,IDIR))
                dEj = abs(Grid%edge(i,j,Grid%ks,JDIR) - Grid%edge(i,j,Grid%ke+1,JDIR))

                if ( (dEi+dEj)> TINY ) then
                    write(*,*) 'edge: i/j, dEi/dEj = ',i,j,dEi,dEj
                endif
            enddo
        enddo


    end subroutine ChkMetricLFM

end module gamdebug