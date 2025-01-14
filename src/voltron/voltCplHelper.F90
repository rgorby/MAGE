module voltCplHelper

    use kdefs
    use voltTypes
    use voltCplTypes
    use imagtubes

    implicit none

    contains

    subroutine genVoltTubes(vApp)
        !! Update Tube_Ts on Voltron's grid
        type(voltApp_T), intent(inout) :: vApp

        integer :: i,j
        real(rp) :: seedR, eqR, mhd_Rin
        type(magLine_T) :: magLine

        associate(sh=>vApp%shGrid, ebApp=>vApp%ebTrcApp, Gr=>vApp%gApp%Grid)
            mhd_Rin = norm2(Gr%xyz(Gr%is,Gr%js,Gr%ks,:))
            seedR = sh%radius  ! Ionosphere radius in Rp
            ! Do field line tracing, populate fromV%ijTubes
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,eqR,magLine)
            do i=sh%isg,sh%ieg+1
                do j=sh%jsg,sh%jeg+1
                    call CleanLine(magLine)

                    eqR = DipColat2L(sh%thRp(i))  ! Function assumes colat coming from 1 Rp, make sure we use the right theta value
                    if (eqR < mhd_Rin) then
                        call DipoleTube(vApp, sh%th(i), sh%ph(j), vApp%State%ijTubes(i,j))
                    else
                        call MHDTube(ebApp, vApp%planet,   & !ebTrcApp, planet
                            sh%th(i), sh%ph(j), seedR, &  ! colat, lon, r
                            vApp%State%ijTubes(i,j), magLine, &  ! IMAGTube_T, magLine_T
                            doShiftO=.true.)
                    endif

                enddo
            enddo

        end associate
    end subroutine genVoltTubes

end module voltCplHelper