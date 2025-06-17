module voltCplHelper

    use kdefs
    use voltTypes
    use voltCplTypes
    use tubehelper

    implicit none

    contains

    subroutine genVoltTubes(vApp)
        !! Update Tube_Ts on Voltron's grid
        class(voltApp_T), intent(inout) :: vApp

        ! calculate field lines
        call calcTubes(vApp)

        ! Now pack into tubeShell
        call tubes2Shell(vApp%shGrid, vApp%State%ijTubes, vApp%State%tubeShell)

    end subroutine genVoltTubes

    subroutine calcTubes(vApp)
        !! Calculate field lines in ijTubes
        class(voltApp_T), intent(inout) :: vApp

        integer :: i,j
        real(rp) :: seedR, eqR, mhd_Rin
        real(rp), dimension(NDIM) :: xyz0
        logical :: doSH,doNH
        type(magLine_T) :: magLine

        !associate(sh=>vApp%shGrid, ebApp=>vApp%ebTrcApp, Gr=>vApp%gApp%Grid)
            !mhd_Rin = norm2(Gr%xyz(Gr%is+2,Gr%js,Gr%ks,:))
        associate(sh=>vApp%shGrid, ebApp=>vApp%ebTrcApp, ebGr=>vApp%ebTrcApp%ebState%ebGr)
            mhd_Rin = norm2(ebGr%xyz(ebGr%is+2,ebGr%js,ebGr%ks,:))
            seedR = sh%radius  ! Ionosphere radius in Rp
            ! Do field line tracing, populate ijTubes
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,eqR,magLine,doSH,doNH,xyz0)
            do j=sh%js,sh%je+1
                do i=sh%is,sh%ie+1
                    !Calculate seed point

                    xyz0 = seedR*[sin(sh%th(i))*cos(sh%ph(j)), &
                                  sin(sh%th(i))*sin(sh%ph(j)), &
                                  cos(sh%th(i))]
                    eqR = DipColat2L(sh%thRp(i))  ! Function assumes colat coming from 1 Rp, make sure we use the right theta value
                    if (eqR .le. mhd_Rin) then
                        !No MHD to tube from
                        call DipoleTube(vApp%planet,xyz0,vApp%State%ijTubes(i,j))
                    else
                        xyz0 = DipoleShift(xyz0, seedR + TINY)
                        if (xyz0(ZDIR) < 0) then
                            doNH = .false.
                            doSH = .true.
                        else
                            doNH = .true.
                            doSH = .false.
                        endif

                        call CleanLine(magLine)
                        !Note: Not using volt time b/c chimp wants time in its units
                        call genLine(ebApp%ebModel,ebApp%ebState,xyz0,ebApp%ebState%eb1%time, magLine,&
                                     doShueO=.false.,doNHO=doNH,doSHO=doSH)
                        call Line2Tube(ebApp,vApp%planet,magLine,vApp%State%ijTubes(i,j))
                    endif

                enddo
            enddo

        end associate

    end subroutine calcTubes

end module voltCplHelper
