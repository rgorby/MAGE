program MIX
    ! base
    use xml_input

    ! gamera
    use types
    use gioH5
    use init

    ! mix
    use mixtypes

    ! interfaces
    use mix_mhd_interface
    use mhd_mix_interface

    implicit none

    ! Gamera file
    character(len=strLen) :: inpXML,inH5

    ! Gamera variables
    type(Model_T) :: Model
    type(Grid_T)  :: Grid
    type(State_T) :: State

    ! Remix variables
    type(mixConductance_T) :: conductance

    ! Input deck
    type(XML_Input_T) :: xmlInp

    integer :: Narg

    ! Input deck
    Narg = command_argument_count()
    if (Narg .eq. 0) then
       write(*,*) 'No input deck specified, defaulting to Input.xml'
       inpXML = "Input.xml"
    else
       call get_command_argument(1,inpXML)
    endif

    write(*,*) 'Reading input deck from ', trim(inpXML)
    inquire(file=inpXML,exist=fExist)
    if (.not. fExist) then
       write(*,*) 'Error opening input deck, exiting ...'
       write(*,*) ''
       stop
    endif

    !Partial inclusion of new XML reader
    xmlInp = New_XML_Input(trim(inpXML),'Gamera',.true.)
    !--------------------------------------------

    inH5 = "/glade/scratch/skareem/octruntemp/msphere_0013_0001_0001_0000_0000_0000_000000000000.Res.00136.h5"

    ! Initalize Gamera model data structure
    call initModel(Model,xmlInp)
    call readH5Grid(Model,Grid,inH5)
    !Turn corners into grid
    call Corners2Grid(Model,Grid)
    call readH5Restart(Model,Grid,State,inH5)

    ! Note, remix reads the xml file on its own
    ! FIXME: should just pass the xmlInp handler to it
    call InitCMI(Model,Grid)
    call init_mix_mhd_interface(ion,hmsphrs,mhdJGrid,mhdPsiGrid,conductance)  ! the arguments are globals in mhd_mix_interface
    call PrepRemixData(Model,Grid,State)

!     !! EVERYTHING ENCLOSED IN !!!!! IS A HACK FOR TESTING IN STANDALONE MODE
!     !! TO BE REMOVED IN THE CMI_REMIX INTERFACE

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     call h5open_f(herr) !Setup H5 Fortran interface
!     call mhd_fromFile("/glade/u/home/skareem/hao/magnetosphere_data/quadexample/msphere_0007_0001_0001_0000_0000_0000_000000000000.h5", 1, mhdg)
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ! Init mix (both hemispheres) and interpolation maps
!     call init_mix_mhd_interface(ion, hmsphrs, mhdg, mhdg)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !    call readGMHDVar("../data/msphere.h5", "Step#1", "Jz", gJz)
!     call readGMHDVar("/glade/u/home/skareem/hao/magnetosphere_data/quadexample/msphere_0007_0001_0001_0000_0000_0000_000000000000.h5","Step#100","Jz",gJz)

!     ni = size(gJz, 1); nj = size(gjz, 2); nk = size(gJz, 3)
!     allocate (mhdvars(1, nj, nk/2, 2, 2)) ! (i,j,k,var,hemisphere)
!     mhdvars(1, :, :, 1, 1) = -gJz(1, 1:nj, 1:nk/2)
!     mhdvars(1, :, :, 1, 2) = -gJz(1, 1:nj, 1:nk/2)
!     mhdvars(1, :, :, 2, 1) = 2*gJz(1, 1:nj, 1:nk/2)
!     mhdvars(1, :, :, 2, 2) = -2*gJz(1, 1:nj, 1:nk/2)
!     call h5close_f(herr) ! Close H5 Fortran interface
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     ! get tilt from MHD eventually
!     call mhd2mix(ion, mhdvars, -9999._rp)

!     call run_mix(remixApp%ion,remixApp%tilt,remixApp%conductance)

!     ! now interpolate from mix to gamera
!     allocate (psi(4, nj + 1, nk/2 + 1, 1, 2)) ! (i,j,k,var,hemisphere)
!     call mix2mhd(ion, psi)
!     write (*, *) minval(psi(1, :, :, 1, 1)), maxval(psi(1, :, :, 1, 1))
!     write (*, *) minval(psi(1, :, :, 1, 2)), maxval(psi(1, :, :, 1, 2))

!     call writeMIX('mixtest.h5', ion, hmsphrs)

end program MIX
