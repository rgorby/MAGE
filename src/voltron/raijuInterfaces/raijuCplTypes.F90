module raijuCplTypes
    !! Types used to couple RAIJU to voltron
    use volttypes
    use imagtubes
    use ebtypes

    use raijutypes


    implicit none

    type raiju_fromV_T
        real(rp) :: tLastUpdate
            !! Time of last update, according to voltron
        type(fLine_T), dimension(:,:), allocatable :: fLines
        type(IMAGTube_T), dimension(:,:), allocatable :: ijTubes
            !! imagTubes with field-line averaged info
        real(rp), dimension(:,:), allocatable :: pot
            !! electrostatic potential from ionosphere [kV]
        procedure(raijuMHD2SpcMap_T), pointer, nopass :: mhd2spcMap => NULL()
            !! This is called by convertToRAIJU to map MHD fluids to RAIJU species

    end type raiju_fromV_T
    

    type raiju_toV_T
        real(rp) :: tLastUpdate
            !! Time of last update, according to RAIJU

    end type raiju_toV_T


    type raiju_cplBase_T
        type(raiju_fromV_T) :: fromV
        type(raiju_toV_T  ) :: toV

        procedure(raijuCpl_v2s_T), pointer ::  convertToRAIJU     => NULL()
            !! Updates fromV object using Voltron information, then puts it into RAIJU
        procedure(raijuCpl_s2v_T), pointer ::  convertToVoltron => NULL()
            !! Updates toV object using RAIJU information, then puts it into Voltron
    end type raiju_cplBase_T



    abstract interface

        subroutine raijuCpl_v2s_T(cplBase, vApp, sApp)
            Import :: raiju_cplBase_T, voltApp_T, raijuApp_T
            class(raiju_cplBase_T), intent(in) :: cplBase
            type(voltApp_T), intent(in   ) :: vApp
            type(raijuApp_T ), intent(inout) :: sApp
        end subroutine raijuCpl_v2s_T


        subroutine raijuCpl_s2v_T(cplBase, vApp, sApp)
            Import :: raiju_cplBase_T, voltApp_T, raijuApp_T
            class(raiju_cplBase_T), intent(inout) :: cplBase
            type(voltApp_T), intent(inout) :: vApp
            type(raijuApp_T) , intent(in   ) :: sApp
        end subroutine raijuCpl_s2v_T


        subroutine raijuMHD2SpcMap_T(Model, Grid, State, ijTubes)
            Import :: raijuModel_T, raijuGrid_T, raijuState_T, IMAGTube_T
            type(raijuModel_T) , intent(in) :: Model
            type(raijuGrid_T)  , intent(in) :: Grid
            type(raijuState_T) , intent(inout) :: State
            type(IMAGTube_T),  dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                                        Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: ijTubes
        end subroutine raijuMHD2SpcMap_T
    end interface


end module raijuCplTypes