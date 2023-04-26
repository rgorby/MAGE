module sifCplTypes
    !! Types used to couple sif to voltron
    use volttypes
    use imagtubes
    use ebtypes

    use siftypes


    implicit none

    type sif_fromV_T
        real(rp) :: tLastUpdate
            !! Time of last update, according to voltron
        type(fLine_T), dimension(:,:), allocatable :: fLines
        type(IMAGTube_T), dimension(:,:), allocatable :: ijTubes
            !! imagTubes with field-line averaged info
        real(rp), dimension(:,:), allocatable :: pot
            !! electrostatic potential from ionosphere [kV]

    end type sif_fromV_T
    

    type sif_toV_T
        real(rp) :: tLastUpdate
            !! Time of last update, according to SIF

    end type sif_toV_T


    type sif_cplBase_T
        type(sif_fromV_T) :: fromV
        type(sif_toV_T  ) :: toV

        procedure(sifCpl_v2s_T), pointer ::  convertToSIF     => NULL()
            !! Updates fromV object using Voltron information, then puts it into SIF
        procedure(sifCpl_s2v_T), pointer ::  convertToVoltron => NULL()
            !! Updates toV object using SIF information, then puts it into Voltron

    end type sif_cplBase_T


   


    abstract interface
        subroutine sifCpl_v2s_T(cplBase, vApp, sApp)
            Import :: sif_cplBase_T, voltApp_T, sifApp_T
            class(sif_cplBase_T), intent(inout) :: cplBase
            type(voltApp_T), intent(in   ) :: vApp
            type(sifApp_T ), intent(inout) :: sApp
        end subroutine sifCpl_v2s_T

        subroutine sifCpl_s2v_T(cplBase, vApp, sApp)
            Import :: sif_cplBase_T, voltApp_T, sifApp_T
            class(sif_cplBase_T), intent(inout) :: cplBase
            type(voltApp_T), intent(inout) :: vApp
            type(sifApp_T) , intent(in   ) :: sApp
        end subroutine sifCpl_s2v_T

    end interface


end module sifCplTypes