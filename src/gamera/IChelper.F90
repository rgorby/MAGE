module ichelper

    ! contain helper functions to mask
    ! the need for mpi comands in the
    ! non-mpi compilations/runs

    use gdefs
    use types
    implicit none

    !VectorField_T
    !Generic vector field
    abstract interface
        subroutine reduced_T(Rio)
            Import :: rp
            real(rp), intent(inout) :: Rio(:,:)
        end subroutine reduced_T
    end interface

    procedure(reduced_T), pointer :: reduced

    abstract interface
        subroutine GlobalGrid_T(Model,lgrid,Grid)
            Import Model_T, Grid_T
            type(Model_T), intent(in) :: Model
            type(Grid_T),  intent(in) :: lGrid 
            type(Grid_T),  intent(inout) :: Grid
        end subroutine GlobalGrid_T
    end interface

    procedure(GlobalGrid_T), pointer :: GlobalGrid

    contains

    subroutine single_reduced(R)
          real(rp), intent(inout) :: R(:,:)
    
          return
    end subroutine single_reduced

    !For some inital conditions we need to know the global grid cordinates
    ! But NO ALLOCATED ARRAYS!!!!!!!!!!!!!!!!!!!!!
    subroutine single_GlobalGrid(Model,lgrid,Grid)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: lGrid
        type(Grid_T),  intent(inout) :: Grid

        !Derived quantities 
        Grid%Ni = lGrid%Nip + 2*Model%nG
        Grid%Nj = lGrid%Njp + 2*Model%nG
        Grid%Nk = lGrid%Nkp + 2*Model%nG

        Grid%Nip = lGrid%Nip
        Grid%Njp = lGrid%Njp
        Grid%Nkp = lGrid%Nkp

        Grid%is = 1; Grid%ie = lGrid%Nip
        Grid%js = 1; Grid%je = lGrid%Njp
        Grid%ks = 1; Grid%ke = lGrid%Nkp

        Grid%isg = Grid%is-Model%nG
        Grid%ieg = Grid%ie+Model%nG

        Grid%jsg = Grid%js-Model%nG
        Grid%jeg = Grid%je+Model%nG

        Grid%ksg = Grid%ks-Model%nG
        Grid%keg = Grid%ke+Model%nG
    end subroutine single_GlobalGrid




     subroutine initIChelper()

        reduced => single_reduced
        GlobalGrid => single_GlobalGrid 

     end subroutine initIChelper

end module ichelper
