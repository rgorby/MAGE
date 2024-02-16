! Heap of helpful functions for ShellGrids and ShellGridVars
module shellUtils

    use shellGrid

    implicit none

    contains


    subroutine wrapJ_SGV(sh, sgVar)
        !! Wrap a ShellGridVar around the periodic j/phi boundary
        !! Assumes all vars within active domain are valid to wrap
        !! Note: anything in the i-direction ghost cells get wrapped too
        type(ShellGrid_T), intent(in) :: sh
        type(ShellGridVar_T), intent(inout) :: sgVar

        associate(Q => sgVar%data)

        ! Theta faces are cell-centered w.r.t. j direction
        if (sgVar%loc == SHCC .or. sgVar%loc == SHFTH) then
            ! Starting ghost cells
            Q(:, sh%jsg:sh%js-1) = Q(:, sh%je-sh%Ngw+1:sh%je)
            ! Ending ghosts cells
            Q(:, sh%je+1:sh%jeg) = Q(:, sh%js:sh%js+sh%Nge-1)
        elseif (sgVar%loc==SHCORNER .or. sgVar%loc == SHFPH) then
            ! Starting ghost cells
            Q(:, sh%jsg:sh%js) = Q(:, sh%je-sh%Ngw+1:sh%je+1)
            ! Ending ghosts cells
            Q(:, sh%je+1:sh%jeg+1) = Q(:, sh%js:sh%js+sh%Nge)
        endif

        end associate

    end subroutine wrapJ_SGV


end module shellUtils