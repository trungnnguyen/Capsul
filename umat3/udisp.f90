! user subroutine for a prescribed boundary conditions
  subroutine disp(u, kstep, kinc, time, node, noel, jdof, coords)
    implicit none

    real(kind = 8),    intent(out) :: u(3)
    integer(kind = 4), intent(in)  :: kstep
    integer(kind = 4), intent(in)  :: kinc
    real(kind = 8),    intent(in)  :: time(2)
    integer(kind = 4), intent(in)  :: node
    integer(kind = 4), intent(in)  :: noel
    integer(kind = 4), intent(in)  :: jdof
    real(kind = 8),    intent(in)  :: coords(3)

    real(kind = 8), parameter :: LenY = 1.0

    if (node == 999999 .and. jdof == 2) then
      u(1) = LenY*(dexp(time(2)) - 1.0d0)
    end if


  end subroutine disp
