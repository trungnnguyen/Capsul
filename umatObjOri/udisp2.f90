! user subroutine for a prescribed boundary conditions
  subroutine disp(u, kstep, kinc, time, node, noel, jdof, coords)
    use init

    implicit none

    real(kind = 8),    intent(out) :: u(3)
    integer(kind = 4), intent(in)  :: kstep
    integer(kind = 4), intent(in)  :: kinc
    real(kind = 8),    intent(in)  :: time(2)
    integer(kind = 4), intent(in)  :: node
    integer(kind = 4), intent(in)  :: noel
    integer(kind = 4), intent(in)  :: jdof
    real(kind = 8),    intent(in)  :: coords(3)

    ! nominal strain rate
    real(kind = 8), parameter :: tensilRate = 5000.0

    real(kind = 8) :: X0(3)

    X0 = getInitCoords(node)
    if (jdof == 3) then  
      u(1) = tensilRate*time(2)*X0(3)
    end if

  end subroutine disp
