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

    real(kind = 8) :: engShear

    real(kind = 8), parameter :: shearRate = 6.667d-4
    real(kind = 8), parameter :: eps = 1.0d-5

    engShear = 2.0d0*(exp(shearRate) - 1.0d0)*time(2)
    if (abs(coords(3)) < eps .or. abs(coords(3) - 1.0d0) < eps) then
      if (jdof == 1 .and. jdof == 3) then
        u(1) = 0.0d0
      else if (jdof == 2) then 
        u(1) = engShear*coords(3)
      end if
    end if

  end subroutine disp
