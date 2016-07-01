  program testTimeCost
    implicit none

    real(kind = 8) :: array1(90), array2(30, 3)
    integer(kind = 4), parameter :: loopCnt = 1e6
    real(kind = 8) :: timeStart, timeFinish

    integer :: i

    array1 = 1.0d0

    call cpu_time(timeStart)
    do i = 1, loopCnt
      array2 = reshape(array1, (/30, 3/))
    end do

    call cpu_time(timeFinish)

    write(*, *) "Start: ", timeStart, " Finish: ", timeFinish, " Time Cost: ", timeFinish - timeStart

    call testStaticVar()
    call testStaticVar()

     
  end program testTimeCost

  subroutine testStaticVar()
    implicit none

    logical :: ltest = .true.
    real    :: dtest = 1.0d0

    write(*, *) "Before:  ltest = ", ltest, " dtest = ", dtest
    ltest = .not. ltest
    dtest = dtest + 1.0d0

    write(*, *) "After :  ltest = ", ltest, " dtest = ", dtest

  end subroutine testStaticVar
