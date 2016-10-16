  module expansion
    use utils, only : RKIND, IKIND, MAX_POLY_ORDER

    implicit none
    private

    public :: ReadExpansion
    public :: fDfGrdTh

    real(kind = RKIND) :: expnCoef(MAX_POLY_ORDER + 1, 3)
    real(kind = RKIND) :: tempRef
 
  contains

    subroutine ReadExpansion(fID)
      use utils, only : ReadLine

      integer(kind = IKIND), intent(in) :: fID

      ! local variables
      integer(kind = IKIND)  :: nOrder
      character(len = 256)   :: line
      integer(kind = IKIND)  :: ii, iStat

      call ReadLine(fID, line)
      read(line, *) tempRef
     
      call ReadLine(fID, line)
      read(line, *) nOrder

      expnCoef = 0.0d0
      do ii = 1, nOrder + 1
        call ReadLine(fID, line)
        read(line, *, iostat=iStat) expnCoef(ii, 1), expnCoef(ii, 2), expnCoef(ii, 3)
        if (iStat /= 0) then
          expnCoef(ii, 2) = expnCoef(ii, 1)  
          expnCoef(ii, 3) = expnCoef(ii, 1)  
        end if
      end do

      write(*, *) "expnCoef = ", expnCoef

    end subroutine ReadExpansion


    function fDfGrdTh(tempCur)
      use utils,   only : UNITMAT
      use algebra, only : Polynomial

      real(kind = RKIND), intent(in)  :: tempCur
      real(kind = RKIND) :: fDfGrdTh(3, 3)  

      real(kind = RKIND) :: alpha(3)
      integer(kind = IKIND) :: ii

      alpha(1) = Polynomial(expnCoef(:, 1), tempCur)   
      alpha(2) = Polynomial(expnCoef(:, 2), tempCur)   
      alpha(3) = Polynomial(expnCoef(:, 3), tempCur)   

      fDfGrdTh = 0.0d0
      do ii = 1, 3
        fDfGrdTh(ii, ii) = 1 + alpha(ii)*(tempCur - tempRef)
      end do
        
    end function fDfGrdTh

  end module expansion

