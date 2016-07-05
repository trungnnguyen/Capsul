  module expansion
    use utils, only : RKIND, IKIND, MAX_POLY_ORDER

    implicit none
    private

    public :: ReadExpansion
    public :: fDfGrdTh

    real(kind = RKIND) :: expnCoef(MAX_POLY_ORDER + 1)
    real(kind = RKIND) :: tempRef
 
  contains

    subroutine ReadExpansion(fID)
      use utils, only : ReadLine

      integer(kind = IKIND), intent(in) :: fID

      ! local variables
      integer(kind = IKIND)  :: nOrder
      character(len = 100)   :: line
      integer(kind = IKIND)  :: ii

      call ReadLine(fID, line)
      read(line, *) tempRef, nOrder
      
      expnCoef = 0.0d0
      do ii = 1, nOrder + 1
        call ReadLine(fID, line)
        read(line, *) expnCoef(ii)
      end do

    end subroutine ReadExpansion


    function fDfGrdTh(tempCur)
      use utils,   only : SMALL, UNITMAT
      use algebra, only : Polynomial

      real(kind = RKIND), intent(in)  :: tempCur
      real(kind = RKIND) :: fDfGrdTh(3, 3)  

      real(kind = RKIND) :: alpha

      alpha = Polynomial(expnCoef, tempCur)   
      fDfGrdTh = (1 + alpha*(tempCur - tempRef)) * UNITMAT
        
    end function fDfGrdTh

  end module expansion

    
