  module init
    use utils, only : IKIND, RKIND, LKIND, MAX_SLIP_SYSTEMS
    
    implicit none

    logical(kind = LKIND) :: initialized = .false.

    ! global variables initialized at the beginning of the run and would not be changed later.
    real(kind = RKIND)    :: oriMatx(3, 3)

    integer(kind = IKIND) :: nSlipSys

    real(kind = RKIND)    :: schmidt(3, 3, MAX_SLIP_SYSTEMS)

    real(kind = RKIND)    :: toler
    real(kind = RKIND)    :: tolerJac
    integer(kind = IKIND) :: maxIters
    integer(kind = IKIND) :: maxItersJac
    real(kind = RKIND)    :: epsInc
 
  contains

    subroutine Initialization(props) 
  
      use utils,     only : CRYS_FILE_NAME, CRYS_FILE_ID, ReadLine
      use algebra,   only : matRot
      use crystal,   only : InitCrysParas
      use expansion, only : ReadExpansion
  
      implicit none
      ! interface arguments
      real(kind = RKIND),    intent(in)  :: props(6)
  
      ! local variables
      real(kind = RKIND)    :: orientMatrix(3, 3)
      real(kind = RKIND)    :: schmidtLoc(3, 3, MAX_SLIP_SYSTEMS)

      character(len = 256)   :: outDir, fileName
      integer(kind = IKIND) :: lenOutDir
      character(len = 256)  :: line
      integer(kind = IKIND) :: ii
  
  
      call GetOutDir(outDir, lenOutDir)   ! abaqus utility function
      fileName = outDir(1:lenOutDir) // "/" // CRYS_FILE_NAME
      open(unit = CRYS_FILE_ID, file = fileName, status = 'old')
  
      call InitCrysParas(CRYS_FILE_ID, nSlipSys, schmidtLoc)

      oriMatx = fOriMatx(props(1:3), props(4:6))

      do ii = 1, nSlipSys
        schmidt(:, :, ii) = matRot(schmidtLoc(:, :, ii), transpose(oriMatx))
        write(*, *) "schmidt ", ii, " = ", schmidt(:, :, ii)
      end do
  

      call ReadExpansion(CRYS_FILE_ID)

      call ReadLine(CRYS_FILE_ID, line)
      read(line, *) toler, tolerJac, maxIters, maxItersJac, epsInc

  
      close(unit = CRYS_FILE_ID)
  
    end subroutine Initialization
  


    function fOriMatx(oriVec1, oriVec2)
      use algebra, only : normalize, cross_product

      real(kind = RKIND), intent(in)   :: oriVec1(3)
      real(kind = RKIND), intent(in)   :: oriVec2(3)
      real(kind = RKIND)  :: fOriMatx(3, 3)


      fOriMatx(1:3, 1) = oriVec1
      fOriMatx(1:3, 2) = oriVec2

      call normalize(fOriMatx(1:3, 1))
      call normalize(fOriMatx(1:3, 2))

      fOriMatx(1:3, 3) = cross_product(fOriMatx(1:3, 1), fOriMatx(1:3, 2))

    end function fOriMatx

 
  end module init

