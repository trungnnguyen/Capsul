  module init
    use utils,            only : IKIND, RKIND, LKIND
    use typeCrysPlasMat,  only : CrysPlasMat 
    
    implicit none

    logical(kind = LKIND)           :: initialized = .false.
    class(CrysPlasMat), pointer     :: pCrysPlasMat
    integer(kind = IKIND)           :: numNode
    real(kind = RKIND), allocatable :: initCoords(:, :)

  contains

    subroutine Initialization(ansyType, tempInit) 
      use utils, only : ReadLine

      integer(kind = IKIND), parameter :: kCrysFileID   = 74
      character(len = *),    parameter :: kCrysFileName = "crystal.prop"   
      integer(kind = IKIND), parameter :: kNodeFileID   = 75
      character(len = *),    parameter :: kNodeFileName = "node.inp"   

      integer(kind = IKIND), intent(in) :: ansyType
      real(kind = RKIND),    intent(in) :: tempInit

      logical(kind = LKIND) :: isTempDep
      character(len = 256)  :: outDir, crysFileName, nodeFileName, line
      integer(kind = IKIND) :: lenOutDir
      integer(kind = IKIND) :: iNode
  
      if (ansyType == 0) then
        isTempDep = .false.
      else if (ansyType == 1) then
        isTempDep = .true.
      end if

      call GetOutDir(outDir, lenOutDir)   ! abaqus utility function
      crysFileName = outDir(1:lenOutDir) // "/" // kCrysFileName
      open(unit = kCrysFileID, file = crysFileName, status = 'old')
      pCrysPlasMat => CrysPlasMat(kCrysFileID, isTempDep, tempInit)
      close(unit = kCrysFileID)

      nodeFileName = outDir(1:lenOutDir) // "/" // kNodeFileName
      open(unit = kNodeFileID, file = nodeFileName, status = 'old')
      call readLine(kNodeFileID, line) 
      read(line, *) numNode
      allocate(initCoords(numNode, 4))
      do iNode = 1, numNode
        call readLine(kNodeFileID, line) 
        read(line, *) initCoords(iNode, 1:4)
      end do
      close(unit = kNodeFileID)
  
    end subroutine Initialization
   


    function GetInitCoords(nodeID) result(coords)
      integer(kind = IKIND), intent(in) :: nodeID
      real(kind = RKIND)    :: coords(3)
      integer(kind = IKIND) :: iNode

      do iNode = 1, numNode
        if (initCoords(iNode, 1) == nodeID) then
          coords = initCoords(iNode, 2:4)
          return
        end if
      end do

      write(*, *) "Can not get the initial coordinates of node ", nodeID

    end function

  end module init

