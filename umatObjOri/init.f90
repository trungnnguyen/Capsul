  module init
    use utils,             only : IKIND, RKIND, LKIND
    use typeCrysPlasMat,   only : CrysPlasMat 
!    use typeCrysPlasMatPt, only : CrysPlasMatPt 
    
    implicit none

    logical(kind = LKIND) :: initialized = .false.
    class(CrysPlasMat),   pointer :: pCrysPlasMat   => null()
!    class(CrysPlasMatPt), pointer :: pCrysPlasMatPt => null() 

  contains

    subroutine Initialization(ansyType, tempInit) 
      use utils, only : CRYS_FILE_NAME, CRYS_FILE_ID

      integer(kind = IKIND), intent(in) :: ansyType
      real(kind = RKIND),    intent(in) :: tempInit

      logical(kind = LKIND) :: isTempDep
      character(len = 256)  :: outDir, fileName
      integer(kind = IKIND) :: lenOutDir
  
      if (ansyType == 0) then
        isTempDep = .false.
      else if (ansyType == 1) then
        isTempDep = .true.
      end if

      call GetOutDir(outDir, lenOutDir)   ! abaqus utility function
      fileName = outDir(1:lenOutDir) // "/" // CRYS_FILE_NAME
      open(unit = CRYS_FILE_ID, file = fileName, status = 'old')
      pCrysPlasMat   => CrysPlasMat(CRYS_FILE_ID, isTempDep, tempInit)
!      pCrysPlasMatPt => CrysPlasMatPt(pCrysPlasMat)
      close(unit = CRYS_FILE_ID)
  
    end subroutine Initialization
   
  end module init

