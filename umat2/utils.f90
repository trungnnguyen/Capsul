module utils
  ! definitions of global constants
  implicit none

  integer(kind = 4), parameter :: IKIND = 4  
  integer(kind = 4), parameter :: RKIND = 8
  integer(kind = 4), parameter :: LKIND = 1
  
  real(kind = RKIND), parameter :: PI    = 3.1415926535897932d0
  real(kind = RKIND), parameter :: SMALL = 1.0d-12

  real(kind = RKIND), parameter :: UNITMAT(3, 3)        &
      = reshape((/1.0d0, 0.0d0, 0.0d0,                  &
                  0.0d0, 1.0d0, 0.0d0,                  &
                  0.0d0, 0.0d0, 1.0d0/), (/3, 3/))

  character(len = *),    parameter :: CRYS_FILE_NAME   = "crystal.prop"
  integer(kind = IKIND), parameter :: CRYS_FILE_ID     = 74 


  integer(kind = IKIND), parameter :: MAX_SLIP_SYSTEMS = 30
  integer(kind = IKIND), parameter :: MAX_POLY_ORDER   = 5

  integer(kind = IKIND), parameter :: UPDSIG = 1
  integer(kind = IKIND), parameter :: UPDJAC = 2
 
contains

  subroutine ReadLine(fID, line)

    integer(kind = IKIND), intent(in)  :: fID
    character(len = 100),  intent(out) :: line

    character(len = 2) :: filter = "#!"
    character(len = 6) :: strFmt = "(A256)"

    do while(.TRUE.)
      read(fID, strFmt) line
      line = adjustl(line)
      if (line(1:1) /= " " .and. index(filter, line(1:1)) == 0) exit
!      if (line(1:1) /= " " .and. line(1:2) /= "**") exit
    end do

  end subroutine ReadLine


  subroutine ToLower(str)

    character(len = *) :: str

    integer(kind = IKIND) :: lenStr
    integer(kind = IKIND) :: ii, iCode

    lenStr = len(str)
    do ii = 1, lenStr
      iCode = iachar(str(ii:ii))
      if (iCode >= 65 .and. iCode <= 90) then
        str(ii:ii) = achar(iCode + 32)
      end if
    end do

  end subroutine ToLower 

end module utils

