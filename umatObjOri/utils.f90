module utils
  ! definitions of global constants
  implicit none

  integer(kind = 4), parameter :: IKIND = 4  
  integer(kind = 4), parameter :: RKIND = 8
  integer(kind = 4), parameter :: LKIND = 1
  
  character(len = *),    parameter :: CRYS_FILE_NAME   = "crystal.prop"
  integer(kind = IKIND), parameter :: CRYS_FILE_ID     = 74 


  integer(kind = IKIND), parameter :: UPDSIG = 1
  integer(kind = IKIND), parameter :: UPDJAC = 2
 
contains

  subroutine ReadLine(fID, line)

    integer(kind = IKIND), intent(in)  :: fID
    character(len = 256),  intent(out) :: line
    integer(kind = IKIND), parameter :: MAX_LINE_LEN = 256
    character(len = 2) :: filter = "#!"
    integer(kind = IKIND) :: ii

    do while(.TRUE.)
      read(fID, '(a)') line
      line = adjustl(line)
      if (line(1:1) /= "" .and. index(filter, line(1:1)) == 0) exit
!      if (line(1:1) /= " " .and. line(1:2) /= "**") exit
    end do

    do ii = 2, MAX_LINE_LEN
      if (index(filter, line(ii:ii)) /= 0) then
        line(ii:MAX_LINE_LEN) = " "
        exit
      end if
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

