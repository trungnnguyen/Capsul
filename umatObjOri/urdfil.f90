
subroutine urdfil(lstop, lovrwrt, kstep, kinc, dtime, time)
  use utils, only : RKIND, IKIND

  implicit none

  integer(kind = RKIND), intent(in) :: lstop
  integer(kind = RKIND), intent(in) :: lovrwrt
  integer(kind = RKIND), intent(in) :: kstep
  integer(kind = RKIND), intent(in) :: kinc
  real(kind = RKIND),    intent(in) :: dtime
  real(kind = RKIND),    intent(in) :: time(2)


  real(kind = RKIND) :: array(513)
  character(len=256) :: outFile 


  real(kind = RKIND) :: aMises, force1, force2
  integer(kind = IKIND) :: key, jrcd, k1, jnode

  outFile = '/home/lijf/sigEps.dat' 
  open(unit=32, file=outFile, status='unknown') 

! FIND CURRENT INCREMENT 
  call posfil(kstep, kinc, array, jrcd)
  write(*, *) "KINC = ", KINC
  do k1 = 1, 999999 
    call dbfile(0, array, jrcd) 
    write(*, *) "k1 = ", k1, "array = ", array(1:10)
    if (jrcd /= 0) exit
    KEY = nint(array(2))

    write(32, *) 'KEY=', key
! RECORD 12 CONTAINS VALUES FOR SINV 
    if(key == 12) then
      aMises = array(3) 
      WRITE(32, *) 'MISES STRESS=', AMISES 
    end if

! RECORD 15 CONTAINS VALUES FOR NODAL FORCES 
    if (key == 15 ) then
      JNODE  = nint(array(3))
      FORCE1 = ARRAY(4) 
      FORCE2 = ARRAY(5) 
      WRITE(32, *) 'NODE NUMBER = ', JNODE 
      WRITE(32, *) 'F_1, F_2    = ', FORCE1, FORCE2 
    end if 

  end do

end subroutine
