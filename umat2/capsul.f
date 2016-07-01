!> main entry of the umat called by abaqus 
  subroutine umat(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,                            &  
                  RPL, DDSDDT, DRPLDE, DRPLDT,                                      &
                  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,   &
                  NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,    &
                  CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

    use utils, only : RKIND, IKIND
    implicit none

    ! Explicit Declaration of arguments called by ABAQUS
    include 'umatArgs.inc'

    if (CMNAME(1:8) == "UMAT_CP1") then
      call umat_cp1(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,                            &  
                    RPL, DDSDDT, DRPLDE, DRPLDT,                                      &
                    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,   &
                    NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,    &
                    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
    else if (CMNAME(1:8) == "UMAT_CP2") then
      call umat_cp2(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,                            &  
                    RPL, DDSDDT, DRPLDE, DRPLDT,                                      &
                    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,   &
                    NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,    &
                    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)


    else
      write(*, *) "Unknown Material Type: ", CMNAME
    end if

!    drpldt = 1.0d1
!    drplde = 2.0d2
!    ddsddt = 3.0d3

    if (npt == 1) then
!      write(*, *) "In umat   incr = ", Kinc, " temp = ", temp, " dtemp = ", dtemp
    end if

  end subroutine umat



  subroutine umatht(u, dudt, dudg, flux, dfdt, dfdg,                           &
                    statev, temp, dtemp, dtemdx, time, dtime, predef, dpred,   &
                    cmname, ntgrd, nstatv, props, nprops, coords, pnewdt,      &
                    noel, npt, layer, kspt, kstep, kinc)

    use utils
    implicit none

    ! Declaration of interface variables
    real(kind = RKIND)    :: U
    real(kind = RKIND)    :: dudt
    real(kind = RKIND)    :: DTIME
    real(kind = RKIND)    :: TEMP
    real(kind = RKIND)    :: DTEMP
    integer(kind = IKIND) :: ntgrd
    integer(kind = IKIND) :: NSTATV
    integer(kind = IKIND) :: NPROPS
    real(kind = RKIND)    :: PNEWDT
    integer(kind = IKIND) :: NOEL
    integer(kind = IKIND) :: NPT
    integer(kind = IKIND) :: LAYER
    integer(kind = IKIND) :: KSPT
    integer(kind = IKIND) :: KSTEP
    integer(kind = IKIND) :: KINC
    character(len = 80)   :: CMNAME

    real(kind = RKIND)    :: dudg(ntgrd, ntgrd)
    real(kind = RKIND)    :: flux(ntgrd)
    real(kind = RKIND)    :: dfdt(ntgrd)
    real(kind = RKIND)    :: dfdg(ntgrd, ntgrd)
    real(kind = RKIND)    :: STATEV(NSTATV)
    real(kind = RKIND)    :: dtemdx(ntgrd)
    real(kind = RKIND)    :: time(2)
    real(kind = RKIND)    :: PREDEF(1)
    real(kind = RKIND)    :: DPRED(1)
    real(kind = RKIND)    :: PROPS(NPROPS)
    real(kind = RKIND)    :: COORDS(3)

    ! local variables
    real(kind = RKIND)    :: cond 
    real(kind = RKIND)    :: specht
    real(kind = RKIND)    :: du

    integer(kind = IKIND) :: i


    if (npt == 1) then
!      write(*, *) "In umatht incr = ", Kinc, " temp = ", temp, " dtemp = ", dtemp
    end if

    cond   = props(1)
    specht = props(2)
    
    dudt = specht
    du   = dudt*dtemp
    u    = u + du

    dfdg = 0.0d0
    do i = 1, ntgrd
      flux(i) = -cond*dtemdx(i)
      dfdg(i, i) = -cond
    end do


  end subroutine


  subroutine uexpan(expan, dexpandt, temp, time, dtime, predef, dpred, statev, cmname, nstatv, noel)
    use utils, only : RKIND, IKIND
    implicit none

    real(kind = RKIND) :: expan(:)
    real(kind = RKIND) :: dexpandt(:)
    real(kind = RKIND) :: temp(2)
    real(kind = RKIND) :: time(2)
    real(kind = RKIND) :: dtime
    real(kind = RKIND) :: predef(:)
    real(kind = RKIND) :: dpred(:)
    integer(kind = IKIND) :: nstatv
    real(kind = RKIND) :: statev(nstatv)
    integer(kind = IKIND) :: noel
    character(len = 80)   :: cmname


  end subroutine uexpan


