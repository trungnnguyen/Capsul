
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


