
    ! Declaration of interface variables
    real(kind = RKIND)    :: SSE
    real(kind = RKIND)    :: SPD
    real(kind = RKIND)    :: SCD
    real(kind = RKIND)    :: RPL
    real(kind = RKIND)    :: DRPLDT
    real(kind = RKIND)    :: DTIME
    real(kind = RKIND)    :: TEMP
    real(kind = RKIND)    :: DTEMP
    integer(kind = IKIND) :: NDI
    integer(kind = IKIND) :: NSHR
    integer(kind = IKIND) :: NTENS
    integer(kind = IKIND) :: NSTATV
    integer(kind = IKIND) :: NPROPS
    real(kind = RKIND)    :: PNEWDT
    integer(kind = IKIND) :: NOEL
    integer(kind = IKIND) :: NPT
    integer(kind = IKIND) :: LAYER
    integer(kind = IKIND) :: KSPT
    integer(kind = IKIND) :: KSTEP
    integer(kind = IKIND) :: KINC
    character(len = *)    :: CMNAME
    character(len = *)    :: CELENT

    real(kind = RKIND)    :: STRESS(NTENS)
    real(kind = RKIND)    :: STATEV(NSTATV)
    real(kind = RKIND)    :: DDSDDE(NTENS, NTENS)
    real(kind = RKIND)    :: DDSDDT(NTENS)
    real(kind = RKIND)    :: DRPLDE(NTENS)
    real(kind = RKIND)    :: STRAN(NTENS)
    real(kind = RKIND)    :: DSTRAN(NTENS)
    real(kind = RKIND)    :: TIME(2)
    real(kind = RKIND)    :: PREDEF(1)
    real(kind = RKIND)    :: DPRED(1)
    real(kind = RKIND)    :: PROPS(NPROPS)
    real(kind = RKIND)    :: COORDS(3)
    real(kind = RKIND)    :: DROT(3, 3)
    real(kind = RKIND)    :: DFGRD0(3, 3)
    real(kind = RKIND)    :: DFGRD1(3, 3)

