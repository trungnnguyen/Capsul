  module crystal
    use utils, only : RKIND, IKIND, LKIND, MAX_POLY_ORDER

    implicit none
    private

    public  :: InitCrysParas
    public  :: fStifLoc
    public  :: fGammaDot
    public  :: fGammaDotDDTauCrit
    public  :: fStressDiverged
    public  :: fTauResl
    public  :: fLp

    integer(kind = IKIND), parameter :: FCC_SLIP_VECTOR(6, 12) =    &
            reshape((/ 1,  1,  1,  1, -1,  0,                       &
                       1,  1,  1,  0, -1,  1,                       &
                       1,  1,  1, -1,  0,  1,                       &
                       1, -1,  1,  1,  1,  0,                       &
                       1, -1,  1,  0,  1,  1,                       &
                       1, -1,  1, -1,  0,  1,                       &
                      -1, -1,  1,  1, -1,  0,                       &
                      -1, -1,  1,  0,  1,  1,                       &
                      -1, -1,  1,  1,  0,  1,                       &
                      -1,  1,  1,  1,  1,  0,                       &
                      -1,  1,  1,  0, -1,  1,                       &
                      -1,  1,  1,  1,  0,  1/), (/6, 12/))


    integer(kind = IKIND) :: nSlipSys
    real(kind = RKIND)    :: elasProps(MAX_POLY_ORDER + 1, 6)
    real(kind = RKIND)    :: gammaDot0, mm, qActive, mmInv, mmBig, gammaDot0MMInv
 
  contains


    subroutine InitCrysParas(fID, nSlipSystems, schmidtLoc)
      use utils,     only : MAX_SLIP_SYSTEMS, ReadLine, ToLower
      use algebra,   only : out_product, normalize
      use hardening, only : InitHardParas

      integer(kind = IKIND), intent(in)   :: fID
      integer(kind = IKIND), intent(out)  :: nSlipSystems
      real(kind = RKIND),    intent(out)  :: schmidtLoc(3, 3, MAX_SLIP_SYSTEMS)

      ! local variables
      integer(kind = IKIND) :: nSlipSets
      integer(kind = IKIND) :: systemSet(MAX_SLIP_SYSTEMS)

      real(kind = RKIND)    :: slipVec(3, MAX_SLIP_SYSTEMS)
      real(kind = RKIND)    :: normVec(3, MAX_SLIP_SYSTEMS)

      real(kind = RKIND)    :: c11, c12, c44, c33, c13, c66

      integer(kind = IKIND)  :: nOrder
      character(len = 100)   :: line
      character(len = 10)    :: crysType
      integer(kind = IKIND)  :: ii, istat

      !--------------------------------------------------------------
      call ReadLine(fID, line)
      read(line, *) crysType
      call ToLower(crysType)

      if (crysType(1:3) == "fcc") then
        nSlipSets = 1
        nSlipSys  = 12
        systemSet(1:nSlipSys) = 1

        normVec(1:3, 1:nSlipSys) = FCC_SLIP_VECTOR(1:3, 1:nSlipSys)
        slipVec(1:3, 1:nSlipSys) = FCC_SLIP_VECTOR(4:6, 1:nSlipSys)
      end if

      nSlipSystems = nSlipSys

      do ii = 1, nSlipSys
        call normalize(slipVec(:, ii))
        call normalize(normVec(:, ii))
        schmidtLoc(:, :, ii) = out_product(slipVec(:, ii), normVec(:, ii))
      end do

      call ReadLine(fID, line)
      read(line, *) nOrder  ! order of polynomials for elastic properties as functions of temperature 
      elasProps = 0.0d0
      do ii = 1, nOrder + 1
        call ReadLine(fID, line)
        read(line, *, iostat=istat) c11, c12, c44, c13, c33, c66
        if (istat == 0) then
          c13 = 0.0d0
          c33 = 0.0d0
          c66 = 0.0d0
        end if
        elasProps(ii, 1:6) = (/c11, c12, c44, c13, c33, c66/)
      end do
      

      call ReadLine(fID, line)
      read(line, *)  gammaDot0, mm, qActive 
      
      mmInv = 1.0d0/mm
      mmBig = 10.0d0**(mm*150.d0)
      gammaDot0MMInv = gammaDot0*mmInv

      call InitHardParas(fID, nSlipSets, nSlipSys, systemSet(1:nSlipSys))

    end subroutine InitCrysParas




    function  fStifLoc(temp)
      use utils,   only : SMALL
      use algebra, only : Polynomial

      real(kind = RKIND), intent(in)  :: temp    ! current temperature
      real(kind = RKIND) :: fStifLoc(3, 3, 3, 3)  

      real(kind = RKIND) :: c11, c12, c44, c33, c13, c66


      c11 = Polynomial(elasProps(:, 1), temp)   
      c12 = Polynomial(elasProps(:, 2), temp)   
      c44 = Polynomial(elasProps(:, 3), temp)   
      c33 = Polynomial(elasProps(:, 4), temp)   
      c13 = Polynomial(elasProps(:, 5), temp)   
      c66 = Polynomial(elasProps(:, 6), temp)   
        
      fStifLoc = 0.0d0

      fStifLoc(1, 1, 1, 1) = c11
      fStifLoc(2, 2, 2, 2) = c11
     
      fStifLoc(1, 1, 2, 2) = c12
      fStifLoc(2, 2, 1, 1) = c12
    
      fStifLoc(1, 3, 1, 3) = c44
      fStifLoc(3, 1, 3, 1) = c44
      fStifLoc(1, 3, 3, 1) = c44
      fStifLoc(3, 1, 1, 3) = c44
     
      fStifLoc(2, 3, 2, 3) = c44
      fStifLoc(3, 2, 3, 2) = c44
      fStifLoc(2, 3, 3, 2) = c44
      fStifLoc(3, 2, 2, 3) = c44

      if (abs(c33) < SMALL) then
        fStifLoc(3, 3, 3, 3) = c11

        fStifLoc(1, 1, 3, 3) = c12
        fStifLoc(3, 3, 1, 1) = c12
        fStifLoc(2, 2, 3, 3) = c12
        fStifLoc(3, 3, 2, 2) = c12

        fStifLoc(1, 2, 1, 2) = c44
        fStifLoc(2, 1, 2, 1) = c44
        fStifLoc(1, 2, 2, 1) = c44
        fStifLoc(2, 1, 1, 2) = c44
      else
        fStifLoc(3, 3, 3, 3) = c33

        fStifLoc(1, 1, 3, 3) = c13
        fStifLoc(3, 3, 1, 1) = c13
        fStifLoc(2, 2, 3, 3) = c13
        fStifLoc(3, 3, 2, 2) = c13

        fStifLoc(1, 2, 1, 2) = c66
        fStifLoc(2, 1, 2, 1) = c66
        fStifLoc(1, 2, 2, 1) = c66
        fStifLoc(2, 1, 1, 2) = c66
      end if

    end function fStifLoc


    

    function fGammaDot(tauRatio, tauSign, temp)
      use utils, only : SMALL

      real(kind = RKIND), intent(in) :: tauRatio(nSlipSys)
      real(kind = RKIND), intent(in) :: tauSign(nSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: fGammaDot(nSlipSys)

      real(kind = RKIND) :: aux
      integer :: i

    
      if (temp > SMALL) then
        aux = exp(-qActive/temp)
      else
        aux = 1.0d0
      end if

      do i = 1, nSlipSys
        fGammaDot(i) = gammaDot0*tauSign(i)*tauRatio(i)**mmInv*aux
      end do
     
    end function fGammaDot


    ! derivatives of gammaDot with respect to tauCrit in each slip systems
    function fGammaDotDDTauCrit(tauRatio, tauCrit, temp)
      
      real(kind = RKIND), intent(in) :: tauRatio(nSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in) :: temp
      real(kind = RKIND) :: fGammaDotDDTauCrit(nSlipSys)

      integer :: i

      do i = 1, nSlipSys
        fGammaDotDDTauCrit(i) = gammaDot0MMInv*tauRatio(i)**(MMInv - 1)/tauCrit(i)
      end do

    end function fGammaDotDDTauCrit 



    function fStressDiverged(tauRatio)
      real(kind = RKIND), intent(in) :: tauRatio(nSlipSys)
      logical(kind = LKIND) :: fStressDiverged

      fStressDiverged = any(tauRatio > mmBig)

    end function fStressDiverged



    function fTauResl(KirchStress, schmidt)

      real(kind = RKIND), intent(in) :: kirchStress(3, 3)
      real(kind = RKIND), intent(in) :: schmidt(3, 3, nSlipSys)
      real(kind = RKIND) :: fTauResl(nSlipSys)

      integer(kind = IKIND) :: i

      do i = 1, nSlipSys
        fTauResl(i) = sum(kirchStress*schmidt(:, :, i)) 
      end do

    end function fTauResl


    ! Calculate the plastic velocity gradients Lp with the shear strain rate
    ! gammaDot
    function fLp(gammaDot, schmidt)

      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in) :: schmidt(3, 3, nSlipSys)
      real(kind = RKIND) :: fLp(3, 3)
      
      integer(kind = IKIND) :: i

      fLp = 0.0d0
      do i = 1, nSlipSys
        fLp = fLp + gammaDot(i) * schmidt(:, :, i)
      end do

    end function fLp

   
  end module crystal


