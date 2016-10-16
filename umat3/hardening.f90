  module hardening
    use utils, only : RKIND, IKIND, MAX_SLIP_SYSTEMS

    implicit none

    private
    public InitHardParas, fTauCritInit, fTauCritDot, fDTauCritDotDGammaDot
   
  
    integer(kind = IKIND), parameter :: MAX_SLIP_SETS = 6

    integer(kind = IKIND) :: nSlipSys
    character(len = 16)   :: hardType
    real(kind = RKIND)    :: hardParas(MAX_SLIP_SYSTEMS, MAX_SLIP_SYSTEMS + 5) 
    real(kind = RKIND)    :: workArray(MAX_SLIP_SYSTEMS)

  contains

  
    subroutine InitHardParas(fID, nSlipSets, nSlipSystems, systemSet)
      use utils, only : ReadLine, toLower
  
      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: nSlipSets
      integer(kind = IKIND), intent(in) :: nSlipSystems
      integer(kind = IKIND), intent(in) :: systemSet(nSlipSystems)

      character(len = 256)  :: line
      integer(kind = IKIND) :: ii
      

      nSlipSys  = nSlipSystems
      hardParas = 0.0d0
      do ii = 1, nSlipSys
        hardParas(ii, ii) = 1.0d0
      end do

      call ReadLine(fID, line)
      read(line, *) hardType
      call toLower(hardType)

      if (hardType(1:8) == "assaneed") then
        call InitANHardParas(fID, nSlipSets, systemSet)
      else if (hardType(1:8) == "vocetome") then
        call InitVTHardParas(fID, nSlipSets, systemSet)
      else if (hardType(1:8) == "vocekock") then
        call InitVKHardParas(fID, nSlipSets, systemSet)
      end if

    end subroutine InitHardParas

    

    subroutine InitANHardParas(fID, nSlipSets, systemSet)
      use utils, only : ReadLine

      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: nSlipSets
      integer(kind = IKIND), intent(in) :: systemSet(nSlipSys)

      character(len = 256) :: line

      real(kind = RKIND) :: qsys(MAX_SLIP_SETS, MAX_SLIP_SETS)
      real(kind = RKIND) :: tau0_set(MAX_SLIP_SETS)
      real(kind = RKIND) :: taus_set(MAX_SLIP_SETS)
      real(kind = RKIND) :: h0_set(MAX_SLIP_SETS)
      real(kind = RKIND) :: tau0(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: taus(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: h0(MAX_SLIP_SYSTEMS)

      integer(kind = IKIND) :: ii, jj, iSet1, iSet2


      do ii = 1, nSlipSets
        call readLine(fID, line)
        read(line, *) (qsys(ii, jj), jj = 1, nSlipSets)
      end do
   
      do ii = 1, nSlipSys
        do jj = 1, nSlipSys
          if(ii /= jj) then
            iSet1 = systemSet(ii)
            iSet2 = systemSet(jj)
            hardParas(ii, jj) = qsys(iSet1, iSet2)
          end if
        end do
      end do
      
      do ii = 1, nSlipSets
        call readLine(fID, line)
        read(line, *) tau0_set(ii), taus_set(ii), h0_set(ii)
      end do

      do ii = 1, nSlipSys
        iSet1 = systemSet(ii)
        tau0(ii) = tau0_set(iSet1)
        taus(ii) = taus_set(iSet1)
        h0       = h0_set(iSet1)
      end do

      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 1) = tau0(1:nSlipSys)
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 2) = taus(1:nSlipSys)
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 3) = h0(1:nSlipSys)

      workArray(1:nSlipSys) = h0(1:nSlipSys) / (taus(1:nSlipSys) - tau0(1:nSlipSys))
     
    end subroutine InitANHardParas



    subroutine InitVTHardParas(fID, nSlipSets, systemSet)
      use utils, only : ReadLine

      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: nSlipSets
      integer(kind = IKIND), intent(in) :: systemSet(nSlipSys)

      character(len = 256) :: line

      real(kind = RKIND) :: qsys(MAX_SLIP_SETS, MAX_SLIP_SETS)
      real(kind = RKIND) :: tau0_set(MAX_SLIP_SETS)
      real(kind = RKIND) :: taus_set(MAX_SLIP_SETS)
      real(kind = RKIND) :: h0_set(MAX_SLIP_SETS)
      real(kind = RKIND) :: h1_set(MAX_SLIP_SETS)
      real(kind = RKIND) :: tau0(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: taus(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: h0(MAX_SLIP_SYSTEMS)
      real(kind = RKIND) :: h1(MAX_SLIP_SYSTEMS)
      integer(kind = IKIND) :: ii, jj, iSet1, iSet2


      do ii = 1, nSlipSets
        call readLine(fID, line)
        read(line, *) (qsys(ii, jj), jj = 1, nSlipSets)
      end do

      do ii = 1, nSlipSys
        do jj = 1, nSlipSys
          if(ii /= jj) then
            iSet1 = systemSet(ii)
            iSet2 = systemSet(jj)
            hardParas(ii, jj) = qsys(iSet1, iSet2)
          end if
        end do
      end do

      do ii = 1, nSlipSets
        call readLine(fID, line)
        read(line, *) tau0_set(ii), taus_set(ii), h0_set(ii), h1_set(ii)
      end do

      if (any(taus_set < tau0_set)) then
        write(*, *) 'ERROR IN HARDENING', taus_set(1:nSlipSys)
        stop
      end if

      do ii = 1, nSlipSys
        iSet1    = systemSet(ii)
        tau0(ii) = tau0_set(iSet1)
        taus(ii) = taus_set(iSet1)
        h0       = h0_set(iSet1)
        h1       = h1_set(iSet1)
      end do

      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 1) = tau0(1:nSlipSys)
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 2) = taus(1:nSlipSys) - tau0(1:nSlipSys)
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 3) = h0(1:nSlipSys)
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 4) = h1(1:nSlipSys)

      workArray(1:nSlipSys) = h0(1:nSlipSys) / taus(1:nSlipSys)

    end subroutine InitVTHardParas



    subroutine InitVKHardParas(fID, nSlipSets, systemSet)
      use utils, only : ReadLine

      integer(kind = IKIND), intent(in) :: fID
      integer(kind = IKIND), intent(in) :: nSlipSets
      integer(kind = IKIND), intent(in) :: systemSet(nSlipSys)

      real(kind = RKIND)   :: h0, tau0, taus0, gammas0, ma
      character(len = 256) :: line

     
      call readLine(fID, line)
      read(line, *)  tau0, taus0, h0, gammas0, ma
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 1) = tau0
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 2) = taus0
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 3) = h0
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 4) = gammas0
      hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 5) = ma

    end subroutine InitVKHardParas



    function fTauCritInit()
      real(kind = RKIND) :: fTauCritInit(nSlipSys)
      real(kind = RKIND) :: tau0(nSlipSys)

      tau0 = hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 1)
      fTauCritInit = tau0

    end function fTauCritInit


    function fTauCritDot(gamma, gammaDot, tauCrit, temp)
      
      real(kind = RKIND), intent(in)  :: gamma(nSlipSys)
      real(kind = RKIND), intent(in)  :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in)  :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in)  :: temp

      real(kind = RKIND) :: fTauCritDot(nSlipSys)

      real(kind = RKIND) :: hh(nSlipSys, nSlipSys)
      real(kind = RKIND) :: absGammaDot(nSlipSys)
      integer(kind = IKIND) :: ii


      hh          = fHH(gamma, gammaDot, tauCrit, temp)
      absGammaDot = abs(gammaDot)
      do ii = 1, nSlipSys
        fTauCritDot(ii) = dot_product(hh(ii, :), absGammaDot)
      end do

    end function fTauCritDot 



    function fHH(gamma, gammaDot, tauCrit, temp)

      real(kind = RKIND), intent(in)  :: gamma(nSlipSys)
      real(kind = RKIND), intent(in)  :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in)  :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in)  :: temp

      real(kind = RKIND) :: fHH(nSlipSys, nSlipSys)

      if (hardType(1:8) == "assaneed") then
        fHH = fHHAN(gamma, gammaDot, temp)
      else if (hardType(1:8) == "vocetome") then
        fHH = fHHVT(gamma, gammaDot, temp)
      else if (hardType(1:8) == "vocekock") then
        fHH = fHHVK(gammaDot, tauCrit, temp)
      end if

    end function fHH


    
    function fHHAN(gamma, gammaDot, temp)
      use algebra, only : sech

      real(kind = RKIND), intent(in) :: gamma(nSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fHHAN(nSlipSys, nSlipSys)

      real(kind = RKIND) :: qq(nSlipSys, nSlipSys)
      real(kind = RKIND) :: h0(nSlipSys)

      real(kind = RKIND) :: aux(nSlipSys)
      real(kind = RKIND) :: hii
      integer(kind = IKIND) :: ii, jj

      qq  = hardParas(1:nSlipSys, 1:nSlipSys)
      h0  = hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 3)
      aux = sech(workArray(1:nSlipSys) * sum(gamma))
      do ii = 1, nSlipSys
        hii = h0(ii)*aux(ii)*aux(ii) 
        do jj = 1, nSlipSys
          fHHAN(ii, jj) = qq(ii, jj)*hii
        end do
      end do

    end function fHHAN



    function fHHVT(gamma, gammaDot, temp)
      
      real(kind = RKIND), intent(in) :: gamma(nSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fHHVT(nSlipSys, nSlipSys)

      real(kind = RKIND) :: qq(nSlipSys, nSlipSys)
      real(kind = RKIND) :: taus(nSlipSys)
      real(kind = RKIND) :: h1(nSlipSys)
      real(kind = RKIND) :: gammaTot
      real(kind = RKIND) :: aux(nSlipSys)
      real(kind = RKIND) :: hii

      integer(kind = IKIND) :: ii, jj


      qq   = hardParas(1:nSlipSys, 1:nSlipSys)
      taus = hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 2)
      h1   = hardParas(1:nSlipSys, MAX_SLIP_SYSTEMS + 4)

      gammaTot = sum(gamma)
      aux = dexp(workArray(ii)*(-gammaTot))
      do ii = 1, nSlipSys
        hii = h1(ii)*(1.0d0 - aux(ii)) + (taus(ii) + h1(ii)*gammaTot) * workArray(ii)*aux(ii)
        do jj = 1, nSlipSys
          fHHVT(ii, jj) = qq(ii, jj)*hii
        end do
      end do

    end function fHHVT


    function fHHVK(gammaDot, tauCrit, temp)

      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fHHVK(nSlipSys, nSlipSys)

      real(kind = RKIND) :: tau0  
      real(kind = RKIND) :: taus0  
      real(kind = RKIND) :: h0  
      real(kind = RKIND) :: gammas0  
      real(kind = RKIND) :: ma  

      real(kind = RKIND) :: gammaDotSum
      real(kind = RKIND) :: taus
      real(kind = RKIND) :: aux
      real(kind = RKIND) :: mm
      integer(kind = IKIND) :: ii


      tau0    = hardParas(1, MAX_SLIP_SYSTEMS + 1)
      taus0   = hardParas(1, MAX_SLIP_SYSTEMS + 2)
      h0      = hardParas(1, MAX_SLIP_SYSTEMS + 3)
      gammas0 = hardParas(1, MAX_SLIP_SYSTEMS + 4)
      ma      = hardParas(1, MAX_SLIP_SYSTEMS + 5)

      mm = ma

      gammaDotSum = sum(abs(gammaDot))
      taus = taus0 * (gammaDotSum/gammas0) ** (mm)
      aux = h0/(taus - tau0)
      do ii = 1, nSlipSys
        fHHVK(1:nSlipSys, ii) = aux*(taus - tauCrit(ii))
      end do

    end function fHHVK

   

    !> \[ \frac{\partial \dot{\tau_i}}{\partial \dot{\gamma_j}} \]
    function fDTauCritDotDGammaDot(gamma, gammaDot, tauCrit, temp)
      real(kind = RKIND), intent(in)  :: gamma(nSlipSys)
      real(kind = RKIND), intent(in)  :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in)  :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in)  :: temp

      real(kind = RKIND) :: fDTauCritDotDGammaDot(nSlipSys, nSlipSys)
      real(kind = RKIND) :: gammaDotSign(nSlipSys)
      real(kind = RKIND) :: absGammaDot(nSlipSys)
      real(kind = RKIND) :: hh(nSlipSys, nSlipSys)
      real(kind = RKIND) :: gg(nSlipSys, nSlipSys)
      integer(kind = RKIND) :: ii, jj

      gammaDotSign = sign(1.0d0, gammaDot)
      absGammaDot  = abs(gammaDot)
      hh = fHH(gamma, gammaDot, tauCrit, temp)
      gg = fGG(gamma, gammaDot, tauCrit, temp)

      do ii = 1, nSlipSys
        do jj = 1, nSlipSys
          fDTauCritDotDGammaDot(ii, jj) = hh(ii, jj)*gammaDotSign(jj) + gg(ii, jj)*absGammaDot(jj)
        end do
      end do

    end function fDTauCritDotDGammaDot



    function fGG(gamma, gammaDot, tauCrit, temp)

      real(kind = RKIND), intent(in)  :: gamma(nSlipSys)
      real(kind = RKIND), intent(in)  :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in)  :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in)  :: temp

      real(kind = RKIND) :: fGG(nSlipSys, nSlipSys)

      if (hardType(1:8) == "assaneed") then
        fGG = fGGAN(gamma, gammaDot, temp)
      else if (hardType(1:8) == "vocetome") then
        fGG = fGGVT(gamma, gammaDot, temp)
      else if (hardType(1:8) == "vocekock") then
        fGG = fGGVK(gammaDot, tauCrit, temp)
      end if

    end function fGG


      
    function fGGAN(gamma, gammaDot, temp)
      
      real(kind = RKIND), intent(in) :: gamma(nSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND)  :: fGGAN(nSlipSys, nSlipSys)


      fGGAN = 0.0d0

    end function fGGAN



    function fGGVT(gamma, gammaDot, temp)
      
      real(kind = RKIND), intent(in) :: gamma(nSlipSys)
      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys, 1)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fGGVT(nSlipSys, nSlipSys)


      fGGVT = 0.0d0

    end function fGGVT

  

    function fGGVK(gammaDot, tauCrit, temp)

      real(kind = RKIND), intent(in) :: gammaDot(nSlipSys)
      real(kind = RKIND), intent(in) :: tauCrit(nSlipSys)
      real(kind = RKIND), intent(in) :: temp

      real(kind = RKIND) :: fGGVK(nSlipSys, nSlipSys) ! dTauCrit(i)/dGamma(j)

      real(kind = RKIND) :: tau0, taus0, h0, gammas0, ma
      real(kind = RKIND) :: gammaDotSign(nSlipSys)
      real(kind = RKIND) :: mm, taus, gammaDotSum
      real(kind = RKIND) :: aux1, aux2, aux3 
      integer(kind = IKIND) :: i, j


      tau0    = hardParas(1, MAX_SLIP_SYSTEMS + 1)
      taus0   = hardParas(1, MAX_SLIP_SYSTEMS + 2)
      h0      = hardParas(1, MAX_SLIP_SYSTEMS + 3)
      gammas0 = hardParas(1, MAX_SLIP_SYSTEMS + 4)
      ma      = hardParas(1, MAX_SLIP_SYSTEMS + 5)

!      mm = temp/ma
      mm = ma
      
      gammaDotSign = sign(1.0d0, gammaDot)
      gammaDotSum = sum(abs(gammaDot))
      if (gammaDotSum < 1.0d-12) then
        fGGVK = 0.0d0
        return 
      end if

      taus = taus0 * (gammaDotSum/gammas0) ** mm
      aux1 = taus - tau0
      !aux2 = -h0*mm*taus/(aux1 * aux1 * gammaDotSum)
      aux2 = h0*mm*taus0*(gammaDotSum/gammas0)**(mm-1)/(aux1 * aux1 * gammas0)
      do i = 1, nSlipSys
        aux3 = aux2*(tauCrit(i) - tau0)
        do j = 1, nSlipSys
          fGGVK(i, j) = aux3*gammaDotSign(j)
        end do
      end do

    end function fGGVK


  end module hardening


