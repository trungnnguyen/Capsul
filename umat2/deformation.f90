
  module deformation
    use utils

    implicit none

  contains


    function fGreenStrain(Ftot)
      implicit none

      real(kind = RKIND), intent(in) :: Ftot(3, 3)
      real(kind = RKIND) :: fGreenStrain(3, 3)

      real(kind = RKIND) :: CC(3, 3)


      CC  = matmul(transpose(Ftot), Ftot)
      fGreenStrain = 0.5d0*(CC - UNITMAT) 

    end function fGreenStrain
    
    
 ! Calculate the skirchoff_rot stress: \f$\sigma_{ij}=C_{ijkl}\epsilon_{kl}\f$
    function fElasStress(stiff, strain)

      real(kind = RKIND), intent(in)  :: strain(3, 3)
      real(kind = RKIND), intent(in)  :: stiff(3, 3, 3, 3)
      real(kind = RKIND) :: fElasStress(3, 3)

      integer(kind = IKIND) :: i, j, k, l

      do i = 1, 3
        do j = 1, 3
          fElasStress(i, j) = 0.0d0
          do k = 1, 3
            do l = 1, 3
              fElasStress(i, j) = fElasStress(i, j) + stiff(i, j, k, l)*strain(k, l)
            end do
          end do
        end do
      end do

    end function fElasStress



    function fDfGrdTh(thExpdCoef, dtemp)
      
      real(kind = RKIND), intent(in) :: thExpdCoef
      real(kind = RKIND), intent(in) :: dtemp
      real(kind = RKIND) :: fDfGrdTh(3, 3)

      fDfGrdTh = thExpdCoef*dtemp*UNITMAT

    end function fDfGrdTh




    function fPrtlGreenStrain(dfGrdElsPrd)
      implicit none

      real(kind = RKIND), intent(in) :: dfGrdElsPrd(3, 3)
      real(kind = RKIND) :: fPrtlGreenStrain(3, 3, 3, 3)

      integer(kind = IKIND) :: i, j, k, l

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              fPrtlGreenStrain(i, j, k, l) = 0.5d0*(UNITMAT(i, l)*dfGrdElsPrd(k, j)  &
                                                  + UNITMAT(j, l)*dfGrdElsPrd(k, i))
            end do
          end do
        end do
      end do

    end function fPrtlGreenStrain

  end module deformation

