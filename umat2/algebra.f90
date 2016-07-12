  module algebra
    use utils, only : RKIND, IKIND

    implicit none

    contains

  
    pure elemental function sech(x)

      real(kind = RKIND), intent(in)  :: x
      real(kind = RKIND) :: sech

      real(kind = RKIND) :: tmp

      if(dabs(x) > 50.0d0) then
        sech = 0.0d0
      else
        tmp   = dexp(x)
        sech = 2.0d0/(tmp + 1.0/tmp)
      endif

    end function sech


    function polynomial(coef, x)
      real(kind = RKIND), intent(in) :: coef(:)
      real(kind = RKIND), intent(in) :: x

      real(kind = RKIND) :: polynomial

      integer(kind = IKIND) :: order, i

      order = size(coef) 
      polynomial = coef(order)
      do i = order, 2, -1
        polynomial = polynomial*x + coef(i-1)
      end do

    end function polynomial


    !> Euclidean norm of the given vector
    pure function vecNorm(vec1)
      
      real(kind = RKIND), intent(in)  :: vec1(:)  !> the input vector
      real(kind = RKIND) :: vecnorm

      vecNorm = dsqrt(dot_product(vec1, vec1))

    end function vecNorm


    !> Normalization of the given vector
    subroutine normalize(vec1)

      real(kind = RKIND) :: vec1(:)
      real(kind = RKIND) :: vnorm

      vnorm = vecnorm(vec1)
      vec1  = vec1/vnorm

    end subroutine
      

    !> Cross product of two vectors in 3D Euclidean space
    pure function cross_product(vec1, vec2)

      real(kind = RKIND), intent(in)  :: vec1(3)
      real(kind = RKIND), intent(in)  :: vec2(3)
      real(kind = RKIND) :: cross_product(3)

      cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
      cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
      cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

    end function cross_product


    ! Tensorial product of two vectors in 3D Euclidean space
    pure function out_product(vec1, vec2)

      real(kind = RKIND), intent(in)  :: vec1(3)
      real(kind = RKIND), intent(in)  :: vec2(3)
      real(kind = RKIND) :: out_product(3, 3)
      
      integer(kind = IKIND) :: i, j

      do i = 1, 3
        do j = 1, 3
          out_product(j, i) = vec1(j)*vec2(i)
        end do
      end do

    end function out_product


    !> Norm of the given matrix
    pure function matNorm(mat)

      real(kind = RKIND), intent(in) :: mat(:, :)
      real(kind = RKIND) :: matNorm

      matNorm = dsqrt(sum(mat*mat))

    end function matNorm
  

    !> Rotate a matrix: \f$ A^R = R^TAR \f$
    function matRot(mat, rot)

      real(kind = RKIND), intent(in) :: mat(3, 3)
      real(kind = RKIND), intent(in) :: rot(3, 3)
      real(kind = RKIND) :: matRot(3, 3)

      real(kind = RKIND) :: rotTr(3, 3)

      rotTr = transpose(rot)

      matRot = matmul(matmul(rotTr, mat), rot)

    end function matRot

    

    !> Tensorial product of two 2nd order tensor in 3D Eculidean space 
    function tenMul(ten1, ten2)
      real(kind = RKIND), intent(in) :: ten1(3, 3), ten2(3, 3)
      real(kind = RKIND) :: tenMul(3, 3, 3, 3)

      integer(kind = IKIND) :: i, j, k, l

      do l = 1, 3
        do k = 1, 3
          do j = 1, 3
            do i = 1, 3
              tenMul(i, j, k, l) = ten1(i, j) * ten2(k, l)
            end do
          end do
        end do
      end do

    end function tenMul


    function ten4Rot(ten4, rot)
      
      real(kind = RKIND) :: ten4(3, 3, 3, 3)
      real(kind = RKIND) :: rot(3, 3)
      real(kind = RKIND) :: ten4Rot(3, 3, 3, 3)

      real(kind = RKIND) :: ten4Aux(3, 3, 3, 3)
      integer(kind = IKIND) :: i, j, k, l, m, n, p, q


      ten4Rot = 0.0d0

      !ten4Aux = tenMul(Rot, Rot)

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              do m = 1, 3
                do n = 1, 3
                  do p = 1, 3
                    do q = 1, 3
                      ten4Rot(i, j, k, l) = ten4Rot(i, j, k, l)                      &
                                          + rot(m, i)*rot(n, j)*rot(p, k)*rot(q, l)  &
                                       !   + ten4Aux(m, i, n, j)*ten4Aux(p, k, q, l)  &
                                          * ten4(m, n, p, q)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

    end function ten4Rot



    !> Contraction of two given 4th order tensors
    function tenCtrct44(ten1, ten2)

      real(kind = RKIND) :: ten1(3, 3, 3, 3)
      real(kind = RKIND) :: ten2(3, 3, 3, 3)
      real(kind = RKIND) :: tenCtrct44(3, 3, 3, 3)
      integer(kind = IKIND) :: i, j, k, l

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            do l = 1, 3
              tenCtrct44(i, j, k, l) = sum(ten1(i, j, :, :)*ten2(:, :, k, l))
            end do
          end do
        end do
      end do

    end function tenCtrct44



    subroutine upper(mat)

      real(kind = RKIND)    :: mat(:, :)

      integer(kind = IKIND) :: row, col
      integer(kind = IKIND) :: i, j
      real(kind = RKIND)    :: aux

      row = size(mat, 1)
      col = size(mat, 2)

      do i = 1, row - 1
        do j = i + 1, row
          aux = mat(j, i)/mat(i, i)
          mat(j, i:row)  = mat(j, i:row)  - mat(i, i:row)*aux
          mat(j, row+1:) = mat(j, row+1:) - mat(i, row+1:)*aux
        end do
      end do

    end subroutine upper
    

    subroutine lower(mat)

      real(kind = RKIND)    :: mat(:, :)

      integer(kind = IKIND) :: row, col
      integer(kind = IKIND) :: i, j
      real(kind = RKIND)    :: aux

      row = size(mat, 1)
      col = size(mat, 2)

      do i = row, 2, -1
        do j = i - 1, 1, -1
          aux = mat(j, i)/mat(i, i)
          mat(j, 1:row) = mat(j, 1:row) - mat(i, 1:row)*aux
          mat(j, row+1:) = mat(j, row+1:) - mat(i, row+1:)*aux
        end do
      end do

    end subroutine lower


    !> Determinant of the given matrix 
    function determinant(mat)

      real(kind = RKIND), intent(in) :: mat(:, :) 
      real(kind = RKIND) :: determinant

      real(kind = RKIND), allocatable :: ma(:, :)
      integer(kind = IKIND) :: row, col
      integer(kind = IKIND) :: i


      row = size(mat, 1)
      col = size(mat, 2)

      if (row /= col) then
          write(*, *) "Error::determinant: invalid mat size"
          stop
      end if

      if (row == 3) then
        determinant = mat(1, 1)*(mat(2, 2)*mat(3, 3) - mat(2, 3)*mat(3, 2))  &
                    - mat(1, 2)*(mat(2, 1)*mat(3, 3) - mat(2, 3)*mat(3, 1))  &
                    + mat(1, 3)*(mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1))
      else if (row == 2) then
        determinant = mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1)
      else if (row == 1) then
        determinant = mat(1, 1) 
      else
        allocate(ma(row, col))
        ma = mat
        call upper(ma)
        determinant = 1.0d0
        do i = 1, row
          determinant = determinant * ma(i, i)
        end do

        if (allocated(ma)) deallocate(ma)
      end if

    end function determinant


    !> 
    function matInvars(mat)

      real(kind = RKIND), intent(in) :: mat(3, 3)
      real(kind = RKIND) :: matInvars(3)

      real(kind = RKIND) :: mat2(3, 3), tr2


      matInvars(1) = mat(1, 1) + mat(2, 2) + mat(3, 3)

      mat2 = matmul(mat, mat)
      tr2  = mat2(1, 1) + mat2(2, 2) + mat2(3, 3)
      matInvars(2) = 0.5d0*(matInvars(1)*matInvars(1) - tr2)

      matInvars(3) = determinant(mat)

    end function matInvars


    !> Calculate the inverse of a 3X3 matrix
    function matInv(mat)

      real(kind = RKIND) :: mat(3, 3)
      real(kind = RKIND) :: matInv(3, 3)

      integer(kind = IKIND) :: row1, col1, row2, col2
      real(kind = RKIND) :: aux(25), det
      integer(kind = IKIND) :: i

      aux(1)  =  mat(1, 1)
      aux(2)  =  mat(1, 2)
      aux(3)  =  mat(1, 3)
      aux(4)  =  mat(2, 1)
      aux(5)  =  mat(2, 2)
      aux(6)  =  mat(2, 3)
      aux(7)  =  mat(3, 1)
      aux(8)  =  mat(3, 2)
      aux(9)  =  mat(3, 3)
      aux(10) = -aux(2)  * aux(4)
      aux(11) =  aux(3)  * aux(4)
      aux(12) =  aux(1)  * aux(5)
      aux(13) = -aux(3)  * aux(5)
      aux(14) = -aux(1)  * aux(6)
      aux(15) =  aux(2)  * aux(6)
      aux(16) =  aux(13) * aux(7)
      aux(17) =  aux(15) * aux(7)
      aux(18) =  aux(2)  * aux(7)
      aux(19) = -aux(3)  * aux(7)
      aux(20) = -aux(5)  * aux(7)
      aux(7)  =  aux(6)  * aux(7)
      aux(21) = -aux(1)  * aux(8)
      aux(22) =  aux(11) * aux(8)
      aux(23) =  aux(14) * aux(8)
      aux(3)  =  aux(3)  * aux(8)
      aux(24) =  aux(4)  * aux(8)
      aux(6)  = -aux(6)  * aux(8)
      aux(1)  =  aux(1)  * aux(9)
      aux(8)  =  aux(10) * aux(9)
      aux(25) =  aux(12) * aux(9)
      aux(2)  = -aux(2)  * aux(9)
      aux(4)  = -aux(4)  * aux(9)
      aux(5)  =  aux(5)  * aux(9)
      aux(9)  =  aux(10) + aux(12)
      aux(10) =  aux(11) + aux(14)
      aux(11) =  aux(13) + aux(15)
      aux(12) =  aux(18) + aux(21)
      aux(13) =  aux(20) + aux(24)
      aux(1)  =  aux(1)  + aux(19)
      aux(8)  =  aux(16) + aux(17) + aux(22) + aux(23) + aux(25) + aux(8)
      aux(2)  =  aux(2)  + aux(3)
      aux(3)  =  aux(4)  + aux(7)
      aux(4)  =  aux(5)  + aux(6)
      aux(5)  =  1.0d0/aux(8)
      aux(6)  =  aux(5)  * aux(9)
      aux(7)  =  aux(10) * aux(5)
      aux(8)  =  aux(11) * aux(5)
      aux(9)  =  aux(12) * aux(5)
      aux(10) =  aux(13) * aux(5)
      aux(1)  =  aux(1)  * aux(5)
      aux(2)  =  aux(2)  * aux(5)
      aux(3)  =  aux(3)  * aux(5)
      aux(4)  =  aux(4)  * aux(5)

      matInv(1, 1) = aux(4)
      matInv(2, 1) = aux(3)
      matInv(3, 1) = aux(10)
      matInv(1, 2) = aux(2)
      matInv(2, 2) = aux(1)
      matInv(3, 2) = aux(9)
      matInv(1, 3) = aux(8)
      matInv(2, 3) = aux(7)
      matInv(3, 3) = aux(6)

    end function matInv


    subroutine polarDcmp(FF, UU, RR)
      use utils, only : UNITMAT, PI
      
      real(kind = RKIND), intent(in)  :: FF(3, 3)
      real(kind = RKIND), intent(out) :: UU(3, 3)
      real(kind = RKIND), intent(out) :: RR(3, 3)

      real(kind = RKIND) :: UUInv(3, 3)
      real(kind = RKIND) :: CC(3, 3), CC2(3, 3)
      real(kind = RKIND) :: iu1, iu2, iu3
      real(kind = RKIND) :: x(3), invars(3), lambda(3)
      real(kind = RKIND) :: b, c, m, n, t, D, tmp1

      integer(kind = IKIND) :: ii
 

      CC = matmul(transpose(FF), FF)

      invars = matInvars(CC)

      tmp1 = invars(1)*invars(1)
      b = invars(2) - tmp1/3.0d0
      c = -2.0d0/27.0d0*tmp1*invars(1) + invars(1)*invars(2)/3.0d0 - invars(3)

      if(dabs(b) < 1.0d-15) then
        x = -c**(1.0d0/3.0d0)
      else
        m  = 2.0d0*dsqrt(dabs(-b)/3.0d0)
        n  = 3.0d0*c/(m*b)
        t  = datan2(dsqrt(dabs(1.0d0 - n*n)), n)/3.0d0
        tmp1 = 2.0d0/3.0d0*PI
        do ii = 1, 3
          x(ii) = m*dcos(t + tmp1*(ii-1))
        end do
      end if

      lambda = dsqrt(x + invars(1)/3.0d0)

!     Invariants of U
      iu1 = lambda(1) + lambda(2) + lambda(3)
      iu2 = lambda(1) * lambda(2) + lambda(1)*lambda(3) + lambda(2)*lambda(3)
      iu3 = lambda(1) * lambda(2) * lambda(3)

      D = iu1*iu2 - iu3
      if (dabs(D) <= 1d-15) write(*,*) 'ERROR D in subroutine polarDecomp', D

      CC2 = matmul(CC, CC)

      UU = (-CC2 + (iu1*iu1 - iu2)*CC + (iu1*iu3)*UNITMAT)/D
      UUInv = (CC - iu1*UU + iu2*UNITMAT)/iu3
      
      RR = matmul(FF, UUInv)

    end subroutine polarDcmp


    function Ten3333ToA99(Ten3333)

      real(kind = RKIND), intent(in) :: Ten3333(3, 3, 3, 3)
      real(kind = RKIND) :: Ten3333ToA99(9, 9)

      integer(kind = IKIND) :: idxii(9), idxjj(9), ii, jj

      idxii = (/1, 2, 3, 1, 2, 3, 1, 2, 3/)
      idxjj = (/1, 1, 1, 2, 2, 2, 3, 3, 3/)

      do ii = 1, 9
        do jj = 1, 9
          Ten3333ToA99(ii, jj) = Ten3333(idxii(ii), idxjj(ii), idxii(jj), idxjj(jj))
        end do
      end do

    end function Ten3333ToA99


    function Ten3333ToA66(Ten3333)

      real(kind = RKIND), intent(in) :: Ten3333(3, 3, 3, 3)
      real(kind = RKIND) :: Ten3333ToA66(6, 6)

      integer(kind = IKIND) :: idxii(6), idxjj(6), ii, jj

      idxii = (/1, 2, 3, 1, 1, 2/)
      idxjj = (/1, 2, 3, 2, 3, 3/)

      do ii = 1, 6
        do jj = 1, 6
          Ten3333ToA66(ii, jj) = Ten3333(idxii(ii), idxjj(ii), idxii(jj), idxjj(jj))
        end do
      end do

    end function Ten3333ToA66



    subroutine GaussJordan(A, S, ANS, iflag)
      use utils, only : SMALL

      real(kind = RKIND),    intent(in)  :: A(:, :)
      real(kind = RKIND),    intent(in)  :: S(:)
      real(kind = RKIND),    intent(out) :: ANS(:)
      integer(kind = IKIND), intent(out) :: iflag

      real(kind = RKIND), allocatable :: B(:, :)
      real(kind = RKIND)     :: maxV
      integer(kind = IKIND)  :: N
      integer(kind = IKIND)  :: i

      
      iflag = 0
      maxV = maxval(dabs(A))
      if (maxV < SMALL) then
        iflag = 1
        return
      end if

      N = size(A, 1)
      allocate(B(N, N+1))

      B(:, 1:N) = A
      B(:, N+1) = S

      call upper(B)
      call lower(B)

      forall(i = 1:N) ANS(i) = B(i, N+1)/B(i, i)

      if (allocated(B)) deallocate(B)

    end subroutine GaussJordan


    subroutine luDcmp(A, N, indx, D, code)

      integer(kind = IKIND), parameter :: NMAX = 100
      real(kind = RKIND),    parameter :: tiny = 1.0d-16

      integer(kind = IKIND), intent(in)    :: N
      real(kind = RKIND),    intent(inout) :: A(N, N)
      integer(kind = IKIND), intent(out) :: indx(N)
      integer(kind = IKIND), intent(out) :: D
      integer(kind = IKIND), intent(out) :: code

      real(kind = RKIND)    :: amax, dum, sum, vv(nmax)
      integer(kind = IKIND) :: i, j, k, imax 


      D    = 1
      code = 0

      do i = 1, N
        amax = 0.0d0
        do j = 1, N
          if (dabs(A(i, j)) > amax) amax = dabs(A(i, j))
        end do
        if(amax < TINY) then
          code = 1
          return
        end if

        vv(i) = 1.0d0 / amax
      end do                    ! i loop
      
      do j = 1, N
        do i = 1, j-1
          sum = A(i, j)
          do k = 1, i-1
            sum = sum - A(i, k)*A(k, j) 
          end do       ! k loop
          A(i, j) = sum
        end do          ! i loop

        amax = 0.0d0
        do i = j, N
          sum = A(i, j)
          do k = 1, j-1
            sum = sum - A(i, k)*A(k, j) 
          end do       ! k loop
          A(i, j) = sum
          dum = vv(i)*DABS(sum)
          if (dum >= amax) then
            imax = i
            amax = dum
          end if
        end do          ! i loop  

        if (j /= imax) then
          do k = 1, N
            dum = A(imax, k)
            A(imax, k) = A(j, k)
            A(j, k) = dum
          end do      ! k loop
          D = -D
          vv(imax) = vv(j)
        end if
 
        indx(j) = imax
        if (DABS(A(j, j)) < TINY) A(j, j) = TINY
 
        if(j /= N) then
          dum = 1.0d0 / A(j, j)
          do i = j+1, N
            A(i, j) = A(i, j)*dum
          end do  ! i loop
        end if 
      end do        ! j loop
 
    end subroutine ludcmp


    subroutine lubksb(A, N, indx, B)
      integer(kind = IKIND), intent(in)  :: N
      integer(kind = IKIND), intent(in)  :: indx(N)
      real(kind = RKIND),    intent(in)  :: A(N, N)
      real(kind = RKIND),    intent(inout) :: B(N)
      
      real(kind = RKIND) :: sum
      integer(kind = IKIND) :: i, j, ii, ll


      ii = 0
      do i = 1, N
        ll = indx(i)
        sum = B(ll)
        B(ll) = B(i)
        if (ii /= 0) then
          do j = ii, i-1
            sum = sum - A(i, j)*B(j)
          end do              
        else if (sum /= 0.0d0) then
          ii = i
        end if
        B(i) = sum
      end do        
      
      do i = N, 1, -1
        sum = B(i)
        if (i < N) then
          do j = i+1, N
            sum = sum - A(i, j)*B(j)
          end do              
        end if
        B(i) = sum / A(i, i)
      end do                  
      
    end subroutine lubksb

    !> Solve the linear equations \f$ Ax=b \f$
    subroutine gauss(AA, bb, xx, nmax, n, iflag)

      integer(kind = IKIND), intent(in)  :: nmax
      integer(kind = IKIND), intent(in)  :: n
      real(kind = RKIND),    intent(in)  :: AA(nmax, nmax) 
      real(kind = RKIND),    intent(in)  :: bb(nmax)
      real(kind = RKIND),    intent(out) :: xx(nmax)
      integer(kind = IKIND), intent(out) :: iflag

      real(KIND = RKIND)    :: A(n,n), b(n)
      integer(kind = IKIND) :: i, j, d, indx(n)


      do i = 1, n
        b(i) = bb(i)
        do j = 1, n
          A(i, j) = AA(i, j)
        end do
      end do

      call ludcmp(A, N, indx, d, iflag)
      call lubksb(A, N, indx, b)

      do i = 1, n
        xx(i) = b(i)
      end do

    end subroutine gauss




    subroutine output_mat(mat)

      real(kind = RKIND) :: mat(:, :)

      integer(kind = IKIND) :: row, col
      integer(kind = IKIND) :: i
      
      character (len = 20)  :: myFormat = '(??(1x, f6.3))' 

      row = size(mat, 1)
      col = size(mat, 2)

      write(myFormat(2:3), '(I2)') col

      do i = 1, row
        write(*, fmt=myFormat) mat(i, :)
      end do
    end subroutine output_mat


  end module algebra

