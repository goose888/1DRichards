module tridiagonalmod
! Description: Perform LU decomposition to solve tridiagonal equation set
  implicit none
  public :: tridiagonal
!
contains
!
  subroutine tridiagonal(siz, A, f, x)
  ! Description: solve tridiagonal equation set A * x = f by LU decomposition
  ! USES:
    implicit none
  !
  ! Arguments:
    integer, intent(in) :: siz                       ! Size of matrix
    real, dimension(:,:), intent(in) :: A            ! Coefficient matrix
    real, dimension(:,:), intent(in) :: f            ! Right side vector
    real, dimension(:,:), intent(out) :: x
  !
  ! Local variables:
    integer :: i, j         ! indices
    real, dimension(siz, siz) :: L          ! Lower matrix
    real, dimension(siz, siz) :: U          ! Upper matrix
    real, dimension(siz, 1) :: y            ! temporary vector
  !
    ! LU decoposition
    do i = 1,siz
      do j = 1,siz
        if(i == (j+1)) then
          L(i, j) = A(i ,j)
          U(i, j) = 0.0
        else if(i == j) then
          if (i == 1) then
            L(i, j) = A(i, j)
            U(i, j) = 1.0
          else
            L(i, j) = A(i ,j) - L(i, j-1)*U(i-1, j)
            U(i, j) = 1.0
          end if
        else if((i+1) == j) then
          L(i, j) = 0.0
          U(i, j) = A(i, j) / L(i, j-1)
        else
          L(i, j) = 0.0
          U(i, j) = 0.0
        end if          
      end do
    end do
  !
    ! Now A = LU, solve Ly = f and Ux = y
    do i = 1,siz
      if(i == 1) then
        y(i, 1) = f(i, 1) / L(i, i)
!        write(6,*) 'L: ', L(i, i)
      else
        y(i, 1) = (f(i, 1) - L(i, i-1)*y(i-1, 1)) / L(i, i)
!        write(6,*) 'L-1: ', L(i, i-1)
!        write(6,*) 'L: ', L(i, i)
      end if
    end do
!    write(6,*) 'y: ', y
!    pause
    do i = siz,1,-1
      if(i == siz) then
        x(i, 1) = y(i, 1)
      else
        x(i, 1) = y(i, 1) - U(i, i+1)*x(i+1, 1)
!        write(6,*) 'U+1: ', U(i, i+1)
      end if
    end do
!    write(6,*) 'x: ', x
!    pause
    
  end subroutine tridiagonal
  
end module tridiagonalmod
