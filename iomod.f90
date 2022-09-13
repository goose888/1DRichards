module iomod
  ! Description: I/O module
  implicit none
  public :: wrtout
!
contains
!
  subroutine wrtout(theta, phi, siz, t, it)
    ! Description: Output of result
    implicit none
  !
    ! Arguments:
    real, dimension(:, :), intent(in) :: theta
    real, dimension(:, :), intent(in) :: phi
    integer, intent(in) :: siz
    integer, intent(in) :: t
    integer, intent(in) :: it
  !
    ! Local variables:
    integer :: i              ! indice
    logical :: alive
  !
    write(6, *) 'Step ', t, ' finished!'
    write(6, *) 'Start writting soil water states into LOG file!'
  
    inquire(file = 'richards.out', exist = alive)
    if(.not. alive) then
      open(unit=10, file='richards.out', status='new')
      write(10, *) '************* Output of Richards Equation ****************'
    else
      open(unit=10, file='richards.out', status='old', position='append')
    end if
    write(10, *) 'nstep = ', t
    write(10, *) 'iteration = ', it
    do i = 1,siz
      write(10, *) 'Layer ', i, ', theta = ', theta(i, 1), ', phi = ', phi(i, 1)
    end do
    write(10, *) ' '
    close(10)
  
    write(6, *) 'Successfully write soil water states into "richards.out"!'
    write(6, *) ' '
  
  end subroutine wrtout
  
end module iomod
