module Reading

  implicit none

contains
  ! Fills y vector with ILI+ datapoints.
  ! Parameter t is used to determine the number of points used
  ! STARTS FROM WEEK 46
  subroutine readData(y,t)
    real*8,intent(out) :: y(:)
    integer,intent(in) :: t
    integer :: week,i
    real*8 :: ILI,ILI_prop,ILI_pos,ILI_plus
    open(7,file='data_week.txt')

    do i = 1,(t+5)
      read(7,*) week,ILI,ILI_prop,ILI_pos,ILI_plus
      if (i>5)then
      y(i-5) = ILI_plus
    endif
    enddo

    close(7)

end subroutine

 ! Reads from control.txt
!  subroutine readData(y,t)
!    real*8,intent(out) :: y(:)
!    integer,intent(in) :: t
!    integer :: week,i
!    real*8 :: ILI_plus
!
!    open(unit=7,file='Control.txt')
!
!    do i = 1,t
!      read(7,*) ILI_plus
!      y(i) = ILI_plus
!    enddo
!    close(7)
!  end subroutine

end module
