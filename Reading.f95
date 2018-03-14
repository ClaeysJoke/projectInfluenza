module Reading


  implicit none

contains
  ! Fills y vector with ILI+ datapoints.
  ! Parameter t is used to determine the number of points used
  subroutine ReadData(y,t)
    real*8,intent(out) :: y(:)
    integer,intent(in) :: t
    integer :: week,i
    real*8 :: ILI,ILI_prop,ILI_pos,ILI_plus
    open(1,'data_week.txt')

    do i = 1,t
      read(1,*) week,ILI,ILI_prop,ILI_pos,ILI_plus
      y(i) = ILI_plus
    enddo

    close(1)

end subroutine

end module
