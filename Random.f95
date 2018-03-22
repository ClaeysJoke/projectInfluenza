module Random
  use Normal
implicit none
contains

  RECURSIVE FUNCTION rand_gamma(shape, scale) RESULT(ans)
    implicit none
    real (kind=selected_real_kind(8) ) :: ans
    real (kind=selected_real_kind(8) ), intent(in) :: shape, scale !shape and scale
        real (kind=selected_real_kind(8) ) :: u,w,d,c,x,xsq,g,v
       IF (shape <= 0.0) THEN

          WRITE(*,*) "Shape PARAMETER must be positive"
        END IF
        IF (scale <= 0.0) THEN

          WRITE(*,*) "Scale PARAMETER must be positive"
        END IF

        IF (shape >= 1.0) THEN
          d = shape - 1.0/3.0
          c = 1.0/(9.0*d)**0.5
          DO while (.true.)
              x = rand_normal(0.0_8, 1.0_8)
      !print*,x
              v = 1.0 + c*x
              DO while (v <= 0.0)
                  x = rand_normal(0.0_8, 1.0_8)
                  v = 1.0 + c*x
              END DO

              v = v*v*v
        call init_random_seed()
              CALL RANDOM_NUMBER(u)
              xsq = x*x
              IF ((u < 1.0 -.0331*xsq*xsq) .OR.  &
                (log(u) < 0.5*xsq + d*(1.0 - v + log(v))) )then
                  ans=scale*d*v
                  RETURN
              END IF

          END DO
        ELSE
          g = rand_gamma(shape+1.0, 1.0_8)
    call init_random_seed()
          CALL RANDOM_NUMBER(w)
          ans=scale*g*(w)**(1.0/shape)
          RETURN
        END IF
    END FUNCTION

    FUNCTION rand_beta(a, b) RESULT(ans)
  implicit none
  real*8, intent(in) :: a,b
  real*8 :: ans
        real*8 :: u,v
      IF ((a <= 0.0) .OR. (b <= 0.0)) THEN

        WRITE(*,*) "Beta PARAMETERs must be positive"
      END IF

       u = rand_gamma(a, 1.0_8)
       v = rand_gamma(b, 1.0_8)
       ans = u / (u + v)
  END FUNCTION

  subroutine rand_Dirichlet(m,param,theta)
implicit none
!Subroutine that draws from m-dimensional Dirichlet
!using gamma distribution Gamm(alpha,1)
integer, intent(in) :: m
real*8, intent(in) :: param(m)
real*8, intent(out) ::theta(m)
integer	:: i
real*8 :: y,s
s=0.0
do i=1,m
y= rand_gamma(1.0D0,param(i))
theta(i)=y
s=s+y
end do
theta=(1.0/s)*theta
end subroutine


end module
