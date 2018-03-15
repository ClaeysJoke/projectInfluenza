module Prior
  use Sampling
  use Normal

contains

  subroutine SamplePrior(X,y)
    real*8 :: X(:),y(:)
    integer :: i
    real*8 :: shape,rate
    real*8 :: alpha,beta
    real*8,dimension(3) :: theta

    shape = 2.0D0
    rate = 0.0001D0
    alpha = 1.62D0
    beta = 7084.1D0


    ! Sample Kappa
    print *, "============ PRIOR ================="
    print *, "Sampling Prior kappa"
    X(1) = rand_gamma(shape,1.0D0/rate)
    ! Sample Lambda
    print *, "Sampling Prior Lambda"
    X(2) = rand_gamma(shape,1.0D0/rate)
    ! "Sample" S0
    print *, "Sampling Prior S0"
    X(3) = 0.9D0
    ! Sample I0
    print *, "Sampling Prior I0"
    X(4) = rand_beta(alpha,beta)
    ! "Sample" R0
    print *, "Sampling Prior R0"
    X(5) = 1.0D0 - X(3) - X(4)

    ! "Sample" PI and PT
    print *, "Sampling Prior PI and PT"
    X(6) = 0.0144D0
    X(7) = 17.9D0


    ! "Sample" rho
    print *, "Sampling Prior Rho"
    X(8) = 1.0D0
    X(8) = getRho(X)
    ! "Sample" beta
    print *, "Sampling Prior Beta"
    X(9) = getBeta(X)
    ! "Sample" thetas, by putting them equal to theta0
    print *, "Sampling thetas"
    theta(1) = X(3)
    theta(2) = X(4)
    theta(3) = X(5)

    do i = 1,size(y)
      X(10+3*(i-1)) = 0.9D0
      X(11+3*(i-1)) = y(i)
      X(12+3*(i-1)) = 1 - y(i) - 0.9D0
    enddo

  end subroutine

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




end module
