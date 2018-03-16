module Prior
  use Sampling
  use Normal
  use ODE

contains

  subroutine SamplePrior(X,y)
    real*8 :: X(:),y(:)
    integer :: i
    real*8 :: shape,rate
    real*8 :: alpha,beta

    shape = 2.0D0
    rate = 0.0001D0
    alpha = 1.62D0
    beta = 7084.1D0


    ! Sample Kappa
    print *, "============ PRIOR ================="

    X(1) = rand_gamma(shape,1.0D0/rate)
    ! Sample Lambda

    X(2) = rand_gamma(shape,1.0D0/rate)
    ! "Sample" S0

    X(3) = 0.9D0
    ! Sample I0
    ! randGamma is still giving shitty values :'(

    X(4) = 0.005D0


    ! "Sample" R0

    X(5) = 1.0D0 - X(3) - X(4)

    ! "Sample" PI and PT

    X(6) = 0.0144D0
    X(7) = 17.9D0


    ! "Sample" rho

    X(8) = 1.0D0
    X(8) = getRho(X)
    ! "Sample" beta

    X(9) = getBeta(X)
    ! "Sample" thetas, by simulating them with the ODE solver from theta0


    do i = 1,size(y)
      select case(i)
      case(1)
        X(10:12) = rkvec(X(3:5),X(9),X(8)*X(9))
      case default
        X(10+3*(i-1):12+3*(i-1)) = rkvec(X(10+3*(i-2):12+3*(i-2)),X(9),X(8)*X(9))
      end select

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
