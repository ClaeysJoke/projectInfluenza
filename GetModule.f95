module GetModule

contains


  real*8 function getBeta(X)
    real*8 :: X(:)
    real*8 :: sigma = 0.0421D0
    real*8,dimension(17) :: tau,log_reg
    real*8 :: PT,I0,rho
    tau(1) = -49.7540D0
    tau(2) = -0.9577D0
    tau(3) = -0.0065D0
    tau(4) = -9.4896D0
    tau(5) = -0.3761D0
    tau(6) = -590.0001D0
    tau(7) = -2537.6102D0
    tau(8) = -4756.1828D0
    tau(9) = -3265.2458D0
    tau(10) = -102.2665D0
    tau(11) = -4.0162D0
    tau(12) = -430.9596D0
    tau(13) = -16.7104D0
    tau(14) = -798.3443D0
    tau(15) = -30.6638D0
    tau(16) = -543.8857D0
    tau(17) = -20.7459D0

    PT = X(7)
    I0 = X(4)
    rho = X(8)

    log_reg(1) = 1
    log_reg(2) = log(PT)
    log_reg(3) = log(PT)**2
    log_reg(4) = log(I0)
    log_reg(5) = log(I0)**2
    log_reg(6) = log(rho)
    log_reg(7) = log(rho)**2
    log_reg(8) = log(rho)**3
    log_reg(9) = log(rho)**4
    log_reg(10) = log(I0)*log(rho)
    log_reg(11) = log(I0)**2*log(rho)
    log_reg(12) = log(I0)*log(rho)**2
    log_reg(13) = log(I0)**2*log(rho)**2
    log_reg(14) = log(I0)*log(rho)**3
    log_reg(15) = log(I0)**2*log(rho)**3
    log_reg(16) = log(I0)*log(rho)**4
    log_reg(17) = log(I0)**2*log(rho)**4

    getBeta = exp(DOT_PRODUCT(log_reg,tau)+0.5D0*sigma**2)

  end function

  real*8 function getRho(X)
        real*8 :: a,b
        real*8 :: I0,S0,PI
        real*8, intent(in) :: X(:)
        real*8 :: fa, fb, temp, sol
        integer :: n,i
        integer :: m = 100


  	!print*,i
  	a = 0.01
  	b = 100
  	I0 = X(4)
  	PI = X(6)
    PT = X(7)
  	!print*,i
  	S0 = X(3)
          fa = evaluate_g_inverse(a,I0,S0,PI)
          fb = evaluate_g_inverse(b,I0,S0,PI)
          if (abs(fa) >  abs(fb)) then
            temp = a
            a = b
            b = temp
            temp = fa
            fa = fb
            fb = temp
          end if
        !print *,"    n        x(n)         f(x(n))"
        !print *," 1 ", b, fb
        !print *," 0 ", a, fa
        do n = 2,m
           if (abs(fa) >  abs(fb)) then
              temp = a
              a = b
              b = temp
              temp = fa
              fa = fb
              fb = temp
           end if
           temp = (b - a)/(fb - fa)
    	 b = a
  	 fb = fa
           a = a - fa*temp
           fa = evaluate_g_inverse(a,I0,S0,PI)
           !print *,n,a,fa
        end do
  	getRho = a

  end function


  function evaluate_g_inverse(x,I0,S0,PI) result(fx)
  	real*8,intent(in) :: x,I0,S0,PI
  	fx = I0 + S0 - PI - x*(log(S0) + 1 - log(x))
  end function


  ! Halley's method, rate of convergence 3
  real*8 function getRhoIter(X)
    real*8 :: X(:)
    real*8 :: I,PI
    real*8 :: g,der_g,der2_g
    real*8 :: rho

    I = X(4)
    PI = X(6)

    ! Initial value close to 0
    rho = 0.1D0

    g = PI + 0.1D0 - rho*(log(0.1D0-I)+1.0D0-log(rho))

    do while (abs(g)<1000*EPSILON(g))
      der_g = -log(0.1D0-I) + log(rho)
      der2_g = 1.0D0/rho

      rho = rho - g*der_g/(der_g**2-0.5D0*g*der2_g)

      if (rho<tiny(rho)) then
        rho = tiny(rho)
      endif

      g = PI + 0.1D0 - rho*(log(0.1D0-I)+1.0D0-log(rho))

    enddo

    getRhoIter = rho


  end function


end module
