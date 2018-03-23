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

  real*8 function getRho(S0,I0,PI)
        real*8 :: x1,x2
        real*8,intent(in) :: I0,S0,PI
        !real*8, intent(in) :: X(:)
        real*8 :: f1, f2, temp, sol
        integer :: n
        integer :: m = 10


  	!print*,i
  	x2 = 0.6D0
  	x1 = 0.8D0

          f2 = evaluate_g_inverse(x2,I0,S0,PI)
          !print *, "Function value left bound",fa
          f1 = evaluate_g_inverse(x1,I0,S0,PI)
          !print *, "Function value right bound",fb
        !  print *, "===== Iteration 1 ======"
        !   print *,"n-2 point",x2
        !   print *, "Function value n-2",f2
        !   print *,"n-1 point",x1
        !   print *, "Function value n-1",f1

        !print *,"    n        x(n)         f(x(n))"
        !print *," 1 ", b, fb
        !print *," 0 ", a, fa
        do n = 2,m
        !  print *, "===== Iteration",n," ======"
        !   print *,"n-2 point",x2
        !   print *, "Function value n-2",f2
        !   print *,"n-1 point",x1
        !   print *, "Function value n-1",f1
           temp = (x1 - x2)/(f1 - f2)
    	 x2 = x1
  	 f2 = f1
           x1 = x1 - f1*temp
           if (x1<tiny(x1)) then
             x1 = tiny(x1)
           endif
           f1 = evaluate_g_inverse(x1,I0,S0,PI)
           if (abs(x1-x2)<1000*epsilon(x1)) then
             exit
           endif
           !print *,n,a,fa
        end do
  	getRho = x1

  end function


  real*8 function evaluate_g_inverse(x,I0,S0,PI)
  	real*8,intent(in) :: x,I0,S0,PI

  	evaluate_g_inverse = I0 + S0 - PI - x*(log(S0) + 1 - log(x))
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

    g = -PI + I + 0.9D0 - rho*(log(0.9D0)+1.0D0-log(rho))

    do while (abs(g)>0.01D0)
      der_g = -log(0.9D0) + log(rho)
      der2_g = 1.0D0/rho

      rho = rho - g*der_g/(der_g**2-0.5D0*g*der2_g)

      if (rho<tiny(rho)) then
        rho = tiny(rho)
      endif

      g = PI + I + 0.9D0 - rho*(log(0.9D0)+1.0D0-log(rho))

    enddo

    getRhoIter = rho


  end function


end module
