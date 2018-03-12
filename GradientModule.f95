module GradientModule
  use Pdf
  use PsiModule

contains
  ! Computes the gradient of the potential function V w.r.t. the parameter
  ! vector and the hidden variables
  subroutine Gradient(X,y,gradX)
    real*8,intent(in) :: X(:),y(:)
    real*8,intent(out) :: gradX(:)

    gradX(1) = gradKappa(X,y)
    gradX(2) = gradLambda(X,y)
    gradX(3) = gradS(X,y)
    gradX(4) = gradI(X,y)
    gradX(5) = gradR(X,y)
    gradX(6) = gradPI(X,y)
    gradX(7) = gradPT(X,y)
    gradX(8) = gradRho(X,y)
    gradX(9) = gradBeta(X,y)

    call gradTheta(X,y,gradX)

  end subroutine

  real*8 function gradKappa(X,y)
    real*8 :: X(:),y(:)
    real*8 :: alpha,beta, kappa
    integer :: t
    real*8,dimension(3) :: alpha_dirichlet
    real*8 :: g,h,der_g,der_h
    real*8 :: gammaKappa,gammaS,gammaI,gammaR
    real*8 :: digammaKappa,digammaS,digammaI,digammaR
    integer :: ifault
    real*8 :: pdf
    real*8 :: thetaS,thetaI,thetaR
    real*8 :: alphaS,alphaI,alphaR

    ! Support variable kappa
    kappa = X(1)
    gammaKappa = GAMMA(kappa)
    digammaKappa = digamma(kappa,ifault)
    ! Parameters initializing
    alpha = 2.0D0
    beta = 0.0001D0

    ! Gradient due to Gamma distribution
    gradKappa = -beta**alpha*((alpha-1)*kappa**(alpha-2)*exp(-beta*kappa)  &
            - beta*kappa**(alpha-2)*exp(-beta*kappa))/(pdfGamma(kappa)*GAMMA(alpha))

    ! Gradient due to Dirichlet distribution of noise added during week-to-week update

    ! Looping over all t
    do t = 1,size(y)
    select case(t)
    case(1)
      alpha_dirichlet = rkvec(X(3:5),X(9),X(8)*X(9))
      thetaS = X(10)
      thetaI = X(11)
      thetaR = X(12)
      alphaS = alpha_dirichlet(1)
      alphaI = alpha_dirichlet(2)
      alphaR = alpha_dirichlet(3)
      gammaS = GAMMA(kappa*alphaS)
      gammaI = GAMMA(kappa*alphaI)
      gammaR = GAMMA(kappa*alphaR)
      digammaS = digamma(kappa*alphaS)
      digammaI = digamma(kappa*alphaI)
      digammaR = digamma(kappa*alphaR)

      pdf = pdfDirichlet(X(10:12),kappa*alpha_dirichlet)

      g = thetaS**(kappa*alphaS)*thetaI**(kappa*alphaI)*thetaR**(kappa*alphaR)*gammaKappa
      h = gammaS*gammaI*gammaR

      der_g = (alphaS*thetaS**(kappa*alphaS)*log(thetaS)*thetaI**(kappa*alphaI)*thetaR**(kappa*alphaR) &
      + thetaS**(kappa*alphaS)*alphaI*thetaI**(kappa*alphaI)*log(thetaI)*thetaR**(kappa*alphaR) &
      + thetaS**(kappa*alphaS)*thetaI**(kappa*alphaI)*alphaR*thetaR**(kappa*alphaR)*log(thetaR) &
      + thetaS**(kappa*alphaS)*thetaI**(kappa*alphaI)*thetaR**(kappa*alphaR)*digammaKappa)*gammaKappa
      der_h = (alphaS*digammaS+alphaI*digammaI+alphaR*digammaR)&
      *gammaS*gammaI*gammaR

      gradKappa = gradKappa - (der_g*h-g*der_h)/(pdf*h**2)

    case default
      alpha_dirichlet = rkvec(X(10+3*(t-1):9+3*t),X(9),X(8)*X(9))
      thetaS = X(10+3*(t-1))
      thetaI = X(11+3*(t-1))
      thetaR = X(12+3*(t-1))
      alphaS = alpha_dirichlet(1)
      alphaI = alpha_dirichlet(2)
      alphaR = alpha_dirichlet(3)
      gammaS = GAMMA(kappa*alphaS)
      gammaI = GAMMA(kappa*alphaI)
      gammaR = GAMMA(kappa*alphaR)
      digammaS = digamma(kappa*alphaS)
      digammaI = digamma(kappa*alphaI)
      digammaR = digamma(kappa*alphaR)

      pdf = pdfDirichlet(X(10+3*t:12+3*t),kappa*alpha_dirichlet)

      g = thetaS**(kappa*alphaS)*thetaI**(kappa*alphaI)*thetaR**(kappa*alphaR)*gammaKappa
      h = gammaS*gammaI*gammaR

      der_g = (alphaS*thetaS**(kappa*alphaS)*log(thetaS)*thetaI**(kappa*alphaI)*thetaR**(kappa*alphaR) &
      + thetaS**(kappa*alphaS)*alphaI*thetaI**(kappa*alphaI)*log(thetaI)*thetaR**(kappa*alphaR) &
      + thetaS**(kappa*alphaS)*thetaI**(kappa*alphaI)*alphaR*thetaR**(kappa*alphaR)*log(thetaR) &
      + thetaS**(kappa*alphaS)*thetaI**(kappa*alphaI)*thetaR**(kappa*alphaR)*digammaKappa)*gammaKappa
      der_h = (alphaS*digammaS+alphaI*digammaI+alphaR*digammaR)&
      *gammaS*gammaI*gammaR

      gradKappa = gradKappa - (der_g*h-g*der_h)/(pdf*h**2)
    end select
    enddo
  end function

  real*8 function gradLambda(X,y)
    real*8 :: X(:),y(:)

  end function

  real*8 function gradS(X,y)
    real*8 :: X(:),y(:)

  end function

  real*8 function gradI(X,y)
    real*8 :: X(:),y(:)

  end function

  real*8 function gradR(X,y)
    real*8 :: X(:),y(:)

  end function

  real*8 function gradPI(X,y)
    real*8 :: X(:),y(:)

  end function

  real*8 function gradPT(X,y)
    real*8 :: X(:),y(:)

  end function

  real*8 function gradRho(X,y)
    real*8 :: X(:),y(:)

  end function

  real*8 function gradBeta(X,y)
    real*8 :: X(:),y(:)

  end function

  subroutine gradTheta(X,y,gradX)
    real*8,intent(in) :: X(:),y(:)
    real*8,intent(out) :: gradX(:)

  end subroutine






end module
