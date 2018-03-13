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
    real*8 :: lambda
    real*8 :: alpha,beta
    integer :: t
    real*8 :: thetaI
    real*8 :: g,h,der_g,der_h
    real*8 :: gammaLambda,gammaThetaI,gammaNegThetaI
    real*8 :: digammaLambda
    real*8 :: ONE = 1.0D0


    lambda = X(2)
    gammaLambda = GAMMA(lambda)
    digammaLambda = digamma(lambda)
    ! Parameters initializing
    alpha = 2.0D0
    beta = 0.0001D0
    ! Gradient due to Gamma distribution
    gradLambda = -beta**alpha*((alpha-1)*lambda**(alpha-2)*exp(-beta*lambda)  &
            - beta*lambda**(alpha-2)*exp(-beta*lambda))/(pdfGamma(lambda)*GAMMA(alpha))

    ! Loop over all y Beta distributions
    do t = 1,size(y)
      thetaI = X(11 + 3*(i-1))
      pdf = pdfBeta(y(i),X(2)*thetaI,X(2)*(ONE-thetaI))

      gammaThetaI = GAMMA(lambda*thetaI)
      gammaNegThetaI = GAMMA(lambda*(ONE-thetaI))

      g = y(t)**(lambda*thetaI-ONE)*(ONE-y(t))**(lambda*thetaI-ONE)
      h = gammaThetaI*gammaNegThetaI/gammaLambda

      der_g = (thetaI*log(y(t))+thetaI*log(ONE-y(t))+digammaLambda)&
      * gammaLambda*y(t)**(lambda*thetaI-ONE)*(ONE-y(t))**(lambda*thetaI-ONE)
      der_h = (thetaI*digamma(lambda*thetaI)+(ONE-thetaI)*digamma(ONE-thetaI))&
      * gammaThetaI*gammaNegThetaI

      gradLambda = gradLambda - (der_g*h- der_h*g)/(pdf*h**2)
    end do

  end function

  real*8 function gradS(X,y)
    real*8 :: X(:),y(:)
    gradS = 0.0D0
  end function

  real*8 function gradI(X,y)
    real*8 :: X(:),y(:)
    real*8 :: alpha,beta
    real*8 :: pdf
    real*8 :: I
    real*8 :: ONE = 1.0D0
    real*8 :: TWO = 2.0D0
    real*8 :: h = 0.001D0
    real*8,dimension(2) :: z
    real*8,dimension(3) :: dirichletWeights
    real*8 :: thetaPlus,thetaMinus

    I = X(4)

    ! Beta distribution for I0
    alpha = 1.62D0
    beta = 7084.1D0

    pdf = pdfBeta(I,alpha,beta)

    gradI = -((alpha-ONE)*I**(alpha-TWO)*(ONE-I)**(beta-ONE) &
    - (beta-ONE)*I**(alpha-ONE)*(ONE-I)**(beta-TWO) &
    )*GAMMA(alpha+beta)/(GAMMA(alpha)*GAMMA(beta))/pdf

    ! Derivative w.r.t truncated normal distribution, using finite central difference
    z(1) = X(6)
    z(2) = X(7)

    pdf = pdfTruncatedNormal(z,I)

    gradI = gradI - (pdfTruncatedNormal(z,I+0.5D0*h) - pdfTruncatedNormal(z,I-0.5D0*h)) &
    /(h*pdf)

    ! Derivative w.r.t update step theta0 => theta1
    pdf = pdfDirichlet(X(10:12),X(1)*rkvec(X(3:5),X(9),X(9)*X(8)))
      ! thetaPlus
    dirichletWeights(1) = X(3)
    dirichletWeights(2) = X(4) + 0.5D0
    dirichletWeights(3) = X(5)

    thetaPlus = pdfDirichlet(X(10:12),X(1)*rkvec(dirichletWeights,X(9),X(9)*X(8)))

      ! thetaMinus
    dirichletWeights(1) = X(3)
    dirichletWeights(2) = X(4) - 0.5D0
    dirichletWeights(3) = X(5)

    thetaMinus = pdfDirichlet(X(10:12),X(1)*rkvec(dirichletWeights,X(9),X(9)*X(8)))

      ! Combining
    gradI = gradI - (thetaPlus-thetaMinus)/(h*pdf)
  end function

  real*8 function gradR(X,y)
    real*8 :: X(:),y(:)
    gradR = 0.0D0
  end function

  real*8 function gradPI(X,y)
    real*8 :: X(:),y(:)
    real*8,dimension(2) :: z
    real*8 :: thetaI
    real*8 :: pdf
    real*8 :: PIPlus,PIMinus
    real*8 :: h = 0.001D0

    thetaI = X(4)

    ! Original pdf
    z(1) = X(6)
    z(2) = X(7)
    pdf = pdfTruncatedNormal(z,thetaI)

    ! Central finite difference
    z(1) = z(1) + h*0.5D0
    PIPLus = pdfTruncatedNormal(z,thetaI)
    z(1) = z(1) - h
    PIMinus = pdfTruncatedNormal(z,thetaI)

    gradPI = -(PIPLus-PIMinus)/(h*pdf)
  end function

  real*8 function gradPT(X,y)
    real*8 :: X(:),y(:)
    real*8,dimension(2) :: z
    real*8 :: thetaI
    real*8 :: pdf
    real*8 :: PTPlus,PTMinus
    real*8 :: h = 0.001D0

    thetaI = X(4)

    ! Original pdf
    z(1) = X(6)
    z(2) = X(7)
    pdf = pdfTruncatedNormal(z,thetaI)

    ! Central finite difference
    z(2) = z(2) + h*0.5D0
    PTPLus = pdfTruncatedNormal(z,thetaI)
    z(2) = z(2) - h
    PTMinus = pdfTruncatedNormal(z,thetaI)

    gradPT = -(PTPLus-PTMinus)/(h*pdf)
  end function

  real*8 function gradRho(X,y)
    real*8 :: X(:),y(:)
    gradRho = 0.0D0
  end function

  real*8 function gradBeta(X,y)
    real*8 :: X(:),y(:)
    gradBeta = 0.0D0
  end function

  subroutine gradTheta(X,y,gradX)
    real*8,intent(in) :: X(:),y(:)
    real*8,intent(out) :: gradX(:)

  end subroutine






end module
