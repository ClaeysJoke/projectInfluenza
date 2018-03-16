module Pdf
  USE PRECISION_MODEL
  use MVSTAT


  implicit none
  real*8,parameter :: PI = 4*ATAN(1.0_8)

contains
  ! Calculates the Gamma probablitiy density function of the scalar x.
  ! Gamma distribution has the following parameters:
  !   shape : alpha = 2.0
  !   rate : beta = 0.0001
  real*8 function pdfGamma(x)
    real*8 :: x
    real*8 :: alpha, beta

    alpha = 2.0D0
    beta = 0.0001D0

    if (x > 0) then
      pdfGamma = (beta**alpha)*(x**(alpha-1))*EXP(-beta*x)/DGAMMA(alpha)
    else
      pdfGamma = 0.0D0
    endif

  end function

  real*8 function logPdfGamma(x)
    real*8 :: x
    real*8 :: gamma

    gamma = pdfGamma(x)

    if (gamma > 0) then
      logPdfGamma = log(gamma)
    else
      logPdfGamma = -huge(0.0D0)
    endif

  end function

  ! Calculates the Beta probablitiy density function of the scalar x.
  ! Beta distribution has given parameters alpha and beta
  real*8 function pdfBeta(x,alpha,beta)
    real*8 :: x,alpha,beta
    real*8 :: betaFunction
    real*8 :: ONE = 1.0D0
    real*8 :: ZERO = 0.0D0


    if ((x>ZERO) .and. (x<ONE)) then
      betaFunction = DGAMMA(alpha)*DGAMMA(beta)/DGAMMA(alpha+beta)
      pdfBeta = x**(alpha-ONE)*(ONE-x)**(beta-ONE)/betaFunction
    else
      pdfBeta = ZERO
    endif
  end function

  ! Calculates the log beta probability density function of the scalar x.
  ! Beta distribution has given parameters alpha and beta
  real*8 function logPdfBeta(x,alpha,beta)
  real*8 :: x,alpha,beta
  real*8 :: ONE = 1.0D0
  real*8 :: ZERO = 0.0D0

  if ((x>ZERO) .and. (x<ONE)) then
  logPdfBeta = (alpha-ONE)*log(x) + (beta-ONE)*log(ONE-x) + &
    dlgama(alpha+beta) - dlgama(alpha) - dlgama(beta)
  else
    logPdfBeta = -huge(ONE)
  endif

  end function

  ! Calculates the probability density function of the vector X of size 2
  ! of the truncated normal distribution.
  ! Original mean and covariance matrix are hard coded.
  ! Lower left bound is [bound; 1]
  ! Upper right bound is [1;,35]
  real*8 function pdfTruncatedNormal(X,bound)
  real*8, dimension(2) :: X
  real*8 :: bound

  real*8,dimension(2) :: lower
  real*8 :: partition
  real*8 :: abs_error, rel_error
  real*8 :: error
  integer :: nevals,inform

  real*8,dimension(2) :: mu
  real*8,dimension(2) :: upper
  real*8,dimension(2,2) :: COVRNC,CONSTR
  integer,dimension(2) :: INFIN
  real*8 :: s_X,s_Y,rho
  mu(1) = 0.0144D0
  mu(2) = 17.0D0
  upper(1) = 1.0D0
  upper(2) = 35.0D0
  COVRNC(1,1) = 0.000036D0
  COVRNC(2,1) = -0.0187D0
  COVRNC(1,2) = -0.0187D0
  COVRNC(2,2) = 16.09
  CONSTR = 0
  CONSTR(1,1) = 1.0
  CONSTR(2,2) = 1.0
  INFIN(1) = 2
  INFIN(2) = 2
  s_X = sqrt(COVRNC(1,1))
  s_Y = sqrt(COVRNC(2,2))
  rho = COVRNC(1,2)/(s_X*s_Y)


  abs_error = 1D-7
  rel_error = 1D-7

  lower(1) = bound
  lower(2) = -1.0D0

  ! Calculating partition function for truncated distribution
  CALL MVDIST(2,COVRNC,0,2,lower,CONSTR,upper,INFIN,mu,10000,abs_error,rel_error, &
    error, partition, nevals,inform)
    if (inform/=0) THEN
      print *,"Integration error"
      pdfTruncatedNormal = 0
      return
    ENDIF

    if ((X(1)<lower(1)).or.(X(1)>upper(1)).or.(X(2)<lower(2)).or.(X(2)>upper(2))) THEN
      pdfTruncatedNormal = 0.0D0
else

  pdfTruncatedNormal = (1.0D0/partition)*EXP((-0.5D0*(1-rho**2))*(((X(1)-mu(1))**2)/COVRNC(1,1) + ((X(2)-mu(2))**2)/COVRNC(2,2) &
    - 2.0D0*rho*(X(1)-mu(1))*(X(2)-mu(2))/(s_X*s_Y)))/(2.0D0*pi*s_X*s_Y*(sqrt((1.0D0-rho**2.0D0))))
	end if
end function

 real*8 function logTruncatedNormal(x,bound)
   real*8,dimension(2) :: x
   real*8 :: bound
   real*8 :: truncated

   truncated = pdfTruncatedNormal(x,bound)

   if (truncated > 0.0D0) then
     logTruncatedNormal = log(truncated)
   else
     logTruncatedNormal = -huge(0.0D0)
   endif

 end function


  ! Calculates the probability density function of a Dirichlet distribution
  ! for vector X of size 3 with parameter vector alpha (size 3)
  real*8 function pdfDirichlet(X,alpha)
    real*8, dimension(3) :: X,alpha
    real*8 :: beta
    real*8 :: one, zero
    one = 1.0D0
    zero = 0.0D0

    if ((zero<X(1)).and.(X(1)<one).and.(zero<X(2)).and.(X(2)<one).and. &
    (zero<X(3)).and.(X(3)<one).and.((one-X(1)-X(2)-X(3))<(100.0D0*EPSILON(X(1))))) then

      beta = PRODUCT(DGAMMA(alpha))/DGAMMA(SUM(alpha))

      pdfDirichlet = X(1)**(alpha(1)-1)*X(2)**(alpha(2)-1)*X(3)**(alpha(3)-1)/beta

    else
      pdfDirichlet = 0
    endif

  end function

  ! Calculates the log probability density function of a Dirichlet distribution
  ! for vector X of size 3 with parameter vector alpha (size 3)
  real*8 function logPdfDirichlet(X,alpha)
    real*8, dimension(3) :: X,alpha
    real*8 :: beta
    real*8 :: one, zero
    one = 1.0D0
    zero = 0.0D0

    if ((zero<X(1)).and.(X(1)<one).and.(zero<X(2)).and.(X(2)<one).and. &
    (zero<X(3)).and.(X(3)<one).and.((one-X(1)-X(2)-X(3))<(100.0D0*EPSILON(X(1))))) then

      logPdfDirichlet = (alpha(1)-one)*log(X(1)) + (alpha(2)-one)*log(X(2)) + (alpha(3)-one)*log(X(3)) &
        + dlgama(sum(alpha)) - dlgama(alpha(1)) - dlgama(alpha(2)) - dlgama(alpha(3))

    else
      logPdfDirichlet = -huge(one)
    endif

  end function
end module
