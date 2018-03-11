module Pdf
  implicit none



contains
  real*8 function pdfGamma(x)
    real*8 :: x
    real*8 :: alpha, beta

    alpha = 2.0D0
    beta = 0.0001D0

    if (x > 0) then
      pdfGamma = (beta**alpha)*(x**(alpha-1))*EXP(-beta*x)/DGAMMA(alpha)
    else
      pdfGamma = 0
    endif

  end function

  real*8 function pdfBeta(x,alpha,beta)
    real*8 :: x,alpha,beta
    real*8 :: betaFunction

    if ((x>0) .and. (x<1)) then
      betaFunction = DGAMMA(alpha)*DGAMMA(beta)/DGAMMA(alpha+beta)
      pdfBeta = x**(alpha-1)*(1-x)**(beta-1)/betaFunction
    else
      pdfBeta = 0
    endif
end function

real*8 function pdfTruncatedNormal(X,bound)
  real*8, dimension(2) :: X
  real*8 :: bound

  ????????????????????
end function



real*8 function pdfDirichlet(X,alpha)
  real*8, dimension(3) :: X,alpha
  real*8 :: beta

  if ((0 < X(1)) .and. (X(1) < 1) .and. (0 < X(2)) .and. (X(2) < 1) &
  (0 < X(3)) .and. (X(3) < 1) .and. (1-X(1)-X(2)-X(3) < 100*EPSILON(X(1)))) then
    beta = PRODUCT(GAMMA(alpha))/GAMMA(SUM(alpha))

    pdfDirichlet = X(1)**(alpha(1)-1)*X(2)**(alpha(2)-1)*X(3)**(alpha(3)-1)/beta

  else
    pdfDirichlet = 0
  endif

end function
end module
