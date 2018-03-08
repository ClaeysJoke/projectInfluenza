module MCMC
  use gradient ! Separate module to compute potention gradient
  implicit none

  save

  private

  integer, parameter      :: KREAL = kind(0.d0)
  integer, parameter      :: burnIn = 1000 !Really needed?
  integer, parameter      :: sampleRate = 200 !Use only every 200th sample
  real(KREAL), parameter  :: dt !Time step Euler-Maruyama


contains

  subroutine newSamplePosterior(kappa, lambda, theta0, &
    PI, PT, rho, beta, theta, y)
    !Create one single new sample, for use in bigger algorithm
    real(KREAL), intent(inout)  :: kappa, lamba, PI, PT, rho, beta
    real(KREAL), intent(inout)  :: theta0(3), theta(:,:)
    real(KREAL), intent(in)     :: y(:)
    real(KREAL)                 :: dKappa, dLambda, dPI, dPT, dRho, dBeta
    real(KREAL), dimension(3)   :: dTheta0
    real(KREAL), allocatable    :: dTheta(:,:)
    integer                     :: thetaSize
    real(KREAL)                 :: gamma, newGamma
    real(KREAL)                 :: newKappa,newLambda,newPI,newPT,newRho,newBeta
    real(KREAL), dimension(3)   :: newTheta0
    real(KREAL), allocatable    :: newTheta(:,:)

    gamma = beta*rho

    thetaSize = size(theta)
    allocate(dTheta(3,thetaSize))
    allocate(newTheta(3,thetaSize))

    !Compute gradient "vector"
    dKappa = gradKappa(kappa, theta0, rho, beta, theta)
    dLambda = gradLambda(lamba, theta, y)
    dTheta0 = gradTheta0(theta0, rho, beta, theta(:,1))
    dPI = gradPI(PI, PT, rho)
    dPT = gradPT(PI, PT, rho)
    drho = gradRho(kappa, beta, gamma, theta)
    dBeta = gradBeta(......)
    dTheta = gradTheta(..........)

    !Generate new sample
    newKappa = kappa - dKappa*dt + sqrt(2*dt)*StandardNormalSample()
    newLambda = lamba - dLambda*dt + sqrt(2*dt)*StandardNormalSample()
    newTheta0 = theta0
    newTheta0(2) = theta0(2) - dTheta0(2)*dt + sqrt(2*dt)*StandardNormalSample()
    newTheta0(3) = 1 - newTheta0(1) - newTheta0(2)
    newPI = PI - dPI*dt + sqrt(2*dt)*StandardNormalSample()
    newPT = PT - dPT*dt + sqrt(2*dt)*StandardNormalSample()
    newRho = rho - dRho*dt + sqrt(2*dt)*StandardNormalSample()
    newBeta = beta - dBeta*dt + sqrt(2*dt)*StandardNormalSample()
    newTheta???????

    newGamma = newBeta*newRho

    





    deallocate(dTheta)
    deallocate(newTheta)

  end subroutine






end module
