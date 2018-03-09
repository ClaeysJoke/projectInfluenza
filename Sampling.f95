module Sampling
  use Normal
  implicit none
  real*8,parameter :: PI = 4*ATAN(1.0_8)

contains
  subroutine GenerateSample(X,gradX,pdfX,y,dt)
    real*8, intent(inout) :: X(:), gradX(:), pdfX
    real*8, intent(in) :: y(:)
    real*8, intent(in) :: dt
    real*8, allocatable :: trialX(:),gradTrialX(:)

    real*8 :: fromTrial, toTrial
    real*8 :: U
    real*8 :: pdfTrialX

    allocate(trialX(size(X)),gradTrialX(size(gradX)))

    do !Loop until acceptable sample is generated
    ! Generate trial
      call GenerateTrial(X,gradX,dt,trialX)
    ! Calculate trial gradient
      call Gradient(trialX,y,gradTrialX)

    ! Sample Uniform distribution for accept-reject
      call RANDOM_NUMBER(U)

    ! Call functions for accept reject
      call Posterior(trialX, y, pdfTrialX)
      fromTrial =  Transition(trialX,gradTrialX,dt,X)
      toTrial = Transition(X,gradX,dt,trialX)

    ! Accept reject
      !if (U <= pdfTrialX*fromTrial/(pdfX*toTrial)) then
        X = trialX
        gradX = gradTrialX
        pdfX = pdfTrialX
        exit
      !endif
    enddo

    deallocate(trialX,gradTrialX)
  end subroutine

  subroutine GenerateTrial(X,gradX,dt,trialX)
    real*8, intent(out) :: trialX(:)
    real*8, intent(in) :: X(:), gradX(:), dt
    integer :: i
    real*8, allocatable :: xi(:)

    allocate(xi(size(X)))

    do i = 1,size(xi)
      xi(i) = rand_normal(0.0D0,1.0D0)
    enddo

    trialX = X - gradX*dt + SQRT(2.0*dt)*xi

    deallocate(xi)
  end subroutine

  subroutine Gradient(X,y,gradX)
    real*8, intent(in) :: X(:),y(:)
    real*8, intent(out) :: gradX(:)
    gradX = -X
  end subroutine

  subroutine Posterior(X,y,pdfX)
    real*8, intent(in) :: X(:),y(:)
    real*8, intent(out) :: pdfX
    pdfX = 0.5
  end subroutine

  real*8 function Transition(X, gradX, dt, destination)
    real*8 :: X(:), gradX(:), dt, destination(:)
    real*8, allocatable :: mu(:)


    allocate(mu(size(X)))
    mu = destination - (X-gradX*dt)

    Transition = EXP(-0.5*(sqrt(2*dt)**(-1)*DOT_PRODUCT(destination-mu,destination-mu))) &
    /(sqrt((2*PI)**size(X)*sqrt(2*dt)**size(X)))

    deallocate(mu)
  end function

end module
