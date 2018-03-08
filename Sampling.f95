module Sampling
  implicit none
  integer, parameter :: KREAL = (0.d0)

contains
  subroutine GenerateSample(X,gradX,pdfX,y,dt)
    real(KREAL), intent(inout) :: X(:), gradX(:), pdfX
    real(KREAL), intent(in) :: y(:)
    real(KREAL), intent(in) :: dt

    real(KREAL) :: fromTrial, toTrial
    real(KREAL) :: U
    real(KREAL) :: pdfTrialX

    ! Generate trial
    call GenerateTrial(X,gradX,trialX)
    ! Calculate trial gradient
    call Gradient(trialX,y,gradTrialX)

    ! Sample Uniform distribution for accept-reject
    call RANDOM_NUMBER(U)

    ! Call functions for accept reject
    call Posterior(trialX, y, pdfTrialX)
    fromTrial =  Transition(trialX,gradTrialX,dt,X)
    toTrial = Transition(X,gradX,dt,trialX)

    ! Accept reject
    if (U <= pdfTrialX*fromTrial/(pdfX*toTrial)) then
      X = trialX
      gradX = gradTrialX
      pdfX = pdfTrialX
      exit
    endif
    enddo
  end subroutine

  subroutine GenerateTrial(X,gradX,trialX)
    real(KREAL), intent(out) :: trialX(:)
    real(KREAL), intent(in) :: X(:), gradX(:)

  end subroutine

  subroutine Gradient(X,y,gradX)
    real(KREAL), intent(in) :: X,y
    real(KREAL), intent(out) :: gradX

  end subroutine

  subroutine Posterior(X,y,pdfX)
    real(KREAL), intent(in) :: X,y
    real(KREAL), intent(out) :: pdfX

  end subroutine

  function Transition(X, gradX, dt, destination)
    real(KREAL) :: X, gradX, dt, destination

  end function

end module
