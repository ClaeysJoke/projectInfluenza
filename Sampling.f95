module Sampling
  integer, parameter :: KREAL = (0.d0)

contains
  subroutine GenerateSample(X,gradX,pdfX,y)
    real(KREAL), intent(inout) :: X(:), gradX(:), pdfX
    real(KREAL), intent(in) :: y(:)
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
end module
