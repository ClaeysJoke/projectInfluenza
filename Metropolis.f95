program Metropolis

use Sampling

integer :: t !Timewindow in weeks
integer, parameter :: KREAL = (0.d0)

real(KREAL), dimension(9) :: phi
real(KREAL), allocatable :: theta(:)
real(KREAL), allocatable :: y
real(KREAL), allocatable :: X, gradX
real(KREAL), allocatable :: pdfX
integer :: count
integer :: burnIn, sampleRate, nbSamples
real(KREAL) :: dt
integer :: nbComponents


allocate(y(t))
allocate(theta(3*t))
nbComponents = 9 + 3*t
allocate(X(nbComponents), gradX(nbComponents))

! Read y from data
call readData(y,t)

! Initialize parameter vectors
call SamplePrior(X)
call Gradient(X,y,gradX)
call Posterior(X,y,pdfX)

! MCMC Burn-in
do count = 1, burnIn
  ! Generate new sample during burn-in
  call GenerateSample(X,gradX,pdfX,y,dt)
enddo

! Generate real samples
do count = 1, nbSamples
  !Generate samples
  call GenerateSample(X,gradX,pdfX,y,dt)
  if (MODULO(count,nbSamples)==0)then
    !Write to file, print out, ....
    WRITE
  endif
enddo

deallocate(y,theta,X,gradX)

end program
