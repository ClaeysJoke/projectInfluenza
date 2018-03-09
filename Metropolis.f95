program Metropolis
  use Sampling
  implicit none

  integer :: t !Timewindow in weeks

  real*8, dimension(9) :: phi
  real*8, allocatable :: theta(:)
  real*8, allocatable :: y(:)
  real*8, allocatable :: X(:), gradX(:)
  real*8 :: pdfX
  integer :: count
  integer :: burnIn, sampleRate, nbSamples
  real*8 :: dt
  integer :: nbComponents

  burnIn = 5
  sampleRate = 1
  nbSamples = 50
  t = 5

  allocate(y(t))
  allocate(theta(3*t))
  nbComponents = 9 + 3*t
  allocate(X(nbComponents), gradX(nbComponents))

  ! Read y from data
  !call readData(y,t)

  y = 0

  ! Initialize parameter vectors
  !call SamplePrior(X)
  X = 0.5
  call Gradient(X,y,gradX)
  call Posterior(X,y,pdfX)

  ! MCMC Burn-in
  do count = 1, burnIn
    ! Generate new sample during burn-in
    call GenerateSample(X,gradX,pdfX,y,dt)
    print *, "Burn-in"
  enddo
  count = 1
  ! Generate real samples
  do count = 1, nbSamples
    !Generate samples
    call GenerateSample(X,gradX,pdfX,y,dt)
    if (MODULO(count,sampleRate)==0)then
      !Write to file, print out, ....
      print *, count
    endif
  enddo

  deallocate(y,theta,X,gradX)

end program
