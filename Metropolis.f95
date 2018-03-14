program Metropolis
  use Sampling
  use Reading
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
  integer :: nbTries

  open(1,"BurnIn.txt")
  open(2,"Samples.txt")

  burnIn = 200
  sampleRate = 1
  nbSamples = 500
  t = 5
  dt = 1

  allocate(y(t))
  allocate(theta(3*t))
  nbComponents = 9 + 3*t
  allocate(X(nbComponents), gradX(nbComponents))

  ! Read y from data
  call readData(y,t)

  ! Initialize parameter vectors
  call SamplePrior(X,y)
  call Gradient(X,y,gradX)
  call Posterior(X,y,pdfX)

  ! MCMC Burn-in
  do count = 1, burnIn
    ! Generate new sample during burn-in
    call GenerateSample(X,gradX,pdfX,y,dt,nbTries)
    write(1,*) count,pdfX,nbTries
    print *, "Burn-in sample",count
  enddo
  count = 1
  ! Generate real samples
  do count = 1, nbSamples
    !Generate samples
    call GenerateSample(X,gradX,pdfX,y,dt,nbTries)
    if (MODULO(count,sampleRate)==0)then
      write(2,*) count,pdfX,nbTries
      print *, "Sampling number",count
    endif
  enddo

  deallocate(y,theta,X,gradX)

end program
