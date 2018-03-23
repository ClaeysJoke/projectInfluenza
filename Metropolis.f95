program Metropolis
  use Sampling
  use Reading
  use Prior
  implicit none

  integer :: t !Timewindow in weeks

  real*8, dimension(9) :: phi
  real*8, allocatable :: theta(:)
  real*8, allocatable :: y(:)
  real*8, allocatable :: X(:)
  real*8 :: pdfX
  integer :: count
  integer :: burnIn, sampleRate, nbSamples
  integer :: nbComponents
  integer :: nbTries
  integer :: i

  open(unit=8,file='/Users/Nick/Documents/2eMaster/Project/Code/projectInfluenza/BurnIn.txt')
  open(unit=9,file='/Users/Nick/Documents/2eMaster/Project/Code/projectInfluenza/Samples.txt')

  burnIn = 5000
  sampleRate = 1
  nbSamples = 10000
  t = 10

  allocate(y(t))
  allocate(theta(3*t))
  nbComponents = 9 + 3*t
  allocate(X(nbComponents))

  ! Read y from data
  call readData(y,t)


  ! Initialize parameter vectors
  print *, "Sampling Prior"
  call SamplePrior(X,y)
  print *, "Computing initial posterior distribution"
  call LogPosterior(X,y,pdfX)
  print *, "Posterior equals ", pdfX
  !print *, "Computing initial gradient"
  !call GradientFinite(X,y,pdfX,gradX)

  !do i = 1,size(X)
!    print *, "gradX(",i,") = ",gradX(i)
!  enddo





  ! MCMC Burn-in
 print *, "================== Burn-in start ====================="
  do count = 1, burnIn
    ! Generate new sample during burn-in
    call GenerateSample(X,pdfX,y,nbTries)
    !write(8,*) count,pdfX,nbTries
    write(8,*) X(1),X(2),pdfX,nbTries
    print *, "Burn-in sample",count
  enddo
  count = 1
  ! Generate real samples
  do count = 1, nbSamples
    !Generate samples
    call GenerateSample(X,pdfX,y,nbTries)
    if (MODULO(count,sampleRate)==0)then
      !write(9,*) count,pdfX,nbTries
      write(9,*) X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9)
      print *, "Sampling number",count
    endif
  enddo

  deallocate(y,theta,X)

end program
