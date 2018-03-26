! MAIN SOURCE FOR CDF'S: http://www.aip.de/groups/soe/local/numres/bookfpdf/f6-2.pdf

!********************************************
!		PART ONE: CDF's
!
!********************************************
function cdfStudent(t,N)
integer ::N
real(kind=selected_real_kind(8)) :: x,t
real(kind=selected_real_kind(8)) :: cdfStudent,cdfBeta
x=N/((t**2)+N)
cdfStudent=1-0.5*cdfBeta(x,real((N/2),8),0.5_8)
end function

FUNCTION cdfBeta(x,a,b)
REAL(kind=selected_real_kind(8)):: cdfBeta,a,b,x
REAL(kind=selected_real_kind(8)):: bt,betacf,gammln
if(x.lt.0..or.x.gt.1.)stop 'bad argument x in betai'
if(x.eq.0..or.x.eq.1.)then
bt=0.
else
bt=exp(DLGAMA(a+b)-DLGAMA(a)-DLGAMA(b)+a*log(x)+b*log(1.-x))
endif
if(x.lt.(a+1.)/(a+b+2.))then 
cdfBeta=bt*betacf(a,b,x)/a
return
else
cdfBeta=1.-bt*betacf(b,a,1.-x)/b
endif
END FUNCTION

FUNCTION betacf(a,b,x)
INTEGER MAXIT
REAL(kind=selected_real_kind(8)) :: betacf,a,b,x,EPS,FPMIN
PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
INTEGER m,m2
REAL aa,c,d,del,h,qab,qam,qap
qab=a+b 
qap=a+1.
qam=a-1.
c=1.
d=1.-qab*x/qap
if(abs(d).lt.FPMIN)d=FPMIN
d=1./d
h=d
do m=1,MAXIT
m2=2*m
aa=m*(b-m)*x/((qam+m2)*(a+m2))
d=1.+aa*d
if(abs(d).lt.FPMIN)d=FPMIN
c=1.+aa/c
if(abs(c).lt.FPMIN)c=FPMIN
d=1./d
h=h*d*c
aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
d=1.+aa*d
if(abs(d).lt.FPMIN)d=FPMIN
c=1.+aa/c
if(abs(c).lt.FPMIN)c=FPMIN
d=1./d
del=d*c
h=h*del
if(abs(del-1.).lt.EPS)goto 1
enddo
stop 'a or b too big, or MAXIT too small in betacf'
1 betacf=h
return
END FUNCTION
function cdfNormal(x,mu,sigma)
real (kind=selected_real_kind(8)) :: cdfNormal, mu, sigma,x
real (kind=selected_real_kind(8)),parameter :: pi=4*ATAN(real(1.0,8))
x=(x-mu)/sigma
x=x/dsqrt(2.0_8)
cdfNormal=0.5*(1+ERF(x))
end function

real*8 function pdfGamma(x,alpha,beta)!shape and rate
    real*8 :: x
    real*8 :: alpha, beta

    !alpha = 2.0D0
    !beta = 0.0001D0

    if (x > 0) then
      pdfGamma = (beta**alpha)*(x**(alpha-1))*EXP(-beta*x)/DGAMMA(real(alpha,8))
    else
      pdfGamma = 0.
    endif

  end function

  real*8 function pdfBeta(x,alpha,beta)
    real*8 :: x,alpha,beta
    real*8 :: betaFunction

    if ((x>0) .and. (x<1)) then
      betaFunction = DGAMMA(real(alpha,8))*DGAMMA(real(beta,8))/DGAMMA(real(alpha+beta,8))
      pdfBeta = x**(alpha-1)*(1-x)**(beta-1)/betaFunction
    else
      pdfBeta = 0
    endif
end function

function cdfGamma(x, alpha, theta)
real(kind=selected_real_kind(8))::x,alpha,theta,gammp
cdfGamma=gammp(alpha,x/theta)
end function
FUNCTION gammp(a,x)
real(kind=selected_real_kind(8)) a,gammp,x
real(kind=selected_real_kind(8)) gammcf,gamser,gln
if(x.lt.0..or.a.le.0.)stop 'bad arguments in gammp'
if(x.lt.a+1.)then
call gser(gamser,a,x,gln)
gammp=gamser
else
call gcf(gammcf,a,x,gln)
gammp=1.-gammcf
endif
return
END function

SUBROUTINE gser(gamser,a,x,gln)
INTEGER ITMAX
real(kind=selected_real_kind(8)) a,gamser,gln,x,EPS
PARAMETER (ITMAX=1000,EPS=3.e-7)
INTEGER n
real(kind=selected_real_kind(8)) ap,del,sum,gammln
gln=DLGAMA(a)
if(x.le.0.)then
if(x.lt.0.)stop 'x < 0 in gser'
gamser=0.
return
endif
ap=a
sum=1./a
del=sum
do n=1,ITMAX
ap=ap+1.
del=del*x/ap
sum=sum+del
if(abs(del).lt.abs(sum)*EPS)goto 1
enddo
stop 'a too large, ITMAX too small in gser'
1 gamser=sum*exp(-x+a*log(x)-gln)
return
end subroutine

SUBROUTINE gcf(gammcf,a,x,gln)
INTEGER ITMAX
real(kind=selected_real_kind(8)) a,gammcf,gln,x,EPS,FPMIN
PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
INTEGER i
real(kind=selected_real_kind(8)) an,b,c,d,del,h,gammln
gln=DLGAMA(a)
b=x+1.-a
c=1./FPMIN
d=1./b
h=d

do i=1,ITMAX
an=-i*(i-a)
b=b+2.
d=an*d+b
if(abs(d).lt.FPMIN)d=FPMIN
c=b+an/c
if(abs(c).lt.FPMIN)c=FPMIN
d=1./d
del=d*c
h=h*del
if(abs(del-1.).lt.EPS)goto 1
enddo

stop 'a too large, ITMAX too small in gcf'
1 gammcf=exp(-x+a*log(x)-gln)*h
return
END subroutine


!********************************************
!		PART TWO: MAIN TESTS
!
!********************************************
!program TEST_DISTRIBUTIONS
!use drawinvbay
!use DIFF_SOLVER
!use DIFF
!use DBSSM_MOD
!variables for sampler tests
!integer :: DRAWS, BURNS,NUMBEROFBINS,PRUNEDDRAWS,clock_1,clock_2,count, count_rate, count_max,i,j
!real(kind=selected_real_kind(8)) :: r,shape,scale,mean,PERCENTILE,WIDTH,temp,EXPECT,chi2,p
!real(kind=selected_real_kind(8)) :: mu, sigma,s2,var,sw,varw,s
!real(kind=selected_real_kind(8)) :: gammp,pdfGamma,cdfNormal,cdfBeta,cdfStudent
!real(kind=selected_real_kind(8)) :: a,b
!integer :: N_NORMAL,N_GAMMA,N_BETA
!real(kind=selected_real_kind(8)),allocatable :: BINS(:)
!real(kind=selected_real_kind(8)),allocatable :: z(:),w(:)

!variables for RK test & DBSSM test
!real(kind=selected_real_kind(8)),allocatable :: theta(:,:),o(:),y(:)
!real(kind=selected_real_kind(8)) :: theta0(3),beta, kappa,lambda, u(3),u0(3),dt
!real(kind=selected_real_kind(8)) :: gamm, t0,dttest,y_and_theta(4,100)
!real(kind=selected_real_kind(8)) :: y_and_theta_best(4,100)
!real(kind=selected_real_kind(8)) :: S0,I0,R0,rho,minmiss,currentmiss
!integer :: K,m
!K=100
!allocate(y(K))
!allocate(theta(3,K))
!open(unit=8,file="TEST_OUTPUT.txt")
!N_NORMAL=0
!N_GAMMA=0
!N_BETA=0
!do j=1,100
!DRAWS=100
!NUMBEROFBINS=20
!allocate(y(DRAWS+BURNS))
!allocate(BINS(NUMBEROFBINS))
!allocate(z(DRAWS))
!allocate(w(DRAWS))
!K_alpha=1.94947/dsqrt(real(DRAWS,8))
!temp=0.0
!K=0.0
!r=0.0
!********************************************
!		NORMAL TEST
!
!********************************************

!mu=1.0
!sigma=1.0
!s2=sigma**2
!chi2=0.0
!var=0.0
!p=0.0
!open(unit=4,file="PVAL.txt")
!TIME THE PROCEDURE
!do i=1,DRAWS+BURNS
!if(i.eq.BURNS) then
!call SYSTEM_CLOCK(clock_1,count_rate,count_max)
!end if
!temp=rand_normal(mu,s2)
!end do
!call SYSTEM_CLOCK(clock_2,count_rate,count_max)
!i=1
!do i=1,DRAWS
!z(i)=rand_normal(mu,sigma)
!w(i)=rand_normal(mu,sigma)
!end do
!s=sum(z)/DRAWS
!sw=sum(w)/DRAWS
!var=(1.0/(DRAWS-1))*sum((z-s)**2)
!varw=(1.0/(DRAWS-1))*sum((z-s)**2)
!r=sum((z-s)*(w-sw))/((DRAWS-1)*sqrt(var*varw))
!p=2*(1.-cdfStudent(abs(r*dsqrt(real((DRAWS-2.0)/(1-r**2),8))),DRAWS-2))
!WIDTH=maxval(z)-minval(z)
!WIDTH=WIDTH/NUMBEROFBINS
!do i=1,NUMBEROFBINS
!BINS(i)=COUNT(((i-1)*width+minval(z) .le. z).and.(z < i*WIDTH+minval(z)))
!EXPECT=(cdfNormal((minval(z)+i*WIDTH),mu,sigma)-cdfNormal((minval(z)+(i-1)*WIDTH),mu,sigma))*DRAWS
!temp=abs(Bins(i)-EXPECT)
! chi2=chi2+(temp**2)/(EXPECT)
!end do
!chi2=chi2/DRAWS

!if (p<0.05)then
!N_NORMAL=N_NORMAL+1
!end if

!write(8,*) "*******NORMAL DISTRIBUTION TEST*******"
!write(8,*) "TIME: ", DRAWS, " DRAWS IN ", real(clock_2-clock_1,8)/count_rate,"s." 
!write(8,*) "REL. ERR. ON MEAN: ", abs((mu-s)/mu)
!write(8,*) "REL. ERR. ON VARIANCE: ", abs((var-s2)/s2)
!write(8,*) "CHI SQ. VALUE: ",chi2,"/// P-VALUE: ",1.0-cdfNormal(chi2,real(DRAWS,8),sqrt(2*real(DRAWS,8)))
!write(8,*)"PEARSON CORR. COEFF.: ",r,"/// P-VALUE: ", p
!write(8,*) "*************************************"

!********************************************
!		  GAMMA TEST
!
!********************************************
!do while (r .ge. 0)
!r=rand_gamma(0.000001_8,1.0_8)
!print*, r
!end do
!DRAWS=100
!BURNS=50
!shape=2.0
!scale=1000.0
!var=0.0
!s=0.0
!r=0.0
!K=0.0
!PERCENTILE=0.0
!z=0.0
!central=0.0
!mean=shape*scale
!chi2=0.0
!p=0.0
!s2=(shape*(scale**2))
!clock_1=0
!clock_2=0
!CALCULATE 99.9% PERCENTILE 
!CALCULATE SAMPLE MEAN AND TIME SAMPLING
!CALCULATE LEFT, CENTRAL AND RIGHT MEAN SQUARED ERROR
!i=1
!do i=1,DRAWS+BURNS
!if (i .eq. BURNS) then
!call SYSTEM_CLOCK(clock_1,count_rate,count_max) !start the clock after burn-in
!end if
!temp=rand_gamma(shape,scale)
!end do
!call SYSTEM_CLOCK(clock_2,count_rate,count_max)
!p=0.0
!i=1
!do while(i.le. DRAWS)
!z(i)=rand_gamma(shape,scale)
!w(i)=rand_gamma(shape,scale)
!if((z(i).ne.z(i)).or.(w(i).ne.w(i)))then !exclude NaN exceptions
!i=i-1
!end if
!i=i+1
!end do
!s=sum(z)/DRAWS
!sw=sum(w)/DRAWS
!var=(1.0/(DRAWS-1))*sum((z-s)**2)
!varw=(1.0/(DRAWS-1))*sum((w-sw)**2)
!r=sum((z-s)*(w-sw))/((DRAWS-1)*sqrt(var*varw))
!p=2*(1.-cdfStudent(abs(r*dsqrt(real((DRAWS-2.0)/(1-r**2),8))),DRAWS-2))
!PERCENTILE=maxval(z)
!WIDTH=PERCENTILE/NUMBEROFBINS
!do i=1,NUMBEROFBINS
! BINS(i)=COUNT(((i-1)*WIDTH .le. z).and.(z < i*WIDTH))
! EXPECT=(cdfGAMMA((i)*WIDTH,shape,scale)-cdfGamma((i-1)*WIDTH,shape,scale))*DRAWS
! temp=abs(Bins(i)-EXPECT)
! chi2=chi2+(temp**2)/(EXPECT)
!end do
!chi2=chi2/DRAWS
!if (p<0.05)then
!N_GAMMA=N_GAMMA+1
!end if
!write(8,*)"*******GAMMA DISTRIBUTION TEST*******"
!write(8,*) "TIME: ", DRAWS, "DRAWS IN", real(clock_2-clock_1,8)/count_rate, "seconds"
!write(8,*) "REL. ERROR ON MEAN ", abs(mean-s)/mean
!write(8,*) "REL. ERR. ON VARIANCE: ", abs((var-s2)/s2)
!write(8,*) "CHI SQ. VALUE: ",chi2,"/// P-VALUE: ",1.0-cdfNormal(chi2,real(DRAWS,8),sqrt(2*real(DRAWS,8)))
!write(8,*)"PEARSON CORR. COEFF.: ",r,"/// P_VALUE ",p
!write(8,*) "*************************************"

!********************************************
!		  BETA TEST
!
!********************************************
!DRAWS=100
!BURNS=50
!a=2.0
!b=2.0
!s=0.0
!p=0.0
!mean=a/(a+b)
!s2=(a*b)/((a+b+1.0)*(a+b)**2)
!chi2=0.0
!clock_1=0
!clock_2=0
!CALCULATE 99.9% PERCENTILE 
!CALCULATE SAMPLE MEAN AND TIME SAMPLING
!CALCULATE LEFT, CENTRAL AND RIGHT MEAN SQUARED ERROR
!i=1
!do i=1,DRAWS+BURNS
!if (i .eq. BURNS) then
!call SYSTEM_CLOCK(clock_1,count_rate,count_max) !start the clock after burn-in
!end if
!temp=rand_beta(a,b)
!end do
!call SYSTEM_CLOCK(clock_2,count_rate,count_max)
!i=1
!do while(i.le. DRAWS)
!z(i)=rand_beta(a,b)
!w(i)=rand_beta(a,b)
!if((z(i).ne.z(i)).or.(w(i).ne.w(i)))then !exclude NaN exceptions
!i=i-1
!end if
!i=i+1
!end do
!s=sum(z)/DRAWS
!sw=sum(w)/DRAWS
!var=(1.0/(DRAWS-1))*sum((z-s)**2)
!varw=(1.0/(DRAWS-1))*sum((w-sw)**2)
!r=sum((z-s)*(w-sw))/((DRAWS-1)*sqrt(var*varw))
!p=2*(1.-cdfStudent(abs(r*dsqrt(real((DRAWS-2.0)/(1-r**2),8))),DRAWS-2))
!PERCENTILE=1.0
!WIDTH=PERCENTILE/NUMBEROFBINS
!do i=1,NUMBEROFBINS
! BINS(i)=COUNT(((i-1)*WIDTH .le. z).and.(z < i*WIDTH))
! EXPECT=(cdfBeta((i)*WIDTH,a,b)-cdfBeta((i-1)*WIDTH,a,b))*DRAWS
! temp=abs(Bins(i)-EXPECT)
! chi2=chi2+(temp**2)/(EXPECT)
!end do
!chi2=chi2/DRAWS
!if (p<0.05)then
!N_BETA=N_BETA+1
!end if
!write(8,*)"*******BETA DISTRIBUTION TEST*******"
!write(8,*) "TIME: ", DRAWS, "DRAWS IN", real(clock_2-clock_1,8)/count_rate, "seconds"
!write(8,*) "REL. ERROR ON MEAN ", abs(mean-s)/mean
!write(8,*) "REL. ERR. ON VARIANCE: ", abs((var-s2)/s2)
!write(8,*) "CHI SQ. VALUE: ",chi2,"/// P-VALUE: ",1.0-cdfNormal(chi2,real(DRAWS,8),sqrt(2*real(DRAWS,8)))
!write(8,*)"PEARSON CORR. COEFF.: ",r,"/// P_VALUE ",p
!write(8,*) "*************************************"
!print*, "TEST COMPLETE"
!deallocate(z)
!deallocate(w)
!deallocate(BINS)
!deallocate(y)

!end do
!write(4,*) "NORMAL CORRELATED %",real(N_NORMAL,8)/(100*DRAWS)
!write(4,*) "GAMMA CORRELATED %",real(N_GAMMA,8)/(100*DRAWS)
!write(4,*) "BETA CORRELATED %",real(N_BETA,8)/(100*DRAWS)


!********************************************
!		   LINPROG TEST
!
!********************************************


!********************************************
!		    RK TEST
!
!********************************************
!t0=0.0
!m=3
!u0=(/1.0,1.0,1.0/)
!beta=0.0
!gamm=0.0
!dt=1
!u=u0
!do i=1,10
!call rk4Vec (m, t0, u, dt, fvec, u,beta,gamm)
!end do
!print*, dsqrt(DOT_PRODUCT(u-dexp(10.0_8),u-dexp(10.0_8)))
!dt=0.5
!u=u0
!do i=1,20
!call rk4Vec (m, t0, u, dt, fvec, u,beta,gamm)
!end do
!print*, dsqrt(DOT_PRODUCT(u-dexp(10.0_8),u-dexp(10.0_8)))
!dt=0.25
!u=u0
!do i=1,40
!call rk4Vec (m, t0, u, dt, fvec, u,beta,gamm)
!end do
!print*, dsqrt(DOT_PRODUCT(u-dexp(10.0_8),u-dexp(10.0_8)))
!dt=0.125
!u=u0
!do i=1,80
!call rk4Vec (m, t0, u, dt, fvec, u,beta,gamm)
!end do
!print*, dsqrt(DOT_PRODUCT(u-dexp(10.0_8),u-dexp(10.0_8)))
!dt=dt/2
!u=u0
!do i=1,160
!call rk4Vec (m, t0, u, dt, fvec, u,beta,gamm)
!end do
!print*, dsqrt(DOT_PRODUCT(u-dexp(10.0_8),u-dexp(10.0_8)))
!dt=dt/2
!u=u0
!do i=1,320
!call rk4Vec (m, t0, u, dt, fvec, u,beta,gamm)
!end do
!print*, dsqrt(DOT_PRODUCT(u-dexp(10.0_8),u-dexp(10.0_8)))

!********************************************
!		    DBSSM TEST
!
!********************************************
!DRAWS=1000
!do i=1,DRAWS
!if (i .eq. 200) then
!call SYSTEM_CLOCK(clock_1,count_rate,count_max)
!end if


!open(unit=1,file='I0.txt')
!open(unit=2,file='R0.txt')
!open(unit=3,file='kappa.txt')
!open(unit=4,file='lambda.txt')
!open(unit=5,file='betas.txt')
!open(unit=7,file='rho.txt') !gamm
!minmiss = huge(0.0)
!call SYSTEM_CLOCK(clock_1,count_rate,count_max)
!do j=1,50
!print*,j
!t0=0.0
!I0=0.0995
!S0 = 0.9
!R0=1-(I0+S0)
!theta0(1) = S0
!theta0(2) = I0
!theta0(3) = R0
!beta=10.0
!kappa=60000.0
!lambda=60000.0
!rho=0.9
!gamm = rho*beta
!dt=.1
!y_and_theta=0.0
!print*,theta0

!print*, y_and_theta(:,j)
!do i = 1,10
!	call DBSSM(K,dt, t0, beta,gamm,kappa,lambda,theta0,theta,y)
!	print*, theta
!	print *, y
!	print*, y_and_theta
	!currentmiss = calc_likelihood(y_and_theta)

	!if (currentmiss .LT. minmiss) then
	!	minmiss = currentmiss
	!	y_and_theta_best(:,:) = y_and_theta(:,:)
	!endif
!print*, i, "DONE"
!enddo

!print *,minmiss
!print*, beta,gamm,R0
!print*, find_u(10.0_8,0.001_8,R0,I0,beta,gamm,0.01_8)
!print*,DEXP(-1*(beta/gamm)*R0)
!enddo

!call SYSTEM_CLOCK(clock_2,count_rate,count_max)

!print*, real(clock_2-clock_1,8)/count_rate


!(K,dt, t0, beta,gamm,kappa,lambda,theta0,theta,y)
!end do
!call SYSTEM_CLOCK(clock_2,count_rate,count_max)
!write(8,*) "DBSSM TIME: ", DRAWS,"CALLS IN", real(clock_2-clock_1,8)/count_rate, " seconds"
!end program

program main
use drawinvbay
use DIFF_SOLVER
use DIFF
use DBSSM_MOD
implicit none
real*8 :: S0,dt,I0,R0,kappa,lambda,beta,gamm,t0,rho,minmiss,currentmiss,mean
real*8, dimension(3) :: theta0
real*8, dimension(:,:), allocatable :: y_and_theta, y_and_theta_best
real*8, dimension(:,:), allocatable::y_mean
integer :: K,i,j,NumberOfSteps !NumberOfSteps is the number of steps in one week
NumberOfSteps=2
dt=.5
K =40 !K is number of weeks
t0 = 0.0 !voorlopig
allocate(y_and_theta(4,NumberOfSteps*K))
allocate(y_and_theta_best(4,NumberOfSteps*K))
allocate(y_mean(4,NumberOfSteps*K))
open(unit=1,file='I0.txt')
open(unit=2,file='R0.txt')
open(unit=3,file='kappa.txt')
open(unit=4,file='lambda.txt')
open(unit=5,file='betas.txt')
open(unit=7,file='rho.txt') !gamm
open(unit=11,file='results1.txt')
open(unit=12,file='results2.txt')
open(unit=13,file='results3.txt')
open(unit=14,file='results4.txt')
open(unit=9,file='thetaS.txt')
minmiss = huge(0.0)
currentmiss=0.0
do j=1,1
print*,j

read(1,*) I0
print *, I0
read(2,*) R0
print*, R0
S0 = 0.9
theta0(1) = S0
theta0(2) = I0
theta0(3) = R0
print *, "THETA0 is: ",theta0
print *, "SUM IS: ",SUM(theta0)
read(3,*) kappa
read(4,*) lambda
read(5,*) beta
read(7,*) rho
gamm = .6!rho*beta
beta=2.5
kappa=20000
lambda=20000
y_mean=0.0
do i = 1,50000
	print*,i
	y_and_theta = DBSSM_FUNC(NumberOfSteps*K,dt, t0, beta,gamm,kappa,lambda,theta0)
	y_mean=y_mean+y_and_theta
	currentmiss = calc_likelihood(y_and_theta,NumberOfSteps,K)
   print *,"CURRENTMISS: ", currentmiss
	if (currentmiss .LT. minmiss) then
		minmiss = currentmiss
		y_and_theta_best(:,:) = y_and_theta(:,:)
	endif

enddo
y_mean=y_mean/50000
print *,"MINMISS: ",minmiss	

enddo
!print*,'part1'
!write(9,*) y_and_theta_best(1,:)
!print*,'part2'
!print*,y_and_theta_best(2,:)
!print*,'part3'
!print*,y_and_theta_best(3,:)
!print*,'part4'
write(11,*) y_mean(1,:)
write(12,*) y_mean(2,:)
write(13,*) y_mean(3,:)
write(14,*) y_mean(4,:)

!print *, minmiss

	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	close(7)
	close(8)
	close(9)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

function calc_likelihood(y_and_theta,NumberOfSteps,K) result (miss)
	implicit none
	integer :: K
	real*8 :: dt
	real*8 :: y_and_theta(4,K)
	real*8 :: miss
	integer :: i,NumberOfSteps
	real*8, dimension(17) :: y_week
	y_week = (/ 21.752,13.566,20.835,10.612,15.355,21.55,26.911,33.207,41.233,49.664,58.386,67.528,85.89,118.46,193.18,290.95,391.85 /)
	miss = 0.0
	do i=1,16 !enkel eerste 16 punten vergelijken
		miss = miss + (abs(y_and_theta(4,i+NumberOfSteps*(i-1)) - y_week(i+1))/y_week(i+1))
	enddo
end function

!maximale fout minimaliseren, dus deze methode berekend maximale miss
function calc_likelihood2(y_and_theta) result (maxmiss)
	implicit none
	real*8, dimension(4,K) :: y_and_theta
	real*8 :: maxmiss,miss
	integer :: i
	real*8, dimension(17) :: y_week
	y_week = (/ 21.752,13.566,20.835,10.612,15.355,21.55,26.911,33.207,41.233,49.664,58.386,67.528,85.89,118.46,193.18,290.95,391.85 /)
	maxmiss = 0.0
	do i=1,16 !enkel eerste 16 punten vergelijken
		miss = abs(y_and_theta(4,i) - y_week(i+1))/y_week(i+1)
		if (miss > maxmiss) then
			maxmiss = miss
		endif
	enddo
end function



end program

