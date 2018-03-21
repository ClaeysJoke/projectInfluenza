!source: http://www.aip.de/groups/soe/local/numres/bookfpdf/f6-2.pdf
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
if(x.lt.0..or.x.gt.1.)pause 'bad argument x in betai'
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
pause 'a or b too big, or MAXIT too small in betacf'
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
if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
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
pause 'a too large, ITMAX too small in gser'
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

pause 'a too large, ITMAX too small in gcf'
1 gammcf=exp(-x+a*log(x)-gln)*h
return
END subroutine

program TEST_DISTRIBUTIONS
use drawinvbay
use DIFF_SOLVER
use DIFF
integer :: DRAWS, BURNS,NUMBEROFBINS,PRUNEDDRAWS,clock_1,clock_2,count, count_rate, count_max,i,j
real(kind=selected_real_kind(8)) :: r,shape,scale,mean,PERCENTILE,WIDTH,temp,EXPECT,chi2,K_alpha,K,K_alpha2,p
real(kind=selected_real_kind(8)) :: mu, sigma,s2,var,sw,varw,s
real(kind=selected_real_kind(8)) :: gammp,pdfGamma,cdfNormal,cdfBeta,cdfStudent
real(kind=selected_real_kind(8)) :: a,b

real(kind=selected_real_kind(8)),allocatable :: y(:),BINS(:)
real(kind=selected_real_kind(8)),allocatable :: z(:),w(:)

DRAWS=100
BURNS=50
NUMBEROFBINS=20
allocate(y(DRAWS+BURNS))
allocate(BINS(NUMBEROFBINS))
allocate(z(DRAWS))
allocate(w(DRAWS))
K_alpha=1.94947/dsqrt(real(DRAWS,8))
temp=0.0
K=0.0
r=0.0
!********************************************
!				NORMAL TEST
!
!********************************************

mu=1.0
sigma=1.0
s2=sigma**2
chi2=0.0
var=0.0
open(unit=5,file="Gaussian.txt")
!TIME THE PROCEDURE
do i=1,DRAWS+BURNS
if(i.eq.BURNS) then
call SYSTEM_CLOCK(clock_1,count_rate,count_max)
end if
temp=rand_normal(mu,s2)
end do
call SYSTEM_CLOCK(clock_2,count_rate,count_max)
do j=1,10
i=1
do i=1,DRAWS
z(i)=rand_normal(mu,sigma)
w(i)=rand_normal(mu,sigma)
end do
s=sum(z)/DRAWS
sw=sum(w)/DRAWS
var=(1.0/(DRAWS-1))*sum((z-s)**2)
varw=(1.0/(DRAWS-1))*sum((z-s)**2)
r=sum((z-s)*(w-sw))/((DRAWS-1)*sqrt(var*varw))
p=p+2*(1.-cdfStudent(abs(r*dsqrt(real((DRAWS-2.0)/(1-r**2),8))),DRAWS-2))
end do
p=p/10
WIDTH=maxval(z)-minval(z)
WIDTH=WIDTH/NUMBEROFBINS
do i=1,NUMBEROFBINS
BINS(i)=COUNT(((i-1)*width+minval(z) .le. z).and.(z < i*WIDTH+minval(z)))
EXPECT=(cdfNormal((minval(z)+i*WIDTH),mu,sigma)-cdfNormal((minval(z)+(i-1)*WIDTH),mu,sigma))*DRAWS
temp=abs(Bins(i)-EXPECT)
 chi2=chi2+(temp**2)/(EXPECT)
end do
chi2=chi2/DRAWS

open(unit=8,file="TEST_OUTPUT.txt")
print*, cdfNormal(0.0_8,0.0_8,1.0_8)
write(8,*) "*******NORMAL DISTRIBUTION TEST*******"
write(8,*) "TIME: ", DRAWS, " DRAWS IN ", real(clock_2-clock_1,8)/count_rate,"s." 
write(8,*) "REL. ERR. ON MEAN: ", abs((mu-s)/mu)
write(8,*) "REL. ERR. ON VARIANCE: ", abs((var-s2)/s2)
write(8,*) "CHI SQ. VALUE: ",chi2,"/// P-VALUE: ",1.0-cdfNormal(chi2,real(DRAWS,8),sqrt(2*real(DRAWS,8)))
write(8,*)"PEARSON CORR. COEFF.: ",r,"/// P-VALUE: ", p
write(8,*) "*************************************"

!********************************************
!				GAMMA TEST
!
!********************************************
DRAWS=100
BURNS=50
shape=2.0
scale=1000.0
var=0.0
s=0.0
r=0.0
K=0.0
PERCENTILE=0.0
z=0.0
central=0.0
mean=shape*scale
chi2=0.0
s2=(shape*(scale**2))
clock_1=0
clock_2=0
!CALCULATE 99.9% PERCENTILE 
!CALCULATE SAMPLE MEAN AND TIME SAMPLING
!CALCULATE LEFT, CENTRAL AND RIGHT MEAN SQUARED ERROR
i=1
do i=1,DRAWS+BURNS
if (i .eq. BURNS) then
call SYSTEM_CLOCK(clock_1,count_rate,count_max) !start the clock after burn-in
end if
temp=rand_gamma(shape,scale)
end do
call SYSTEM_CLOCK(clock_2,count_rate,count_max)
p=0.0
do j=1,10
i=1
do while(i.le. DRAWS)
z(i)=rand_gamma(shape,scale)
w(i)=rand_gamma(shape,scale)
if((z(i).ne.z(i)).or.(w(i).ne.w(i)))then !exclude NaN exceptions
i=i-1
end if
i=i+1
end do
s=sum(z)/DRAWS
sw=sum(w)/DRAWS
var=(1.0/(DRAWS-1))*sum((z-s)**2)
varw=(1.0/(DRAWS-1))*sum((w-sw)**2)
r=sum((z-s)*(w-sw))/((DRAWS-1)*sqrt(var*varw))
p=p+2*(1.-cdfStudent(abs(r*dsqrt(real((DRAWS-2.0)/(1-r**2),8))),DRAWS-2))
end do
p=p/10
PERCENTILE=maxval(z)
WIDTH=PERCENTILE/NUMBEROFBINS
do i=1,NUMBEROFBINS
 BINS(i)=COUNT(((i-1)*WIDTH .le. z).and.(z < i*WIDTH))
 EXPECT=(cdfGAMMA((i)*WIDTH,shape,scale)-cdfGamma((i-1)*WIDTH,shape,scale))*DRAWS
 temp=abs(Bins(i)-EXPECT)
 K=max(K,temp/DRAWS)
 chi2=chi2+(temp**2)/(EXPECT)
end do
K=dsqrt(real(DRAWS,8))*K
chi2=chi2/DRAWS
write(8,*)"*******GAMMA DISTRIBUTION TEST*******"
write(8,*) "TIME: ", DRAWS, "DRAWS IN", real(clock_2-clock_1,8)/count_rate, "seconds"
write(8,*) "REL. ERROR ON MEAN ", abs(mean-s)/mean
write(8,*) "REL. ERR. ON VARIANCE: ", abs((var-s2)/s2)
write(8,*) "CHI SQ. VALUE: ",chi2,"/// P-VALUE: ",1.0-cdfNormal(chi2,real(DRAWS,8),sqrt(2*real(DRAWS,8)))
write(8,*)"PEARSON CORR. COEFF.: ",r,"/// P_VALUE ",p
write(8,*) "*************************************"

!********************************************
!				BETA TEST
!
!********************************************
DRAWS=100
BURNS=50
a=2.0
b=2.0
s=0.0
K=0.0
mean=a/(a+b)
s2=(a*b)/((a+b+1.0)*(a+b)**2)
chi2=0.0
clock_1=0
clock_2=0
!CALCULATE 99.9% PERCENTILE 
!CALCULATE SAMPLE MEAN AND TIME SAMPLING
!CALCULATE LEFT, CENTRAL AND RIGHT MEAN SQUARED ERROR
i=1
do i=1,DRAWS+BURNS
if (i .eq. BURNS) then
call SYSTEM_CLOCK(clock_1,count_rate,count_max) !start the clock after burn-in
end if
temp=rand_beta(a,b)
end do
call SYSTEM_CLOCK(clock_2,count_rate,count_max)
do j=1,10
i=1
do while(i.le. DRAWS)
z(i)=rand_beta(a,b)
w(i)=rand_beta(a,b)
if((z(i).ne.z(i)).or.(w(i).ne.w(i)))then !exclude NaN exceptions
i=i-1
end if
i=i+1
end do
s=sum(z)/DRAWS
sw=sum(w)/DRAWS
var=(1.0/(DRAWS-1))*sum((z-s)**2)
varw=(1.0/(DRAWS-1))*sum((w-sw)**2)
r=sum((z-s)*(w-sw))/((DRAWS-1)*sqrt(var*varw))
p=p+2*(1.-cdfStudent(abs(r*dsqrt(real((DRAWS-2.0)/(1-r**2),8))),DRAWS-2))
end do
p=p/10
PERCENTILE=1.0
WIDTH=PERCENTILE/NUMBEROFBINS
do i=1,NUMBEROFBINS
 BINS(i)=COUNT(((i-1)*WIDTH .le. z).and.(z < i*WIDTH))
 EXPECT=(cdfBeta((i)*WIDTH,a,b)-cdfBeta((i-1)*WIDTH,a,b))*DRAWS
 temp=abs(Bins(i)-EXPECT)
 K=max(K,temp/DRAWS)
 chi2=chi2+(temp**2)/(EXPECT)
end do
K=dsqrt(real(DRAWS,8))*K
chi2=chi2/DRAWS
write(8,*)"*******BETA DISTRIBUTION TEST*******"
write(8,*) "TIME: ", DRAWS, "DRAWS IN", real(clock_2-clock_1,8)/count_rate, "seconds"
write(8,*) "REL. ERROR ON MEAN ", abs(mean-s)/mean
write(8,*) "REL. ERR. ON VARIANCE: ", abs((var-s2)/s2)
write(8,*) "CHI SQ. VALUE: ",chi2,"/// P-VALUE: ",1.0-cdfNormal(chi2,real(DRAWS,8),sqrt(2*real(DRAWS,8)))
write(8,*)"PEARSON CORR. COEFF.: ",r,"/// P_VALUE ",p
write(8,*) "*************************************"
print*, "TEST COMPLETE"
deallocate(z)
deallocate(w)
deallocate(BINS)
deallocate(y)
end program
