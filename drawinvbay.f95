PROGRAM drawinvbay
use normaldistribution
  real (kind = selected_real_kind(8) ) :: kappa, lambda, S0, I0,R0
  real,dimension(2) :: mu,lower,upper
  real,dimension(2,2) :: sigma, SigmaInv
  real,dimension(4,2) :: A
  real,dimension(4) :: b
  integer :: N
  
  mu(1) = 0.0144
  mu(2) = 17.9
  
  lower(2)=1
  upper(1)=1
  upper(2)=35
  
  A = reshape((/1, 0, -1, 0, 0, 1, 0, -1/), shape(A))
  sigma = reshape((/0.000036, -0.0187, -0.0187, 16.09/), shape(sigma))
  SigmaInv = reshape((/70093.66, 81.46373,81.46373, 0.1568285/), shape(SigmaInv))
  
  kappa = rand_gamma(2.0_8,0.0001_8)
  lambda = rand_gamma(2.0_8,0.0001_8)
  S0 = 0.9
  I0 = rand_beta(1.62,7084.10) !these values by fitting a beta distribution to historical ILI+ for t=0
  R0 = 1-S0-I0
  
  lower(1) = I0
  b(1) = upper(1)
  b(2) = upper(2)
  b(3) = -lower(1)
  b(4) = -lower(2)
  
  N=1 !how many samples do you want
  
  rows_A = SIZE(A,1)
  columns_A = SIZE(A,2)
  rows_b = SIZE(b,1)
  z = rand_truncated_normal(mu,sigma,SigmaInv,N,A,b)
  
 
contains

!OK
RECURSIVE FUNCTION rand_gamma(shape, scale) RESULT(ans)
	implicit none
	real*8 :: ans
	real*8, intent(in) :: shape, scale
      	real*8 :: u,w,d,c,x,xsq,g,v
      IF (shape <= 0.0) THEN

        WRITE(*,*) "Shape PARAMETER must be positive"
      END IF
      IF (scale <= 0.0) THEN

        WRITE(*,*) "Scale PARAMETER must be positive"
      END IF
      
      IF (shape >= 1.0) THEN
        d = shape - 1.0/3.0
        c = 1.0/(9.0*d)**0.5
        DO while (.true.)
            x = rand_normal(0.0_8, 1.0_8)
            v = 1.0 + c*x
            DO while (v <= 0.0)
                x = rand_normal(0.0_8, 1.0_8)
                v = 1.0 + c*x
            END DO

            v = v*v*v
            CALL RANDOM_NUMBER(u)
            xsq = x*x
            IF ((u < 1.0 -.0331*xsq*xsq) .OR.  &
              (log(u) < 0.5*xsq + d*(1.0 - v + log(v))) )then
                ans=scale*d*v
                RETURN
            END IF

        END DO
      ELSE
        g = rand_gamma(shape+1.0, 1.0_8)
        CALL RANDOM_NUMBER(w)
        ans=scale*g*(w)**(1.0/shape)
        RETURN
      END IF
END FUNCTION

!OK
FUNCTION rand_beta(a, b) RESULT(ans)
	implicit none
	real*8, intent(in) :: a,b
	real*8 :: ans
      	real*8 :: u,v
      IF ((a <= 0.0) .OR. (b <= 0.0)) THEN

        WRITE(*,*) "Beta PARAMETERs must be positive"
      END IF
      
       u = rand_gamma(a, 1.0_8)
       v = rand_gamma(b, 1.0_8)
       ans = u / (u + v)
END FUNCTION

!OK
FUNCTION rand_truncated_normal(mu,sigma,SigmaInv,N,Atemp,btemp) RESULT(r) !MULTIVARIATE
	implicit none
	real*8,dimension(2,2),intent(in) :: sigma,SigmaInv
	real*8,dimension(4,2),intent(in) :: Atemp
	real*8,dimension(4,1),intent(in) :: btemp
	real*8,dimension(2,4) :: A
	real*8,dimension(1,4) :: b
	real*8,dimension(1,2),intent(in) :: mu
	real*8 :: r
	integer,intent(in) :: N
	real*8 :: trials,passes,s,maxsample,rhoThr,defaultRhoThr,rho,ngibbs,q
	integer :: nar,p,m
	real*8,dimension(2) :: range_lb_ub
	real*8,dimension(:),allocatable::c,c_temp,Ai
	real*8,dimension(:,:),allocatable::c1,c2
	real*8,dimension(:),allocatable::Sigmai_i,SigmaInv_i,x_i,mu_i
	REAL*8,DIMENSION(:,:),ALLOCATABLE::xf,x,Xbig,sigmaInv_i_i,A_i,Sigma_i_iInv
	real*8,dimension(4,1) :: Atheta,Btheta
	real*8 :: ftheta
	integer :: discard
	integer,parameter :: CMAX=10, VMAX=10
	integer :: NC,NV,XERR
	integer :: i
	real*8 :: lb,ub,mui,s2i
	real*8 :: theta
	real*8 TS(CMAX,VMAX)
	defaultRhoThr = 0.00029
	rhoThr = defaultRhoThr

	A = transpose(Atemp) !2,4
	b = transpose(btemp) !1,4
	

	p = SIZE(mu) !2
	m = SIZE(A,2) !4
	
	ALLOCATE(xf(p,1))
	ALLOCATE(Xbig(N,p))
	ALLOCATE(x(1,p)) !1,2
	ALLOCATE(Sigmai_i(p-1))
	ALLOCATE(SigmaInv_i_i(p-1,p-1))
	ALLOCATE(Sigma_i_iInv(p-1,p-1))
	ALLOCATE(SigmaInv_i(p-1))
	ALLOCATE(x_i(p-1))
	ALLOCATE(mu_i(p-1))
	ALLOCATE(A_i(p-1,m))
	ALLOCATE(c(m))
	ALLOCATE(c_temp(m))
	ALLOCATE(Ai(m))
	Xbig = 0.0d0
	nar = 0
	ngibbs = 0.0
	rho = 1.0

	IF (nar.LT.N) THEN !nog waarden gevraagd
		if (nar.GT.0) THEN
			x(1,:)=Xbig(nar,:) !1,p
			discard = 0
		ELSE IF ( ALL(MATMUL(transpose(A),transpose(mu)).LE.transpose(b)) ) THEN
			x(1,:) = mu(1,:) !1,2
			discard = 1
		ELSE
			xf = chebycenter(transpose(A),transpose(b),1.0_8)
				IF (.NOT. (ALL(MATMUL(transpose(A),xf) .LE. transpose(b))) ) THEN
					print *,'error: failed to find a feasible point'
				END IF
			Atheta = -MATMUL(transpose(A),xf-transpose(mu))
			btheta = transpose(b) - MATMUL(transpose(A),xf)
			
			allocate(c1(4+2,1))
			c1(1:4,1) = Atheta
			c1(4+1,1) = 1
			c1(4+2,1) = -1
			Atheta = transpose(c1)
			
			allocate(c2(SIZE(transpose(btheta))+2))
			c2(1:SIZE(transpose(btheta))) = transpose(btheta)
			c2(SIZE(transpose(btheta))+1) = 1
			c2(SIZE(transpose(btheta))+2) = 0
		  	btheta = transpose(c2)
			
			deallocate(c1)
			deallocate(c2)
			ftheta = -1
			!now we need to solve al linear optimization problem, therefore we need the simplex method
			theta = linprog1(ftheta,Atheta,btheta)
			
			x = mu + MATMUL(1-theta,transpose(xf)-mu)
			discard = 1
		END IF
	n = nar
	DO WHILE (n.LT.N)
		DO i=1,p
			IF (.NOT. (i .EQ. p)) THEN !ith column of sigma, remove ith row
			Sigmai_i(1:(i-1)) = sigma(1:(i-1),i)
			Sigmai_i(i:(p-1)) = sigma((i+1):p,i) 
			ELSE
			Sigmai_i(1:(p-1)) = sigma(1:(p-1),p)
			END IF
			
			!inverse of sigma in which ith row and ith column is removed
			SigmaInv_i_i(1:(i-1),1:(i-1)) = SigmaInv(1:(i-1),1:(i-1))
			SigmaInv_i_i(1:(i-1),i:(p-1)) = SigmaInv(1:(i-1),(i+1):p)
			SigmaInv_i_i(i:(p-1),1:(i-1)) = SigmaInv((i+1):p,1:(i-1))
			SigmaInv_i_i(i:(p-1),i:(p-1)) = SigmaInv((i+1):p,(i+1):p)
			SigmaInv_i(1:(i-1)) = SigmaInv(1:(i-1),i)
			SigmaInv_i(i:(p-1)) = SigmaInv((i+1):p,i)
			Sigma_i_iInv = SigmaInv_i_i - MATMUL(SigmaInv_i,transpose(SigmaInv_i))/SigmaInv(i,i)
		
			x_i(1:(i-1))=x(1:(i-1))
			x_i(i:(p-1))=x((i+1):p)

			mu_i(1:(i-1))=mu(1:(i-1))
			mu_i(i:(p-1))=mu((i-1):p)
		
		mui = mu(i) + MATMUL(transpose(Sigmai_i),MATMUL(Sigma_i_iInv,(transpose(x_i) - transpose(mu_i))))
		s2i = sigma(i,i) - MATMUL( transpose(Sigmai_i),MATMUL(Sigma_i_iInv,(transpose(x_i) - transpose(mu_i))) )

		A_i(1:(i-1),:) = A(1:(i-1),:)
		A_i(i:(p-1),:) = A((i+1):p,:)

		Ai = A(i,:)
		c_temp = b - MATMUL(x_i,A_i)
		
		DO q=1,m
			c(q) = c_temp(q)/Ai(1,q)
		ENDDO
		
		lb = max(c(Ai<0))
		IF (SIZE(lb)<1) THEN
			lb = log(0)
		ENDIF
		ub = min(c(Ai>0))
		IF (SIZE(ub)<1) THEN
			ub = -log(0)
		ENDIF
		range_lb_ub(1)=lb-mui
		range_lb_ub(2)=ub-mui
		x(i) = TruncatedGaussian(-sqrt(s2i),range_lb_ub) !-mui want je wilt een verdeling rond nul, zodat je kan sampelen als een normale verdeling
		x(i) = x(i) + mui !+mui om het gemiddelde juist in te stellen
		END DO
	IF (discard .LE. 0) THEN
		n = n + 1
		Xbig(n,:) = x
		ngibbs = ngibbs+1
	END IF
	discard = discard - 1
	END DO 
	END IF
END FUNCTION

!OK
FUNCTION rand_normal(mean,stdev) RESULT(c)
	implicit none
	real*8, intent(in) :: mean, stdev
	real*8 :: c
      	real*8 :: theta,r
	real*8, dimension(2) :: temp
	real*8,parameter :: PI_8 = 4*ATAN(1.0_8)

        CALL RANDOM_NUMBER(temp)
        r=(-2.0*log(temp(1)))**0.5
        theta = 2.0*PI_8*temp(2)
        c= mean+stdev*r*sin(theta)
      
END FUNCTION

!TODO
FUNCTION linprog2(f,A,b,lb,ub) RESULT(theta)
	implicit none
	real*8,dimension(:,:),intent(in) :: A,b
	real*8,dimension(1,2),intent(in) :: lb,ub
	real*8,dimension(:,:) :: f
	real*8,dimension(:),allocatable::theta
	
	allocate(theta(size(f,1)))
	theta=0.0d0
END FUNCTION

!TODO
FUNCTION linprog1(f,A,b) RESULT(theta)
	implicit none
	real*8,dimension(:,:),intent(in) :: A,b
	real*8,dimension(:,:) :: f
	real*8,dimension(:),allocatable::theta
	
	allocate(theta(size(f,1)))
	theta=0.0d0
END FUNCTION

!computes center of largest hypersphere enclosed by the polytope Ax <= b
FUNCTION chebycenter(A,b,r0) RESULT (c)
	implicit none
	real*8,dimension(4,2),intent(in)::A
	real*8,dimension(4,1),intent(in)::b
	real*8,intent(in)::r0
	integer :: p
	REAL*8,DIMENSION(:,:),ALLOCATABLE :: A1
	REAL*8,DIMENSION(:,:),ALLOCATABLE :: f
	REAL*8,DIMENSION(:,:),ALLOCATABLE :: lb
	REAL*8,DIMENSION(:,:),ALLOCATABLE :: ub
	REAL*8,DIMENSION(:,:),ALLOCATABLE :: c
	REAL*8,DIMENSION(:),ALLOCATABLE :: ctemp	
	real*8,DIMENSION(4,1) :: an
	REAL*8 :: r
	n = size(A,1)
	p = size(A,2)
	ALLOCATE(A1(n,p+1))
	A1 = 0.0d0
	ALLOCATE(f(p+1,1))
	f = 0.0d0
	ALLOCATE(lb(p+1,1))
	lb = 1.0d0
	ALLOCATE(ub(p+1,1))
	ub = 1.0d0
	ALLOCATE(c(p,1))
	ALLOCATE(ctemp(p+1))
	an(:,1) = SQRT(SUM(A**2,2)) !should be pointwise
	A1(:,1:p) = A
	A1(:,p+1) = an(:,1)
	f(p+1,1) = -1
	
	lb = lb*huge(0.0)
	ub = -ub*huge(0.0)
	lb(p+1,1) = huge(0.0)
	ub(p+1,1) = r0
	
	ctemp=linprog2(f,A1,b,lb,ub)
	r = ctemp(p+1)
	c(:,1) = ctemp(1:p)
END FUNCTION

FUNCTION TruncatedGaussian(sigma ,range) RESULT(X)
	implicit none
	real*8 :: sigma,range
	real*8 :: X
	call TRUNCATED_NORMAL_AB_SAMPLE(0.0_8,sigma,range(1),range(2),8,X) !hier seed=8 genomen?
!use truncated_normal library: TRUNCATED_NORMAL_AB_SAMPLE(0,sigma,range(1),range(2),8,X)
END FUNCTION

!FUNCTION inv(A) RESULT(Ainv)
!	real(dp), dimension(:,:), intent(in) :: A
! 	real(dp), dimension(size(A,1),size(A,2)) :: Ainv
!
!  	real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
!  	integer, dimension(size(A,1)) :: ipiv   ! pivot indices
!  	integer :: n, info
!
!  	! External procedures defined in LAPACK
!  	external DGETRF
!  	external DGETRI
!
!  	! Store A in Ainv to prevent it from being overwritten by LAPACK
!  	Ainv = A
!  	n = size(A,1)
!
!  	! DGETRF computes an LU factorization of a general M-by-N matrix A
!  	! using partial pivoting with row interchanges.
!  	call DGETRF(n, n, Ainv, n, ipiv, info)
!
!  	if (info /= 0) then
!     		stop 'Matrix is numerically singular!'
!  	end if
!
!  	! DGETRI computes the inverse of a matrix using the LU factorization
!  	! computed by DGETRF.
!  	call DGETRI(n, Ainv, n, ipiv, work, n, info)
!
!  	if (info /= 0) then
!     		stop 'Matrix inversion failed!'
!  	end if
!END FUNCTION

end program
