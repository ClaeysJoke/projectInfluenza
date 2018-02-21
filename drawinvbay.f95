PROGRAM drawinvbay
use normaldistribution
  real (selected_real_kind(8) ) :: kappa, lambda, S0, I0,R0
  
  kappa = call rand_gamma(2, 0.0001)
  lambda = call rand_gamma(2,0.0001)
  S0 = 0.9
  I0 = call rand_beta(1.62,7084.10) !these values by fitting a beta distribution to historical ILI+ for t=0
  R0 = 1-S0-I0
END

!OK
RECURSIVE FUNCTION rand_gamma(shape, SCALE) RESULT(ans)
      DOUBLE PRECISION SHAPE,scale,u,w,d,c,x,xsq,g
      IF (shape <= 0.0d0) THEN

        WRITE(*,*) "Shape PARAMETER must be positive"
      END IF
      IF (scale <= 0.0d0) THEN

        WRITE(*,*) "Scale PARAMETER must be positive"
      END IF
      
      IF (shape >= 1.0d0) THEN
        d = SHAPE - 1.0d0/3.0d0
        c = 1.0d0/(9.0d0*d)**0.5
        DO while (.true.)
            x = rand_normal(0.0d0, 1.0d0)
            v = 1.0 + c*x
            DO while (v <= 0.0d0)
                x = rand_normal(0.0d0, 1.0d0)
                v = 1.0d0 + c*x
            END DO

            v = v*v*v
            CALL RANDOM_NUMBER(u)
            xsq = x*x
            IF ((u < 1.0d0 -.0331d0*xsq*xsq) .OR.  &
              (log(u) < 0.5d0*xsq + d*(1.0d0 - v + log(v))) )then
                ans=scale*d*v
                RETURN
            END IF

        END DO
      ELSE
        g = rand_gamma(shape+1.0d0, 1.0d0)
        CALL RANDOM_NUMBER(w)
        ans=scale*g*(w)**(1.0d0/shape)
        RETURN
      END IF
END FUNCTION

!OK
FUNCTION rand_beta(a, b) RESULT(ans)
      DOUBLE PRECISION a,b,ans,u,v
      IF ((a <= 0.0d0) .OR. (b <= 0.0d0)) THEN

        WRITE(*,*) "Beta PARAMETERs must be positive"
      END IF
      
       u = rand_gamma(a, 1.0d0)
       v = rand_gamma(b, 1.0d0)
       ans = u / (u + v)
END FUNCTION

FUNCTION rand_truncated_normal(mu,sigma,SigmaInv,N,A,b) RESULT r !MULTIVARIATE
	DOUBLE PRECISION trials,passes,s,R,p,m,mu,maxsample,sigma,N,A,b,rhoThr,defaultRhoThr,rho,nar,ngibbs
	REAL,DIMENSION(:,:),ALLOCATABLE::X
	real,dimension(:),allocatable::c
	real,dimension(:),allocatable::c2
	integer,parameter :: CMAX=10, VMAX=10
	integer :: NC,NV,XERR
	real*8 TS(CMAX,VMAX)
	defaultRhoThr = 0.00029
	rhoThr = defaultRhoThr

	A = transpose(A)
	B = transpose(B)
	mu = transpose(mu)

	p = SIZE(mu)
	m = SHAPE(A)(2)
	
	ALLOCATE(X(N,p))
	X = 0.0d
	nar = 0
	ngibbs = 0
	rho = 1

	IF nar.LT.N THEN !nog waarden gevraagd
		if nar.GT.0 THEN
			x=X(nar,:)
			discard = 0
		ELSE IF ALL(MATMUL(transpose(A),transpose(mu)).LE.transpose(b)) THEN
		x = transpose(mu)
		discard = 1
		ELSE
		xf = call chebycenter(transpose(A),transpose(b),1)
			IF .NOT. ALL(MATMUL(transpose(A),xf) .LE. transpose(b)) THEN
			print *,'error: failed to find a feasible point'
			ENDIF
		Atheta = -MATMUL(transpose(A),xf-transpose(mu))
		btheta = transpose(b) - MATMUL(transpose(A),xf)
		allocate(c(SIZE(transpose(Atheta))+2))
		c(1:SIZE(transpose(Atheta))) = transpose(Atheta)
		c(SIZE(transpose(Atheta))+1) = 1
		c(SIZE(transpose(Atheta))+2) = -1
		Atheta = transpose(c)
		allocate(c2(SIZE(transpose(btheta))+2))
		c2(1:SIZE(transpose(btheta))) = transpose(btheta)
		c2(SIZE(transpose(btheta))+1) = 1
		c2(SIZE(transpose(btheta))+2) = 0
		btheta = transpose(c2)
		deallocate(c1)
		deallocate(c2)
		ftheta = -1
		!now we need to solve al linear optimization problem, therefore we need the simplex method
		linprog
		theta = call Results(NV,NC,TS,XERR)
		x = transpose(mu) + MATMUL(1-theta,xf-transpose(mu))
		discard = 1
		END IF
	n = nar
	DO WHILE n.LT.N
		DO i=1,p
			Sigmai_i = sigma([1:(i-1) (i+1):p],i)
			Sigma_i_iInv = SigmaInv([1:(i-1) (i+1):p],[1:(i-1) (i+1):p]) - MATMUL(SigmaInv([1:(i-1) (i+1):p],i),transpose(SigmaInv([1:(i-1) (i+1):p],i))) / SigmaInv(i,i)
		x_i = x([1:(i-1) (i+1):p])
		mu_i = mu([1:(i-1) (i+1):p])
		mui = mu(i) + MATMUL(transpose(Sigmai_i),MATMUL(Sigma_i_iInv,(transpose(x_i) - transpose(mu_i))))
		s2i = sigma(i,i) - MATMUL( transpose(Sigmai_i),MATMUL(Sigma_i_iInv,(transpose(x_i) - transpose(mu_i))) )

		A_i = A([1:(i-1) (i+1):p],:)
		Ai = A(i,:)
		c = (b-MATMUL(x_i,A_i))/Ai
		lb = max(c(Ai<0))
		IF (SIZE(lb)<1) THEN
			lb = log(0)
		ENDIF
		ub = min(c(Ai>0))
		IF (SIZE(ub)<1) THEN
			ub = -log(0)
		ENDIF
		x(i) = call TruncatedGaussian(-sqrt(s2i),[lb ub]-mui) !-mui want je wilt een verdeling rond nul, zodat je kan sampelen als een normale verdeling
		x(i) = x(i) + mui !+mui om het gemiddelde juist in te stellen
	END DO
	IF (discard .LE. 0) THEN
		n = n + 1
		X(n,:) = x
		ngibbs = ngibbs+1
	END IF
	discard = discard - 1
END DO 
END FUNCTION

!OK
FUNCTION rand_normal(mean,stdev) RESULT(c)
       DOUBLE PRECISION :: mean,stdev,c,temp(2)
      IF(stdev <= 0.0d0) THEN

        WRITE(*,*) "Standard Deviation must be +ve"
      ELSE
        CALL RANDOM_NUMBER(temp)
        r=(-2.0d0*log(temp(1)))**0.5
        theta = 2.0d0*PI*temp(2)
        c= mean+stdev*r*sin(theta)
      END IF
END FUNCTION

FUNCTION linprog(f,A,b,lb,ub) RESULT(theta)
	
END FUNCTION

!computes center of largest hypersphere enclosed by the polytope Ax <= b
FUNCTION chebycenter(A,b,r0) RESULT (c,r)
	REAL,DIMENSION(:,:),ALLOCATABLE :: A1
	REAL,DIMENSION(:,:),ALLOCATABLE :: f
	REAL,DIMENSION(:,:),ALLOCATABLE :: lb
	REAL,DIMENSION(:,:),ALLOCATABLE :: ub
	n = shape(A)(1)
	p = shape(A)(2)
	ALLOCATE(A1(n,p+1))
	ALLOCATE(f(p+1,1))
	ALLOCATE(lb(p+1,1))
	lb = 1.0d
	ALLOCATE(ub(p+1,1))
	ub = 1.0d
	an = SQRT(SUM(A**2,2)) !should be pointwise
	A1 = 0.0d
	A1(:,1:p) = A
	A1(:,p+1) = an
	f(p+1) = -1
	
	lb = lb*log(0)
	ub = -ub*log(0)
	lb(p+1) = log(0)
	ub(p+1) = r0
	
	c=linprog(f,A1,b,lb,ub)
	r = c(p+1)
	c = c(1:p)
END FUNCTION

FUNCTION TruncatedGaussian(sigma ,range) RESULT(X)
	call TRUNCATED_NORMAL_AB_SAMPLE(0,sigma,range(1),range(2),8,X) !hier seed=8 genomen?
!use truncated_normal library: TRUNCATED_NORMAL_AB_SAMPLE(0,sigma,range(1),range(2),8,X)
END FUNCTION

FUNCTION inv(A) RESULT(Ainv)
	real(dp), dimension(:,:), intent(in) :: A
  	real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  	real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  	integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  	integer :: n, info

  	! External procedures defined in LAPACK
  	external DGETRF
  	external DGETRI

  	! Store A in Ainv to prevent it from being overwritten by LAPACK
  	Ainv = A
  	n = size(A,1)

  	! DGETRF computes an LU factorization of a general M-by-N matrix A
  	! using partial pivoting with row interchanges.
  	call DGETRF(n, n, Ainv, n, ipiv, info)

  	if (info /= 0) then
     		stop 'Matrix is numerically singular!'
  	end if

  	! DGETRI computes the inverse of a matrix using the LU factorization
  	! computed by DGETRF.
  	call DGETRI(n, Ainv, n, ipiv, work, n, info)

  	if (info /= 0) then
     		stop 'Matrix inversion failed!'
  	end if
END FUNCTION


