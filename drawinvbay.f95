PROGRAM drawinvbay

  real (kind = selected_real_kind(8) ) :: kappa, lambda, S0, I0,R0
  real*8,dimension(2) :: mu,lower,upper
  real*8,dimension(2,2) :: sigma, SigmaInv
  real*8,dimension(4,2) :: A
  real*8,dimension(4) :: b
  integer :: N
  real*8, dimension(:), allocatable :: array_kappa
  integer :: nr_values, i 
  real*8 :: normal

  nr_values = 10
  allocate(array_kappa(nr_values))
	
  mu(1) = 0.0144
  mu(2) = 17.9
  
  lower(2)=1
  upper(1)=1
  upper(2)=35
  
  A = reshape((/1, 0, -1, 0, 0, 1, 0, -1/), shape(A))
  sigma = reshape((/0.000036, -0.0187, -0.0187, 16.09/), shape(sigma))
  SigmaInv = reshape((/70093.66, 81.46373,81.46373, 0.1568285/), shape(SigmaInv))
  
  !kappa = rand_gamma(2.0_8,0.0001_8)
  !print*,kappa


  !lambda = rand_gamma(2.0_8,0.0001_8)
  !print*,lambda
  !S0 = 0.9

  !I0 = rand_beta(1.62_8,7084.10_8) !these values by fitting a beta distribution to historical ILI+ for t=0
  !print*,I0
  !R0 = 1-S0-I0
  
  !lower(1) = I0
  !b(1) = upper(1)
  !b(2) = upper(2)
  !b(3) = -lower(1)
  !b(4) = -lower(2)
  
  !N=1 !how many samples do you want
  
  !rows_A = SIZE(A,1)
  !columns_A = SIZE(A,2)
  !rows_b = SIZE(b,1)
  !z = rand_truncated_normal(mu,sigma,SigmaInv,N,A,b)
  
  !rho = calculate_rho(10,I0,S0,PI)

  normal = rand_normal(5.0_8,1.0_8)
  print*,normal
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
		!print*,x
            v = 1.0 + c*x
            DO while (v <= 0.0)
                x = rand_normal(0.0_8, 1.0_8)
                v = 1.0 + c*x
            END DO

            v = v*v*v
	    call init_random_seed()
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
	call init_random_seed()
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
FUNCTION rand_truncated_normal(mu,sigma,SigmaInv,Nbig,Atemp,btemp) RESULT(r) !MULTIVARIATE
	implicit none
	integer :: q
	real*8,dimension(2,2),intent(in) :: sigma,SigmaInv
	real*8,dimension(4,2),intent(in) :: Atemp
	real*8,dimension(4,1),intent(in) :: btemp
	real*8,dimension(2,4) :: A
	real*8,dimension(1,4) :: b
	real*8,dimension(1,2),intent(in) :: mu
	real*8 :: r
	integer,intent(in) :: Nbig
	real*8 :: trials,passes,s,maxsample,rhoThr,defaultRhoThr,rho,ngibbs
	integer :: nar,p,m
	real*8,dimension(2) :: range_lb_ub
	real*8,dimension(:),allocatable::c
	real*8,dimension(:,:),allocatable::c1,c2,x_i,mu_i,c_temp,Ai
	real*8,dimension(:,:),allocatable::Sigmai_i
	REAL*8,DIMENSION(:,:),ALLOCATABLE::xf,x,Xbig,sigmaInv_i_i,A_i,Sigma_i_iInv,SigmaInv_i
	real*8,dimension(4,1) :: Atheta,Btheta
	real*8,dimension(1,1) :: ftheta
	integer :: discard
	integer,parameter :: CMAX=10, VMAX=10
	integer :: NC,NV,XERR
	integer :: i
	real*8 :: lb,ub,maxc,minc
	real*8,dimension(2,1) :: theta
	real*8,dimension(1,1) :: mui,s2i
	real*8 TS(CMAX,VMAX)
	integer :: hasmax,hasmin
	defaultRhoThr = 0.00029
	rhoThr = defaultRhoThr

	A = transpose(Atemp) !2,4
	b = transpose(btemp) !1,4
	

	p = SIZE(mu) !2
	m = SIZE(A,2) !4
	
	ALLOCATE(xf(p,1))
	ALLOCATE(Xbig(N,p))
	ALLOCATE(x(1,p)) !1,2
	ALLOCATE(Sigmai_i(1,p-1))
	ALLOCATE(SigmaInv_i_i(p-1,p-1))
	ALLOCATE(Sigma_i_iInv(p-1,p-1))
	ALLOCATE(SigmaInv_i(p-1,1))
	ALLOCATE(x_i(1,p-1))
	ALLOCATE(mu_i(1,p-1))
	ALLOCATE(A_i(p-1,m))
	ALLOCATE(c(m))
	ALLOCATE(c_temp(1,m))
	ALLOCATE(Ai(1,m))
	Xbig = 0.0d0
	nar = 0
	ngibbs = 0.0
	rho = 1.0

	IF (nar.LT.Nbig) THEN !nog waarden gevraagd
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
			c1(1:4,1) = Atheta(:,1)
			c1(4+1,1) = 1
			c1(4+2,1) = -1
			Atheta = transpose(c1)
			
			allocate(c2(4+2,1))
			c2(1:4,1) = Btheta(:,1)
			c2(4+1,1) = 1
			c2(4+2,1) = 0
		  	btheta = transpose(c2)
			
			deallocate(c1)
			deallocate(c2)
			ftheta(1,1) = -1
			!now we need to solve al linear optimization problem, therefore we need the simplex method
			theta = linprog1(ftheta,Atheta,btheta)
			
			x = mu + MATMUL(1-theta,xf-transpose(mu))
			discard = 1
		END IF
	n = nar
	DO WHILE (n.LT.Nbig)
		DO i=1,p
			IF (.NOT. (i .EQ. p)) THEN !ith column of sigma, remove ith row
			Sigmai_i(1,1:(i-1)) = sigma(1:(i-1),i)
			Sigmai_i(1,i:(p-1)) = sigma((i+1):p,i) 
			ELSE
			Sigmai_i(1,1:(p-1)) = sigma(1:(p-1),p)
			END IF
			
			!inverse of sigma in which ith row and ith column is removed
			SigmaInv_i_i(1:(i-1),1:(i-1)) = SigmaInv(1:(i-1),1:(i-1))
			SigmaInv_i_i(1:(i-1),i:(p-1)) = SigmaInv(1:(i-1),(i+1):p)
			SigmaInv_i_i(i:(p-1),1:(i-1)) = SigmaInv((i+1):p,1:(i-1))
			SigmaInv_i_i(i:(p-1),i:(p-1)) = SigmaInv((i+1):p,(i+1):p)
			SigmaInv_i(1:(i-1),1) = SigmaInv(1:(i-1),i)
			SigmaInv_i(i:(p-1),1) = SigmaInv((i+1):p,i)
			Sigma_i_iInv = SigmaInv_i_i - MATMUL(SigmaInv_i,transpose(SigmaInv_i))/SigmaInv(i,i)
		
			x_i(1,1:(i-1))=x(1,1:(i-1))
			x_i(1,i:(p-1))=x(1,(i+1):p)

			mu_i(1,1:(i-1))=mu(1,1:(i-1))
			mu_i(1,i:(p-1))=mu(1,(i-1):p)
		
		mui(:,:) = mu(1,i) + MATMUL(Sigmai_i,MATMUL(Sigma_i_iInv,(transpose(x_i) - transpose(mu_i))))

		s2i(:,:) = sigma(i,i) - MATMUL( transpose(Sigmai_i),MATMUL(Sigma_i_iInv,(transpose(x_i) - transpose(mu_i))) )

		A_i(1:(i-1),:) = A(1:(i-1),:)
		A_i(i:(p-1),:) = A((i+1):p,:)

		Ai(1,:) = A(i,:)
		c_temp(1:1,:) = b - MATMUL(x_i,A_i)
		
		DO q=1,m
			c(q) = c_temp(1,q)/Ai(1,q)
		ENDDO
		
		maxc=-huge(0.0)
		hasmax = 0
		DO q=1,SIZE(c)
			IF ((c(q)>maxc) .AND. (Ai(1,q)<0)) THEN
				maxc = c(q)
				hasmax = 1
			ENDIF
		ENDDO

		lb = maxc

		IF (hasmax .EQ. 0) THEN
			lb = huge(0.0)
		ENDIF

		minc=huge(0.0)
		hasmin = 0
		DO q=1,SIZE(c)
			IF ((c(q)<minc) .AND. (Ai(1,q)>0)) THEN
				minc = c(q)
				hasmin = 1
			ENDIF
		ENDDO

		ub = minc		
		
		IF (hasmin .EQ. 0) THEN
			ub = -huge(0.0)
		ENDIF
		range_lb_ub(1)=lb-mui(1,1)
		range_lb_ub(2)=ub-mui(1,1)
		x(1,i) = TruncatedGaussian(-sqrt(s2i(1,1)),range_lb_ub) !-mui want je wilt een verdeling rond nul, zodat je kan sampelen als een normale verdeling
		x(1,i) = x(1,i) + mui(1,1) !+mui om het gemiddelde juist in te stellen
		END DO
	IF (discard .LE. 0) THEN
		n = n + 1
		Xbig(n,:) = x(1,:)
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
	

	call init_random_seed()
	
        CALL RANDOM_NUMBER(temp)
	!print*,temp
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
	real*8,dimension(:,:),allocatable::theta
	
	allocate(theta(size(f,1),1))
	theta=0.0d0
END FUNCTION

!TODO
FUNCTION linprog1(f,A,b) RESULT(theta)
	implicit none
	real*8,dimension(:,:),intent(in) :: A,b
	real*8,dimension(:,:) :: f
	real*8,dimension(:,:),allocatable::theta

	allocate(theta(size(f,1),1))
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
	REAL*8,DIMENSION(:,:),ALLOCATABLE :: ctemp	
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
	ALLOCATE(ctemp(p+1,1))
	an(:,1) = SQRT(SUM(A**2,2)) !should be pointwise
	A1(:,1:p) = A
	A1(:,p+1) = an(:,1)
	f(p+1,1) = -1
	
	lb = lb*huge(0.0)
	ub = -ub*huge(0.0)
	lb(p+1,1) = huge(0.0)
	ub(p+1,1) = r0
	
	ctemp=linprog2(f,A1,b,lb,ub)
	r = ctemp(p+1,1)
	c(:,1) = ctemp(1:p,1)
END FUNCTION

RECURSIVE FUNCTION TruncatedGaussian(sigma,range) RESULT(X)
	implicit none
	real*8 :: sigma
	real*8,dimension(2) :: range
	real*8 :: X
	X = rand_normal(0.0_8,sigma)
	IF (x<range(1) .OR. x>range(2)) THEN
		X = TruncatedGaussian(sigma,range)
	ENDIF
	
END FUNCTION

function calculate_rho(m,I0,S0,PI) result(sol)
      real*8 :: a,b
      real*8, intent(in) :: I0,S0,PI
      integer, intent(in) :: m
      real :: fa, fb, temp, sol
      integer :: n
      interface
      function f(x)
      real, intent(in) :: x
      end function f
      end interface
	a = -10
	b = 10

      fa = evaluate_g_inverse(a,I0,S0,PI)
      fb = evaluate_g_inverse(b,I0,S0,PI)
      if (abs(fa) >  abs(fb)) then
         temp = a
         a = b
         b = temp
         temp = fa
         fa = fb
         fb = temp
      end if
      print *,"    n        x(n)         f(x(n))"
      print *," 1 ", b, fb
      print *," 0 ", a, fa	
      do n = 2,m
         if (abs(fa) >  abs(fb)) then
            temp = a
            a = b
            b = temp
            temp = fa
            fa = fb
            fb = temp
         end if
         temp = (b - a)/(fb - fa)
  	 b = a
	 fb = fa
         a = a - fa*temp
         fa = evaluate_g_inverse(a,I0,S0,PI)
         !print *,n,a,fa
      end do  
	sol = a 
      end function calculate_rho

function evaluate_g_inverse(x,I0,S0,PI) result(fx)
	real*8,intent(in) :: x,I0,S0,PI
	fx = I0 + S0 - PI - x*(log(S0) + 1 - log(x))
end function


 SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
 END SUBROUTINE

        
!FUNCTION TruncatedGaussian(sigma ,range) RESULT(X)
!	implicit none
!	real*8 :: sigma
!	real*8,dimension(2) :: range
!	real*8 :: X
!	call TRUNCATED_NORMAL_AB_SAMPLE(0.0_8,sigma,range(1),range(2),8,X) !hier seed=8 !genomen?
!!use truncated_normal library: TRUNCATED_NORMAL_AB_SAMPLE(0,sigma,range(1),range(2),8,X)
!END FUNCTION

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
