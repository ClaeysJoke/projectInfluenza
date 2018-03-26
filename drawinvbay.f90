MODULE drawinvbay
 implicit none
  real*8,dimension(2) :: lower,upper
  real*8,dimension(2,2) :: SigmaInv
  integer :: N
  real*8, dimension(:), allocatable :: array_kappa
  integer :: nr_values
  real*8 :: normal, PI, PT
  real*8, dimension(17) :: X_beta, tau
  !mean = 0.0

!	open(unit = 5, file='I0.txt')
!	open(unit = 7, file='kappa.txt')
!	open(unit = 8, file='lambda.txt')
!	open(unit = 9, file='R0.txt')


!  nr_values = 10
!  allocate(array_kappa(nr_values))
	
!  mu(1) = 0.0144
!  mu(2) = 17.9
  
!  lower(2)=1
!  upper(1)=1
!  upper(2)=35
  
!  A = reshape((/1, 0, -1, 0, 0, 1, 0, -1/), shape(A))
!  sigma = reshape((/0.000036, -0.0187, -0.0187, 16.09/), shape(sigma))
 ! SigmaInv = reshape((/70093.66, 81.46373,81.46373, 0.1568285/), shape(SigmaInv))
  
!  do i=1,5000
!  kappa = rand_gamma(real(2.0,8),real(10000.0,8))
!  write(7,*) kappa
  !print*,kappa
!  enddo
!  close(7)

!  do i=1,5000
!  lambda = rand_gamma(real(2.0,8),real(10000.0,8))
!  write(8,*) lambda
  !print*,lambda
!  mean = mean + lambda
!  enddo
!  close(8)

!  mean = mean/5000.0
!  print*,mean

!  S0 = 0.9

!  do i=1,5000
!  I0 = rand_beta(1.62_8,7084.10_8) !these values by fitting a beta distribution to historical ILI+ for t=0
!  !print*,I0
!  write(5,*) I0
!  enddo
!  close(5)

!  open(unit=5,file='I0.txt')
!  do i=1,5000
!  read(5,*) I0
!  R0 = 1-0.9-I0
!  write(9,*) R0
!  enddo
!  close(5)
!  close(9)
  
  !lower(1) = I0
  !b(1) = upper(1)
  !b(2) = upper(2)
  !b(3) = -lower(1)
  !b(4) = -lower(2)

  
!  N=1 !how many samples do you want
  
  !vanaf hier niet meer loopen, gwn in 1 keer, geen init_seed meer nodig

  !rows_A = SIZE(A,1)
  !columns_A = SIZE(A,2)
  !rows_b = SIZE(b,1)

!  call accept_reject
  !print *,z(1,1)
  !print *,z(1,2)

  !z = rand_truncated_normal(mu,sigma,SigmaInv,N,A,b)

  !PI = z(1)
  !PT = z(2)  

!  call calculate_rho(10)




!  call calc_betas !writes the betavalues to beta, requires truncated_normal.txt, rho.txt, I0.txt


contains

subroutine calc_betas
	implicit none
	real*8, dimension(17) :: X_beta, tau, beta
	real*8 :: rho, I0
	integer :: i
	tau(1) = -49.7540
	tau(2) = -0.9577
	tau(3) = -0.0065
	tau(4) = -9.4896
	tau(5) = -0.3761
	tau(6) = -590.0001
	tau(7) = -2537.6102
	tau(8) = -4756.1828
	tau(9) = -3265.2458
	tau(10) = -102.2665
	tau(11) = -4.0162
	tau(12) = -430.9596
	tau(13) = -16.7104
	tau(14) = -798.3443
	tau(15) = -30.6638
	tau(16) = -543.8857
	tau(17) = -20.7459
	open(unit = 2, file='truncated_normal.txt')
	open(unit = 3, file='betas.txt')
	open(unit = 4, file='rho.txt')
	open(unit = 5, file='I0.txt')
		do i=1,5000
			READ(2,*) PI ,PT
			read(4,*) rho
			read(5,*) I0
			X_beta(1) = 1
  			X_beta(2) = log(PT)
  			X_beta(3) = (log(PT))**2
  			X_beta(4) = log(I0)
  			X_beta(5) = (log(I0))**2
  			X_beta(6) = log(rho)
  			X_beta(7) = (log(rho))**2
  			X_beta(8) = (log(rho))**3
  			X_beta(9) = (log(rho))**4
  			X_beta(10) = log(I0)*log(rho)
  			X_beta(11) = ((log(I0))**2)*log(rho)
  			X_beta(12) = log(I0)*((log(rho))**2)
  			X_beta(13) = ((log(I0))**2)*((log(rho))**2)
  			X_beta(14) = log(I0)*((log(rho))**3)
  			X_beta(15) = ((log(I0))**2)*((log(rho))**3)
  			X_beta(16) = log(I0)*((log(rho))**4)
  			X_beta(17) = ((log(I0))**2)*((log(rho))**4)
			beta = EXP(dot_product(X_beta,tau) + 0.5*(0.0421**2))
			write(3,*) beta
		enddo
	close(2)
	close(3)
	close(4)
	close(5)
end subroutine

RECURSIVE FUNCTION rand_gamma(shape, scale) RESULT(ans)
	implicit none
	real (kind=selected_real_kind(8) ) :: ans
	real (kind=selected_real_kind(8) ), intent(in) :: shape, scale !shape and scale
     	real (kind=selected_real_kind(8) ) :: u,w,d,c,x,xsq,g,v
	call init_random_seed()
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
				if (ans.ne. ans) then
				print*, "NaN exception"
				ans=rand_gamma(shape,scale)
				end if
                RETURN
            END IF

        END DO
      ELSE
        g = rand_gamma(shape+1.0, 1.0_8)
	call init_random_seed()
        CALL RANDOM_NUMBER(w)
        ans=scale*g*(w)**(1.0/shape)
		if (ans.ne. ans) then
		print*, "NaN exception"
		ans=rand_gamma(shape,scale)
		end if
        RETURN
      END IF
END FUNCTION

!OK
FUNCTION rand_beta(a, b) RESULT(ans)
	implicit none
	real*8, intent(in) :: a,b
	real*8 :: ans
      	real*8 :: u,v
		call init_random_seed()
      IF ((a <= 0.0) .OR. (b <= 0.0)) THEN

        WRITE(*,*) "Beta PARAMETERs must be positive"
      END IF
      
       u = rand_gamma(a, 1.0_8)
       v = rand_gamma(b, 1.0_8)
       ans = u / (u + v)
END FUNCTION

subroutine accept_reject
	implicit none
	real*8, dimension(1,2) :: z
	real*8 :: lb1,PI,PT
	integer :: succeed,i,j
	succeed = 0
	i=1
	j=1
	OPEN(unit=1,FILE='normal_02_15000.txt') !you become this by running normal_distribution.f95
	open(unit = 2, file='truncated_normal.txt')
	open(unit=5, file='I0.txt')
    		do while (i .LT. 5001) !opgelet j mag niet terug naar nul worden gezet, want je leest de j-de lijn van de 15000 lijnen
			
			read(5, *) lb1
			succeed = 0
			do while ((j .LT. 15001) .AND. (succeed .EQ. 0))
    				READ(1,*) PI ,PT
				!print *,PI
    				IF ((PI.LT.1) .AND. (PI.GT.lb1) .AND. (PT.LT.35) .AND. (PT.GT.1)) THEN
					write(2,*) PI, PT
					succeed = 1
				ENDIF
			j = j + 1
			enddo
		i = i+1
    		enddo
   	CLOSE(1)
	close(2)
	close(5)
END subroutine

!OK
!FUNCTION rand_truncated_normal(mu,sigma,SigmaInv,Nbig,Atemp,btemp) RESULT(r) !MULTIVARIATE
!	implicit none
!	integer :: q
!	real*8,dimension(2,2),intent(in) :: sigma,SigmaInv
!	real*8,dimension(4,2),intent(in) :: Atemp
!	real*8,dimension(4,1),intent(in) :: btemp
!	real*8,dimension(2,4) :: A
!	real*8,dimension(1,4) :: b
!	real*8,dimension(1,2),intent(in) :: mu
!	real*8 :: r
!	integer,intent(in) :: Nbig
!	real*8 :: trials,passes,s,maxsample,rhoThr,defaultRhoThr,rho,ngibbs
!	integer :: nar,p,m
!	real*8,dimension(2) :: range_lb_ub
!	real*8,dimension(:),allocatable::c
!	real*8,dimension(:,:),allocatable::c1,c2,x_i,mu_i,c_temp,Ai
!	real*8,dimension(:,:),allocatable::Sigmai_i
!	REAL*8,DIMENSION(:,:),ALLOCATABLE::xf,x,Xbig,sigmaInv_i_i,A_i,Sigma_i_iInv,SigmaInv_i
!	real*8,dimension(4,1) :: Atheta,Btheta
!	real*8,dimension(1,1) :: ftheta
!	integer :: discard
!	integer,parameter :: CMAX=10, VMAX=10
!	integer :: NC,NV,XERR
!	integer :: i
!	real*8 :: lb,ub,maxc,minc
!	real*8,dimension(2,1) :: theta
!	real*8,dimension(1,1) :: mui,s2i
!	real*8 TS(CMAX,VMAX)
!	integer :: hasmax,hasmin
!	defaultRhoThr = 0.00029
!	rhoThr = defaultRhoThr
!
!	A = transpose(Atemp) !2,4
!	b = transpose(btemp) !1,4
!	
!
!	p = SIZE(mu) !2
!	m = SIZE(A,2) !4
!	
!	ALLOCATE(xf(p,1))
!	ALLOCATE(Xbig(N,p))
!	ALLOCATE(x(1,p)) !1,2
!	ALLOCATE(Sigmai_i(1,p-1))
!	ALLOCATE(SigmaInv_i_i(p-1,p-1))
!	ALLOCATE(Sigma_i_iInv(p-1,p-1))
!	ALLOCATE(SigmaInv_i(p-1,1))
!	ALLOCATE(x_i(1,p-1))
!	ALLOCATE(mu_i(1,p-1))
!	ALLOCATE(A_i(p-1,m))
!	ALLOCATE(c(m))
!	ALLOCATE(c_temp(1,m))
!	ALLOCATE(Ai(1,m))
!	Xbig = 0.0d0
!	nar = 0
!	ngibbs = 0.0
!	rho = 1.0
!
!	IF (nar.LT.Nbig) THEN !nog waarden gevraagd
!		if (nar.GT.0) THEN
!			x(1,:)=Xbig(nar,:) !1,p
!			discard = 0
!		ELSE IF ( ALL(MATMUL(transpose(A),transpose(mu)).LE.transpose(b)) ) THEN
!			x(1,:) = mu(1,:) !1,2
!			discard = 1
!		ELSE
!			xf = chebycenter(transpose(A),transpose(b),1.0_8)
!				IF (.NOT. (ALL(MATMUL(transpose(A),xf) .LE. transpose(b))) ) THEN
!					print *,'error: failed to find a feasible point'
!				END IF
!			Atheta = -MATMUL(transpose(A),xf-transpose(mu))
!			btheta = transpose(b) - MATMUL(transpose(A),xf)
!			
!			allocate(c1(4+2,1))
!			c1(1:4,1) = Atheta(:,1)
!			c1(4+1,1) = 1
!			c1(4+2,1) = -1
!			Atheta = transpose(c1)
!			
!			allocate(c2(4+2,1))
!			c2(1:4,1) = Btheta(:,1)
!			c2(4+1,1) = 1
!			c2(4+2,1) = 0
!		  	btheta = transpose(c2)
!			
!			deallocate(c1)
!			deallocate(c2)
!			ftheta(1,1) = -1
!			!now we need to solve al linear optimization problem, therefore we need the simplex method
!			theta = linprog1(ftheta,Atheta,btheta)
!			
!			!x = mu + MATMUL(1-theta,xf-transpose(mu))
!			discard = 1
!		END IF
!	n = nar
!	DO WHILE (n.LT.Nbig)
!		DO i=1,p
!			IF (.NOT. (i .EQ. p)) THEN !ith column of sigma, remove ith row
!			Sigmai_i(1,1:(i-1)) = sigma(1:(i-1),i)
!			Sigmai_i(1,i:(p-1)) = sigma((i+1):p,i) 
!			ELSE
!			Sigmai_i(1,1:(p-1)) = sigma(1:(p-1),p)
!			END IF
!			
!			!inverse of sigma in which ith row and ith column is removed
!			SigmaInv_i_i(1:(i-1),1:(i-1)) = SigmaInv(1:(i-1),1:(i-1))
!			SigmaInv_i_i(1:(i-1),i:(p-1)) = SigmaInv(1:(i-1),(i+1):p)
!			SigmaInv_i_i(i:(p-1),1:(i-1)) = SigmaInv((i+1):p,1:(i-1))
!			SigmaInv_i_i(i:(p-1),i:(p-1)) = SigmaInv((i+1):p,(i+1):p)
!			SigmaInv_i(1:(i-1),1) = SigmaInv(1:(i-1),i)
!			SigmaInv_i(i:(p-1),1) = SigmaInv((i+1):p,i)
!			Sigma_i_iInv = SigmaInv_i_i - MATMUL(SigmaInv_i,transpose(SigmaInv_i))/SigmaInv(i,i)
!		
!			x_i(1,1:(i-1))=x(1,1:(i-1))
!			x_i(1,i:(p-1))=x(1,(i+1):p)
!
!			mu_i(1,1:(i-1))=mu(1,1:(i-1))
!			mu_i(1,i:(p-1))=mu(1,(i-1):p)
!		
!		mui(:,:) = mu(1,i) + MATMUL(Sigmai_i,MATMUL(Sigma_i_iInv,(transpose(x_i) - transpose(mu_i))))
!
!		s2i(:,:) = sigma(i,i) - MATMUL( transpose(Sigmai_i),MATMUL(Sigma_i_iInv,(transpose(x_i) - transpose(mu_i))) )
!
!		A_i(1:(i-1),:) = A(1:(i-1),:)
!		A_i(i:(p-1),:) = A((i+1):p,:)
!
!		Ai(1,:) = A(i,:)
!		c_temp(1:1,:) = b - MATMUL(x_i,A_i)
!		
!		DO q=1,m
!			c(q) = c_temp(1,q)/Ai(1,q)
!		ENDDO
!		
!		maxc=-huge(0.0)
!		hasmax = 0
!		DO q=1,SIZE(c)
!			IF ((c(q)>maxc) .AND. (Ai(1,q)<0)) THEN
!				maxc = c(q)
!				hasmax = 1
!			ENDIF
!		ENDDO
!
!		lb = maxc
!
!		IF (hasmax .EQ. 0) THEN
!			lb = huge(0.0)
!		ENDIF
!
!		minc=huge(0.0)
!		hasmin = 0
!		DO q=1,SIZE(c)
!			IF ((c(q)<minc) .AND. (Ai(1,q)>0)) THEN
!				minc = c(q)
!				hasmin = 1
!			ENDIF
!		ENDDO
!
!		ub = minc		
!		
!		IF (hasmin .EQ. 0) THEN
!			ub = -huge(0.0)
!		ENDIF
!		range_lb_ub(1)=lb-mui(1,1)
!		range_lb_ub(2)=ub-mui(1,1)
!		x(1,i) = TruncatedGaussian(-sqrt(s2i(1,1)),range_lb_ub) !-mui want je wilt een verdeling rond nul, zodat je kan sampelen als een normale verdeling
!		x(1,i) = x(1,i) + mui(1,1) !+mui om het gemiddelde juist in te stellen
!		END DO
!	IF (discard .LE. 0) THEN
!		n = n + 1
!		Xbig(n,:) = x(1,:)
!		ngibbs = ngibbs+1
!	END IF
!	discard = discard - 1
!	END DO 
!	END IF
!END FUNCTION
!!OK
function rand_normal ( a, b)
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rand_normal
  real ( kind = 8 ), parameter :: r4_pi = 3.141592653589793
  real ( kind = 8 ) r4_uniform_01
  real ( kind = 8 ) x
  call init_random_seed()
  r1 = rand()
  r2 = rand()
  x = sqrt ( - 2.0E+00 * log ( r1 ) ) * sin ( 2.0E+00 * r4_pi * r2 )

  rand_normal = a + b * x

  return
end function

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
!FUNCTION chebycenter(A,b,r0) RESULT (c)
!	implicit none
!	real*8,dimension(4,2),intent(in)::A
!	real*8,dimension(4,1),intent(in)::b
!	real*8,intent(in)::r0
!	integer :: p,n
!	REAL*8,DIMENSION(:,:),ALLOCATABLE :: A1
!	REAL*8,DIMENSION(:,:),ALLOCATABLE :: f
!	REAL*8,DIMENSION(:,:),ALLOCATABLE :: lb
!	REAL*8,DIMENSION(:,:),ALLOCATABLE :: ub
!	REAL*8,DIMENSION(:,:),ALLOCATABLE :: c
!	REAL*8,DIMENSION(:,:),ALLOCATABLE :: ctemp	
!	real*8,DIMENSION(4,1) :: an
!	REAL*8 :: r
!	n = size(A,1)
!	p = size(A,2)
!	ALLOCATE(A1(n,p+1))
!	A1 = 0.0d0
!	ALLOCATE(f(p+1,1))
!	f = 0.0d0
!	ALLOCATE(lb(p+1,1))
!	lb = 1.0d0
!	ALLOCATE(ub(p+1,1))
!	ub = 1.0d0
!	ALLOCATE(c(p,1))
!	ALLOCATE(ctemp(p+1,1))
!	an(:,1) = SQRT(SUM(A**2,2)) !should be pointwise
!	A1(:,1:p) = A
!	A1(:,p+1) = an(:,1)
!	f(p+1,1) = -1
!	
!	lb = lb*huge(0.0)
!	ub = -ub*huge(0.0)
!	lb(p+1,1) = huge(0.0)
!	ub(p+1,1) = r0
	
!	ctemp=linprog2(f,A1,b,lb,ub)
!	r = ctemp(p+1,1)
!	c(:,1) = ctemp(1:p,1)
!END FUNCTION

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

subroutine calculate_rho(m)
      real*8 :: a,b
      real*8 :: I0,S0,PI
      integer, intent(in) :: m
      real*8 :: fa, fb, temp, sol
      integer :: n,i
      interface
      function f(x)
      real*8, intent(in) :: x
      end function f
      end interface

	open(unit=1,file='I0.txt')
	open(unit = 2, file='truncated_normal.txt')
	open(unit = 9, file='rho.txt')
do i=1,5000
	!print*,i
	a = 0.2
	b = 20
	read(1,*) I0
	read(2,*) PI,PT
	!print*,i
	S0 = 0.9
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
      !print *,"    n        x(n)         f(x(n))"
      !print *," 1 ", b, fb
      !print *," 0 ", a, fa	
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
	!print*,i
	write(9,*) sol
enddo
	close(1)
	close(2)
	close(9)
end subroutine

function evaluate_g_inverse(x,I0,S0,PI) result(fx)
	real*8,intent(in) :: x,I0,S0,PI
	real*8 :: fx
	fx = I0 + S0 - PI - x*(log(S0) + 1 - log(x))
end function


SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock,clock2
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
            real	:: y
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
	    clock2=clock
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
	    
	    do while(clock .eq. clock2)
		CALL SYSTEM_CLOCK(COUNT=clock2)
	    end do
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

end module
