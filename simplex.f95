module simplex

!min cx s.t. Ax=b; x>=0
FUNCTION simplex_method(A,b,c,irule) RESULT (istatus,X,eta,iB,iN,xB)

REAL,DIMENSION(:,:),ALLOCATABLE :: X

[istatus,iB,iN,xB] = simplex_init(A,b,c)

IF ((istatus .EQ. 4) .OR. (istatus .EQ. 16)) THEN
	print *,'not feasible or initialization procedure had failed)
ELSE
	m = shape(A)(1)
	n = shape(A)(2)
	Binv = call inv(A(:,iB))
	
	DO WHILE (istatus .EQ. 0)
		[istatus,iB,iN,xB] = call simplex_step(A,b,c,iB,iN,xB,irule)
	END DO
	
	IF (istatus .EQ. -1) THEN
		X = allocate(X(n,1))
		X = 0.0d
		X(iB) = xB
		eta = c(iB)*xB
		istatus = 0
	ELSE
		istatus = 32
		print *,program is unbounded
	END IF
END IF

END FUNCTION


FUNCTION simplex_init(
	integer :: j
	REAL,DIMENSION(:,:),ALLOCATABLE :: E,E1,A1,c1,Binv1
	REAL,DIMENSION(:),ALLOCATABLE :: iB1,iN1
	integer :: m
	integer :: n
	
	m = shape(A)(1)
	n = shape(A)(2)
	
	ALLOCATE(E(m,m))
	ALLOCATE(E1(m,m))
	E = 0.0d
	E1 = 0.0d
	forall(j=1:m) 
		E(j,j)=1.0d
		E1(j,j)=1.0d
	end forall
	
	DO index=1,m
		IF b(index)<0 THEN
			E1(index,index) = -1
		END IF
	ENDDO
	
	ALLOCATE(A1(m,n+m))
	A1(:,1:n) = MATMUL(E,A)
	A1(:,(n+1):(n+m)) = E1
	
	b1 = MATMUL(E,b)
	
	ALLOCATE(c1(1,n+m))
	c1=1.0d
	forall(j=1:n)
		c1(j) = 0.0d
	end forall
	
	ALLOCATE(iB1(m))
	forall(j=1:m)
		iB1(j)=n+j
	end forall
	
	ALLOCATE(iN1(n))
	forall(j=1:n)
		iN1(j)=j
	end forall
	
	ALLOCATE(Binv1(m,m))
	Binv1 = 0.0d
	forall(j=1:m)
		Binv1(j,j)=1.0d
	end forall
	
	xB1 = b1
	irule = 1
	
	!OPMERKING: in matlab hier lege iB,iN,xB; nog kijken wat hiermee gedaan wordt
	istatus = 0
	
	DO WHILE .NOT.(istatus .EQ. -1)
		[istatus,iB1,iN1,xB1]=simplex_step(A1,b1,c1,iB1,iN1,xB1,irule)
	ENDDO
	
	!https://github.com/kkaushik7/Simplex/blob/master/src/simplex_init.m
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

end module
