

!************************************************************
!			LINPROG				    *
!							    *
!************************************************************
! The following code tries to solve the linear programming  *
! problem						    *
!			minimize (f^t)*x		    *
!			s.t. A*x<=b			    *
! when f,A,b are given					    *
! or the problem					    *
!			minimize (f^t)*x		    *
!			s.t. A*x<=b && lb<=x<=ub	    *
! when f,A,b,lb,ub are given.				    *
! NO EQUALITY CONSTRAINTS				    *
! It uses the simplex algorithm				    *
!							    *
! In what follows suffix _1 will refer to the first problem,*
! suffix _2 refers to the second			    *
! Both are converted to the standard form		    *
! 			minimize (g^t)*y		    *
!			s.t. C*y=d && y>=0		    *
!							    *
!************************************************************




!************************************************************
!			PART ONE: SETUP
!
!************************************************************


module SETUP

implicit none

!Subroutine to check validity of the problem
!The flag is true if no error message is returned,
!meaning the problem is admissible
 interface Check
	 module procedure Check_1, Check_2
 end interface
 
!Subroutine converting to the form
!	minimize (f^t)*x
!	s.t. A*x=b, 0<=x
interface Convert
	module procedure Convert_1,Convert_2
end interface

contains


subroutine Check_1(f, A, b, nineqcstr, nvars)
 integer							:: nineqcstr, nvars
 real(kind=selected_real_kind(8)), dimension(nineqcstr,nvars) 	:: A
 real(kind=selected_real_kind(8))				:: f(nvars)
 real(kind=selected_real_kind(8))				:: b(nineqcstr)
 
 if (nvars .eq. 0) then
  stop "Badly defined problem"
 else if (size(A).eq. 0) then
  stop "faulty input matrix A"
 else if ((size(A,1).ne. nineqcstr) .or. (size(A,2).ne. nvars)) then
  stop "Argument mismatch: size of A does not correspond to number of inequalities or number of variables."
 else if(size(b).eq. 0) then
  stop "empty input bound b"
 end if
end subroutine

subroutine Check_2(f, A, b, nineqcstr, nvars, lb, ub)

 integer							:: nineqcstr, nvars,i
 real(kind=selected_real_kind(8)), dimension(nineqcstr,nvars) 	:: A
 real(kind=selected_real_kind(8)), dimension(nvars) 		:: f,lb,ub
 real(kind=selected_real_kind(8)), dimension(nineqcstr) 	:: b
 logical 							:: msg

 do i=1,nvars
 msg=(lb(i)<=ub(i))
 end do

 if (nvars .eq. 0) then
  stop "Badly defined problem"
 else if (size(A).eq. 0) then
  stop "faulty input matrix A"
 else if(size(b).eq. 0) then
  stop "empty input bound b"
 else if (.not.msg) then
  stop "boundaries infeasable"
 end if

end subroutine

!Convert given A0, f0, b0 to (A0|-A0|Id), (f0^t|-f0^t|0)^t, b0
!
!
!
subroutine Convert_1(f0, A0, b0, f, A, b, nineqcstr, nvars)
integer								:: nineqcstr,nvars,i
real(kind=selected_real_kind(8)), dimension(nineqcstr,nvars) 	:: A0
real(kind=selected_real_kind(8)), dimension(nvars)		:: f0
real(kind=selected_real_kind(8)), dimension(nineqcstr)		:: b0
real(kind=selected_real_kind(8)), allocatable 			:: A(:,:)
real(kind=selected_real_kind(8)), allocatable			:: f(:)
real(kind=selected_real_kind(8)), allocatable			:: b(:)
allocate(A(nineqcstr,(2*nvars)+nineqcstr))
allocate(f(2*nvars+nineqcstr))
allocate(b(nineqcstr))
A(:,:)=0.0
A(:,1:nvars)=A0(:,:)
A(:,(nvars+1):(2*nvars))=-1*A0(:,:)
do i=1,nineqcstr
A(i,(2*nvars+i))=1.0
end do
f(1:nvars)=f0(1:nvars)
f(nvars+1:2*nvars)=-1*f0(1:nvars)
f(2*nvars+1:2*nvars+nineqcstr)=0
b=b0
end subroutine

!Same as above, only with preamble of converting x to y=x-lb,
!and without doubling A0 to (A0 -A0) 
!f remains unchanged in this first step
!	  (A0 Id 0 )
!	A=(Id 0  Id)
subroutine Convert_2(f0, A0,b0, f, A, b, nineqcstr, nvars,lb,ub)
integer								:: nineqcstr,nvars,i
real(kind=selected_real_kind(8)), dimension(nineqcstr,nvars) 	:: A0
real(kind=selected_real_kind(8))				:: f0(nvars),lb(nvars),ub(nvars)
real(kind=selected_real_kind(8))				:: b0(nineqcstr)
real(kind=selected_real_kind(8)), allocatable			:: b(:)
real(kind=selected_real_kind(8)), allocatable			:: A(:,:)
real(kind=selected_real_kind(8)), allocatable			:: f(:)
allocate(f(nvars+nvars+nineqcstr))
allocate(b(nvars+nineqcstr))
allocate(A(nineqcstr+nvars,nvars+nineqcstr+nvars))
f(:)=0.0
f(1:nvars)=f0(1:nvars)
!FIRST STEP
A(:,:)=0.0
A(1:nineqcstr,1:nvars)=A0
do i=1,nvars
A(nineqcstr+i,i)=1.0
end do
b(1:nineqcstr)=b0-MATMUL(A0,lb)
b(nineqcstr+1:nineqcstr+nvars)=ub-lb

do i=1,nineqcstr+nvars
A(i,nvars+i)=1.0
end do
end subroutine
end module






!************************************************************
!		PART TWO: SOLVING STANDARD PROBLEM
!
!************************************************************


subroutine solve(f,A,b, Arows, Acolumns, SOL, nvars)

integer 			:: nvars
integer				:: Arows, Acolumns,j,n,m,i,n0
integer				:: nonbasic(nvars), basic(Arows), L(Arows-1)
real(kind=selected_real_kind(8)):: f(Acolumns), b(Arows), A(Arows,Acolumns), AT(Acolumns,Arows),SOL(nvars),t
real(kind=selected_real_kind(8)):: Z(Acolumns)

SOL(:)=0.0


AT=TRANSPOSE(A)

do i=nvars+1,Acolumns
basic(i-nvars)=i
end do
Z=MATMUL(AT,f(basic))
j=0
n0=0
i=0
do while(MAXVAL(f-Z)>0)
j=j+1
m=MAXLOC(f-Z,1 )
n=MINLOC(b/A(:,m),1, b/A(:,m) .ge. 0 )

if (j>100) then
 stop "LOOP FORMS, perterb problem slightly"
end if
 n0=n

!PIVOT OPERATION
b(n)=(1/A(n,m))*b(n)
A(n,:)=(1/A(n,m))*A(n,:)
do i=1,Arows
if (i .ne. n) then
b(i)=b(i)-A(i,m)*b(n)
A(i,:)=A(i,:)-A(i,m)*A(n,:)
end if
end do

AT=TRANSPOSE(A)
basic(n)=m

Z=MATMUL(AT,f(basic))
end do
do i=1,nvars
do j = 1, Arows
    if (basic(j) .eq. i) then
       Sol(i) = b(j)
      exit
    endif
end do
end do
end subroutine





!************************************************************
!			     MAIN
!
!************************************************************

program main
use SETUP
 real(kind=selected_real_kind(8)), dimension(6,2)	:: A0
 real(kind=selected_real_kind(8))			:: f0(2),lb(2),ub(2)
 real(kind=selected_real_kind(8)),allocatable 		:: f(:), A(:,:), b(:)
 real(kind=selected_real_kind(8)) 			:: b0(6)
 real(kind=selected_real_kind(8)) 			:: SOL(2),Time(100)
 real(kind=selected_real_kind(8))			:: TOLFUN,T
 INTEGER						:: T1,T2, count_rate, count_max,S1,S2
 integer						:: Bounded
 !DEFINE THE PROBLEM
 Bounded=1 !Set to 0 for unbounded, 1 for bounded problem
 lb(1)=-1
 lb(2)=-1
 ub(1)=4
 ub(2)=5
 !f(:)=0.0
 f0(1)=1.0
 f0(2)=.333
 !f(3)=-1.0
 !f(4)=-3.333
 b0(1)=2.0
 b0(2)=1.0
 b0(3)=2.0
 b0(4)=1.0
 b0(5)=-1.0
 b0(6)=2.0
 nvars=2
 nineqcstr=6
 A0(:,:)=0.0
 A0(1,1)=1.0
 A0(1,2)=1.0
 A0(2,1)=1.0
 A0(2,2)=0.25
 A0(3,1)=1.0
 A0(3,2)=-1.0
 A0(4,1)=-0.25
 A0(4,2)=-1.0
 A0(5,1)=-1.0
 A0(5,2)=-1.0
 A0(6,1)=-1.0
 A0(6,2)=1

 !A(1:6,3:4)=-A(1:6,1:2)
 !do i=1,6
 !A(i,4+i)=1.0
 !end do
 Time(:)=0
 S1=0
 S2=0
 do i=1,100
 
 call SYSTEM_CLOCK(S1,count_rate,count_max)
 call Convert(f0,A0,b0,f,A,b,nineqcstr,nvars,lb,ub)
 A=A+0.000000000001 !PERTURBATION TO AVOID SINGULARITY (Turn off unless program recommends it)
 call solve(f,A,b, 8,10,SOL,2)
 SOL=SOl+Bounded*lb
 print*, SOL
 call SYSTEM_CLOCK(S2,count_rate,count_max)
 deallocate(A)
 deallocate(f)
 deallocate(b)
 Time(i)=real(S2-S1,kind=selected_real_kind(8))/count_rate
 end do

 print *, SUM(Time) /100
 

end
