MODULE Normal
implicit none

contains
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

 END MODULE
