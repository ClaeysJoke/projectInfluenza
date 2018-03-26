program Test


  implicit none

  print *, testFunction(-huge(0.0D0),-huge(0.0D0))

  print *, (testFunction(-huge(0.0D0),-huge(0.0D0))<-huge(0.0D0))




contains
  real*8 function testFunction(x,y)
    real*8 x,y

    testFunction = x+y

  end function

end program
