program gibbs

	implicit none

	integer :: n

	real (kind = selected_real_kind(8) ) :: kappaprev,lambdaprev,Sprev,Iprev,Rprev,zprev,rhoprev,betaprev

	real (kind = selected_real_kind(8) ) :: kappacurr,lambdacurr,Scurr,Icurr,Rcurr,zcurr,rhocurr,betacurr

	real (kind = selected_real_kind(8) ) :: thetaprev,yprev,thetacurr,ycurr
	!Initiele waarden voor bovenstaande parameters (prev)

	integer :: i

	do i=0,n
		!kappacurr = according to all previous values
		!lambdacurr = according to kappacurr en Sprev, Iprev ...
		!Scurr = according to kappacurr, lambdacurr, Iprev, Rprev ...
		!Icurr =
		!Rcurr =
		!zcurr =
		!rhocurr =
		!betacurr =

		!thetacurr =
		!ycurr =

	kappaprev = kappacurr
	lambdaprev = lambdacurr
	Sprev = Scurr
	Iprev = Icurr
	Rprev = Rcurr
	zprev = zcurr
	rhoprev = rhocurr
	betaprev = betacurr
	thetaprev = thetacurr
	yprev = ycurr
	enddo

end program
