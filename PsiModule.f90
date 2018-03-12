MODULE PsiModule
  IMPLICIT NONE
CONTAINS
  REAL*8 FUNCTION digamma ( x, ifault )

    !*****************************************************************************80
    !
    !! DIGAMMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    20 March 2016
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Jose Bernardo.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jose Bernardo,
    !    Algorithm AS 103:
    !    Psi ( Digamma ) Function,
    !    Applied Statistics,
    !    Volume 25, Number 3, 1976, pages 315-317.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the digamma function.
    !    0 < X.
    !
    !    Output, integer IFAULT, error flag.
    !    0, no error.
    !    1, X <= 0.
    !
    !    Output, real ( kind = 8 ) DIGAMMA, the value of the digamma function at X.
    !
    IMPLICIT NONE

    REAL ( kind = 8 ), PARAMETER :: c = 8.5D+00
    REAL ( kind = 8 ), PARAMETER :: euler_mascheroni = 0.57721566490153286060D+00
    REAL ( kind = 8 ) digamma
    INTEGER ( kind = 4 ) ifault
    REAL ( kind = 8 ) r
    REAL ( kind = 8 ) x
    REAL ( kind = 8 ) x2
    !
    !  Check the input.
    !
    IF ( x <= 0.0D+00 ) THEN
       digamma = 0.0D+00
       ifault = 1
       RETURN
    END IF
    !
    !  Initialize.
    !
    ifault = 0
    !
    !  Approximation for small argument.
    !
    IF ( x <= 0.000001D+00 ) THEN
       digamma = - euler_mascheroni - 1.0D+00 / x + 1.6449340668482264365D+00 * x
       RETURN
    END IF
    !
    !  Reduce to DIGAMA(X + N).
    !
    digamma = 0.0D+00
    x2 = x

    DO WHILE ( x2 < c )
       digamma = digamma - 1.0D+00 / x2
       x2 = x2 + 1.0D+00
    END DO
    !
    !  Use Stirling's (actually de Moivre's) expansion.
    !
    r = 1.0D+00 / x2

    digamma = digamma + LOG ( x2 ) - 0.5D+00 * r

    r = r * r

    digamma = digamma &
         - r * ( 1.0D+00 / 12.0D+00 &
         - r * ( 1.0D+00 / 120.0D+00 &
         - r * ( 1.0D+00 / 252.0D+00 &
         - r * ( 1.0D+00 / 240.0D+00 &
         - r * ( 1.0D+00 / 132.0D+00 ) ) ) ) )

    RETURN
  END FUNCTION digamma
  SUBROUTINE psi_values ( n_data, x, fx )

    !*****************************************************************************80
    !
    !! PSI_VALUES returns some values of the Psi or Digamma function.
    !
    !  Discussion:
    !
    !    In Mathematica, the function can be evaluated by:
    !
    !      PolyGamma[x]
    !
    !    or
    !
    !      PolyGamma[0,x]
    !
    !    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
    !
    !    PSI(1) = -Euler's constant.
    !
    !    PSI(X+1) = PSI(X) + 1 / X.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 August 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    National Bureau of Standards, 1964,
    !    ISBN: 0-486-61272-4,
    !    LC: QA47.A34.
    !
    !    Stephen Wolfram,
    !    The Mathematica Book,
    !    Fourth Edition,
    !    Cambridge University Press, 1999,
    !    ISBN: 0-521-64314-7,
    !    LC: QA76.95.W65.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = 8 ) X, the argument of the function.
    !
    !    Output, real ( kind = 8 ) FX, the value of the function.
    !
    IMPLICIT NONE

    INTEGER ( kind = 4 ), PARAMETER :: n_max = 11

    REAL ( kind = 8 ) fx
    REAL ( kind = 8 ), SAVE, DIMENSION ( n_max ) :: fx_vec = (/ &
         -0.5772156649015329D+00, &
         -0.4237549404110768D+00, &
         -0.2890398965921883D+00, &
         -0.1691908888667997D+00, &
         -0.6138454458511615D-01, &
         0.3648997397857652D-01, &
         0.1260474527734763D+00, &
         0.2085478748734940D+00, &
         0.2849914332938615D+00, &
         0.3561841611640597D+00, &
         0.4227843350984671D+00 /)
    INTEGER ( kind = 4 ) n_data
    REAL ( kind = 8 ) x
    REAL ( kind = 8 ), SAVE, DIMENSION ( n_max ) :: x_vec = (/ &
         1.0D+00, &
         1.1D+00, &
         1.2D+00, &
         1.3D+00, &
         1.4D+00, &
         1.5D+00, &
         1.6D+00, &
         1.7D+00, &
         1.8D+00, &
         1.9D+00, &
         2.0D+00 /)

    IF ( n_data < 0 ) THEN
       n_data = 0
    END IF

    n_data = n_data + 1

    IF ( n_max < n_data ) THEN
       n_data = 0
       x = 0.0D+00
       fx = 0.0D+00
    ELSE
       x = x_vec(n_data)
       fx = fx_vec(n_data)
    END IF

    RETURN
  END SUBROUTINE psi_values
  SUBROUTINE timestamp ( )

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    IMPLICIT NONE

    CHARACTER ( len = 8 ) ampm
    INTEGER ( kind = 4 ) d
    INTEGER ( kind = 4 ) h
    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) mm
    CHARACTER ( len = 9 ), PARAMETER, DIMENSION(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    INTEGER ( kind = 4 ) n
    INTEGER ( kind = 4 ) s
    INTEGER ( kind = 4 ) values(8)
    INTEGER ( kind = 4 ) y

    CALL DATE_AND_TIME ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    IF ( h < 12 ) THEN
       ampm = 'AM'
    ELSE IF ( h == 12 ) THEN
       IF ( n == 0 .AND. s == 0 ) THEN
          ampm = 'Noon'
       ELSE
          ampm = 'PM'
       END IF
    ELSE
       h = h - 12
       IF ( h < 12 ) THEN
          ampm = 'PM'
       ELSE IF ( h == 12 ) THEN
          IF ( n == 0 .AND. s == 0 ) THEN
             ampm = 'Midnight'
          ELSE
             ampm = 'AM'
          END IF
       END IF
    END IF

    WRITE ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
         d, TRIM ( month(m) ), y, h, ':', n, ':', s, '.', mm, TRIM ( ampm )

    RETURN
  END SUBROUTINE timestamp





END MODULE PsiModule
