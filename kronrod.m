function [ x, w1, w2 ] = kronrod ( n )
%*****************************************************************************80
%
%% KRONROD adds N+1 points to an N-point Gaussian rule.
%
%  Discussion:
%
%    This subroutine calculates the abscissas and weights of the 2N+1
%    point Gauss Kronrod quadrature formula which is obtained from the
%    N point Gauss quadrature formula by the optimal addition of N+1 points.
%
%    The optimally added points are called Kronrod abscissas.  The
%    abscissas and weights for both the Gauss and Gauss Kronrod rules
%    are calculated for integration over the interval [-1,+1].
%
%    Since the quadrature formula is symmetric with respect to the origin,
%    only the nonnegative abscissas are calculated.
%
%    Note that the code published in Mathematics of Computation
%    omitted the definition of the variable which is here called COEF2.
%
%  Storage:
%
%    Given N, let M = ( N + 1 ) / 2.
%
%    The Gauss-Kronrod rule will include 2*N+1 points.  However, by symmetry,
%    only N + 1 of them need to be listed.
%
%    The arrays X, W1 and W2 contain the nonnegative abscissas in decreasing
%    order, and the weights of each abscissa in the Gauss-Kronrod and
%    Gauss rules respectively.  This means that about half the entries
%    in W2 are zero.
%
%    For instance, if N = 3, the output is:
%
%    I      X               W1              W2
%
%    1    0.960491        0.104656         0.000000
%    2    0.774597        0.268488         0.555556
%    3    0.434244        0.401397         0.000000
%    4    0.000000        0.450917         0.888889
%
%    and if N = 4, (notice that 0 is now a Kronrod abscissa)
%    the output is
%
%    I      X               W1              W2
%
%    1    0.976560        0.062977        0.000000
%    2    0.861136        0.170054        0.347855
%    3    0.640286        0.266798        0.000000
%    4    0.339981        0.326949        0.652145
%    5    0.000000        0.346443        0.000000
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    03 August 2010
%
%  Author:
%
%    Original FORTRAN77 version by Robert Piessens, Maria Branders.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Robert Piessens, Maria Branders,
%    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
%    of Gauss and Lobatto,
%    Mathematics of Computation,
%    Volume 28, Number 125, January 1974, pages 135-139.
%
%  Parameters:
%
%    Input, integer N, the order of the Gauss rule.
%
%    Input, real TOL, the requested absolute accuracy of the
%    abscissas.
%
%    Output, real X(N+1), the abscissas.
%
%    Output, real W1(N+1), the weights for the Gauss-Kronrod rule.
%
%    Output, real W2(N+1), the weights for the Gauss rule.
%
% (Modified slightly by Alex Townsend, October 2016.) 

    tol = 1e-15; 
  m = floor ( ( n + 1 ) / 2 );
  even = ( 2 * m == n );

  d = 2.0;
  an = 0.0;
  for k = 1 : n
    an = an + 1.0;
    d = d * an / ( an + 0.5 );
  end
%
%  Calculation of the Chebyshev coefficients of the orthogonal polynomial.
%
  tau(1) = ( an + 2.0 ) / ( an + an + 3.0 );
  b(m) = tau(1) - 1.0;
  ak = an;

  for l = 1 : m - 1

    ak = ak + 2.0;
    tau(l+1) = ( ( ak - 1.0 ) * ak ...
      - an * ( an + 1.0 ) ) * ( ak + 2.0 ) * tau(l) ...
      / ( ak * ( ( ak + 3.0 ) * ( ak + 2.0 ) ...
      - an * ( an + 1.0 ) ) );
    b(m-l) = tau(l+1);

    for ll = 1 : l
      b(m-l) = b(m-l) + tau(ll) * b(m-l+ll);
    end

  end

  b(m+1) = 1.0;
%
%  Calculation of approximate values for the abscissas.
%
  bb = sin ( 1.570796 / ( an + an + 1.0 ) );
  x1 = sqrt ( 1.0 - bb * bb );
  s = 2.0 * bb * x1;
  c = sqrt ( 1.0 - s * s );
  coef = 1.0 - ( 1.0 - 1.0 / an ) / ( 8.0 * an * an );
  xx = coef * x1;
%
%  Coefficient needed for weights.
%
%  COEF2 = 2^(2*n+1) * n! * n! / (2n+1)!
%
  coef2 = 2.0 / ( 2 * n + 1 );
  for i = 1 : n
    coef2 = coef2 * 4.0 * i / ( n + i );
  end
%
%  Calculation of the K-th abscissa (a Kronrod abscissa) and the
%  corresponding weight.
%
  for k = 1 : 2 : n

    [ xx, w1(k) ] = abwe1 ( n, m, tol, coef2, even, b, xx );
    w2(k) = 0.0;

    x(k) = xx;
    y = x1;
    x1 = y * c - bb * s;
    bb = y * s + bb * c;

    if ( k == n )
      xx = 0.0;
    else
      xx = coef * x1;
    end
%
%  Calculation of the K+1 abscissa (a Gaussian abscissa) and the
%  corresponding weights.
%
    [ xx, w1(k+1), w2(k+1) ] = abwe2 ( n, m, tol, coef2, even, b, xx );

    x(k+1) = xx;
    y = x1;
    x1 = y * c - bb * s;
    bb = y * s + bb * c;
    xx = coef * x1;

  end
%
%  If N is even, we have one more Kronrod abscissa to compute,
%  namely the origin.
%
  if ( even )
    xx = 0.0;
    [ xx, w1(n+1) ] = abwe1 ( n, m, tol, coef2, even, b, xx );
    w2(n+1) = 0.0;
    x(n+1) = xx;
  end
     x = [ -x x(end-1:-1:1) ]';
     w1 = [ w1 w1(end-1:-1:1) ];
     w2 = [ w2 w2(end-1:-1:1) ];
  return
end

function [ x, w ] = abwe1 ( n, m, tol, coef2, even, b, x )

%*****************************************************************************80
%
%% ABWE1 calculates a Kronrod abscissa and weight.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    03 August 2010
%
%  Author:
%
%    Original FORTRAN77 version by Robert Piessens, Maria Branders.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Robert Piessens, Maria Branders,
%    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
%    of Gauss and Lobatto,
%    Mathematics of Computation,
%    Volume 28, Number 125, January 1974, pages 135-139.
%
%  Parameters:
%
%    Input, integer N, the order of the Gauss rule.
%
%    Input, integer M, the value of ( N + 1 ) / 2.
%
%    Input, real TOL, the requested absolute accuracy of the
%    abscissas.
%
%    Input, real COEF2, a value needed to compute weights.
%
%    Input, logical EVEN, is TRUE if N is even.
%
%    Input, real B(M+1), the Chebyshev coefficients.
%
%    Input, real X, an estimate for the abscissa.
%
%    Output, real X, the abscissa.
%
%    Output, real W, the weight.
%
  if ( x == 0.0 )
    ka = 1;
  else
    ka = 0;
  end
%
%  Iterative process for the computation of a Kronrod abscissa.
%
  for iter = 1 : 50

    b1 = 0.0;
    b2 = b(m+1);
    yy = 4.0 * x * x - 2.0;
    d1 = 0.0;

    if ( even )
      ai = m + m + 1;
      d2 = ai * b(m+1);
      dif = 2.0;
    else
      ai = m + 1;
      d2 = 0.0;
      dif = 1.0;
    end

    for k = 1 : m
      ai = ai - dif;
      i = m - k + 1;
      b0 = b1;
      b1 = b2;
      d0 = d1;
      d1 = d2;
      b2 = yy * b1 - b0 + b(i);
      if ( ~even )
        i = i + 1;
      end
      d2 = yy * d1 - d0 + ai * b(i);
    end

    if ( even )
      f = x * ( b2 - b1 );
      fd = d2 + d1;
    else
      f = 0.5 * ( b2 - b0 );
      fd = 4.0 * x * d2;
    end
%
%  Newton correction.
%
    delta = f / fd;
    x = x - delta;

    if ( ka == 1 )
      break
    end

    if ( abs ( delta ) <= tol )
      ka = 1;
    end

  end
%
%  Catch non-convergence.
%
  if ( ka ~= 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'ABWE1 - Fatal error!\n' );
    fprintf ( 1, '  Iteration limit reached.\n' );
    fprintf ( 1, '  Last DELTA was %e\n', delta );
    error ( 'ABWE1 - Fatal error!' );
  end
%
%  Computation of the weight.
%
  d0 = 1.0;
  d1 = x;
  ai = 0.0;
  for k = 2 : n
    ai = ai + 1.0;
    d2 = ( ( ai + ai + 1.0 ) * x * d1 - ai * d0 ) / ( ai + 1.0 );
    d0 = d1;
    d1 = d2;
  end

  w = coef2 / ( fd * d2 );

  return
end

function [ x, w1, w2 ] = abwe2 ( n, m, tol, coef2, even, b, x )

%*****************************************************************************80
%
%% ABWE2 calculates a Gaussian abscissa and two weights.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    30 April 2013
%
%  Author:
%
%    Original FORTRAN77 version by Robert Piessens, Maria Branders.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Robert Piessens, Maria Branders,
%    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
%    of Gauss and Lobatto,
%    Mathematics of Computation,
%    Volume 28, Number 125, January 1974, pages 135-139.
%
%  Parameters:
%
%    Input, integer N, the order of the Gauss rule.
%
%    Input, integer M, the value of ( N + 1 ) / 2.
%
%    Input, real TOL, the requested absolute accuracy of the
%    abscissas.
%
%    Input, real COEF2, a value needed to compute weights.
%
%    Input, logical EVEN, is TRUE if N is even.
%
%    Input, real B(M+1), the Chebyshev coefficients.
%
%    Input, real X, an estimate for the abscissa.
%
%    Output, real X, the abscissa.
%
%    Output, real W1, the Gauss-Kronrod weight.
%
%    Output, real W2, the Gauss weight.
%
  if ( x == 0.0 )
    ka = 1;
  else
    ka = 0;
  end
%
%  Iterative process for the computation of a Gaussian abscissa.
%
  for iter = 1: 50

    p0 = 1.0;
    p1 = x;
    pd0 = 0.0;
    pd1 = 1.0;
%
%  When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
%
    if ( n <= 1 )
      if ( eps < abs ( x ) )
        p2 = ( 3.0 * x * x - 1.0 ) / 2.0;
        pd2 = 3.0 * x;
      else
        p2 = 3.0 * x;
        pd2 = 3.0;
      end
    end

    ai = 0.0;
    for k = 2 : n
      ai = ai + 1.0;
      p2 = ( ( ai + ai + 1.0 ) * x * p1 - ai * p0 ) / ( ai + 1.0 );
      pd2 = ( ( ai + ai + 1.0 ) * ( p1 + x * pd1 ) - ai * pd0 ) / ( ai + 1.0 );
      p0 = p1;
      p1 = p2;
      pd0 = pd1;
      pd1 = pd2;
    end
%
%  Newton correction.
%
    delta = p2 / pd2;
    x = x - delta;

    if ( ka == 1 )
      break
    end

    if ( abs ( delta ) <= tol )
      ka = 1;
    end

  end
%
%  Catch non-convergence.
%
  if ( ka ~= 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'ABWE2 - Fatal error!\n' );
    fprintf ( 1, '  Iteration limit reached.\n' );
    fprintf ( 1, '  Last DELTA was %e\n', delta );
    error ( 'ABWE2 - Fatal error!' );
  end
%
%  Computation of the weight.
%
  an = n;

  w2 = 2.0 / ( an * pd2 * p0 );

  p1 = 0.0;
  p2 = b(m+1);
  yy = 4.0 * x * x - 2.0;
  for k = 1 : m
    i = m - k + 1;
    p0 = p1;
    p1 = p2;
    p2 = yy * p1 - p0 + b(i);
  end

  if ( even )
    w1 = w2 + coef2 / ( pd2 * x * ( p2 - p1 ) );
  else
    w1 = w2 + 2.0 * coef2 / ( pd2 * ( p2 - p0 ) );
  end

  return
end
