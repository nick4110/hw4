function val = ajt253_myquad( k )
%<NET_ID>_MYQUAD   My quadrature command 
% 
% This function takes in an integer k between 1 and 30. The command is 
% designed to output the value for the kth integral from hw4. 
%
% Originally written by Alex Townsend, September 2016. 

% Need this for one of the functions: 
n = 50;
[jj, kk] = meshgrid(1:n);

% The list of functions from the homework 4: 
func = {
    @(x) exp(x)
    @(x) double( x >= 0.3 )
    @(x) x.*sin(1./x)
    @(x) log(1+x)
    @(x) 1 ./ cosh( 10 * (x - 0.2) * 2 ) + ...
         1 ./ cosh( 100 * (x - 0.4) * 4 ) + ...
         1 ./ cosh( 1000 * (x - 0.6) * 8 )
    @(x) 1 ./ (x.^2 + 1.005)
    @(x) cos(x) + 1e-2*randn(size(x))
    @(x) 25 * exp(-25*x)
    @(x) (exp(10*x)-1)./exp((101/10)*x)
    @(x) x.^2+x+1
    @(x) floor(exp(x))
    @(x) sqrt( x.^3 )
    @(x) 1 ./ sqrt(x)
    @(x) (x < 1) .* (x + 1) + ...
         (1 <= x & x <= 3) .* (3 - x) + ...
         (x > 3) * 2
    @(x) 1 ./ (1 + x.^4)
    @(x) sin(100 * pi * x) ./ (pi * x)
    @(x) log(x)
    @(x) cos(10000*acos(x))
    @(x) cos(100*acos(x))
    @(x) det( 1./(jj+kk-1) - x*eye(50) )    % Becareful you have to evaluate 
                                            % this function one value of "x"
                                            % at a time!  
    @(x) (x<=0).*sin(pi*x)+(x>=0).*cos(pi*x)
    @(x) cos(pi*x/2)/1e9 
    @(x) exp(-x)
    @(x) log(2*x)./log(x)
    @(x) (-1).^(floor(1./x))./(floor(1./x))
    @(x) sin(5*sin(x))
    @(x) exp(x.^2)
    @(x) sin(x)./log(x)
    @(x) 1./sqrt(1-x.^2)
    @(x) log(2*x)-log(2) 
    };

% Corresponding intervals from homework 4: 
interval = [ 0 1 ; 
              0 1 ; 
              0 1 ; 
              0 1e-12;
              0 1 ; 
              -1 1 ;
              -1 1 ; 
              0 10 ; 
              0 100 ; 
              0 1 ; 
              0 3; 
              0 1; 
              0 1; 
              0 5 ;
              0 1 ;
              0 1;
              0 1; 
              -1 1;
              -1 1;
              -1 1;
              -1 1;
              -1 1;
              0 inf;
              0 1/2;
              0 1;
              0 pi;
              -1 1;
              0 1;
              -.5 .5;
              0 1];

% Put into variables f, a, b: 
f = func{k}; 
a = interval(k,1);
b = interval(k,2); 

% Select which numerical integration scheme I will use based on k (below is just
% an example, function 21 may be unsuitable for Gauss quadrature.). There are also
% other integration methods that you should employ. You may also not want to use
% the trapezoidal rule.  This is just an example script:

if ( k == 1 || k == 21 )
    
    % DO GAUSSQUADRATURE 
    val = GaussQuadrature(f, a, b);
    
elseif (k == 20 || k == 3 ) 
    
    % DO COMPOSITE TRAPEZOIDAL RULE 
    val = compositeTrapezoidal(f, a, b);
    
elseif (k == 22 || k == 4 )
    % DO LEAST SQUARES FIT 
    val = leastSquaresFit( f, a, b);
else
    % I am going to miss these integrals: 
    val = NaN; 
end

end

function val = GaussQuadrature( f, a, b ) 
% GAUSSQUADRATURE(F, A, B) 

% This command computes the integral of F on [A, B] using Gaussian
% quadrature. 

end

function val = compositeTrapezoidal( f, a, b ) 
% COMPOSITETRAPEZOIDAL(F, A, B) 

% This command computes the integral of F on [A, B] using the composite
% trapezoidal rule. 

end

function val = leastSquaresFit( f, a, b ) 
% LEASTSQUARESFIT(F, A, B) 

% This command computes the integral of F on [A, B] using best least 
% squares fits. 

end
