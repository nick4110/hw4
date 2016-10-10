function val = ndb63_myquad( k )
%ndb63_MYQUAD   My quadrature command
%
% This function takes in an integer k between 1 and 30. The command is
% designed to output the value for the kth integral from hw4.

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
    %@(x) (exp(10*x)-1)./exp((101/10)*x)
    @(x) (exp(x)).^(-1/10) - (exp(x)).^(9/10)./(exp(11*x))
    @(x) x.^2+x+1
    @(x) floor(exp(x))
    @(x) x.^(3/2)
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
    0 1
    ];

% Put into variables f, a, b:
f = func{k};
a = interval(k,1);
b = interval(k,2);

% Select which numerical integration scheme I will use based on k (below is just
% an example, function 21 may be unsuitable for Gauss quadrature.). There are also
% other integration methods that you should employ. You may also not want to use
% the trapezoidal rule.  This is just an example script:

if ( k == 1  || k == 6  ||k==10   || k==15    || k==27 || k==29 )
    tic
    val = GaussQuadrature(200,f, a, b);
    toc
elseif (  k==8 )
    tic
    val = GaussQuadrature(50,f, a, b);
    toc
elseif (  k==9 || k==22 || k==26)
    tic
    val = GaussQuadrature(188,f, a, b);
    toc
elseif (  k==12 )
    tic
    val = GaussQuadrature(750,f, a, b);
    toc
elseif (  k==19 )
    tic
    val = GaussQuadrature(300,f, a, b);
    toc
elseif (  k== 16 )
    tic
    val = GaussQuadrature(261,f, a, b);
    toc
elseif ( k==23)
    tic
    val = GaussQuadrature(195, f, a, 40);
    toc
elseif ( k==23)
    tic
    val = AdaptiveSimpsonwithendinfcheck(f, a, b, 10^(-1));
    toc
elseif ( k==2 || k==4 || k==14)
    tic
    val = AdaptiveSimpson(f, a, b, 10^(-16));
    toc
else
    fprintf('User has not yet defined a quadrature');
    val = NaN;
end

end

function val = GaussQuadrature(p, f, a, b )
[x,c]=legpts(p);
d=b-a;
c=c*d/2;
x=(d*x+b+a)/2;
val = c*f(x);
end

%{
function val = ClenshawCurtis(p, f, a, b )
format long;
[x,c]=chebpts(p);
d=b-a;
c=c*d/2;
x=(d*x+b+a)/2;
val = c*f(x);
end
%}

function s0 = AdaptiveSimpson(f, a, b, tol)

format long;
if (a==inf) 
    a=2^(52);
end
if (a==-inf) 
    a=-2^(52);
end
if (b==inf) 
    b=2^(52);
end
if (b==-inf) 
    b=-2^(52);
end
m=(a+b)/2;
fm=f(m);
fa=f(a);
fb=f(b);
midam=(m+a)/2;
midmb=(m+b)/2;
fmidam=f(midam);
fmidmb=f(midmb);



if (isnan(fm)==1 || isfinite(fm)==0)
    fm=f(m+2^(-53));
end
if (isnan(fa)==1 || isfinite(fa)==0)
    fa=f(a+2^(-53));
end
if (isnan(fm)==1 || isfinite(fb)==0)
    fb=f(b-2^(-53));
end
if (isnan(fmidam)==1 || isfinite(fmidam)==0)
    fmidam=f(midam+2^(-53));
end
if (isnan(fmidmb)==1 || isfinite(fmidmb)==0)
    fmidmb=f(midmb+2^(-53));
end

%fm
%fa
%fb
%fmidam
%fmidmb

s0= (b-a)/6 * (fa+4*fm+fb);
s1= (m-a)/6 * (fa+4*fmidam+fm);
s2= (b-m)/6 * (fm+4*fmidmb+fb);
errorestimate = abs(1/15* (s1+s2-s0));
if (errorestimate < tol)
    s0=s1+s2;
else
    s0=AdaptiveSimpson(f,a,m,tol/2) + AdaptiveSimpson(f,m,b,tol/2);
end
end

function s0 = AdaptiveSimpsonwithendinfcheck(f, a, b, tol)

format long;
if (a==inf) 
    a=2^(52);
end
if (a==-inf) 
    a=-2^(52);
end
if (b==-inf) 
    b=-2^(52);
end
if (b==inf) 
    b=2^(52);
end

m=(a+b)/2;
fm=f(m);
fa=f(a);
fb=f(b);
midam=(m+a)/2;
midmb=(m+b)/2;
fmidam=f(midam);
fmidmb=f(midmb);

s0= (b-a)/6 * (fa+4*fm+fb);
s1= (m-a)/6 * (fa+4*fmidam+fm);
s2= (b-m)/6 * (fm+4*fmidmb+fb);
errorestimate = abs(1/15* (s1+s2-s0));
if (errorestimate < tol)
    s0=s1+s2;
else
    s0=AdaptiveSimpsonwithendinfcheck(f,a,m,tol/2) + AdaptiveSimpsonwithendinfcheck(f,m,b,tol/2);
end
end

%{
function s0 = AdaptiveGK(f, a, b, tol)

format long;
m=(a+b)/2;

[s0,errbnd]= quadgk(f,a,b,10^(-1);
s1= quadgk(f,a,m);
s2= quadgk(f,m,b);
errorestimate = abs(errbnd);
if (errorestimate < tol)
    s0=s1+s2;
else
    s0=AdaptiveGK(f,a,m,tol/2) + AdaptiveGK(f,m,b,tol/2);
end
end
%}

%{
function s0 = RectangleRule(f, a, b, numrect)
format long;
sum = 0;
N=numrect;
h=(b-a)/N;
for j=0:N-1
    sum=sum+f(a+j*h);
end
s0=h*sum;
end
%}

%{
function val = GaussKronrod(p, f, a, b )
format long;
[x,c]=GKIntP(p);
d=b-a;
c=c*d/2;
x=(d*x+b+a)/2;
g=f(x);
val = sum(c.*g);
end
%}

%{
function val = LeastSquaresFunc(f,a,b,p)
format long;
nodes=chebpts(p);
y=f(nodes);
A=zeros(p,p);
for i=1:p
    for j=1:p
        A(i,j)=nodes(i,1)^(j-1);
    end
end
weights=(A'*A) \ (A'*y);

clenprec=p;
% Run Clenshaw Curtis
[x,c]=chebpts(clenprec);
d=b-a;
c=c*d/2;
x=(d*x+b+a)/2;

fx = zeros(clenprec,1);
x
weights
for i = 1:clenprec
    fx(i,1)=x(i,1)^(i-1)*weights(i,1);
end   
fx
val=c*(fx);
end
%}


