%% CS 4210 / MATH 4250
% Lecture 12: Gauss quadrature
%
% Alex Townsend, 29th September 2016

%% Legendre polynomials: 
P = legpoly(0:5); 
plot(P, 'linewidth', 2)
set(gca,'fontsize',16)
xlabel('x','fontsize',16)

%% Orthogonal:
% The Legendre polynomials are orthogonal: 
P10 = legpoly(10); 
P101 = legpoly(101); 
sum(P10.*P10)        % = 2/(2*10+1)
sum(P101.*P101)      % = 2/(2*101+1)
sum(P10.*P101)       % = 0

%% The Gauss weights: 
% The Gauss-Legendre weights are
% 
%   w_k = \int_{-1}^1 l_k(x) dx, 
%
% where l_k(x) is the kth Lagrange polynomial. 

n = 10; 
[x, w] = legpts( n ); 
k = 4;   % Pick a k.  
lk = @(y) prod( y - x([1:k-1 k+1:end]) )./prod( x(k) - x([1:k-1 k+1:end]) );
lk = chebfun(lk,'vectorize');

% These two should be the same: 
sum(lk)
w(k)

%% For analytical functions, occasionally Gauss is a little better than CC:

f = @(x) x.*sin(2*exp(2*sin(2*exp(2*x))));

clf, exact = sum(chebfun(f));
errc = []; errg = [];
nn = round(2.^(1:.01:9));
for n = nn
  [x, w] = legpts( n ); 
  Igauss = w*f(x);
  errg = [ errg abs(Igauss-exact)];
  [s, w] = chebpts( n );
  Iclenshawcurtis = w*f(s);
  errc = [errc abs(Iclenshawcurtis-exact)];
end
errc(errc==0) = 1e-16;
errg(errg==0) = 1e-16;

semilogy(nn,errc,'.k-','linewidth',1,'markersize',16), grid on, hold on,
semilogy(nn,errg,'.-','linewidth',1,'markersize',16), grid on
xlabel('n','fontsize',12), ylabel('Error',FS,12)
loglog(nn,errc,'.-r','linewidth',1,'markersize',16)
title('Clenshaw-Curtis versus Gauss','fontsize',14)
set(gca,'fontsize',16)
legend('CC','Gauss')

%% Differential functions they are about the same

f = @(x) abs(x-.3); 

clf, exact = sum(chebfun(f,'splitting','on'));
errc = []; errg = [];
nn = round(2.^(1:.5:21));
for n = nn
  [x, w] = legpts( n ); 
  Igauss = w*f(x);
  errg = [ errg abs(Igauss-exact)];
  [s, w] = chebpts( n );
  Iclenshawcurtis = w*f(s);
  errc = [errc abs(Iclenshawcurtis-exact)];
end
errc(errc==0) = 1e-16;
errg(errg==0) = 1e-16;

loglog(nn,errc,'.-','linewidth',1,'markersize',16), grid on, hold on,
loglog(nn,errg,'.-','linewidth',1,'markersize',16)
xlabel('n','fontsize',12), ylabel('Error',FS,12)
title('Clenshaw-Curtis versus Gauss','fontsize',16)
set(gca,'fontsize',16)
legend('CC','Gauss')
