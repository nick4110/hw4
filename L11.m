%% CS 4210 / MATH 4250
% Lecture 11: Numerical integration / quadrature
%
% Alex Townsend, 27th September 2016

%% Clenshaw-Curtis quadrature:

n = 100; 
f = @(x) cos(10*x); 
[x, w] = chebpts(n+1);
quadrature = w*f(x)

%%
p = chebfun.interp1( x, f(x) );
interp = sum( p )

%%
exact = quadgk(f, -1, 1)

%% The convergence of Chebyshev coefficients: 
clf
x = chebfun('x');
f = abs(x-.3);
fc = chebfun(@(x) f(x),1e5);
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
plotcoeffs(fc,'loglog','.',MS,8), axis([1 1e5 1e-12 1])
xlabel('n',FS,12), ylabel('Chebyshev coefficient',FS,12)
nn = round(2.^(1:.5:16));
hold on, loglog(nn,.01*nn.^(-2),'--k',LW,2)
text(4e2,.5e-9,'n^{-2}',FS,18)

%% The convergence of quadrature rules: 

clf, exact = sum(f);
errc = [];
nn = round(2.^(1:.5:20));
for n = nn
  [s, w] = chebpts( n );
  Iclenshawcurtis = w*f(s);
  errc = [errc abs(Iclenshawcurtis-exact)];
end
loglog(nn,errc,'.-r',LW,1,MS,16), grid on
axis([1 1e6 1e-15 1]), hold on
xlabel('Npts',FS,12), ylabel('Error',FS,12)
loglog(nn,errc,'.-r',LW,1,MS,16)
title('Clenshaw-Curtis quadrature',FS,14)
loglog(nn,.01*nn.^(-2),'--k',LW,2)
text(7e2,.5e-9,'n^{-2}',FS,18)
