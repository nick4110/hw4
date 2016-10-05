function [xgk,wgk]=GKIntP(nK)
% Gauss-Kronrod integration points in the interval [-1,1].
%
% Description
%     [#xgk#,#wgk#]=GKIntP(#nK#)
%     generates the #nK#-point Gauss-Kronrod quadrature rule. It calculates
%     the abscissas (#xgk#) and weights (#wgk#) for the Gauss-Kronrod
%     quadrature.
%
% Input arguments
%     #nK# (scalar) is the number of the integration points to be
%     calculated.
%
% Output arguments
%     #xgk# ([2*floor(#nK#/2)+1 x 1]) contains the coordinates of the
%     integration points.
%     #wgk# ([2*floor(#nK#/2)+1 x 1]) contains the weights of the
%     integration points.
%
% Parents (calling functions)
%     GKQuad > GKIntP
%
% Children (called functions)
%     GKIntP > JKelements
%     GKIntP > JRecCoeff
%

% __________________________________________________________________________
%% Copyright
%
%  (c) 2016 by George Papazafeiropoulos
%  Captain, Infrastructure Engineer, Hellenic Air Force
%  Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%
% Email: gpapazafeiropoulos@yahoo.gr
%
% Website: http://users.ntua.gr/gpapazaf/
%

N=floor(nK/2);
a=0;
b=0;
% calculate the recurrence coefficients for monic Jacobi polynomials in the
% interval [-1,1]
cd=JRecCoeff(2*N,a,b);
% transform cd (defined on [-1 1]) into ab (defined on [0 1])
ab(1:2*N,1)=(1+cd(1:2*N,1))./2;
ab(1,2)=cd(1,2)/2^(a+b+1);
ab(2:2*N,2)=cd(2:2*N,2)./4;
% calculate the alpha- and beta-elements in the Jacobi-Kronrod matrix
ab0=JKelements(N,ab);
% check if Gauss-Kronrod integration is possible
if (sum(ab0(:,2)>0)~=2*N+1)
    error('Gauss-Kronrod does not exist')
end
% form the Jacobi-Kronrod matrix
k=1:2*N;
J=diag(ab0(1:2*N+1,1),0)+diag(sqrt(ab0(k+1,2)),1)+diag(sqrt(ab0(k+1,2)),-1);
% calculate the integration points
[V,D]=eig(J);
d=diag(D);
e=ab0(1,2).*(V(1,:).^2);
[x,i]=sort(d);
% transform abscissas and weights defined on [0 1] back into [-1 1]
xgk=2*x-1;
wgk=2*e(i)';

