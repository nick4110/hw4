%% CS 4210 / MATH 4250
% Lecture 13: Adaptive quadrature
%
% Alex Townsend, 4th October 2016

%% 
% Subdivision for nonsmooth functions: 

f = @(x) abs(x).^3;   
subdivisionDiagram(f,1)    % We were lucky. 
set(gca,'fontsize',16)

%%
% Subdivision for nonsmooth functions: 

f = @(x) abs(x-pi/6).^3;
subdivisionDiagram(f,3)   % Standard luck. 
set(gca,'fontsize',16)

%%
% Subdivision for oscillatory functions: 
f = @(x) sin(500*x);
subdivisionDiagram(f,3)
set(gca,'fontsize',16)

%%
% Subdivision for functions with singularity just off the interval: 

f = @(x) sqrt(x+1.0001);
subdivisionDiagram(f,3)

%% Gauss quadrature and Gauss-Kronrod:
% Gauss-Kronrod is used in industry-level codes:

f = @(x) cos(10*x); 
exact = quadgk(f,-1,1);
exact

%%
n = 7; 
[x, w] = legpts( n );
G7 = w*f(x)

%%
[x, w1, w2] = kronrod( n );
K15 = w1*f(x)
% G7 = w2*f(x)

%%
errorEstimate = (200*abs(G7-K15)).^(1.5)

%% Refine the quadrature rule: 

n = 10; 
[x, w] = legpts( n );
G10 = w*f(x)

%%
[x, w1, w2] = kronrod( n );
K21 = w1*f(x)
G10 = w2*f(x)

%%
errorEstimate = (200*abs(G10-K21)).^(1.5)

%% QUADGK:

f = @(x) abs(x-pi/6);
[I, error_bound] = quadgk(f, -1, 1, 'RelTol', 1e-13)
exact = sum(chebfun(f,[-1 pi/6 1]))
abs( exact - I )

