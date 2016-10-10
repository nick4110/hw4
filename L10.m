%% CS 4210 / MATH 4250
% Lecture 10: Best approximation
%
% Alex Townsend, 22nd September 2016

%% Comparison of polynomial approximants:  
x = chebfun('x'); 
f = abs(x-.25); 
n = 20;
p_best = remez(f, n); 
p_interp = chebfun(f, n+1); 
p_lsq = polyfit(f, n); 
xk = linspace(-1,1,n+1);
cfs = cubic_spline(xk, f(xk));
plot(f-p_best, 'linewidth', 2 ), hold on
plot(f-p_interp, 'linewidth', 2 )
plot(f-p_lsq, 'linewidth', 2 )
x = linspace(-1,1,1000);
plot(x,f(x)-evaluate_spline(cfs,xk,x), 'linewidth', 2 )
set(gca, 'fontsize', 16)
legend('p_best','p_interp','p_lsq','spline')

%% Noisy functions

f = @(x) cos(x); 
f_noise = @(x) f(x) + 1e-2*randn(1,size(x,2)); 
x = linspace(-1,1,100); 
plot(x, f(x), 'linewidth', 2 ), hold on
plot(x, f_noise(x), 'linewidth', 2)
set(gca, 'fontsize', 16)
hold off

%% Polynomial interpolation
h = chebfun(f_noise);
plot(h, 'linewidth', 2)
set(gca, 'fontsize', 16)

%% Piecewise polynomials 
h = chebfun(f_noise, 'splitting', 'on');
plot( h, 'linewidth', 2 )
set(gca, 'fontsize', 16)

%% Splines
n = 1000; 
xk = linspace(-1,1,n+1);
cfs = cubic_spline(xk, f(xk));
plot(x,f(x)-evaluate_spline(cfs,xk,x), 'linewidth', 2 )

%% Height and weight data in pounds and inches. 
data = LoadBMIData;   % Data is of children in Hong Kong. 
                      % Adult BMI data may differ!!! 
height = data(:,2); 
weight = data(:,3);

BMI = 703*weight./height.^2; 
[mx, idx1] = max(BMI); 
[mn, idx2] = min(BMI); 

plot(height, weight,'k.','markersize',3), hold on
plot(height(idx1),weight(idx1),'.r','markersize',20)
plot(height(idx2),weight(idx2),'.b','markersize',20)
xlabel('Height in inches'), axis([60 75 60 180])
ylabel('Weight in lb') 
set(gca,'fontsize',16)

%% BMI curves (Healthy range)
t = linspace(60,75);
plot(t, 16/703*t.^2,'r-','linewidth',2)
plot(t, 23/703*t.^2,'r-','linewidth',2)
shg
%% Height is a normal curve:
hold off
histogram(height, 'Normalization', 'pdf') 
xlabel('Height in inches')
set(gca,'fontsize',16)

%% Medical BMI = weight/(C*height^2)
A = height.^2; 
x = (A'*A) \ (A'*weight); 

plot(height, weight,'k.','markersize',3), hold on
plot(height(idx1),weight(idx1),'.r','markersize',20)
plot(height(idx2),weight(idx2),'.b','markersize',20)
t = linspace(60,75);
plot(t, 16/703*t.^2,'r-','linewidth',2)
plot(t, 23/703*t.^2,'r-','linewidth',2)
xlabel('Height in inches') 
ylabel('Weight in lb') 
set(gca,'fontsize',16), axis([60 75 60 180])
t = linspace(min(height), max(height)); 
plot(t, x*t.^2,'b-','linewidth',3)
lsq_error = norm( weight - x.*height.^2, 2 ) 

%% What about: weight/( A + B*height + C*height^2 )
A = [ones(numel(height),1) height height.^2]; 
% [Q, R] = qr( A, 0 ); 
% y = Q'*weight; 
% x = R \ y;
x = (A'*A) \ (A'*weight); 

plot(height, weight,'k.','markersize',3), hold on
plot(height(idx1), weight(idx1),'.r','markersize',20)
plot(height(idx2), weight(idx2),'.b','markersize',20)
t = linspace(60,75);
plot(t, 16/703*t.^2,'r-','linewidth',2)
plot(t, 23/703*t.^2,'r-','linewidth',2)
xlabel('Height in inches') 
ylabel('Weight in lb') 
set(gca,'fontsize',16), axis([60 75 60 180])
t = linspace(min(height), max(height)); 
plot(t, x(1) + x(2)*t + x(3)*t.^2,'g-','linewidth',3)
lsq_error = norm( weight - A * x, 2 ) 

%% Or weight/( A + B*height + C*height^2+D*height^3 )
A = [ones(numel(height),1) height height.^2 height.^3]; 
% [Q, R] = qr( A, 0 ); 
% y = Q'*weight; 
% x = R \ y;
x = (A'*A) \ (A'*weight); 


plot(height, weight,'k.','markersize',3), hold on
plot(height(idx1), weight(idx1),'.r','markersize',20)
plot(height(idx2), weight(idx2),'.b','markersize',20)
t = linspace(60,75);
plot(t, 16/703*t.^2,'r-','linewidth',2)
plot(t, 23/703*t.^2,'r-','linewidth',2)
xlabel('Height in inches') 
ylabel('Weight in lb') 
set(gca,'fontsize',16), axis([60 75 60 180])
t = linspace(min(height), max(height)); 
plot(t, x(1) + x(2)*t + x(3)*t.^2 + x(4)*t.^3,'m-','linewidth',3)

lsq_error = norm( weight - A * x, 2 )  


%% BMI 
A = [height.^1.65]; 
% [Q, R] = qr( A, 0 ); 
% y = Q'*weight; 
% x = R \ y;
x = (A'*A) \ (A'*weight); 

plot(height, weight,'k.','markersize',3), hold on
plot(height(idx1), weight(idx1),'.r','markersize',20)
plot(height(idx2), weight(idx2),'.b','markersize',20)
t = linspace(60,75);
plot(t, 16/703*t.^2,'r-','linewidth',2)
plot(t, 23/703*t.^2,'r-','linewidth',2)
xlabel('Height in inches') 
ylabel('Weight in lb') 
set(gca,'fontsize',16), axis([60 75 60 180])
t = linspace(min(height), max(height)); 
plot(t, x(1)*t.^1.65,'m-','linewidth',3)

norm( weight - A * x, 2 ) 

%% Townsend's BMI formula for Hong Kong children :p 
%  BMI = weight/(A + Bheight^(1/2) + Cheight + Dheight^1.5 + Eheight^2)
% Make up your own formula: 

A = bsxfun(@power, height, 0:.5:2);
% [Q, R] = qr( A, 0 ); 
% y = Q'*weight; 
% x = R \ y;
x = (A'*A) \ (A'*weight); 

plot(height, weight,'k.','markersize',3), hold on
plot(height(idx1), weight(idx1),'.r','markersize',20)
plot(height(idx2), weight(idx2),'.b','markersize',20)
t = linspace(60,75);
plot(t, 16/703*t.^2,'r-','linewidth',2)
plot(t, 23/703*t.^2,'r-','linewidth',2)
xlabel('Height in inches') 
ylabel('Weight in lb') 
set(gca,'fontsize',16), axis([60 75 60 180])
t = linspace(min(height), max(height))'; 
plot(t, bsxfun(@power, t, 0:.5:2)*x ,'m-','linewidth',3)
lsq_error = norm( weight - A * x, 2 ) 
