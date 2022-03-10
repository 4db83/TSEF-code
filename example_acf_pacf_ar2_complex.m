% Script: example_acf_pacf_ar2_complex.m
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear;clc;clf;
% addpath(genpath('PATH-TO-FOLDER/db.toolbox'))

% parameter values
a1  = 1.4; a2 = -0.85;
xg  = linspace(0,1.5,1000)';     % grid for plotting

% roots of lag polynomial
aL = [1 -a1 -a2];               % lag polynomial
lag_roots  = roots(fliplr(aL));   
real_roots = real(lag_roots);
imag_roots = imag(lag_roots);
fprintf(' Lag polynomial roots are: %2.4f+/-%2.4f \n', [real_roots(1) imag_roots(1)])
fprintf(' Modulus of roots is : %2.4f \n\n', abs(lag_roots(1)));
% characteristic roots
Phi = [a1 a2;1 0];
char_roots = eig(Phi);            
real_roots = real(char_roots);
imag_roots = imag(char_roots);
fprintf(' Lag polynomial roots are: %2.4f+/-%2.4f \n', [real_roots(1) imag_roots(1)])
fprintf(' Modulus of eigenvalues is: %2.4f \n', abs(char_roots(1)));
% theoretical acf and pacf 
acf_t  = acf0( aL,1,50);
pacf_t = pacf0(aL,1,50);

% simulate from AR(2) with complex roots and plot
seed(123);
T = 2e2;                % sampel size
B = 2e2;                % burn-in period.
c = 0.3;                % constant
u = randn(B+T,1);       % WN errors.
x = zeros(B+T,1);       % space allocation, burn-in plus T
for t = 3:(B+T)
    x(t) = c+a1*x(t-1)+a2*x(t-2)+u(t);
end
% x = filter(1,aL,c+u);   % inbuild matlab command (faster), no loop needed.
x = x(B+1:end);         % drop the burn-in period.

% clear plotting area and set default linewidth
set(groot,'defaultLineLineWidth',1.5); 
fns = 15;         % font size for plots
stp = -1.2;       % subtitle position adjustment
dims = [.3 .20];  % subfigure dims
tspc = .34;       % topspace

% plots
subplot(2,2,1);
  bar(acf_t(2:end), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.15 .6 dims], fns);
setyticklabels([-1:.2:1], 1); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(a) Theoretical ACF', stp)

subplot(2,2,2);
  bar(pacf_t(2:end), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.55 .6 dims], fns);
setyticklabels([-1:.2:1], 1); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(b) Theoretical PACF', stp)

subplot(2,2,3:4);
  plot(x)
hline(0)
box on; grid on;
setplot([.15 tspc .70 .2], fns);
setyticklabels([-8:2:10], 0); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(c) Time series plot', stp)

% computing the ACF from the MA(inf) representation of the AR model.
N   = 1e3;              % approximation order
b		= arma2ma(aL,1,N);	% MA(inf) beta coefficients.
g0  = b'*b;							% variance.
gg  = zeros(N,1);
for i = 1:N
    gg(i) = b(1:end-i)'*b(1+i:end);	% acfs.
end
acf_MA_inf = gg./g0;

subplot(2,2,1);
  bar(acf_MA_inf(1:50), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.15 .6 dims], fns);
setyticklabels([-1:.2:1], 1); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(a) Theoretical ACF', stp)

% uncomment to print to pdf 
% print2pdf('example_acf_pacf_ar2_complex','../graphics')










%EOF 