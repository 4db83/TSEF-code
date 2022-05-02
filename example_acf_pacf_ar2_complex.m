% Script: example_acf_pacf_ar2_complex.m
% NOTE: To be able to run this code, you need the contents of the db.toolbox available from:
% https://github.com/4db83/db.toolbox/archive/refs/heads/main.zip. Unzip the contents locally
% to the same directory as this script, and then uncomment the following line below:
% addpath(genpath('./db.toolbox-main'))
clear;clc;clf;

% lag polynomial
a1 = 1.4; a2 = -0.85;
aL = [1 -a1 -a2];
% companion matrix Phi
Phi = [a1 a2;
        1 0];

% theoretical ACF and PACF 
acf_t  = acf0 (aL,1,50);
pacf_t = pacf0(aL,1,50);

% lag polynomial and characteristic/factored roots 
lag_roots     = roots(fliplr(aL));   
r_lag_roots   = real(lag_roots);
i_lag_roots   = imag(lag_roots);
% characteristic roots (eigenvalues of Phi)
char_roots    = roots(aL);
% char_roots    = eig(Phi); % or eig(Phi)
r_char_roots  = real(char_roots);
i_char_roots  = imag(char_roots);

fprintf(' Lag Polynomial roots are:           %2.4f +/-%2.4fi \n', [r_lag_roots(1) i_lag_roots(1)])
fprintf(' Modulus of roots is:                %2.4f \n', abs(lag_roots(1)));
fprintf(' Factored Polynomial roots are are:  %2.4f +/-%2.4fi \n', [r_char_roots(1) i_char_roots(1)])
fprintf(' Modulus of eigenvalues is:          %2.4f \n', abs(char_roots(1)));

% SIMULATE FROM AR(2) WITH COMPLEX ROOTS AND PLOT
seed(123);
T = 2e2;                % sample size
B = 2e2;                % burn-in period.
c = 0.3;                % constant
u = randn(B+T,1);       % WN errors.
x = zeros(B+T,1);       % space allocation, burn-in plus T
for t = 3:(B+T)
    x(t) = c+a1*x(t-1)+a2*x(t-2)+u(t);
end
% inbuild matlab command (faster), no loop needed, generates zero mean series, 
% so need to add the unconditional mean back into x.
% x = c/sum(aL) + filter(1,aL,u); 
x = x(B+1:end);                 % drop the burn-in period.

% PLOTS
% clear plotting area and set default linewidth
set(groot,'defaultLineLineWidth',1.5); 
xg  = linspace(0,1.5,1e3)';     % grid for plotting
fns = 15;         % font size for plots
stp = -1.2;       % subtitle position adjustment
dims = [.3 .20];  % subfigure dims
tspc = .34;       % topspace

% 2x2 grid of plots
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

% COMPUTING THE ACF FROM THE MA(INF) REPRESENTATION OF THE AR MODEL.
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