% Script: example_acf_non_invertible_ma2.m
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear;clc;clf;
% addpath(genpath('PATH-TO-FOLDER/db.toolbox'))

% parameter values
b1 = -3.5 ; b2 = -2;
xg = linspace(0,1.5,1000)';     % grid for plotting
bL = [1 +b1 +b2];                 % lag polynomial
% lag roots of non-invertible model
lag_roots = roots(fliplr(bL));
% factored roots (eigenvalues) to make invertibel model 
fct_roots = 1./lag_roots; % or roots(bL) but this reverses the order,
% bLplus = [1 -(fact_roots(1)+fact_roots(2)) fact_roots(1)*fact_roots(2)]
b1_plus = -( fct_roots(1) + 1/fct_roots(2) );
b2_plus =  ( fct_roots(1)*1/fct_roots(2) );
bL_plus = [1 +b1_plus +b2_plus]; 
% theoretical acf and pacf 
acf_t  = acf0( 1,bL,50);
pacf_t = pacf0(1,bL,50);
% theoretical acf and pacf of invertible process
acf_t_plus  = acf0( 1,bL_plus,50);
pacf_t_plus = pacf0(1,bL_plus,50);

fprintf(' Polynomial roots are: %2.2f %2.2f \n', lag_roots);
fprintf('   Factored roots are: %2.2f %2.2f \n', fct_roots);

% clear plotting area and set default linewidth
set(groot,'defaultLineLineWidth',1.5); 
fns = 15;         % font size for plots
stp = -1.2;       % subtitle position adjustment
dims = [.3 .20];  % subfigure dims
tspc = .34;       % topspace
ytcks = [-.2:.1:.3]; % yticks 

% plots
subplot(2,2,1);
  bar(acf_t(2:end), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.15 .6 dims], fns);
setyticklabels(ytcks, 1); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(a) ACF $\beta(z)=(1-3.5z-2z^2)$', stp,[],1)

subplot(2,2,2);
  bar(pacf_t(2:end), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.55 .6 dims], fns);
setyticklabels(ytcks, 1); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(b) PACF $\beta(z)=(1-3.5z-2z^2)$', stp,[],1)

subplot(2,2,3); % ACF
  bar(acf_t_plus(2:end), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.15 tspc dims], fns);
setyticklabels(ytcks, 1); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(c) ACF $\beta(z)=(1+0.25z-0.125z^2)$', stp,[],1)

subplot(2,2,4);
  bar(pacf_t_plus(2:end), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.55 tspc dims], fns);
setyticklabels(ytcks, 1); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(d) PACF $\beta(z)=(1+0.25z-0.125z^2)$', stp,[],1)

% uncomment to print to pdf 
% print2pdf('example_acf_non_invertible_ma2','../graphics')


%EOF

