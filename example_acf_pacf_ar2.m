% Script: example_acf_pacf_ar2.m
% NOTE: To be able to run this code, you need the contents of the db.toolbox available from:
% https://github.com/4db83/db.toolbox/archive/refs/heads/main.zip. Unzip the contents locally
% to the same directory as this script, and then uncomment the following line below:
% addpath(genpath('./db.toolbox-main'))
clear;clc;clf;

% parameter values
a1  = 1.5 ; a2 =-0.56;
xg  = linspace(0,1.5,1000)';     % grid for plotting
aL  = [1 -a1 -a2];               % lag polynomial
Phi = [a1 a2;1 0];
% theoretical acf and pacf 
acf_t  = acf0(aL,1,50);
pacf_t = pacf0(aL,1,50);
% clear plotting area and set default linewidth to 2
set(groot,'defaultLineLineWidth',1.5); 
fns = 15;         % font size for plots
stp = -1.2;       % subtitle position adjustment
dims = [.3 .20];  % subfigure dims
tspc = .34;       % topspace
% roots 
lag_roots = roots(fliplr(aL));  % roots of lag polynomial
char_roots = eig(Phi);          % characteristic roots

subplot(2,2,1); hold on;
  plot(xg,polyval(fliplr(aL),xg))
  hline(0,'-k'); vline(lag_roots(1),'--r'); vline(lag_roots(2),'--r')
hold off;
box on; grid on;
xlim([1 1.5]); ylim([-0.01 0.01]); 
setplot([.15 .6 dims], fns, 3);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
setoutsideTicks
tickshrink(.8)
subtitle('(a) Roots of lag polynomial', stp)

subplot(2,2,2); hold on;
  plot(xg,polyval(aL,xg));
  hline(0,'-k'); vline(char_roots(1),'--r'); vline(char_roots(2),'--r')
hold off;
box on; grid on;
ylim([-0.01 0.01]); xlim([.5 .9]);
setplot([.55 .6 dims], fns, 3);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
setoutsideTicks
tickshrink(.8)
subtitle('(b) Characteristic roots', stp)

subplot(2,2,3); % ACF
  bar(acf_t(2:end), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.15 tspc dims], fns);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(c) Theoretical ACF', stp)

subplot(2,2,4);
  bar(pacf_t(2:end), 'FaceColor', [.7 .8 1]);
box on; grid on;
setplot([.55 tspc dims], fns);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
subtitle('(d) Theoretical PACF', stp)

% print roots to screen
fprintf(' Lag polynomial roots are: %2.2f %2.2f \n', lag_roots)
fprintf(' Characteristic roots are: %2.2f %2.2f \n', char_roots)

% uncomment to print to pdf 
% print2pdf('example_acf_pacf_ar2','../graphics')













%EOF 