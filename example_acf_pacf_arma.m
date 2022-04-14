% Script: example_acf_pacf_arma.m
% NOTE: To be able to run this code, you need the contents of the db.toolbox available from:
% https://github.com/4db83/db.toolbox/archive/refs/heads/main.zip. Unzip the contents locally
% to the same directory as this script, and then uncomment the following line below:
% addpath(genpath('./db.toolbox-main'))
clear;clc;clf;

% parameter values
aL = [1 -1.3 0.8 0.1];
bL = [1  0.4  -0.2];

% theoretical acf and pacf 
N = 50;
acf_t  = acf0(aL,bL,N);
pacf_t = pacf0(aL,bL,N);

% clear plotting area and set default linewidth to 2
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

% uncomment to print to pdf 
% print2pdf('example_acf_pacf_arma','../graphics')
