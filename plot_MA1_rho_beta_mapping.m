% Script: plot_MA1_rho_beta_mapping.m
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear;clc;clf;
% addpath(genpath('PATH-TO-FOLDER/db.toolbox'))

% create inline function
f = inline('b./(1+b.^2)','b');
x = linspace(-3,3,100);

% plot
set(groot,'defaultLineLineWidth',1.5); 
fns = 18;         % font size for plots

hold on;
  plot(x,f(x))
  hline(0,'k-'); % vline(0,'k-')
  vline(-1,'k--');vline(1,'k--')
  hline(0.4,'r-',[],[],'r');vline(0.5,'r-');vline(2,'r-')
hold off;		
box on; grid on;
setplot([.84 .28], fns);
setyticklabels([-.5:.1:.5], 1); 
set(gca,'GridLineStyle',':','GridAlpha',1/3);
tickshrink(.8)
setoutsideTicks
% setplot([.9 .6],13);

xlabel('$\beta_1$','Interpreter','Latex')
ylabel('$\rho(1)$','Interpreter','Latex')
moveylabel(-.05);  
movexlabel(-.05);  

% uncomment to print to pdf 
% print2pdf('plot_MA1_rho_beta_mapping','../graphics')
