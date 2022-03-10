% Script: plot_loss_functions.m
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear;clc;
% addpath(genpath('d:/matlab.tools/db.toolbox'))
% parameter values
e  = linspace(-1,1,1e3)';

%	Quadratic Error Loss: 
QEL = @(e) (e.^2);
% Absolute Error Loss: 
AEL = @(e) (abs(e));
% LinEx Error Loss: 
b = 3; a = 1;
LINEX = @(e) ( b*(exp(a*e) - a*e-1) );
% LinLin Error Loss: 
b = 3; a = 1; I = e > 0;
LINLIN = @(e) ( b*abs(e).*I + a*abs(e).*(1-I)  );

% clear plotting area and set default linewidth to 2
clf; 
set(groot,'defaultLineLineWidth',2); 
Fns = 18; % font size
hold on; LG = [];
  LG(1) = plot(e,QEL(e));
  LG(2) = plot(e,AEL(e));
  LG(3) = plot(e,LINEX(e));
  LG(4) = plot(e,LINLIN(e));
vline(0,'k:');
hold off;
box on; grid on;
setplot([.80 .24], Fns, 1);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
% ylim([-.25 1.05])
setoutsideTicks
ylabel('$L(e)$', 'Interpreter','latex')
moveylabel(-.01)
xlabel('$e$', 'Interpreter','latex')
movexlabel(-.05)
add2yaxislabel
tickshrink(.6)
% add legend
legnames = { 'Quadratic',...
             'Absolute',...
             'LinEx',...
             'LinLin'}; 
legendflex(LG, legnames, 'fontsize', Fns - 1, 'anchor',[1 1],'Interpreter','Latex')

% uncomment to print to pdf 
% print2pdf('loss_functions_2022','../graphics')



























%EOF
