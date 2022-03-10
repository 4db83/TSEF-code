% Script: example_acf_pacf_ar2.m
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear;clc;
% addpath(genpath('PATH-TO-FOLDER/db.toolbox'))

% parameter values
a1 = 1.5 ; a2 =-0.56;
xg = linspace(0,1.5,1000)';     % grid for plotting
aL = [1 -a1 -a2];               % lag polynomial
Fns = 13;                       % font size for plots

% lag polynomial roots
lag_roots = roots(fliplr(aL));  % roots of lag polynomial
subplot(2,1,1);
plot(xg,polyval(fliplr(aL),xg),'Linewidth',1);hold on;
hline(0,'-k');vline(lag_roots(1),'--r');vline(lag_roots(2),'--r')
hold off;
ylim([-0.01 0.01]);
xlim([1 1.5]);
setplot([.4 .7], Fns);
setytick(3)
%print2pdf('..\lectures\graphics\lag_roots')
fprintf(' Lag polynomial roots are: %2.2f %2.2f \n', lag_roots)

% characteristic roots
Phi = [a1 a2;1 0];
char_roots = eig(Phi);
subplot(2,1,1);
  plot(xg,polyval(aL,xg),'Linewidth',1);hold on;
hline(0,'-k');vline(char_roots(1),'--r');vline(char_roots(2),'--r');
hold off;
ylim([-0.01 0.01])
xlim([.5 .9]);
setplot([.4 .7], Fns);
setytick(3)
%print2pdf('..\lectures\graphics\char_roots')
fprintf(' Characteristic roots are: %2.2f %2.2f \n', char_roots)

%% plot theoretical acf and pacf
acf_out     = acf0(aL,1,50);
pacf_out    = pacf0(aL,1,50);
bar(acf_out(2:end)) ;xlim([0 50]);          % ACF
setplot([.4 .7], Fns);
%print2pdf('..\lectures\graphics\acf_ar2')
bar(pacf_out(2:end));xlim([0 50]);          % PACF
setplot([.4 .7], Fns);
%print2pdf('..\lectures\graphics\pacf_ar2')
