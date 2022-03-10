%Script: example_ar2
% uncomment print2pdf to fprint to pdf using the print2pdf function.
clc;clear all;clf
%% parameter values
a1 = 1.5 ; a2 =-0.56;
xg = linspace(0,1.5,1000)';     % grid for plotting
aL = [1 -a1 -a2];               % lag polynomial
FN = 13;                        % font size for plots

%% lag polynomial roots
lag_roots = roots(fliplr(aL));  % roots of lag polynomial
plot(xg,polyval(fliplr(aL),xg),'Linewidth',1);hold on;
hline(0,'-k');vline(lag_roots(1),'--r');vline(lag_roots(2),'--r')
hold off;
ylim([-0.01 0.01]);
xlim([1 1.5]);
setplot([.4 .7], FN);
setytick(3)
%print2pdf('..\lectures\graphics\lag_roots')
fprintf(' Lag polynomial roots are: %2.2f %2.2f \n', lag_roots)

%% characteristic roots
F = [a1 a2;1 0];
char_roots = eig(F);
plot(xg,polyval(aL,xg),'Linewidth',1);hold on;
hline(0,'-k');vline(char_roots(1),'--r');vline(char_roots(2),'--r');
hold off;
ylim([-0.01 0.01])
xlim([.5 .9]);
setplot([.4 .7], FN);
setytick(3)
%print2pdf('..\lectures\graphics\char_roots')
fprintf(' Characteristic roots are: %2.2f %2.2f \n', char_roots)

%% plot theoretical acf and pacf
acf_out     = acf0(aL,1,50);
pacf_out    = pacf0(aL,1,50);
bar(acf_out(2:end)) ;xlim([0 50]);          % ACF
setplot([.4 .7], FN);
%print2pdf('..\lectures\graphics\acf_ar2')
bar(pacf_out(2:end));xlim([0 50]);          % PACF
setplot([.4 .7], FN);
%print2pdf('..\lectures\graphics\pacf_ar2')
