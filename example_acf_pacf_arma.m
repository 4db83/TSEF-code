%Script: example_acf_pacf_arma models
clc;clear all;
% parameter values
aL = [1 -1.3 0.8 0.1];
bL = [1  0.4  -0.2];
FN = 13;                        % font size for plots
% factored roots (need to be less than 1 in mod)
ar_roots = roots(aL)
ma_roots = roots(bL)

%% plot of the theoretical acf and pacf
N = 50;
acf_out     = acf0(aL,bL,N);
pacf_out    = pacf0(aL,bL,N);
bar(acf_out(2:end)) ;xlim([0 N+.5]);          % ACF
ylim([-1 1]);
setplot([.4 .7], FN);
%print2pdf('..\lectures\graphics\acf_arma')
bar(pacf_out(2:end));xlim([0 N+.5]);          % PACF
setplot([.4 .7], FN);
ylim([-1 1]);
%print2pdf('..\lectures\graphics\pacf_arma')
