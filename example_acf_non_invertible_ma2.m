%Script: example_ma2
clc;clear all;clf
% parameter values
b1 = -3.5 ; b2 = -2;
xg = linspace(0,1.5,1000)';     % grid for plotting
bL = [1 +b1 +b2];               % lag polynomial
FN = 13;                        % font size for plots

% lag roots and acf of non-invertible model
lag_roots   = roots(fliplr(bL));
acf_out     = acf0(1,bL,50);
bar(acf_out(2:end)) ;xlim([0 50]);          % ACF
setplot([.4 .7], FN);
setytick
%print2pdf('..\lectures\graphics\acf_ma2_ni')

% factored roots and acf of invertibel model 
fact_roots  = 1./lag_roots;
acf_out_i   = acf0(1,[1 -(fact_roots(1)+fact_roots(2)) fact_roots(1)*fact_roots(2)],50);
bar(acf_out_i(2:end));xlim([0 50]);          % PACF
setplot([.4 .7], FN);
setytick
%print2pdf('..\lectures\graphics\acf_ma2_i')

fprintf(' Polynomial roots are: %2.2f %2.2f \n', lag_roots);
fprintf('   Factored roots are: %2.2f %2.2f \n', fact_roots);

%% plot a simulated series from the MA(2)
seed(123);
T = 2e2;                        % sampel size
B = 1e2;                        % burn-in period.
c = 0.3;                        % constant
u = randn(B+T,1);               % WN errors.
x = zeros(B+T,1);               % space allocation, burn-in plus T
for t = 3:(B+T)
    x(t) = c+b1*u(t-1)+b2*u(t-2)+u(t);
end
%y = filter(bL,1,c+u);          % inbuild matlab command (50% faster), no loop needed.
x = x(B+1:end);                 % drop the burn-in period.
plot(x,'Linewidth',1)
setplot([.7 .4], 10);hold on;
hline(0,'k'); hold off;
