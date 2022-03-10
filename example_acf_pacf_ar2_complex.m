%Script: example_ar2_complex
% uncomment print2pdf to fprint to pdf using the print2pdf function.
clc;clear all;clf
% parameter values
a1 = 1.4 ; a2 = -0.85;
xg = linspace(0,1.5,1000)';     % grid for plotting
aL = [1 -a1 -a2];               % lag polynomial
FN = 13;                        % font size for plots

%% lag polynomial roots
lag_roots = roots(fliplr(aL));  % roots of lag polynomial
real_roots = real(lag_roots);
imag_roots = imag(lag_roots);
fprintf(' Lag polynomial roots are: %2.4f+/-%2.4f \n', [real_roots(1) imag_roots(1)])
fprintf(' Modulus of roots is : %2.4f \n\n', abs(lag_roots(1)));

%% characteristic roots
F = [a1 a2;1 0];
char_roots = eig(F);
real_roots = real(char_roots);
imag_roots = imag(char_roots);
fprintf(' Lag polynomial roots are: %2.4f+/-%2.4f \n', [real_roots(1) imag_roots(1)])
fprintf(' Modulus of eigenvalues is: %2.4f \n', abs(char_roots(1)));

%% plot theoretical acf and pacf
acf_out     = acf0(aL,1,50);
pacf_out    = pacf0(aL,1,50);
bar(acf_out(2:end)) ;xlim([0 50]);          % ACF
setplot([.4 .7], FN);
%print2pdf('..\lectures\graphics\acf_ar2_complex')
bar(pacf_out(2:end));xlim([0 50]);          % PACF
setplot([.4 .7], FN);
%print2pdf('..\lectures\graphics\pacf_ar2_complex')

%% plot a simulated series from the AR(2) with complex roots
seed(123);
T = 2e2;                % sampel size
B = 2e2;                % burn-in period.
c = 0.3;                % constant
u = randn(B+T,1);         % WN errors.
x = zeros(B+T,1);       % space allocation, burn-in plus T
for t = 3:(B+T)
    x(t) = c+a1*x(t-1)+a2*x(t-2)+u(t);
end
%  x = filter(1,aL,c+u);% inbuild matlab command (50% faster), no loop needed.
x = x(B+1:end);         % drop the burn-in period.
plot(x,'Linewidth',1)
setplot([.7 .4], 10);hold on;
hline(0,'k'); hold off;
%print2pdf('..\lectures\graphics\ar2_complex')

%% computing the ACF from the MA(inf) representation of the AR model.
A   = 1e3;              % approximation order
b		= arma2ma(aL,1,A);	% MA(inf) beta coefficients.
g0  = b'*b;							% variance.
gg  = zeros(A,1);
for i = 1:A;
    gg(i) = b(1:end-i)'*b(1+i:end);	% acfs.
end;
pp = gg./g0;
bar(pp(1:50));xlim([0 50])
setplot([.4 .7], FN);

