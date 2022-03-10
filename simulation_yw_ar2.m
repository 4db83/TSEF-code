%% Example: YW estimation of a persistent AR(2) model.
clc;clear all;clf
% matlabpool 4;
T   = 2e2;              % sample size
B   = 2e2;              % burn-in
NS  = 1e5;              % number of simulations

% AR(2) parameter setting
C   = 0.3;
a1  = 1.5;  a2 = -0.75;
aL  = [1 -a1 -a2];
fprintf('   Modulus of a(L) Roots of:     %2.6f   %2.6f \n', abs(roots(fliplr(aL)))');

% simulate AR(2) process
Zi= C/sum(aL);
x = filter(1,aL,C + randn(T+B,NS), Zi*ones(cols(aL)-1,1));
x = x(B+1:end,:);       % drop burn in 

% demean the series.
mx = x-repmat(mean(x),T,1);

% Create "by hand" the rho(1) and rho(2) terms. you can also use autocorr(y,2) function
p1 = @(y) (y(2:T)'*y(1:T-1)/(y'*y));
p2 = @(y) (y(3:T)'*y(1:T-2)/(y'*y));

% space allocation
yw_hat  = zeros(NS,2);
ols_hat = zeros(NS,2);
Ins     = ones(T-2,1);      % create constant

tic;
parfor (i = 1:NS)
    % YW based version
    acf2 = [1; p1(mx(:,i)); p2(mx(:,i))];               
    % acf2 = autocorr(mx(:,i),2);
    yw_hat(i,:) = ([1 acf2(2); acf2(2) 1])\acf2(2:end);
    % OLS version by hand as in inv(X'X)*(X'Y)
    Y = x(3:T,i); X = [Ins x(2:T-1,i) x(1:T-2,i)];
    btmp = X\Y;
    ols_hat(i,:) = btmp(2:end);
end;
toc;
    
% Plots of densities.
gx      = linspace(a1-.7,a1+.3,1000)';              % xgrid for Density estimate and plot
fxyw1   = ksdensity(yw_hat(:,1),gx);
fxols1  = ksdensity(ols_hat(:,1),gx);
plot(gx, fxyw1,'r');hold on;
plot(gx,fxols1,'b');hold off;
legend({'YW';'OLS'},'FontSize',10);

% Print output to screen
fprintf(' ----------------------------------------\n')
fprintf('          Yule Walker         OLS \n')
fprintf('       Mean (Variance)   Mean (Variance)\n')
fprintf(' ----------------------------------------\n')
printx = [mean(yw_hat)' var(yw_hat)' mean(ols_hat)' var(ols_hat)'];
labls  = ['a1';'a2'];
for j = 1:2;
    fprintf(' %s:  % 2.4f (%2.4f)  % 2.4f  (%2.4f)\n', labls(j,:), printx(j,:));
end;
fprintf(' ----------------------------------------\n')
fprintf('       All Done Now!\n')
