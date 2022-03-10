%% Example: YW estimation of a persistent AR(2) model.
clc;clear all;clf
%matlabpool 4;
T   = 2e2;              % sample size
B   = 2e2;              % burn-in
NS  = 1e3;              % number of simulations

% AR(2) parameter setting
C   = 0.3;
b1  = 0.5;
bL  = [1 b1];
fprintf('   Modulus of b(L) Roots of:     %2.6f   \n', abs(1/roots(bL))');

% simulate MA(1) process
x = filter(bL,1,C + randn(T+B,NS),0);
x = x(B+1:end,:);       % drop burn in 

% demean the series.
mx = x-repmat(mean(x),T,1);

% Create "by hand" the rho(1) and rho(2) terms. you can also use autocorr(y,2) function
p1 = @(y) (y(2:T)'*y(1:T-1)/(y'*y));

% space allocation
rho1     = zeros(NS,1);
ols_hat  = zeros(NS,1); 

tic;
options = optimset('Display','off','LargeScale','off');
parfor (j = 1:NS)
    % MoM Estimation of \hat\rho(1) 
    rho1(j,1) = p1(mx(:,j));
    
    % Minimizing the sum of squred errors over a grid
    ols_hat(j,1) = fminunc(@(bhat)((filter(1,[1 bhat],mx(:,j),0))'*(filter(1,[1 bhat],mx(:,j),0))/T),0,options);
end;
toc;

%
% MoM method stats: 
I_05    = rho1<0.5;                         % number of times \rho(1)<0.5
No_05   = mean(rho1>0.5);
fprintf('------------------------------------------------\n');
fprintf('Proportion of complex solutions    : % 2.6f \n',No_05);
fprintf('------------------------------------------------\n');
b1_mom      = (1-sqrt(1-4*rho1.^2))./(2*rho1);
mom_hat     = b1_mom(I_05);

asvar_mom   = @(b1)((1+b1^2+4*b1^4+b1^6+b1^8)/(1-b1^2)^2/T);
fprintf('Asymptotic Variance at true b1     : % 2.6f \n', asvar_mom(mean(mom_hat)));
fprintf('Monto Carlo Variance with T = % 2d : % 2.6f \n', [T var(mom_hat)]);

% Conditional MLE = OLS based results
fprintf('------------------------------------------------\n');
asvar_ols   = @(b1)((1-b1^2)/T);
fprintf('Asymptotic Variance at true b1     : % 2.6f \n', asvar_ols(mean(ols_hat)));
fprintf('Monto Carlo Variance with T = % 2d : % 2.6f \n', [T var(ols_hat)]);
%
% Plots of densities.
gx      = linspace(b1-.5,b1+.5,1000)';              % xgrid for Density estimate and plot
fxyw1   = ksdensity(mom_hat(:,1),gx);
fxols1  = ksdensity(ols_hat(:,1),gx);
plot(gx, fxyw1,'r');hold on;
plot(gx,fxols1,'b');hold off;
legend({'MoM';'MLE'},'FontSize',10);

% Print output to screen
fprintf(' ----------------------------------------\n')
fprintf('            MoM               OLS \n')
fprintf('       Mean (Variance)   Mean (Variance)\n')
fprintf(' ----------------------------------------\n')
printx = [mean(mom_hat)' var(mom_hat)' mean(ols_hat)' var(ols_hat)'];
labls  = ['b1'];
for j = 1:size(labls,1);
    fprintf(' %s:  % 2.4f (%2.4f)  % 2.4f  (%2.4f)\n', labls(j,:), printx(j,:));
end;
fprintf(' ----------------------------------------\n')
fprintf('       All Done Now!\n')

