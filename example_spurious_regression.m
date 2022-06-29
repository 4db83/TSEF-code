% Script: example_spurious_regression.m
% Example of what happens to OLS under spurious regression with I(1) veriables.
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear; clc;
% addpath(genpath('PATH-TO-FOLDER/db.toolbox'))
% some controls
N     = 1e4;
T     = 1e2;
seed  = 123;
c     = 0;

% level series
y = cumsum(c + randn(T,N));
x = cumsum(c + randn(T,N));

% differenced series
dy = y(2:T,:)-y(1:T-1,:);
dx = x(2:T,:)-x(1:T-1,:);

% storage allocation
bhat	= zeros(N,1);
tstat = zeros(N,1);

tic
% main loop
for n = 1:N
	[b_hat, olsout] = fastols(y(:,n), x(:,n),1);
	bhat(n,:)     = b_hat;
	tstat(n,:)		=	olsout.tstat;
	% based on differenced data
	[b_hatd, olsoutd] = fastols(dy(:,n), dx(:,n),1);
	bhatd(n,:)		= b_hatd;
	tstatd(n,:)		=	olsoutd.tstat;
end
toc;

% level series
% p-value for 95% CI (it is not 5%)
fprintf('Pr(|t-stat|>1.96) (should be 0.05): % 2.4f \n', mean(abs(tstat)>1.96))

% differenced series
% p-value for 95% CI (it is not 5%)
fprintf('Pr(|t-stat|>1.96) (should be 0.05): % 2.4f \n', mean(abs(tstatd)>1.96))

%% PLOTS
set(groot,'defaultLineLineWidth',1.75); % sets the default linewidth;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
% some controls
p.dm = @(x,y)([x y .4 .2]); % FIG DIMENSION
p.fs = 22;			            % FONT SIZE                                       
p.aj = -1.22;	              % SUBTITLE POSITION ADJUSTMENT                  
p.df = 21;			            % DATE FREQUENCY                                  
xg = linspace(-4,4,1e3)';   % xgrid

% plots level series
figure(1);clf;
subplot(1,2,1);
hold on;
  h1 = histogram(bhat,133,'Normalization','pdf');
  plot(xg,normpdf(xg,mean(bhat),std(bhat)))
hold off; box on; grid on;
setplot(p.dm(.045,.6),p.fs,[],6/5);
setyticklabels([0:0.1:0.6], 1)
setxticklabels([-5:1:5])
set(gca,'GridLineStyle',':','GridAlpha',1/3)
setoutsideTicks; add2yaxislabel; tickshrink(.9)
subtitle('(a) Sampling distribution of $\hat{\beta}$', p.aj, p.fs)

subplot(1,2,2);
hold on;
  h2 = histogram(tstat,133,'Normalization','pdf');
  plot(xg,normpdf(xg,0,1))
hold off; box on; grid on;
setplot(p.dm(.555,.6),p.fs,[],6/5);
setyticklabels([0:0.05:0.45],2)
setxticklabels([-40:10:40])
set(gca,'GridLineStyle',':','GridAlpha',1/3)
setoutsideTicks; add2yaxislabel; tickshrink(.9)
subtitle('(b) Sampling distribution of $t-$statistic', p.aj, p.fs, 1)
% uncomment to pritn to pdf
% print2pdf('spurious_bhat_tstat','../graphics/');

% plots differenced series
figure(2);clf
subplot(1,2,1);
hold on;
  h1 = histogram(bhatd,133,'Normalization','pdf');
  plot(xg,normpdf(xg,mean(bhatd),std(bhatd)))
hold off; box on; grid on;
setplot(p.dm(.045,.6),p.fs,[],6/5);
setyticklabels([0:1:5], 0)
setxticklabels([-5:1:5])
set(gca,'GridLineStyle',':','GridAlpha',1/3)
setoutsideTicks; add2yaxislabel; tickshrink(.9)
subtitle('(a) Sampling distribution of $\hat{\beta}$', p.aj, p.fs)

subplot(1,2,2);
hold on;
  h2 = histogram(tstatd,133,'Normalization','pdf');
  plot(xg,normpdf(xg,0,1))
hold off; box on; grid on;
setplot(p.dm(.555,.6),p.fs,[],6/5);
setyticklabels([0:0.05:0.45],2)
setxticklabels([-40:10:40])
set(gca,'GridLineStyle',':','GridAlpha',1/3)
setoutsideTicks; add2yaxislabel; tickshrink(.9)
subtitle('(b) Sampling distribution of $t-$statistic', p.aj, p.fs, 1)
% uncomment to pritn to pdf
% print2pdf('spurious_bhat_tstat_diff','../graphics/');
 









%EOF