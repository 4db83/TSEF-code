% example spurious regression with two independent random walks. 
clear all;clc;
% check matlabpool if open
if ~matlabpool('size') > 0; matlabpool; end
% sample size and number of simulations
%%
T     = 1e2;
Nsim  = 1e5;
seed(1234);
c = 0.0;
tic;
% level series
y = cumsum(c+randn(T,Nsim));
x = cumsum(c+randn(T,Nsim));

% differenced series
dy = y(2:T,:)-y(1:T-1,:);
dx = x(2:T,:)-x(1:T-1,:);

% storage allocation
bhat	= zeros(Nsim,1);
tstat = zeros(Nsim,1);

% main loop
parfor n = 1:Nsim
	olsout			= fastols(y(:,n), x(:,n));
	bhat(n,:)		= olsout.bhat;
	tstat(n,:)	=	olsout.tstat;
	% based on differenced data
	olsoutd			= fastols(dy(:,n), dx(:,n));
	bhatd(n,:)	= olsoutd.bhat;
	tstatd(n,:)	=	olsoutd.tstat;
end;
toc;

%% level series
% p-value for 95% CI (it is not 5%)
fprintf('Pr(|t-stat|>1.96) (should be 0.05): % 2.4f \n', mean(abs(tstat)>1.96))
% plots
clf;
% bhat
dims = [.4 .4];
%subplot(1,2,1);
setplot(dims,11);
h1 = histogram(bhat,300,[0 1/T]);
ylim([0 4.1])
xlim([-5 5])
setytick(1);
setyticklabels([0:1:5])
%print2pdf('../lectures/graphics/spurious_bhat');
%
% tstat
setplot(dims,10.50);
setytick(2);
h2 = histogram(tstat,300,1);
xlim([-60 60])
setytick(2);
%print2pdf('../lectures/graphics/spurious_tstat');

%% differenced data
% p-value for 95% CI (it is not 5%)
fprintf('Pr(|t-stat|>1.96) (should be 0.05): % 2.4f \n', mean(abs(tstatd)>1.96))
% plots
clf;
% bhat
dims = [.4 .4];
%subplot(1,2,1);
setplot(dims,11);
h1 = histogram(bhatd,300,[0 1/T]);
ylim([0 4.1])
xlim([-5 5])
setytick(1);
setyticklabels([0:1:5])
%print2pdf('../lectures/graphics/spurious_bhatd');
%
% tstat
setplot(dims,10.50);
setytick(2);
h2 = histogram(tstatd,300,1);
xlim([-60 60])
setytick(2);
%print2pdf('../lectures/graphics/spurious_tstatd');
