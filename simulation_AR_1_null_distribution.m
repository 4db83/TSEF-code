clc;clear all;clf
seed(1234);
% check matlabpool if open
if ~matlabpool('size') > 0; matlabpool; end
tic;
Nsim= 5e4;
TT	=	[50 100 200 500 1000];
A		= [.98 .95 .90 .5 .3];

% choose the parsistence of the series
a		= A(5);
% choose the sample size
T		= TT(5);
% constant in the model
c		= 0.3;
% Nsim by T matrix of WN(0,1);
U		= randn(T,Nsim);
% simulate the AR(1) process
Y		= filter(1,[1 -a],c + U, c/(1-a));

%% create data vectors
yt		= trimr(Y(1:T,:),1,0);		% yt
yt_1	= trimr(Y(1:T,:),0,1);		% y_{t-1}
I			= ones(T-1,1);						% for intercept/constant

% space allocation
bhat	= zeros(Nsim,2);

% main loop
parfor i = 1:Nsim
	ar1	= fastols(yt(:,i),[I yt_1(:,i)]);
	bhat(i,:)	= ar1.bhat';
	se(i,:)		= ar1.se';
	tstat(i,:) = (ar1.bhat'-[c a])./ar1.se';
end
% point estimate & t-statistics
c_hat = bhat(:,1); c_tstat = tstat(:,1);
a_hat = bhat(:,2); a_tstat = tstat(:,2);

%% plots
%plot(bhat)
subplot(2,2,1);histogram(c_hat);setytick(1);
subplot(2,2,2);histogram(a_hat);setytick(0);

subplot(2,2,3);hold on; histogram(c_tstat); 
x = linspace(min(c_tstat),max(c_tstat),1000);
plot(x,normpdf(x,0,1),'r','Linewidth',2); 
ylim([0 0.45]);setytick;
hold off;

subplot(2,2,4);hold on; histogram(a_tstat); 
x = linspace(min(a_tstat),max(a_tstat),1000);
plot(x,normpdf(x,0,1),'r','Linewidth',2); 
ylim([0 0.45]);setytick;
hold off;
toc;

% 95%CIs are
fprintf('-----------------------------------------------------\n');
fprintf(' Sample size is: T = %i, true AR parameter is: %1.1f \n', [T a])
fprintf('-----------------------------------------------------\n');
fprintf(' 95%% CI for c-tstat:		[% 2.4f % 2.4f] \n', percentile(c_tstat,[.025 .975]))
fprintf(' 95%% CI for a-tstat:		[% 2.4f % 2.4f] \n', percentile(a_tstat,[.025 .975]))
fprintf('-----------------------------------------------------\n');
fprintf(' Mean(c_hat): % 2.4f;		Bias: % 2.4f \n', [mean(c_hat) c-mean(c_hat)])
fprintf(' Mean(a_hat): % 2.4f;		Bias: % 2.4f \n', [mean(a_hat) a-mean(a_hat)])

