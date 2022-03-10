% example of ARMA model fit to US Real GDP;
clc;clear all;clf;
addpath('d:/matlab.tools/db.toolbox')
g_ = graph_h;
% check matlabpool if open
%if ~matlabpool('size') > 0; matlabpool; end
% load the data with path to data file
% disp('loading the data ...')
% [data0, dates0] = xlsread('./data/US real GDP.xlsx');
% disp('done')
% use addpath(genpath ...) to add subdirectories as well.
%addpath(genpath('d:/matlab.tools/mfe.toolbox'))

if exist('./data/us_real_gdp.mat') 
	load './data/us_real_gdp.mat';
else
	[data0, dates0] = xlsread('./data/US real GDP.xlsx');
	dates = datenum(dates0(8:end,1),'dd.mm.yyyy');				% and for windows xp
	%dates = datenum(dates0(8:end,1));% ,'dd.mm.yyyy');	% use this for windows 7 above
	Y = data0;
end;
disp('done')

% generate logged series and first differences
y			= log(Y);
dy		= 100*delta(y,1);			% delta function keeps the first entry as NAN!
T			= size(dy,1);

%% plots
% Plot the level series of (log) real GDP
plot(dates,y,g_.lw,1)
setplot([.8 .35],10)
xlim([dates(1) dates(end)]);
set(gca,'XTick',dates([1:24:end]));
setytick(1)
ylim([7.4 10])
hline(0,'k')
datetick('x','yyyy:QQ','keepticks','keeplimits');
%print2pdf('../lectures/graphics/USGDP_level');

% Plot the (log) growth rates of real GDP
plot(dates,dy,g_.lw,1)
setplot([.8 .35],10)
xlim([dates(1) dates(end)]);
set(gca,'XTick',dates([1:24:end]));
setytick(1)
ylim([-3.6 4.51])
setyticklabels([-3.5:1:4.5])
hline(0,'k')
datetick('x','yyyy:QQ','keepticks','keeplimits');
%print2pdf('../lectures/graphics/USGDP_growth');

% drop the first observation because it is NAN
dy		= dy(2:end);				% drop the first entry with tha NAN.

% Plot the ACF/PACF for the growh series.
plotacf(dy);
%print2pdf('../lectures/graphics/acfplot');

%% Estimate the ARMA models
% set upper bounds for P and Q to search over the ARMA model:CHOOSE THESE CAREFULLY.
P = 2;
Q = 2;

% this is a simulated WN series of size 5.
%seed(123);dy = randn(5e3,1);

% space allocation for SBIC and AIC values
BIC_pq	= zeros(P+1,Q+1);
AIC_pq  = zeros(P+1,Q+1);
HQC_pq		= zeros(P+1,Q+1);
pq			= {};						% cell array to display the ARMA orders if needed.

for q = 1:(Q+1);
	parfor (p = 1:(P+1))
		% create the PP and QQ vector entries to be used in estiamte_armax function
		PP = 1:(p-1);
		QQ = 1:(q-1);
		% this stores the ARMA(p,q) orders in a cell array.
		pq{p,q}		= [num2str(max([0 PP])) ',' num2str(max([0 QQ]))];
		% estimate the different ARMA models and store the ICs.
		[~,diagn]	= estimate_armax(dy,1,PP,QQ,[],[],[],[],[],1);
		%[~,diagn]	= estimate_armax(dy,1,PP,QQ);
		BIC_pq(p,q) = diagn.SBIC;
		AIC_pq(p,q) = diagn.AIC;
		HQC_pq(p,q)	= diagn.HQC;
	end;
end;
	
% Store the BIC and AIC matrices
AIC = [[nan (0:Q)];[(0:P)' AIC_pq]];
BIC = [[nan (0:Q)];[(0:P)' BIC_pq]];
HQC = [[nan (0:Q)];[(0:P)' HQC_pq]];

% Display the best ARMA(p,q) orders for AIC, BIC and HQC
% fprintf('---------------------------------------------------------------------------\n');
[p_aic q_aic]=find(min(min(AIC_pq))==AIC_pq);
fprintf('AIC best fitting ARMA model is: ARMA(%d,%d)  \n', [p_aic q_aic]-1)
[p_bic q_bic]=find(min(min(BIC_pq))==BIC_pq);
fprintf('BIC best fitting ARMA model is: ARMA(%d,%d)  \n', [p_bic q_bic]-1)
[p_hqc q_hqc]=find(min(min(HQC_pq))==HQC_pq);
fprintf('HQC best fitting ARMA model is: ARMA(%d,%d)  \n', [p_hqc q_hqc]-1)

%%
% these are the 'best' p,q values for the ARMA(p,q). I use here BIC
p0 = 3;
q0 = 3;

% Estimate the final 'best' model
[armaout,diagout,addon]= estimate_armax(dy,1,1:p0,1:q0);
%print_arma_results(armaout,diagout,addon);

% get the alpha(L) and beta(L) polynomials and compute the roots of the process.
aL = [1 -armaout.pars(2:p0+1)'];
bL = [1 armaout.pars(p0+2:end)'];
disp('Inverse Roots of alpha(L) polynomial')
disp(roots(aL)')
if q0>0;
	fprintf('-------------------------------------------------------------------------------\n');
	disp('Inverse Roots of beta(L) polynomial')
	disp(roots(bL)')
	fprintf('===============================================================================\n');
else
	fprintf('-------------------------------------------------------------------------------\n');
end;

%% plot the fitted and actual series
uhat	= armaout.resids;
dyhat = armaout.yhat;

plot(dates,[nan;dy],g_.lw,1);
hold on;
plot(dates,[nan;dyhat],'r-',g_.lw,1);
setplot([.8 .35],10)
xlim([dates(1) dates(end)]);
set(gca,'XTick',dates([1:24:end]));
setytick(1)
%ylim([7.4 10])
hline(0,'k')
datetick('x','yyyy:QQ','keepticks','keeplimits');
legend('Actual','Fitted')
%print2pdf('../lectures/graphics/fitted_gdp');

%
% plot the residuals
plot(dates,[nan;uhat],g_.lw,1)
setplot([.8 .35],10)
xlim([dates(1) dates(end)]);
set(gca,'XTick',dates([1:24:end]));
setytick(1)
%ylim([7.4 10])
hline(0,'k')
datetick('x','yyyy:QQ','keepticks','keeplimits');
%print2pdf('../lectures/graphics/uhat_gdp');

%%
% plot the ACF/PACF of the residuals
plotacf(uhat);
%print2pdf('../lectures/graphics/uhat_acfplot');

%% plot theoretical ACF/PACF values of fitted model to visually compare to sample ACF/PACF
plotacf0(aL,bL,[-.2 -.2]);
%print2pdf('../lectures/graphics/acf0_armafitted');

%%  Other matlab implementations --------------------------------------------------
% %using Pattons wrapper function
[theta,sig2,vcv,resids,yhat] = estimate_armax2(dy,p0,q0);
%

%% using Econometrics Toolsbox
options = optimset('fmincon');
opts		= optimset(	options,'MaxFunEvals',3000,'Display','off',...
										'Diagnostics','off','HessUpdate','bfgs','Algorithm','sqp',...
										'TolX',1e-14,'TolFun', 1e-14);

arma_model = arima(p0,0,q0);
[MTout,vcv,LogL,arinfot] = estimate(arma_model,dy,'options',opts,'print',false);
print(MTout,vcv)
uhat2 = infer(MTout,dy);
plot([uhat uhat2])



