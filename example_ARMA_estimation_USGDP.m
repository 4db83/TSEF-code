% Script: example_ARMA_estimation_USGDP.m
% Example of ARMA model fitted to US Real GDP;
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear;clc;clf;
% addpath(genpath('PATH-TO-FOLDER/db.toolbox'))
get_new_data = 0;

if get_new_data
  % get data from FRED2
  % Billions of Chained 2012 Dollars, Seasonally Adjusted Annual Rate
  usdata = as_timetable(getFredData('GDPC1', '1947-01-01', '2021-12-31','lin','q'),'gdpc1'); 
  usdata = synchronize(usdata,as_timetable(getFredData('USRECQ', '1947-01-01', '2021-12-31','lin','q'),'NBER'));
  save('./data/real_gdp_US_2021Q4.mat', 'usdata');
  disp('done saving the data')
else
  load './data/real_gdp_US_2021Q4.mat';
end

% generate y = log-gdp and dy = annualized gpd growth.
usdata.y  = log(usdata.gdpc1);
usdata.dy = 400*(usdata.y - lag(usdata.y,1));
% uncomment to write to csv file 
% writetimetable(usdata,'real_gdp_US_2021Q4.csv','Delimiter',',')
ss = timerange('Q1-1947', 'Q4-2019', 'closed');
usdata    = usdata(ss,:);
% head2tail(usdata);
T2 = 151;
usdata(T2,:)

% recession bar color
rec_CLR = .8*ones(3,1); 
fig.ds = 23;
fig.fs = 20;
fig.st = -1.21;

% PLOTTING
set(groot,'defaultLineLineWidth',2); % sets the default linewidth;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
clf;
subplot(2,1,1)
hold on;
  bar( usdata.NBER*10.5, 1, 'FaceColor', rec_CLR); 
  plot(usdata.y,'Color',clr(1));
hold off;
box on; grid on;
setplot([.05 .60 .9 .20],[],[],6/5);
setyticklabels(7.6:.2:10, 1)
setdateticks(usdata.Time, fig.ds, 'yyyy:QQ', fig.fs);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(7.6);hline(10);
setoutsideTicks
add2yaxislabel
tickshrink(.9)
subtitle('(a) Log of US real GDP (level series)', fig.st)
% legendflex(LG,{'T-Bill','FFR'}, 'fontsize', Fns, 'anchor',[1 1],'Interpreter','Latex')

subplot(2,1,2)
hold on;
  bar( usdata.NBER*16, 1, 'FaceColor', rec_CLR); 
  bar(-usdata.NBER*12, 1, 'FaceColor', rec_CLR); 
  plot(usdata.dy,'Color',clr(1));
hold off; vline(T2,'r:'); 
box on; grid on;
setplot([.05 .34 .9 .20],[],[],6/5);
setyticklabels(-12:4:16, 0)
setdateticks(usdata.Time, fig.ds, 'yyyy:QQ', fig.fs);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(-12);hline(16);hline(0);
ylim([-12 16])
setoutsideTicks
add2yaxislabel
tickshrink(.9)
subtitle('(b) First difference of log of US real GDP (annualized growth rate)', fig.st)


%% Estimate the ARMA models
dy = usdata{T2+1:end,'dy'};
II = usdata{T2+1:end,'NBER'};
fprintf('Sample size is %d\n', size(dy,1) )
clf;
plotacf(dy,25)

% set upper bounds for P and Q to search over the ARMA model:CHOOSE THESE CAREFULLY.
P = 3;
Q = 2;

% this is a simulated WN series of size 5.
%seed(123);dy = randn(5e3,1);

% space allocation for SBIC and AIC values
BIC_pq	= zeros(P+1,Q+1);
AIC_pq  = zeros(P+1,Q+1);
HQC_pq	= zeros(P+1,Q+1);
pq			= {};						% cell array to display the ARMA orders if needed.

for q = 1:(Q+1)
	for p = 1:(P+1)
		% create the PP and QQ vector entries to be used in estiamte_armax function
		PP = 1:(p-1);
		QQ = 1:(q-1);
		% this stores the ARMA(p,q) orders in a cell array. Not really needed
    %	pq{p,q}		= [num2str(max([0 PP])) ',' num2str(max([0 QQ]))];
		% estimate the different ARMA models and store the ICs.
%     tmp_	= estimate_armax(dy,1,PP,QQ,II,[],[],[],[],1);
		tmp_	= estimate_armax(dy,1,PP,QQ,[],[],[],[],[],1);
		%[~,diagn]	= estimate_armax(dy,1,PP,QQ);
		BIC_pq(p,q) = tmp_.diagnostics.SBIC;
		AIC_pq(p,q) = tmp_.diagnostics.AIC;
		HQC_pq(p,q)	= tmp_.diagnostics.HQC;
  end
end
	
% Store the BIC and AIC matrices
AIC = [[nan (0:Q)];[(0:P)' AIC_pq]];
BIC = [[nan (0:Q)];[(0:P)' BIC_pq]];
HQC = [[nan (0:Q)];[(0:P)' HQC_pq]];
[AIC BIC HQC]

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
q0 = 1;

% Estimate the final 'best' model
[armaout]= estimate_armax(dy,1,1:p0,1:q0);
print_arma_results(armaout);

% get the alpha(L) and beta(L) polynomials and compute the roots of the process.
aL = [1 -armaout.pars(2:p0+1)'];
bL = [1 armaout.pars(p0+2:end)'];
disp('Inverse Roots of alpha(L) polynomial')
disp(roots(aL)')
if q0>0
	fprintf('-------------------------------------------------------------------------------\n');
	disp('Inverse Roots of beta(L) polynomial')
	disp(roots(bL)')
	fprintf('===============================================================================\n');
else
	fprintf('-------------------------------------------------------------------------------\n');
end

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



