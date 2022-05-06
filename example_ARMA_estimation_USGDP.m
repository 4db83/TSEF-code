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
  % print/export to xlsx if needed
  print2xls(usdata,'/data/real_gdp_US_2021Q4.xlsx')
  disp('done saving the data')
else
  load './data/real_gdp_US_2021Q4.mat';
end

%% generate y = log-gdp and dy = annualized gpd growth.
usdata.y  = log(usdata.gdpc1);
usdata.dy = 400*(usdata.y - lag(usdata.y,1));
% % uncomment to write to csv file 
% writetimetable(usdata,'real_gdp_US_2021Q4.csv','Delimiter',',')
ss = timerange('Q1-1947', 'Q4-2019', 'closed');
usdata    = usdata(ss,:);
% head2tail(usdata);
T2 = 151; % break in volatility

% recession bar color
rec_CLR = .8*ones(3,1); 
fig.ds = 23;
fig.fs = 20;
fig.st = -1.21;
fig.dim = [.85 .2];
fig.pos = @(x) ([.07 x]);

% PLOTTING
set(groot,'defaultLineLineWidth',2); % sets the default linewidth;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
clf;
subplot(2,1,1)
hold on;
  bar( usdata.NBER*10, 1, 'FaceColor', rec_CLR); 
  plot(usdata.y,'Color',clr(1));
hold off;
box on; grid on;
setplot([fig.pos(.60) fig.dim],[],[],6/5);
setyticklabels(7.5:.5:10, 1)
setdateticks(usdata.Time, fig.ds, 'yyyy:QQ', fig.fs);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(7.5);hline(10);
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
setplot([fig.pos(.34) fig.dim],[],[],6/5);
setyticklabels(-12:4:16, 0)
setdateticks(usdata.Time, fig.ds, 'yyyy:QQ', fig.fs);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(-12);hline(16);hline(0);
ylim([-12 16])
setoutsideTicks
add2yaxislabel
tickshrink(.9)
subtitle('(b) First difference of log of US real GDP (annualized growth rate)', fig.st)

% uncomment to print to pdf
% print2pdf('USGDP_level_growth','../graphics')

% PLOT SAMPLE ACF/PACF OF ANNUALIZED GDP 
plotacf(usdata.dy);
% uncomment to print to pdf
% print2pdf('acf_USGDP_growth','../graphics')

% Estimate the ARMA models
% set upper bounds for p* and q* to search over the ARMA model: CHOOSE THESE CAREFULLY.
P = 2;
Q = 1;

% SPACE ALLOCATION FOR BIC AND AIC VALUES
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
		tmp_ = estimate_armax(usdata.dy,1,PP,QQ,[],[],[],[],[],1);
		BIC_pq(p,q) = tmp_.diagnostics.SBIC;
		AIC_pq(p,q) = tmp_.diagnostics.AIC;
		HQC_pq(p,q)	= tmp_.diagnostics.HQC;
  end
end
	
% STORE THE BIC AND AIC MATRICES
AIC = [[nan (0:Q)];[(0:P)' AIC_pq]];
BIC = [[nan (0:Q)];[(0:P)' BIC_pq]];
HQC = [[nan (0:Q)];[(0:P)' HQC_pq]];
ICs = [AIC BIC HQC];

% Display the best ARMA(p,q) orders for AIC, BIC and HQC
% fprintf('---------------------------------------------------------------------------\n');
[p_aic q_aic]=find(min(min(AIC_pq))==AIC_pq);
fprintf('AIC best fitting ARMA model is: ARMA(%d,%d)  \n', [p_aic q_aic]-1)
[p_bic q_bic]=find(min(min(BIC_pq))==BIC_pq);
fprintf('BIC best fitting ARMA model is: ARMA(%d,%d)  \n', [p_bic q_bic]-1)
[p_hqc q_hqc]=find(min(min(HQC_pq))==HQC_pq);
fprintf('HQC best fitting ARMA model is: ARMA(%d,%d)  \n', [p_hqc q_hqc]-1)

%% Estimate the final 'best' models based on IC
arma_aic = estimate_armax(usdata.dy,1,1:(p_aic-1),1:(q_aic-1)); print_arma_results(arma_aic);
arma_bic = estimate_armax(usdata.dy,1,1:(p_bic-1),1:(q_bic-1)); print_arma_results(arma_bic);
arma_hqc = estimate_armax(usdata.dy,1,1:(p_hqc-1),1:(q_hqc-1)); print_arma_results(arma_hqc);

% DO THE PLOTTING NOW
clf;
subplot(2,1,1)
hold on; LG = [];
  bar( usdata.NBER*16, 1, 'FaceColor', rec_CLR); 
  bar(-usdata.NBER*12, 1, 'FaceColor', rec_CLR); 
LG(1) = plot(usdata.dy,'Color',clr(1),'LineWidth',2.75);
LG(2) = plot(addnans(arma_aic.yhat,1),'Color',clr(2),'LineWidth',2.5);
LG(3) = plot(addnans(arma_bic.yhat,1),'Color',clr(3),'LineStyle','--');
% LG(4) = plot(addnans(arma_hqc.yhat,1),'Color',clr(5),'LineStyle',':');
hold off; vline(T2,'r:'); 
box on; grid on;
setplot([fig.pos(.6) fig.dim],[],[],6/5);
setyticklabels(-12:4:16, 0)
setdateticks(usdata.Time, fig.ds, 'yyyy:QQ', fig.fs);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(-12);hline(16);hline(0);
ylim([-12 16])
setoutsideTicks
add2yaxislabel
tickshrink(.9)
subtitle('(a) Actual and fitted values', fig.st)
legendflex(LG,{'GDP growth','AIC-ARMA(2,1)','BIC-ARMA(1,0)'}, 'fontsize', fig.fs - 1, 'anchor',3.*[1 1],'Interpreter','Latex')

subplot(2,1,2)
hold on; LG = [];
  bar( usdata.NBER*16, 1, 'FaceColor', rec_CLR); 
  bar(-usdata.NBER*12, 1, 'FaceColor', rec_CLR); 
% LG(1) = plot(usdata.dy,'Color',clr(1));
LG(1) = plot(addnans(arma_aic.uhat,1),'Color',clr(5),'LineWidth',2.75);
LG(2) = plot(addnans(arma_bic.uhat,1),'Color',clr(2),'LineStyle','--');
% LG(4) = plot(addnans(arma_hqc.yhat,1),'Color',clr(5),'LineStyle',':');
hold off; vline(T2,'r:'); 
box on; grid on;
setplot([fig.pos(.34) fig.dim],[],[],6/5);
setyticklabels(-12:4:16, 0)
setdateticks(usdata.Time, fig.ds, 'yyyy:QQ', fig.fs);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(-12);hline(16);hline(0);
ylim([-12 16])
setoutsideTicks
add2yaxislabel
tickshrink(.9)
subtitle('(b) Residuals', fig.st)
legendflex(LG,{'AIC-ARMA(2,1)','BIC-ARMA(1,0)'}, 'fontsize', fig.fs - 1, 'anchor',3.*[1 1],'Interpreter','Latex')
% uncomment to print to pdf 
% print2pdf('fitted_values_US','../graphics');

% [usdata.dy addnans(arma_aic.yhat,1) usdata.dy-addnans(arma_aic.yhat,1)]

%% PLOT SAMPLE  ACF/PACF OF THE RESIDUAL SERIES OF THE AR(1)
plotacf(arma_bic.uhat);
% uncomment to print to pdf
print2pdf('acf_arma_bic_fit','../graphics')

plotacf(arma_aic.uhat);
% uncomment to print to pdf
print2pdf('acf_arma_aic_fit','../graphics')

% plot theoretical ACF/PACF values of fitted model to visually compare to sample ACF/PACF
clf;clc;
aL_bic = [1 -arma_bic.pars(2:(p_bic-1)+1)'];
bL_bic = [1  arma_bic.pars((p_bic-1)+2:end)'];

% plotacf0(aL_bic,bL_bic);
% uncomment to print to pdf
% print2pdf('acf0_ar1', '../graphics');



















































%EOF 
