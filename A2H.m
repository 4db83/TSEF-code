% Script: simulation_DF_unitroot_critical_values.m
% Simulation of Dickey-Fuller critical values. Cleaned up. does not require fastols anymore, uses
% fstols function defined at the end of this script.
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear; clc; clf;
addpath(genpath('./db.toolbox/'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. Computer Exercises 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set simulation parameters
N     = 1e4;
T     = 170;
C     = ones(T,1);
trnd  = (1:T)';
rho   = 1; 
% set seend of random number generator for reproducibility 
rng(123)

% space allocation for storage of tstats and coeffs
tstat0  = zeros(N,1);   rho0  = zeros(N,1);
tstatC  = zeros(N,1);   rhoC  = zeros(N,1);
tstatCT = zeros(N,1);   rhoCT = zeros(N,1);

% -----------------------------------------------------------------------------------------
% MAIN LOOP TO GENERATE RW SERIES AND ESTIMATE 3 DIFFERENT UNIT ROOT TESTING REGRESSIONS
tic
for jj = 1:N
  % generate pure random walk
  y = cumsum(randn(T+1,1));
  
  % MAKE Y X VARIABLES
  Y = y(2:end); 			% y(t)
  X = y(1:end-1);		  % y(t-1)

	% run the 3 separate regressions
  [bhat,se]      = fastols(Y,X);  
  [Cbhat,Cse]    = fastols(Y,[X C]);
  [CTbhat,CTse]  = fastols(Y,[X C trnd]);
  % store bhat (rho_hat) coeffcients
  rho0(jj) 		= bhat(1);
  rhoC(jj) 		= Cbhat(1);
  rhoCT(jj)		= CTbhat(1);
  % store t-stats now for the null of a unit-root
  tstat0(jj)  = (bhat(1)  -rho)/se(1);
  tstatC(jj)  = (Cbhat(1) -rho)/Cse(1);
  tstatCT(jj) = (CTbhat(1)-rho)/CTse(1);
end 
toc

%% COMPUTE PERCENTILS OF CRITICAL VALUES
pctls   = [1 2.5 5 7.5 10 50 90 92.5 95 97.5 99]/100; 
% DF-critical values
pct_N0  = norminv(pctls);
pct_t0  = quantile(tstat0 , pctls ); 
pct_tC  = quantile(tstatC , pctls ); 
pct_tCT = quantile(tstatCT, pctls ); 
% combine for printing
PCT_t   = [pct_N0; pct_t0; pct_tC; pct_tCT];

% make a pretty output table and print to screen
caseNamaes = {'N(0,1)', 'DF0 (No Constant)', 'DFc (Constant)','DFt (Constant & Trend)'};
sep; disp('                  Percentiles of the simulated distribution and the N(0,1) ')
print2screen(PCT_t,[['T = ',num2str(T),' (rho-1)/se(rho)'], caseNamaes],num2str(pctls'),4)

% compute the of the tstatistic density over xg grid
xg  = linspace(-6,4,1e3)';              % xgrid for tstats
t0  = ksdensity(tstat0, xg); 
tC  = ksdensity(tstatC, xg);
tCT = ksdensity(tstatCT,xg);

% compute the density of T(rho_hat -1).
pxg = linspace(-4,4,1e3)';             % xgrid for T(rho_hat -1)
p0  = ksdensity(sqrt(T)*(rho0 -rho), pxg); 
pC  = ksdensity(sqrt(T)*(rhoC -rho), pxg);
pCT = ksdensity(sqrt(T)*(rhoCT-rho), pxg);
% FOR UNIT ROOT PROCESS
if rho == 1
  pxg = linspace(-35,5,1e3)';           % xgrid for T(rho_hat -1)
  p0  = ksdensity(T*(rho0 -rho), pxg); 
  pC  = ksdensity(T*(rhoC -rho), pxg);
  pCT = ksdensity(T*(rhoCT-rho), pxg);
end

%% -----------------------------------------------------------------------------------------
% PLOTS 
set(groot,'defaultLineLineWidth',1.5);  % sets the default linewidth for all lines in plot                                ;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
fig.dim = [.85 .7];
fig.pos = @(x) ([.07 x]);

% plot DF t-statistic distribution
set(figure(1), 'WindowStyle', 'Docked'); clf;
hold on; LG = [];
LG(1) = plot(xg,t0);
LG(2) = plot(xg,tC);
LG(3) = plot(xg,tCT);
LG(4) = plot(xg,normpdf(xg,0,1),'-k');
hold off; 
box on; grid on;
setplot([fig.pos(.2) fig.dim],1,16,3/4);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
setyticklabels([0:0.1:0.6], 1)
setoutsideTicks
add2yaxislabel
tickshrink(.9)
legNames = {'$\tau_0$ no constant, no trend' ;
       			'$\tau_{\mu}$ with constant' ;
       			'$\tau_{\tau}$ with constant and trend' ;
       			'$N(0,1)$' }; 
addlegend(LG,legNames,1)
addsubtitle('Distribution of the $t-$statistic',-1.15,18)
% print2pdf('DF_tdistr',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A. Computer Exercises 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://github.com/4db83/code-TSEF/raw/main/data/real_gdp_US_2022Q4.mat
% this downloads the data and stores it locally in the current directory
websave('US_data.mat','https://github.com/4db83/code-TSEF/raw/main/data/real_gdp_US_2022Q4.mat');
% LOADS THE DATA YOU HAVE DOWNLOADED, WHICH IS CALLED: USDATA
load('US_data.mat') 
% Name Date column correctly
usdata.Properties.DimensionNames{1} = 'Date';
% or  
% load('./data/real_gdp_US_2022Q4.mat');
% GENERATE Y = LOG-GDP AND DY = ANNUALIZED GPD GROWTH.
usdata.y  = log(usdata.gdpc1);
usdata.dy = 100*(usdata.y - lag(usdata.y,1));
% set the sample size to the one required
ss = timerange('Q1-1947', 'Q4-2019', 'closed');
% ss = timerange('01.01.1947', '31.12.2019', 'closed');
usdata = usdata(ss,:);
Dates  = usdata.Properties.RowTimes;

% ESTIMATE ARMA(1,1) MODEL
arma11 = estimate_armax(usdata.dy,1, 1,1); print_arma_results(arma11);
arma01 = estimate_armax(usdata.dy,1, 0,1); print_arma_results(arma01);
arma10 = estimate_armax(usdata.dy,1, 1,0); print_arma_results(arma10);
% USING MATLAB PROPRIETARY ECONOMETRICS TOOLBOX 
% arma11_ET = estimate(arima(1,0,1),usdata.dy);
% arma01_ET = estimate(arima(0,0,1),usdata.dy);
% arma10_ET = estimate(arima(1,0,0),usdata.dy);

%% BN-decompositions based on the lecture notes examples
% BN Based on an AR(1)
mu_AR1      = arma10.pars(1)/(1 - arma10.pars(2));
trnd_BN_AR1 = usdata.y + arma10.pars(2)/(1 - arma10.pars(2))*(usdata.dy - mu_AR1);
cycl_BN_AR1 = usdata.y - trnd_BN_AR1;

% BN base on an MA(1)
cycl_BN_MA1 = -arma01.pars(2)*[NaN; arma01.uhat];

% BN base on an ARMA(1,1)
mu_ARMA11      = arma11.pars(1)/(1 - arma11.pars(2));
kk = arma11.pars(2)*(usdata.dy - mu_ARMA11) - arma11.pars(3)*[NaN; arma11.uhat];
trnd_BN_ARMA11 = usdata.y + arma11.pars(2)/(1 - arma11.pars(2))*kk;
cycl_BN_ARMA11 = usdata.y - trnd_BN_ARMA11;
% add them to the usdata base
usdata.BN_AR1     = cycl_BN_AR1;
usdata.BN_MA1     = cycl_BN_MA1;
usdata.BN_ARMA11  = cycl_BN_ARMA11;

% % compare to the R estiamtes
% R = (parquetread('R_cycles.parquet')); % R cycles only
% R.Properties.DimensionNames{1} = 'Date'
% M = usdata(:,end-2:end); % matlab cycles only
% all = join(M,R)
% parquetwrite('R_out.parquet',R)

% HP-filter cycle
[cycl,trnd] = hp_filter(usdata.y, 1600);

%plot the HP-filter and BN decomposition cycles/temporary components
set(figure(2), 'WindowStyle', 'Docked'); clf;
hold on; LG = [];
addrecessionBars(usdata.NBER,[-6.5,4.5])
LG(1) = plot(100*cycl,      'Color',clr(1));
LG(2) = plot(cycl_BN_AR1,   'Color',clr(2));
LG(3) = plot(cycl_BN_MA1,   'Color',clr(3));
LG(4) = plot(cycl_BN_ARMA11,'Color',clr(4));
hold off; hline(0)
box on; grid on;
setplot([fig.pos(.2) fig.dim],2,16,3/4);
% set(gca,'GridLineStyle',':','GridAlpha',1/3);
setyticklabels(-6:2:4,0); 
ylim([-6.5 4.5])
setdateticks(Dates, 22, 'yyyy:QQ');	
% setyticklabels([0:0.1:0.6], 1)
setoutsideTicks
add2yaxislabel

legNames = {'HP-Filter Cycle'   ;
            'BN-AR(1) Cycle'    ;
            'BN-MA(1) Cycle'    ; 
            'BN-ARMA(1,1) Cycle';
            };
addlegend(LG,legNames,3)
% print2pdf('HPvsBN',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIAL FAST OLS FUNCTION IN SAME FILE SO THAT YOU CAN SEE WHAT IS COMPUTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast ols
function [beta,se_beta] = fastols(y,X)
  [T, k]  = size(X);
  beta    = X\y;
  uhat		= y-X*beta;
  se_beta = sqrt(diag(inv(X'*X)*(uhat'*uhat)/(T-k)));
end

% HP-filter
function [cycle,trend,m] = hp_filter(data,lambda)
  T			= size(data,1);
  I			= speye(T);													% sparse identiy matrix.
  D			= diff(I,2);												% Difference maker matrix.
  m			= (I + lambda*(D'*D));		          % shows the moving average weight matrix.
  trend	= m\data; 												  % permanent component.
  cycle	= data - trend;							        % transitory component. 
end














































%EOF