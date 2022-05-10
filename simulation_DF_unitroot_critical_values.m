% Script: simulation_DF_unitroot_critical_values.m
% Simulation of Dickey-Fuller critical values
% uncomment print2pdf generate pdf from plot using the print2pdf function.
clear; clc;
% addpath(genpath('PATH-TO-FOLDER/db.toolbox'))
% some controls
N     = 1e4;
T     = 251;
C     = ones(T-1,1);
trnd  = (1:T-1)';
seed  = 1234;

% space allocation for storage of tstats and coeffs
tstat0  = zeros(N,1);   rho0 = zeros(N,1);
tstatC  = zeros(N,1);   rhoC = zeros(N,1);
tstatCT = zeros(N,1);   rhoCT= zeros(N,1);

tic;
% main loop
for jj = 1:N
  % generate pure random walk
  y = cumsum(randn(T,1));
  % make Y X variables
  Y = y(2:T); 			% y(t)
  X = y(1:T-1);		  % y(t-1)
  
	% run the 3 separate regressions. ,1 is for no constant
  [bhat,out]      = fastols(Y, X,1);  
  [Cbhat,Cout]    = fastols(Y,[X C],1);
  [CTbhat,CTout]  = fastols(Y,[X C trnd],1);

  % store t-stats and bhat coeffcients
  tstat0(jj)  = (bhat(1)-1)/out.se(1);
  tstatC(jj)  = (Cbhat(1)-1)/Cout.se(1);
  tstatCT(jj) = (CTbhat(1)-1)/CTout.se(1);
  rho0(jj) 		= bhat(1);
  rhoC(jj) 		= Cbhat(1);
  rhoCT(jj)		= CTbhat(1);
end
toc

% DF-critical values
pct_p0 = percentile(tstat0 , [1 2.5 5]) ; % histogram(tstat0, 100, 'Normalization','pdf')
pct_pC = percentile(tstatC , [1 2.5 5]) ; % histogram(tstatC, 100, 'Normalization','pdf')
pct_CT = percentile(tstatCT, [1 2.5 5]); % histogram(tstatCT,100, 'Normalization','pdf')

PCT = [pct_p0; pct_pC; pct_CT];

fprintf('    0.010     0.025     0.050\n')
disp(PCT)

%% PLOT CONTROLS
clf;
set(groot,'defaultLineLineWidth',1.5); % sets the default linewidth;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
fig.dim = [.85 .25];
fig.pos = @(x) ([.07 x]);
% xgrid for Density estimate and plot
xg  = linspace(-6,5,1e3)';              
% compute the density over xg grid
t0  = ksdensity(tstat0, xg); 
tC  = ksdensity(tstatC, xg);
tCT = ksdensity(tstatCT,xg);
 
LW = 'LineWidth';
clf;
hold on;
LG(1) = plot(xg,t0);
LG(2) = plot(xg,tC);
LG(3) = plot(xg,tCT);
LG(4) = plot(xg,normpdf(xg,0,1),'-k');
hold off; 
box on; grid on;
setplot([fig.pos(.60) fig.dim],16,[],6/5);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
setoutsideTicks
add2yaxislabel
tickshrink(.9)

legNames = {'$\tau_0$' ;
       			'$\tau_{\mu}$' ;
       			'$\tau_{\tau}$' ;
       			'$N(0,1)$' }; 
legendflex(LG,legNames,'Interpreter','Latex')







%EOF