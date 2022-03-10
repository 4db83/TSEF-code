clear all; clc;

N     = 1e6;
B     = 50;
T     = 25 ;
C     = ones(T-1,1);
trnd  = (1:T-1)';

tstat     = zeros(N,1); Bhat0 = zeros(N,1);
tstatC    = zeros(N,1); BhatC = zeros(N,1);
tstatCT   = zeros(N,1); BhatCT= zeros(N,1);

tic;
% big loop
parfor (i = 1:N)
%i = 1;
  u   = randn(B+T,1); % generate shocks.
  y   = cumsum(u);    % RW without drift.

  % drop burn in
  y   = y(B+1:B+T);
  
  % make Y X variables
  Y   = y(2:T);
  X   = y(1:T-1);
  
% run the 3 sperate regressions.
  out   = fastols(Y,X);  
  Cout  = fastols(Y,[X C]);
  CTout = fastols(Y,[X C trnd]);

  %BhatCT(1)./resCT.bstd(1)
  tstat(i)   = out.bhat(1)/out.se(1);
  tstatC(i)  = Cout.bhat(1)/Cout.se(1);
  tstatCT(i) = CTout.bhat(1)/CTout.se(1);
  Bhat0(i) = out.bhat(1);
  BhatC(i) = Cout.bhat(1);
  BhatCT(i)= CTout.bhat(1);
end;
toc
%
percentile(tstatC,0.05);hist(tstatC,100)
percentile(tstatCT,0.05);hist(tstatCT,100);


%%
B0	= T*(Bhat0-1);
BC	= T*(BhatC-1);
BCT = T*(BhatCT-1);

pctB0   = percentile(B0,[.01 .025 .05 .1 .5 .9 .95 .975 .99])';
pctBC   = percentile(BC,[.01 .025 .05 .1 .5 .9 .95 .975 .99])';
pctBCT  = percentile(BCT,[.01 .025 .05 .1 .5 .9 .95 .975 .99])';


PCT = [pctB0;pctBC;pctBCT];

fprintf('    0.010     0.025     0.050\n')
disp(PCT)

%%
% gx      = linspace(-70,10,1000)';              % xgrid for Density estimate and plot
% t0      = ksdensity(T*(Bhat0-1),gx);p0 = percentile(t0,5);
% tC      = ksdensity(T*(BhatC-1),gx);
% tCT     = ksdensity(T*(BhatCT-1),gx);
% 
% %%
% LW = 'LineWidth';
% plot(gx,tCT,'-','Color',[0 .7 0],LW,1);
% hold on;
% plot(gx,tC, '-b',LW,1);
% plot(gx,t0, '-r',LW,1);
% %plot(gx,normpdf(gx,0,1),'-k',LW,1);
% hold off;
% h = legend({
%         '$T(\hat\rho_{\tau}-1) $';
%         '$T(\hat\rho_{\mu}-1) $';
%         '$T(\hat\rho_{0}-1) $';
% %        '$\tau_{0_{}}$';
% %        '$N(0,1)$'
%         },'FontSize',14,'Interpreter','Latex','Location','NorthWest');
% hc = get(h,'Children');
% pos = get(hc(9),'Position');
% deltaX = 10.1; % Need to experiment with this to choose the best value for your case
% pos2 = pos - [deltaX 0 0];
% set(hc(9),'Position',pos2)
% xlim([-35 5])
% FN = 11;                        % font size for plots
% setplot([.8 .5], FN);
% setytick(get(gca,'YTick'))
% %print2pdf('..\lectures\graphics\df_dist2');
% %save('Bhat.mat', 'Bhat0', 'BhatC', 'BhatCT');
