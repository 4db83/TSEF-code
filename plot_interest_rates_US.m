% plot US interest rates
% TB3MS
% 3-Month Treasury Bill Secondary Market Rate (TB3MS)
% Federal Funds Effective Rate (DFF)
% ---------------------------------------------------------------------------------------------
clear; clc;
% % % get data from FREDwebsite if not in ./data/ directory
% usdat = as_timetable(                  getFredData('TB3MS'  , '1934-01-01', '2022-03-01','lin','m'),'Tbill');
% usdat = synchronize(usdat,as_timetable(getFredData('DFF'    , '1947-01-01', '2022-03-01','lin','m'),'ffrate'));
% usdat = synchronize(usdat,as_timetable(getFredData('USRECM' , '1934-01-01', '2022-03-01','lin','m'),'NBER'));
% save('./data/interest_rates_US_2022M1.mat', 'usdat');

% load consumption_income_US data for us
load('./data/interest_rates_US_2022M1.mat')
clc;clf;
% usdat.c = log(usdat.pcecc96);
% usdat.y = log(usdat.dpic96);
tt = timerange('01-Jan-1934', '01-Feb-2022', 'closed');
% tt = timerange('Q1-1947', 'Q4-2021', 'closed');
ss_usdat = usdat(tt,:);
head2tail(ss_usdat );

%% PLOTTING
set(groot,'defaultAxesXTickLabelRotationMode','manual')
set(groot,'defaultLineLineWidth',2); % sets the default linewidth to 1.5;
% recession indicators
NBER_I = ss_usdat.NBER;
% recession bar color
rec_CLR = .83*ones(3,1); 
clf;clc
Fns = 30;
LG = [];
whp = [.5 .2825 .20];	% width/hight of plot
% whp = [.5 .48 .38];	% width/hight of plot
sbl = -1.30;
Dsp = 85;
% --------------------------------------------------------------------------------------------------

subplot(1,3,1)
hold on;
%   plot(ss_usdat.RetailSales*1e-6);
  bar( NBER_I*20, 1, 'FaceColor', rec_CLR); 
%   bar( NBER_I*11.65,1,'FaceColor',rec_CLR,'EdgeColor',rec_CLR,'ShowBaseLine','on' );
LG(1) = plot(ss_usdat.Tbill ,'Color',clr(1));
LG(2) = plot(ss_usdat.ffrate,'Color',clr(2),'LineStyle','--');
% LG(1) = plot(ss_usdat.RetailSales*1e-6,'Color',clr(2));
hold off;
box on; grid on;
setplot([.05 whp],[],0);
% setyticklabels(11.6:.2:13.4,1,Fns)
setdateticks(ss_usdat.Time,Dsp,'yy', Fns);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(0);hline(20);
tickshrink(1.9)
setoutsideTicks
subtitle('(a) Levels', sbl)
legendflex(LG,{'T-Bill','FFR'}, 'fontsize', Fns, 'anchor',[1 1],'Interpreter','Latex')
% add2yaxislabel

subplot(1,3,2)
hold on; LG = [];
%   plot(ss_usdat.RetailSales*1e-6);
  bar( NBER_I*3, 1, 'FaceColor', rec_CLR); 
  bar(-NBER_I*5, 1, 'FaceColor', rec_CLR); 
%   bar( NBER_I*11.65,1,'FaceColor',rec_CLR,'EdgeColor',rec_CLR,'ShowBaseLine','on' );
LG(1) = plot(delta(ss_usdat.Tbill),'Color',clr(1));
% LG(1) = plot(ss_usdat.RetailSales*1e-6,'Color',clr(2));
hold off;
box on; grid on;
setplot([.385 whp],[],0);
% setyticklabels(11.5:.5:13.5,1,Fns)
setdateticks(ss_usdat.Time,Dsp,'yy', Fns);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(-5);hline(3);
hline(0,'k');
tickshrink(1.9)
setoutsideTicks
subtitle('(b) Differences', sbl)
% add2yaxislabel
% ADD LEGEND
legendflex(LG,{'T-Bill'}, 'fontsize', Fns, 'anchor',[1 1],'Interpreter','Latex')

subplot(1,3,3)
dhndl = density(delta(ss_usdat.Tbill),[0.10 53]);
subtitle('(c) Density of differences', sbl)
tickshrink(.66)
box on; ax = gca;
ax.LineWidth = 1.25;
grid on;
set(gca,'GridLineStyle',':','GridAlpha',1/3);
% setoutsideTicks
% setplot([.65 whp],Fns,1);
setyticklabels(0:1:4,0,Fns)
setxticklabels(-4:1:2)
xlim([min(delta(ss_usdat.Tbill)) max(delta(ss_usdat.Tbill))])
set(gca,'Position',[.71 whp],'FontSize',Fns);
% set(gca,'Position',[.21 whp],'FontSize',Fns);

% PRINT TO PDF
% print2pdf('interest_rates_US_2022M1','../graphics/');






































%  EOF 

