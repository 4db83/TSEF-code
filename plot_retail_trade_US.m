% plot US retail trade, seasonally unadjusted.
% NOTE: To be able to run this code, you need the contents of the db.toolbox available from:
% https://github.com/4db83/db.toolbox/archive/refs/heads/main.zip. Unzip the contents locally
% to the same directory as this script, and then uncomment the following line below:
% addpath(genpath('./db.toolbox-main'))
clear;clc;clf;
% sets the default linewidth to 1.5;
set(groot,'defaultLineLineWidth',2); 

% % get data from FREDwebsite if not in ./data/ directory
% usdat = as_timetable(getFredData('RSXFSN', '1947-01-01', '2021-12-31','lin','m'),'RetailSales');
% % usdat = synchronize(usdat,as_timetable(getFredData('USRECQ', '1947-01-01', '2021-12-31','lin','q'),'NBER'));
% usdat = synchronize(usdat,as_timetable(getFredData('USRECM', '1947-01-01', '2021-12-31','lin','m'),'NBER'));
% % combs = synchronize(us_data_merged,usdat);
% save('./data/retail_sales_US_2022M1.mat', 'usdat');

% load consumption_income_US data for us
load('./data/retail_sales_US_2022M1.mat')
% usdat.c = log(usdat.pcecc96);
% usdat.y = log(usdat.dpic96);
tt = timerange('01-Jan-1992', '01-Jan-2022', 'closed');
% tt = timerange('Q1-1947', 'Q4-2021', 'closed');
ss_usdat = usdat(tt,:);
head2tail(ss_usdat );

%% recession indicators
NBER_I = ss_usdat.NBER;
% recession bar color
rec_CLR = .83*ones(3,1); 
clf;clc
Fns = 20;
LG = [];
whp = [.5 .41 .18];	% width/hight of plot
stp = -1.24;

subplot(1,2,1)
hold on;
%   plot(ss_usdat.RetailSales*1e-6);
  bar( NBER_I*700, 1, 'FaceColor', rec_CLR); 
%   bar( NBER_I*11.65,1,'FaceColor',rec_CLR,'EdgeColor',rec_CLR,'ShowBaseLine','on' );
LG(1) = plot(ss_usdat.RetailSales*1e-3,'Color',clr(2));
% LG(1) = plot(ss_usdat.RetailSales*1e-6,'Color',clr(2));
hold off;
box on; grid on;
% setplot([.07 tp-0*dp whp]);
setplot([.06 whp],[],0);
% setplot([.36 .86 .21],.25)
% setplot([.80 .27]);
% setyticklabels(11.6:.2:13.4,1,Fns)
setdateticks(ss_usdat.Time,31,'yy', Fns);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(0);hline(700);
tickshrink(1.5)
setoutsideTicks
subtitle('(a) Raw series (Billions of Dollars)', stp)
% add2yaxislabel

subplot(1,2,2)
hold on;
%   plot(ss_usdat.RetailSales*1e-6);
  bar( NBER_I*13.5, 1, 'FaceColor', rec_CLR); 
%   bar( NBER_I*11.65,1,'FaceColor',rec_CLR,'EdgeColor',rec_CLR,'ShowBaseLine','on' );
LG(1) = plot(log(ss_usdat.RetailSales),'Color',clr(2));
% LG(1) = plot(ss_usdat.RetailSales*1e-6,'Color',clr(2));
hold off;
box on; grid on;
setplot([.55 whp]);
% setplot([.07 tp-1*dp whp]);
% setplot([.36 .86 .21],.25)
% setplot([.80 .27]);
setyticklabels(11.5:.5:13.5,1,Fns)
setdateticks(ss_usdat.Time,31,'yy', Fns);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
hline(11.5);hline(13.5);
tickshrink(1.5)
setoutsideTicks
subtitle('(b) Log-transformed', stp)
% add2yaxislabel
% ADD LEGEND
% legendflex(LG,{'log(US Retail Sales)'}, 'fontsize', Fns + 1, 'anchor',[1 1],'Interpreter','Latex')
% PRINT TO PDF
print2pdf('retail_sales_US_2022M1','../graphics/');






































%  EOF 

% TB3MS
% 3-Month Treasury Bill Secondary Market Rate (TB3MS)
% Federal Funds Effective Rate (DFF)