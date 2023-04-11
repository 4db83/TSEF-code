% plot consumption and disposable income data for the US to replicate Hamilton figure on page 601
clear; clc; clf;
addpath(genpath('./db.toolbox/'))

% get data from FREDwebsite if not in ./data/ directory
usdat = as_timetable(getFredData('pcecc96', '1947-01-01', '2021-12-31','lin','q'),'pcecc96');
% usdat = synchronize(usdat,as_timetable(getFredData('USRECQ', '1947-01-01', '2021-12-31','lin','q'),'NBER'));
usdat = synchronize(usdat,as_timetable(getFredData('dpic96', '1947-01-01', '2021-12-31','lin','q'),'dpic96'));
% combs = synchronize(us_data_merged,usdat);
% save('./data/consumption_income_US_2021Q4.mat', 'usdat');

% load consumption_income_US data for us
load('./data/consumption_income_US_2021Q4.mat')
clc;clf;
usdat.c = log(usdat.pcecc96);
usdat.y = log(usdat.dpic96);
% head2tail(usdat);
ss = timerange('Q1-1947', 'Q4-2021', 'closed');
ss = timerange('Q1-1947', 'Q4-1989', 'closed');
ss = timerange('Q1-1947', 'Q4-2019', 'closed');
% extract variables of interst
y = usdat(ss,'y');
c = usdat(ss,'c');

clf;
set(groot,'defaultLineLineWidth',2);
% font size
Fns = 20;
hold on;
  plot(c.Variables);
  plot(y.Variables);
hold off;
setplot([.57 .26]);
setdateticks(usdat.Time,31,'yyyy', Fns);	
setoutsideTicks; add2yaxislabel
% add legend
addlegend({'Consumption $(c_t)$','Disposible Income $(y_t)$'}, 1, Fns + 1, 'Latex')

% UNCOMMENT TO PRINT TO PDF 
% print2pdf('gpd_conusmption_2021Q4', 1);



%  EOF 













% aout = ols(c,y,'Consumption');
% plot(aout.uhat)
% nexttile
% %%
% t = tiledlayout(2,2,'TileSpacing','Compact');
% 
% % Tile 1
% nexttile([1 2])
% plot(aout.uhat)
% 
% % Tile 2
% nexttile
% plot(rand(1,20))
% title('Sample 2')
% 
% % Tile 3
% nexttile
% plotacf(aout.uhat)
% subplot(2,2,3:4)
