% plot consumption and disposable income data for the US to replicate Hamilton figure on page 601
% ---------------------------------------------------------------------------------------------
clear; clc;
set(groot,'defaultLineLineWidth',2); % sets the default linewidth to 1.5;

% % get data from FREDwebsite if not in ./data/ directory
% usdat = as_timetable(getFredData('pcecc96', '1947-01-01', '2021-12-31','lin','q'),'pcecc96');
% % usdat = synchronize(usdat,as_timetable(getFredData('USRECQ', '1947-01-01', '2021-12-31','lin','q'),'NBER'));
% usdat = synchronize(usdat,as_timetable(getFredData('dpic96', '1947-01-01', '2021-12-31','lin','q'),'dpic96'));
% % combs = synchronize(us_data_merged,usdat);
% save('./data/consumption_income_US_2021Q4.mat', 'usdat');

% load consumption_income_US data for us
load('./data/consumption_income_US_2021Q4.mat')
clc;clf;
usdat.c = log(usdat.pcecc96);
usdat.y = log(usdat.dpic96);
% head2tail(usdat);
tt = timerange('Q1-1947', 'Q4-2021', 'closed');
tt = timerange('Q1-1947', 'Q4-1989', 'closed');
tt = timerange('Q1-1947', 'Q4-2019', 'closed');
ss_usdat = usdat(tt,:);

y = ss_usdat.y;
c = ss_usdat.c;

aout = ols(c,y);
% tiledlayout(2,2); 


plot(aout.uhat)
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


%%
clf;
Fns = 22;
hold on;
  plot(ss_usdat.c);
  plot(ss_usdat.y);
hold off;
grid on; box on;
% setplot([.36 .86 .21],.25)
setplot([.57 .26]);
% setyticklabels(9.5:.5:11.5,1,14)
setdateticks(ss_usdat.Time,31,'yyyy', Fns);	
set(gca,'GridLineStyle',':','GridAlpha',1/3);
setoutsideTicks
% add2yaxislabel
% ADD LEGEND
legendflex({'Consumption $(c_t)$','Disposible Income $(y_t)$'}, 'fontsize', Fns + 1, 'anchor',[1 1],'Interpreter','Latex')
% print2pdf('gpd_conusmption_2021Q4','../graphics/');



%  EOF 