% plot HP-filter cycle and BN based on AR(2)
clear; clc;
% LOAD THE DATA
% % uncomment get new data from FRED;
% beg_date = '1947-01-01';
% end_date = '2022-03-31';
% usdata = as_timetable(getFredData('GDPC1',beg_date,end_date,'lin','q'),'gdpc1'); 
% usdata = synchronize(usdata,as_timetable(getFredData('USRECQ',beg_date,end_date,'lin','q'),'NBER'));
% usdata = synchronize(usdata,as_timetable(getFredData('GDPPOT',beg_date,end_date,'lin','q'),'CBO'));
% save('./data/usdata_2022Q1.mat', 'usdata');
% disp(' ... done getting and saving the data')
disp(' ... loading the data ... ');
load './data/usdata_2022Q1.mat';   % usdata

% make log gdp and gdp growth data
usdata.y  = 100*log(usdata.gdpc1); % multiply by 100 to be expressed in percent, no need to scale later
usdata.dy = delta(usdata.y);
usdata.cbo_gap = 100*(usdata.gdpc1 - usdata.CBO)./usdata.CBO;
usdata.tt = (1:size(usdata,1))';
% truncate sample to required time period
ss = timerange('Q1-1947', 'Q4-2019', 'closed');
usdata = usdata(ss,:);
clc;
[hp_c, hp_t] = hp_filter(usdata.y,1600);
usdata.hp_c = hp_c;
usdata.hp_t = hp_t;
head2tail(usdata)

% estimate AR(2) for dy gdp growth
clc; 
ar_order = 1;
ar2_out = estimate_armax(usdata.dy,1,1:ar_order,0);
print_arma_results(ar2_out)
aL  = ar2_out.aL;
aL1 = sum(aL);
% make the BN permanent compoonent as
x   = usdata.y;
x_p = 1/aL1*[x mlag(x,ar_order)]*(aL)';
x_t = x - x_p;
usdata.x_p = x_p;
usdata.x_t = x_t;


xp = x + (ar2_out.pars(2))/aL1*( delta(x,1)-ar2_out.mu )
xt = x - xp;
usdata.xt = xt;

% PLOT CONTROLS
clf
set(groot,'defaultLineLineWidth',2); % sets the default linewidth;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
% recession bar color
rec_CLR = .8*ones(3,1); 
pl.ds = 22;
pl.fs = 20;
pl.st = -1.21;
pl.ps = @(x) ([.07 x .86 .70]);

% adjust sample for plotting
sss = timerange('Q2-1949', 'Q2-2016', 'closed');
usdata = usdata(sss,:);

% now plot
hold on;
	bar( usdata.NBER*8, 1, 'FaceColor', rec_CLR); 
	bar(-usdata.NBER*8, 1, 'FaceColor', rec_CLR); 
%   plot(usdata.hp_c,     'Color', clr(1));
  plot(usdata.cbo_gap,  'Color', clr(2));
  plot(usdata.x_t,      'Color', clr(3));
  plot(usdata.xt,      'Color', clr(4));
hold off; 
box on; grid on;
setyticklabels(-8:2:8)
setplot(pl.ps(.2),pl.fs,0,6/5)
setdateticks(usdata.Time, pl.ds, 'yyyy:QQ', pl.fs)
hline(-8);hline(8);hline(0); % add lines to top and bottom to kill the bar-plot going over the box lines
setgridlinestyle
setoutsideTicks
add2yaxislabel
tickshrink(.9)


% clc
% 
% tt_cbo = ols(usdata.cbo_gap,usdata.tt);
% tt_hpc = ols(usdata.hp_c,usdata.tt);

% [mean(usdata.hp_c); nanmean(usdata.cbo_gap)]


% clc;
% a1 =  0.34;
% a2 = -0.10;
% 
% A  = [a1 a2; 1 0];
% 
% N = 200;
% aa = zeros(2,2,N);
% for i = 1:N
%   aa(:,:,i) = A^i;
% end
% 
% format long
% 
% sum(aa,3)
% 
% 
% inv(eye(2)-A)*A
