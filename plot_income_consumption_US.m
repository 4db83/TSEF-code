% example of ARMA model fit to US Real GDP;
clc;clear all;clf;
% use addpath(genpath ...) to add subdirectories as well.
%addpath(genpath('d:/matlab.tools/mfe.toolbox'))
addpath('d:/matlab.tools/db.toolbox')
g_ = graph_h;
% check matlabpool if open
if ~matlabpool('size') > 0; matlabpool; end
% load the data with path to data file
disp('loading the data ...')
if exist('./data/Consumption and GDP.mat') 
	load './data/Consumption and GDP.mat';
else
	[data0, dates0] = xlsread('./data/Consumption and GDP.xls');
	dates = datenum(dates0(8:end,1),'dd/mm/yyyy');				% and for windows xp
	%dates = datenum(dates0(8:end,1));% ,'dd.mm.yyyy');	% use this for windows 7 above
	X = data0;
    save 
    % data is already in loggd form
    y = X(:,1);
    c = X(:,3);
    T = size(y,1);
    save './data/Consumption and GDP.mat' dates y c
end;
disp('done')

%% plots
% Plot the level series of (log) real GDP and real consumption
plot(dates,y,g_.lw,1); hold on;
plot(dates,c,'r',g_.lw,1);
legend('log(income)','log(consumption)','Location','NorthWest')
hold off;
setplot([.8 .4],10)
xlim([dates(1) dates(end)]);
set(gca,'XTick',dates([1:24:end]));
setytick(1)
ylim([7 10])
hline(0,'k')
datetick('x','yyyy:QQ','keepticks','keeplimits');
%print2pdf('../lectures/graphics/C_Y_US');

%% Scatterplot the level series of (log) real GDP and real consumption
%scatter(y,c,g_.mk,'.',g_.ms,2)
scatter(y,c,10,'b','fill');
box on;
setplot([.8 .4],9.8);
setytick(1);
xlabel('log(income)');
ylabel('log(consumption)');
%axis tight
ylim([floor(min(c)) 9.5])
xlim([7.5 9.75])
%print2pdf('../lectures/graphics/scatter_C_Y_US');
% % 
% % %% Terms structure data 
% % clearvars -except dates g_;clc;
% % d0 = find(datenum('1962:Q1','yyyy:QQ')==dates);
% % d1 = find(datenum('2010:Q3','yyyy:QQ')==dates);
% % dates = dates(d0:d1);
% % load './data/usmacro.mat'
% % y1 = r1yr; y5 = r5yr; y10 = r10yr;
% % plot(dates,y1,g_.lw,1); hold on;
% % plot(dates,y5,'r',g_.lw,1);
% % plot(dates,y10,'g',g_.lw,1);
% % hold off;
% % setplot([.8 .4],10)
% % xlim([dates(1) dates(end)]);
% % DS = 24;		% date spread
% % set(gca,'XTick',dates([1:DS:end]));
% % setytick(0)
% % datetick('x','yyyy:QQ','keepticks','keeplimits');
% % ylim([0 18]);
% % legend('1 year','5 year','10 year')
% % %print2pdf('../lectures/graphics/us_yields');
% % 
% % %% scatter plot now
% % scatter3(y1,y5,y10,10,'b','fill');
% % ylim([0 18]);xlim([0 18]);zlim([0 20]);
% % setplot([.8 .6],10);
% % setytick(0);
% % xlabel('1 Year');
% % ylabel('5 Year');
% % zlabel('10 Year');
% % %set(gca,'ZTickLabel',[0:2:18])
% % %print2pdf('../lectures/graphics/scatter_us_yields');
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
