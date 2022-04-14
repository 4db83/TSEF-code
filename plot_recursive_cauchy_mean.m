% Script: plot_recursive_cauchy_mean.m
% NOTE: To be able to run this code, you need the contents of the db.toolbox available from:
% https://github.com/4db83/db.toolbox/archive/refs/heads/main.zip. Unzip the contents locally
% to the same directory as this script, and then uncomment the following line below:
% addpath(genpath('./db.toolbox-main'))
clear;clc;clf;

% set simulation seed for reproducibility
seed(25);
n = 5e4;
x = 1 + randraw('cauchy', [], n);
mu = zeros(n,1);

% recursively compute the mean
parfor(i = 30:n)
  mu(i) = mean(x(1:i));
end
disp('done')

%% plot the recursive mean estimate
clf; 
set(groot,'defaultLineLineWidth',2.5); 
fig.fs = 21;			                        % FONT SIZE                                       
fig.dm = @(x)([.1 .8-(x-1)*.215 .8 .23]);   % FIG DIMENSION      

plot(mu);
hline(1);
grid on; box on;
setplot(fig.dm(2), fig.fs, 0);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
% setyticklabels(ytcks, 1); 
setoutsideTicks
add2yaxislabel
% tickshrink(.8)

% uncomment to print to pdf 
% print2pdf('cauchy_mu_hat2','../graphics')


% 
% % ADD GRIDLINES AND SECOND AXIS
% 	set(gca,'GridLineStyle',':','GridAlpha',1/3);
% 	setoutsideTicks; add2yaxislabel; 
%   grid on; box on;

%% now plot the densities
xg = linspace(-14,14,5e3)';
cauchypdf = @(x,mu,sig) (pi*sig*(1+((x-mu)./sig).^2)).^(-1);

fig.fs = 20;			                          % FONT SIZE                                       
fig.dm = @(x)([.1 .8-(x-1)*.215 .8 .25]);   % FIG DIMENSION      

clf;clc
hold on;
  plot(xg,cauchypdf(xg, 0,1  )) 
  plot(xg,cauchypdf(xg, 1,2  ))
  plot(xg,cauchypdf(xg,-1,0.5))
hold off;
grid on; box on;
setplot(fig.dm(2), fig.fs, 1);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
setxticklabels(-12:2:12); 
setoutsideTicks
add2yaxislabel

% uncomment to print to pdf 
print2pdf('cauchy_pdfs','../graphics')


