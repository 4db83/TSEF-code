clear all; clc;
% [GS,p]=graph_h(140,8.60);						% ghostscript and graphics handles.

% fout = 'rec0';
% dir_ ='D:/teaching/TSE2011/lectures/graphics/';


seed(4321); 

n = 5e3;
x = 1 + randn(n,1);
x_bar = nan(n,1);

for i = 10:n
%   x_bar(i) = sqrt(i)*(mean(x(1:i)) - 1);
  x_bar(i) = mean(x(1:i));
end

%%
clf; 
set(groot,'defaultLineLineWidth',2.5); 
p.fs = 22;			% FONT SIZE                                       
p.dm = @(x)([.8-(x-1)*.215 .86 .17]);    % FIG DIMENSION      

  plot(x_bar);
 	
  hline(1,'k-');
%   hline(0,'k-');
	grid on; box on;
  setplot([.8 .23],p.fs)
  setyticklabels(0.6:.1:1.2,1)
  
	set(gca,'GridLineStyle',':','GridAlpha',1/3);
	setoutsideTicks; add2yaxislabel
  

print2pdf('mu_hat','../lectures/graphics/')
% print2pdf('Z_hat','../lectures/graphics/')