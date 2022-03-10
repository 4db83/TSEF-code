clear all; clc;
% check matlabpool if open
% if ~matlabpool('size') > 0; matlabpool; end
% ghostscript and graphics handles.
p = graph_h;						

fout = 'rec0';
dir_ ='../lectures/graphics/';
seed(25);

n = 5e4;
x = 1 + randraw('cauchy', [], n);

parfor(i = 30:n)
  x_bar(i) = mean(x(1:i));
end
disp('done')
%%
% plot(x_bar,'LineWidth',1);hold on;
% setplot([.8 .5], 11)
% hline(1,'k-');
% hold off;
%print2pdf([dir_ fout 'b']);

clf; 
set(groot,'defaultLineLineWidth',2.5); 
p.fs = 22;			% FONT SIZE                                       
p.dm = @(x)([.8-(x-1)*.215 .86 .17]);    % FIG DIMENSION      

  plot(x_bar);
 	
  hline(1,'k-');
%   hline(0,'k-');
	grid on; box on;
  setplot([.8 .23],p.fs,0)
%  setyticklabels(0.6:.1:1.2,1)

  % ADD GRIDLINES AND SECOND AXIS
	set(gca,'GridLineStyle',':','GridAlpha',1/3);
	setoutsideTicks; add2yaxislabel; 
  grid on; box on;
  
%%
print2pdf('cauchy_mu_hat','../lectures/graphics/')

