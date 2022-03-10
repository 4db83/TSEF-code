clear all; clc;
% ghostscript and graphics handles.
p = graph_h;						

fout	= 'ma1';
dir_	='../lectures/graphics/';
plt		= 1% ghostscript and graphics handles.p
;

f = inline('b./(1+b.^2)','b');
x = linspace(-3,3,100);
plot(x,f(x),'LineWidth',1);hold on;
setplot([.9 .6],13);
hline(0,'k-');vline(0,'k-')
vline(-1,'k--');vline(1,'k--')
hline(0.4,'r-');vline(0.5,'r-');vline(2,'r-')
xlabel('$\beta_1$','Interpreter','Latex')
ylabel('$\rho(1)$','Interpreter','Latex')
set(gca,'YTick',-.5:.1:.5)
setylabel(.25);
setxlabel(-1e-16);
hold off;		

if plt == 1
	print2pdf([dir_ fout]);
end;
