% testplot
clf;clc;clear
seed(123)
x = randn(100,4);

% plots
p.dm = @(x,y)([x y .4 .2]);% FIG DIMENSION
p.fs = 15;			% FONT SIZE                                       
p.df = 21;			% DATE FREQUENCY                                  
p.aj = -1.225;	% SUBTITLE POSITION ADJUSTMENT                  
 		

subplot(2,2,1);
  plot(x(:,1));
  hline(0); box on; grid on;
	setplot(p.dm(.055,.6),p.fs,1)
	% MAKE OUTSIDE TICKS FOR TWO AXIS Y LABELS 
	set(gca,'GridLineStyle',':','GridAlpha',1/3)
  setoutsideTicks
  add2yaxislabel
  subtitle('(a) $\alpha x_t$', p.aj, p.fs, 1)

subplot(2,2,2);
  plot(x(:,2));
  hline(0); box on; grid on;
	setplot(p.dm(.55,.6),p.fs,1)
	% MAKE OUTSIDE TICKS FOR TWO AXIS Y LABELS 
	set(gca,'GridLineStyle',':','GridAlpha',1/3)
  setoutsideTicks
  add2yaxislabel
  subtitle('(b) $\alpha x_t$', p.aj, p.fs)

subplot(2,2,3);
  plot(x(:,3));
  hline(0); box on; grid on;
	setplot(p.dm(.055,.3),p.fs,1)
	% MAKE OUTSIDE TICKS FOR TWO AXIS Y LABELS 
	set(gca,'GridLineStyle',':','GridAlpha',1/3)
  setoutsideTicks
  add2yaxislabel
  subtitle('(c) $\alpha x_t$', p.aj, p.fs)
  
subplot(2,2,4);
  plot(x(:,4));
  hline(0); box on; grid on;
	setplot(p.dm(.55,.3),p.fs,1)
	% MAKE OUTSIDE TICKS FOR TWO AXIS Y LABELS 
	set(gca,'GridLineStyle',':','GridAlpha',1/3)
  setoutsideTicks
  add2yaxislabel
  subtitle('(d) $\hat{\alpha} x_t$', p.aj, p.fs, 1)
  