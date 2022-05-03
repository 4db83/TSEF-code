% generating a Browninan Motion/Wiener Process with different values of T
clc;clear all;
addpath('d:/matlab.tools/db.toolbox/');

T = 1e6;
t = (0:1:T)'/T;                       % t is the column vector [0 1/T 2/T ... T/T=1]
W = [0; cumsum(randn(T,1))]/sqrt(T);  % S is running sum of N(0,1/T) variables
%seed(123)

%% plot the path of the Wiener process
plot(t,W);          
hline(0,'k-')
%title({['Wiener process with $T = ' int2str(T) '$']},'Interpreter','Latex')
setytick(1)
setplot([.6 .5],14);
%ylim([-.801 .401])
%print2pdf(['../lectures/graphics/bm_' num2str(T)])


