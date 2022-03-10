%% Example MA(1) projection from page 74 in Brockwell and Davies
clc;clear all;clf;
%% analytic (requires symbolic math toolbox)
syms b
p = b/(1+b^2);

R = [1   p   0   0   0;
     p   1   p   0   0;
     0   p   1   p   0;
     0   0   p   1   p;
     0   0   0   p   1;];

r = [p;
     0;
     0;
     0;
     0;];

t   = 2;
sol = inv(R(1:t,1:t))*r(1:t)
sol2=simplify(sol);
pretty((sol2))

% substitute in numbers.
subs(sol2,b,-0.9)

%% Numerical MA(1) example taken from page of 74 Brockwell and Davies with 
T = 50;
bet     = -0.9;
rho1    = bet/(1+bet^2);
R       = diag(ones(T,1)) + diag(repmat(rho1,1,T-1),1) + diag(repmat(rho1,1,T-1),-1);
r       = zeros(T,1);r(1)=rho1;

for t = 1:4
inv(R(1:t,1:t))*r(1:t)
end

comp  = [inv(R)*r bet*[1; -garchar([],bet,T-1)]];
% Note: Minus here because matlab returns the wrong sign for garchar
%       Also, the garchar representation does not return 1 for Psi_0, so we do not need 
%       to multipliy by thete as in the hand formula we would have 

fprintf(' %2.5f %2.5f \n', comp(1:20,:)') % prints by row