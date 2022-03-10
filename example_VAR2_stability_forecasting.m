% VAR(2) stability and correlations
clear all;clc;
% set up parameters
p = 2; % lag order
k = 2; % number of variables

A1	= [0.5 0.1 ;...
			 0.4 0.3];
		
A2	= [-0.2 0.1 ;...
			 -0.3 0.2]; 

c		= [0.2;0.3];

sig = [1.75 0.25 ;...
			 0.25 3.00];

% set up companion form parameters		
Ik	= eye(k);
Ok	= zeros(k,k);

A		= [A1 A2 ;...
		   Ik Ok];

C		= [c;Ok(:,1)];

Sig = [sig Ok;...
			 Ok Ok];
	 
% lagploynomial (I - A1L - A2L) using symbolic math toolbox
syms z
AL		=  Ik - A1*z - A2*z^2;	% symbolic polynomial
detAL = det(AL);							% det of symbolic polynomial
%%
fprintf('Lag polynomial expression when evaluating det(AL) \n');
fprintf(' %s \n\n', char(detAL));

% convert back to numerical polynomial
lag_roots = roots(sym2poly(detAL));
fprintf('Roots of Lag-polynomial:\n');
disp(lag_roots);

% compute eig(A)
eigA	= eig(A);
fprintf('\nRoots of charactistic-polynomial:\n');
disp(eigA);

%% unconditional mean
mu		= inv(Ik - A1 - A2)*c;

% or from companion from
S			= [Ik Ok];						% selection vector
Mu		= inv(eye(p*k)-A)*C;
muC		= S*Mu;

%% variance/covariance matrix
Gam0	= reshape(inv(eye(16)-kron(A,A))*Sig(:),4,4); 
gam0	= Gam0(1:2,1:2);
gam1	= Gam0(1:2,3:4);

% computing the autocovariance and autocorrelation recursions for
J			= 10;							% highest number of ACVs and ACFs to compute
% ACVFs
Gams	= zeros(k,k,J);
Gams(:,:,1) = gam0;
Gams(:,:,2) = gam1;
% ACFs
d			= sqrt(diag(diag(gam0)));
Rhos	= zeros(k,k,J);
Rhos(:,:,1) = inv(d)*gam0*inv(d);
Rhos(:,:,2) = inv(d)*gam1*inv(d);

for j = (p+1):J
	Gams(:,:,j) = A1*Gams(:,:,j-1) + A2*Gams(:,:,j-2);
	Rhos(:,:,j) = inv(d)*Gams(:,:,j)*inv(d);
end;

fprintf('\n----------------------------------------\n');
fprintf('  		  ACVs				  ACFs \n');
fprintf('----------------------------------------\n');
disp([Gams Rhos])
fprintf('----------------------------------------\n');

%% forecasting 
H		= 20;							% upper forecast horizon cap.
Xt	= [1 2 3 4]';			% (pkx1) vector of Observed X_t at time t to use a base to forecast from.
% constructing the Psi weights
Psi	= zeros(k,k,H);
% create psi function
psi = @(A,h) (S*A^h*S');

% looping through different time horizons
for jj = 1:H
	Psi(:,:,jj) = psi(A,jj-1);
end;

for h =	1:H;							
	Xhat_h(:,h)		= S*(Mu + A^h*Xt);						% h-step ahead forecasts
	SigU_h(:,:,h) = Psi(:,:,h)*sig*Psi(:,:,h)';	% h-step ahead Sig weiths that need to be summed.
end

% plot the forecasts 
hold on;
plot(1:H,Xhat_h(1,:),'b');
plot(1:H,Xhat_h(2,:),'r');
hline(mu(1),':b');
hline(mu(2),':r');
hold off;

% forecast error variance converges to gam0
SigU_H = sum(SigU_h,3);
disp(' Forecast error variance matrix');
disp(SigU_H)
disp(' Gamma(0) covariance matrix');
disp(gam0)



