function stat = ols_mt(Y,X)
NObs    = size(Y,1);

% Perform regressions
beta    = inv(X'*X)*(X'*Y);
k       = size(beta,1);

% get the ordinary statics
yhat    = X*beta;
u       = yhat-Y;
SST     = (Y-mean(Y))'*(Y-mean(Y));
SSR     = u'*u;
R2      = 1-SSR/SST ;
std2hat = u'*u/(NObs-k-1);
stdhat  = sqrt(std2hat);

Xinv    = inv(X'*X);
dXi     = diag(Xinv);
sdXi    = sqrt(dXi);
seb     = stdhat*sdXi;
tstat   = beta./seb;    
     
% Generate Output;
stat.beta   = beta;
stat.sde    = seb;
stat.tstat  = tstat;
stat.R2     = R2;
stat.SSR    = SSR;
stat.SST    = SST;
stat.u      = u;
