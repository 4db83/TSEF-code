clear all; clc;
% addpath('d:/matlab.tools/db.toolbox/')

% A = [0.4	0.1; 0.2	0.5];
% S = [0.25	0.3; 0.3	0.9];
% A = [ 0.9569265097 0.007330814144 ; -24.16479184  0.5478001881];
% S = [14.2270332454324	-194.229041071163; -194.229041071163	22150.4568298334];

% Lutkepohl exmple pages 56 onwards
A = [.5 0 0;.1 .1 .3;0 .2 .3];
S = [2.25 0 0; 0 1 .5; 0 .5 .74];


P = chol(S)';
d = diag(P);
D = diag(d);
%P = P*inv(D);
N = 20;
k = size(A,1);
irf = zeros(N+1,k^2);
phi = zeros(k,k,N+1);

diag_sphi = zeros(N+1,k);
sumirf2 = zeros(N+1,k^2);

for i=1:N+1
  tmp_  = A^(i-1)*P;
  irf(i,:) = (tmp_(:))';
  phi(:,:,i) = tmp_;
  phi2(:,:,i) = phi(:,:,i).^2;
  sum_phi2(:,:,i) = sum(phi2(:,:,1:i),3);
  irf2(i,:) = irf(i,:).^2;
  sumirf2(i,:) = sum(irf2(1:i,:),1);
  phiphi(:,:,i) = phi(:,:,i)*phi(:,:,i)';
  sphi(:,:,i) = sum(phiphi(:,:,1:i),3);
  diag_sphi(i,:) = (diag(sphi(:,:,i)))';
  %fevd(i,1) = sum(ej'*phi(:,:,i)*el)^2/diag_sphi(j)
end
fevd = sumirf2./repmat(diag_sphi,1,k);
reshape(fevd,(N+1)*k,k)


subplot(1,k,1)
plot(irf(:,[1:k]))
legend({'$U_1$','$U_2$','$U_2$'},'Interpreter','Latex')
title('Response of $X_1$ to shocks in','Interpreter','Latex')
hline(0,'k-');
subplot(1,k,2);
plot(irf(:,[k-1:2*k]));
legend({'$U_1$','$U_2$','$U_2$'},'Interpreter','Latex')
title('Response of $X_2$ to shocks in','Interpreter','Latex')
hline(0,'k-');
subplot(1,3,3);
plot(irf(:,[2*k-1:3*k]));
legend({'$U_1$','$U_2$','$U_2$'},'Interpreter','Latex')
title('Response of $X_2$ to shocks in','Interpreter','Latex')
hline(0,'k-');
irf
