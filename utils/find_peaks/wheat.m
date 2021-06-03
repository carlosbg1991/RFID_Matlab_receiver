% Comparison of CLS and ILS multivariate calibration
% for an experimental data set.
format compact
clear
load wheatx;
load wheaty;
A=wheatx';
C=wheaty';
Abar=mean(wheatx');
Cbar=mean(wheaty');
C=center(reduce(wheaty))';
NumStds=length(C)
NumPoints=length(A(:,1))
w=[1:NumPoints];
s=[1:NumStds];

% Show that single-wavelength calibration is not feasible
subplot(221)
plot(wheaty',wheatx','x')
xlabel('Actual Concentration');
ylabel('Signal at peak wavelength');
title('Analytical curves at each wavelength')

% Calibration by Classical Least Squares (K matrix method)
% Using straight MLR (least-squares regression) formulation
% C matrix is augmented with a column of 1's for background correction
Caug=[C' ones(size(s))']'; 
K=A*Caug'*inv(Caug*Caug');

subplot(222)
plot(K)
title('Spectra of components by CLS');
xlabel('Wavelength');
ylabel('Absorbance');

% Calibration by Inverse Least Squares (K matrix method)
% Need use only one of the following more-or-less equivalent methods:
% Using straight MLR (least-squares regression) formulation
% P=C*A'*inv(A*A');
% Using pseudoinverse
% P=C*pinv(A);
% Using right division
% P=C/A;
% Using left division
P=(A'\C')';      

% Prediction from original calibration set
B=A;

% CLS method
KMatrixPred=inv(K'*K)*K'*B;
Y=C(1,:);
Ybar=mean(Y);
Yhat=KMatrixPred(1,:);
% Compute Standard Error of Estimate
KStdErrEst1=sqrt(sum((Yhat-Y) .^2)./(NumStds-1));
% Compute R-squared
KRSquared1=(sum((Yhat-Ybar).^2))./(sum((Y-Ybar).^2));

subplot(223)
plot([min(C) max(C)],[min(C) max(C)],C(1,:),KMatrixPred(1,:),'o')
title('Using Classical Least Squares');
xlabel('Actual Concentration');
ylabel('Measured Concentration');

% P-matrix method
PMatrixPred=P*B;
Y=C(1,:);
Ybar=mean(Y);
Yhat=PMatrixPred(1,:);
% Compute Standard Error of Estimate
PStdErrEst1=sqrt(sum((Yhat-Y) .^2)./(NumStds-1));
% Compute R-squared
PRSquared1=(sum((Yhat-Ybar).^2))./(sum((Y-Ybar).^2));

subplot(224)
plot([min(C) max(C)],[min(C) max(C)],C(1,:),PMatrixPred(1,:),'o')
title('Using Inverse Least Squares');
xlabel('Actual Concentration');
ylabel('Measured Concentration');

subplot

