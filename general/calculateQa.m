%Calculate Qa threshold
function Qa=calculateQa(S)
%k=number relevant PC's
%Lambda = sum of squares of all projections
%n = total number of PCs calculated
lambda=diag(S).^2;       %S matrix calculated by SVD function has square root lambda values along diagonal
%index_zscore=150;       %voltage point sample along individual scan to take z score from through recorded time
%c_a=zscore(Isub(index_zscore,:));
c_alpha=1.96;       %z-score for 95% CI 1.645 for 90%
theta1=sum(lambda(2:length(lambda)));
lambda2=lambda.^2;
theta2=sum(lambda2(2:length(lambda2)));
lambda3=lambda.^3;
theta3=sum(lambda3(2:length(lambda3)));
h0=1-(2*theta1*theta3/theta2^2);
Qa=theta1.*(c_alpha.*sqrt(2.*theta2*h0^2)./theta1+1+theta2.*h0.*(h0-1)./(theta1^2)).^(1/h0);
end
