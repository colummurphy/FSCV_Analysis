function [p,R]=calcSlope(x,y)

%Calculate slope of data points assuming linear fit
% p = p1x1 + p2
%y=data(xx1(1),1:100);
%x=(1:100)./10;


[p,S]=polyfit(x,y,1);   %Polynomial fit of order n = 1 (Linear)
Rsqr=(S.normr/norm(y - mean(y)))^2; %Coefficient of determination
R=sqrt(Rsqr);   %Correlation, how good is fit