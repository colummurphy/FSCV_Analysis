function filtered=gauss(y,sigma)
%y=CV_mat_DA1(:,1);

%sigma=5;
sz=length(y);
x = linspace(-sz/2,sz/2,sz);
gaussFilter=exp(-x.^2/(2*sigma^2));
gaussFilter=gaussFilter/sum(gaussFilter);
yfilt=conv(y,gaussFilter,'same');
filtered=yfilt;
end