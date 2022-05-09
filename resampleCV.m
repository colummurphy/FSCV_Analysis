function resampled=resampleCV(data, currentRate, targetRate)
%data is array of CV's with row being sample for each CV, and column is a
%separate CV, currentRate = # of samples currently per CV, targetRate is #
%of samples per CV targeted
sigma=3; ren=1; reb=1;
%resample assumes edge values are zero so need to pad
%08/14/2017
%fs1 = 10;             % Original sampling frequency in Hz
%t1 = 0:1/fs1:1;       % Time vector
x = data;         % Define a linear sequence
%resampled=[];
padratio=0.08;
if targetRate>=214 && targetRate<=252
    padratio=1;
end
if targetRate==175 && currentRate>=214 && currentRate<=240
    sigma=2;
    padratio=1;
end
currentpaddedrate=padratio*2*currentRate;
targetpadrate=targetRate*2*padratio;
for ii=1:size(data,2)
    xpad = [repmat(x(1,ii), 1, round(padratio*currentRate)), x(:,ii)', repmat(x(end,ii), 1, round(padratio*currentRate))];
    %tpad = [-1/fs1*10 : 1/fs1: 0-1/fs1, t1, 1+1/fs1:1/fs1:1+1/fs1*10];
    
    ypad = resample(xpad,targetpadrate,currentpaddedrate,ren,reb);  % Now resample it
    extracted=ypad(round(padratio*targetRate)+1:(end-round(padratio*targetRate)));  %extract without pads now
 % resampled(:,ii)=extracted;
    resampled(:,ii)=gauss(extracted,sigma);
end

end