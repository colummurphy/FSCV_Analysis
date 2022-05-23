function matchingSamples=findNearestTS(targetTS,allTS,deviation)
%Find sample of all TS that is closest in value to the target TS
%allTS can be an array where each row is a separate trial and targetTS is a
%vector that contains the TS to target for each trial. Deviation is the
%maximum allowable difference in value (in s). If the closest value provides a
%difference > deviation then error. Usually a percentage of the sampling
%period. Eg. 1kHz sampling rate, 1 ms sample period, 50% error = 0.5 ms.



TSdiff=allTS-targetTS;
[~,minSamples]=min(abs(TSdiff'));
ind = sub2ind(size(allTS),1:size(allTS,1),minSamples);  %Convert array of indeces to vector of indeces for the TS (Col) for each trial (row)
minTS=allTS(ind);
TSerror=abs(targetTS-minTS');
if find(any(TSerror>deviation))
    STOP;
    warning('TS for alignment > 50% error')
end

matchingSamples=minSamples';