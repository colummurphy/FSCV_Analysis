function rate=getSampleRate(timestamps)
% determine sampling rate for LFP data, taken from dan's code
% We assume that the first two frames were recorded continuously, and
% calculate the sample period from them:
% Find the sample period, ignoring timestamp intervals that
% indicate a break in recording:
sampledurs = timestamps(2:end)-timestamps(1:end-1);
minsd = min(sampledurs);
mediansd = median(sampledurs);
mediansample = mediansd;
sampledurs = sampledurs(sampledurs < 2*minsd);
if std(sampledurs(sampledurs < minsd + 1.01*mediansample)) > 1e-6
    warning('samples have jitter.');
end
if std(sampledurs) > mediansample || mediansd > minsd + mediansample
    warning('sample durations are of excessively variable duration.');
end
samplePeriod = median(sampledurs);
rate = 1/samplePeriod;

end