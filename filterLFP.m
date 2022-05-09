function filteredData=filterLFP(data,samplingRate,filterBand)
filter_band=filterBand;
order = 4;
buttertype = '';
freqlim=filter_band;
if freqlim(1)==0
    buttertype='low';
    freqlim=freqlim(2);
end
if freqlim(2)==inf
    buttertype='high';
    freqlim=freqlim(1);
end
fwdflag = true;
revflag = true;
lfp_SamplePeriod=1/samplingRate;
% 'butter' requires freqs spec'd with 1.0 corresponding to half the sample rate. 
%freqlim must be 0.0 < freqlim < 1.0, with 1.0 corresponding to half the sample rate.
%freqlim(1) cannot be zero for default argument of bandpass
freqlim = freqlim*lfp_SamplePeriod*2;
if isempty(buttertype)
    [z, p, k] = butter(4, freqlim);
else
    [z, p, k] = butter(4, freqlim, buttertype);
end
[sos,g]=zp2sos(z,p,k);
h2=dfilt.df2sos(sos,g);
if revflag
    result = filter(h2, data(end:-1:1));
    if ~fwdflag
        result = result(end:-1:1);
    end
end
if fwdflag
    if revflag
        result = filter(h2, result(end:-1:1));
    else
        result = filter(h2, data(:));
    end
end

filteredData=result;
end