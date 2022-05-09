function smoothedData=smoothwin(data,winsize)
%smooth with hann function, got from lfp_bandpassPower.m dan gibson code
%04/02/2018
%to smooth beta power samples after filtering, squaring operations
%make window 35 samples wide in paper or half-width of one cycle of 13 Hz
%target .0384 s wide window

        winfunc = hanning(2*winsize+1);

s1 = conv(data, winfunc);
% s1 has an "extra" <smoothing> points on each end, so the first
% "valid" point is (smoothing + 1) and the last is (end - smoothing).
smoothedData = s1(winsize + 1 : end - winsize);
if winsize > 0
    smoothedData = smoothedData / sum(winfunc);
end
end
