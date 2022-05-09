function [ts, ts2] = lfp_findSpikesFSCV(filenum, thresh, lfp_Samples,lfp_SamplePeriod,varargin)
%HELEN MODIFIED 08/2021, remove globals and put inputs as called from 
%extractSpikesFSCV
%09/11/2021 remove frame necessity from index2time code below by inputting
%timestamps directly
%Just output sample IDS not timestamps
%ts and ts2 are in samples now


%ts = lfp_findSpikes(filenum, thresh)
%   Searches CSC channel <filenum> for waveforms that resemble action
% potentials.  It is assumed that the channel is suitably high-pass
% filtered to remove high amplitude low frequency components, but not so
% much as to distort a typical spike waveform (e.g. with a cutoff frequency
% in the 50 - 100 Hz range).
%INPUTS
% filenum: filenum of the CSC channel to analyze.
% thresh: the spike waveform must exceed this threshold in the
%   positive-going direction. 
%OUTPUTS
% ts: column vector of timestamps denoting the peak time of each
%   positive-going spike.
% ts2: column vector of timestamps denoting the peak time of each biphasic
%   spike comprising a valley deeper than -thresh followed by a peak
%   greater than thresh.
%OPTIONS
% 'invert' - inverts the signal in <filenum> before analyzing.
%NOTES
% For each sample at which the signal crosses <thresh>, the interval
% <win> around that sample is analyzed to determine whether it contains a
% plausible spike waveform by applying the following criteria:
%   0. The entire window around the trigger must exist.
%   1. It must not be a local minimum that just touches <thresh>; this
%   makes threshold crossing a three-point criterion, so it is defined to
%   be unsatisfied at the first sample.
%   2. It must not contain any samples of magnitude greater than
%   10*<thresh>.
%   3. The maximum must be in the middle 50% of <win>.
%   4. The last sample must be within <thresh>/2 of the first sample.
% A good default value for <thresh> is 5*SD, where lfp_XLimAll and
% lfp_AlignmentRef have been set appropriately to select the time window to
% be analyzed, and then lfp_findbadtrials with an appropriate 'magfactor'
% is used to eliminate trials with huge signals in the analysis window:
%     lfp_BadTrials = lfp_findbadtrials(filenum, ...
%         'windows', {lfp_AlignmentRef}, lfp_XLimAll, ...
%         'magfactor', 10);
%     sampledata = lfp_getSamples([], filenum, []);
%     SD = std(sampledata(:));
% In severe cases, the lfp_findbadtrials call can be run iteratively until
% no bad trials are found.
%   Some highpass filtering is necessary to achieve reliable triggering on
% potential spikes, but the cutoff frequency should be no higher than 100
% Hz to avoid grossly distorting the spike waveform (50 Hz is better if
% that doesn't disturb the triggering too much).  Some lowpass filtering is
% also desirable so that <peakidx> doesn't get thrown off by high frequency
% noise.


% 08/02/21 from DG:
%lfp_Samples variable can be in frames or a vector, up to you. 
%It's addressed as a vector, which will work with both formats.  
%However, if you're planning to pass it to 'dg_writeCSC' eventually,
%better to format it in frames so you don't have to remember to reshape it before calling 'dg_writeCSC'.

%$Rev: 403 $
%$Date: 2019-08-27 16:52:36 -0400 (Tue, 27 Aug 2019) $
%$Author: dgibson $

%global lfp_Samples lfp_SamplePeriod

pol = 1; % polarity of analysis
argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'invert'
            pol = -1;
        otherwise
            error('lfp_findSpikes:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

win = [-0.0005 0.0015]; % window around trigger in seconds
winidx = round(win(1)/lfp_SamplePeriod) : round(win(end)/lfp_SamplePeriod);
q1idx = round(length(winidx)/4);
q3idx = round(3*length(winidx)/4);

% <xingidx> points at the last sample that was below <thresh>:
xingidx = find( [ false
    reshape(pol * lfp_Samples{filenum}(3 : end) > thresh, [], 1) ...
    & reshape(pol * lfp_Samples{filenum}(2 : end - 1) <= thresh, [], 1) ...
    & reshape(pol * lfp_Samples{filenum}(1 : end - 2) <= thresh, [], 1)
    ] );
peakidx = zeros(size(xingidx));
peak2idx = zeros(size(xingidx));
peakidx2 = 1;
peak2idx2 = 1;
for xingidx2 = 1:length(xingidx)
    if xingidx(xingidx2) + winidx(1) < 1 || ...
            xingidx(xingidx2) + winidx(end) > numel(lfp_Samples{filenum})
        continue
    end
    if any( abs(pol * lfp_Samples{filenum}(xingidx(xingidx2) + winidx)) ...
            > 10 * thresh ) ...
            || pol * lfp_Samples{filenum}(xingidx(xingidx2) + winidx(end)) ...
            - pol * lfp_Samples{filenum}(xingidx(xingidx2) + winidx(1)) ...
            > thresh/2
        continue
    else
        [~, maxidx] = max(pol * lfp_Samples{filenum}(xingidx(xingidx2) + ...
            winidx));
        if maxidx < q1idx || maxidx > q3idx
            continue
        end
        prespikewin = maxidx + (-round(0.0002/lfp_SamplePeriod):0);
        if any(pol * lfp_Samples{filenum}( xingidx(xingidx2) ...
                + winidx(prespikewin) ) < -thresh)
            peak2idx(peak2idx2) = maxidx + xingidx(xingidx2) + winidx(1) ...
                - 1;
            peak2idx2 = peak2idx2 + 1;
        else
            peakidx(peakidx2) = maxidx + xingidx(xingidx2) + winidx(1) - 1;
            peakidx2 = peakidx2 + 1;
        end
    end
end
peakidx(find(peakidx==0, 1) : end) = [];
%ts = reshape(lfp_index2time(peakidx), [], 1);           %THIS REQUIRE FRAMES
%adjusted here 09/11/2021
%just export indexes in sample indexes rather than converting to timestamp
ts=peakidx;
%ts = reshape(lfp_index2time(peakidx), [], 1); 
peak2idx(find(peak2idx==0, 1) : end) = [];
%ts2 = reshape(lfp_index2time(peak2idx), [], 1);
ts2=peak2idx;
   

