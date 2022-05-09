function [morletdata, allP] = morletgram(data,samplingfreq,freqlim,varargin)
%copied from dg_morletgram 04/10/2018 hns
%morletdata = lfp_morletgram(trials, filenums, window)
% For each trial or each filenum, computes a Morlet wavelet scalogram, and
% returns their average.  The calculated scalogram is trimmed at the
% beginning and end by a time interval (<margintime>) sufficient to avoid
% large (~30 dB below real signal) windowing artifacts, and the analysis
% window is lengthened by the same amount so that the final result shows
% the time interval specified by <window>.  Note that if recorded segments
% are not long enough to accommodate that lengthening, then the displayed
% result will be shorter than specified. When using the default value for
% 'param', <margintime> is 1.55 seconds / lfp_FreqLim(1)) at both ends.  If
% lfp_FreqLim is empty, the default is [5 1/(3*lfp_SamplePeriod)].
%INPUTS
% trials, filenums, window: as in lfp_getSamples, except that only ONE of
%   <trials> or <filenums> may be a non-scalar array.  Computes a Morlet
%   wavelet scalogram for each trial or file.  Scalograms are averaged over
%   all trials or files.
%OUTPUTS
% morletdata: a struct suitable for use with lfp_plotMorletgram.  Contains
%       the following fields, the last batch of which contain the values of
%       the like-named local options variables:
%   P - the pseudo-power matrix in scales X samples format.
%   ntrigs - number of triggers that contributed to the matrix.
%   trials - the list of trials actually contributing to the result (i.e.
%       excluding trials that lacked a reference event)
%   trialslabel - as computed when lfp_TrialStyle = 'rule'.
%   timepts - for the P matrix, in seconds.
%   freqs - the pseudo-frequency at each scale (see dg_morletgram).
%   win - the time window relative to lfp_AlignmentRef that was analyzed.
%   filenames - cell array of string value(s) from lfp_FileNames for the
%       channel(s) that was(were) analyzed.
%   sessionnames - value of lfp_SessionNames.
%   align - value of lfp_AlignmentRef.
%   evtmatrix - {foobar: returned by lfp_getSpikes}.
%   trigfuncArgStr - a string representation of the 'trigfunc' function's
%       arguments (which could be arbitrarily large, and that would be
%       bad).
%   countsflag
%   rawflag
%   evtavg2flag
%   evtbounds
%   multitrigflag
%   norefOKflag
%   trigfuncflag
%   trigfuncH
% allP: 3D array of wave-by-wave scalograms, the average of which over the
%   third dimension is returned as <morletdata.P>.
%OPTIONS
% All lfp_getSamples options are available, plus:
%   'fspacing', fspacing - ratio between scales on successive rows, which
%       gets passed in to dg_morletgram.  Use a value of 10^(1/n) where <n>
%       is an integer if you want nice round frequency tick mark labels.
%       Default: 10^(1/64).
%   'lin' - uses a linear frequency scale instead of the default log scale.
%       Also changes default value of <fspacing> to 1 Hz; to override that
%       default, you must use the 'fspacing' option *after* the 'lin'
%       option in your argument list.
%   'norm', b - normalizes pseudo-power values to a pink noise spectrum
%       given by pseudo-power = a * f ^ (-b).  Example: "'norm', 1.5" would
%       normalize to (a / f^1.5).  <a> is fitted to make the average
%       difference over all frequencies in log power between the actual
%       pseudo-spectrum and the pink spectrum equal to zero.
%   'param', paramval - the Morlet wavelet parameter value, which gets
%       passed in to dg_morletgram.  This is approximately the number of
%       cycles in the wavelet. Default = 6.
%NOTES
%   To produce figures, see lfp_plotMorletgram.
%   If this function should ever acquire a 'trigfunc' option, it should
% really go in lfp_CSCboilerplate/lfp_getSamples.

%$Rev: 379 $
%$Date: 2016-03-30 15:21:20 -0400 (Wed, 30 Mar 2016) $
%$Author: dgibson $


lfp_SamplePeriod=1/samplingfreq;
fspacing = 10^(1/64);
paramval = 6;

argnum = 1;
normb = 0;
fscale = 'log';
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'fspacing'
                argnum = argnum + 1;
                fspacing = varargin{argnum};
            case 'lin'
                fscale = 'lin';
                fspacing = 1;
            case 'norm'
                argnum = argnum + 1;
                normb = varargin{argnum};
            case 'param'
                argnum = argnum + 1;
                paramval = varargin{argnum};
            otherwise
                error('lfp_morletgram:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('lfp_morletgram:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if isempty(freqlim)
    freqlim = [5 1/(3*lfp_SamplePeriod)];
end

window=[0 size(data,1)/samplingfreq];
FourierFactor = 4*pi/(paramval+sqrt(2+paramval^2));
margintime = 1.5 / (FourierFactor * freqlim(1));
marginpts = round(margintime/lfp_SamplePeriod);
window = window + [-margintime margintime];
%trials = setdiff(trials, badtrials);
% Put <sampledata> into samples X waves format:
%sampledata = squeeze(sampledata);
sampledata=data;%
%numwaves = size(sampledata,2);
%for k = 1:numwaves
    %[P, freqs] = dg_morletgram(sampledata(:,k), lfp_SamplePeriod, ...
    %    fspacing, freqlim, paramval, fscale);
        [P, freqs] = dg_morletgram(sampledata(:,1), lfp_SamplePeriod, ...
        fspacing, freqlim, paramval);
    allP = P(:, marginpts+1 : end-marginpts);
 %   if k == 1
%        allP = P(:, marginpts+1 : end-marginpts);
%    else
%        allP(:,:,k) = P(:, marginpts+1 : end-marginpts);
%    end
%end

morletdata.P = mean(allP, 3);
if normb ~= 0
    avgP = mean(morletdata.P, 2);
    logpinkspec = reshape(log(freqs .^ -normb), [], 1);
    logdiffspec = log(avgP) - logpinkspec;
    Pscale = exp(mean(logdiffspec));
    pinkspec = reshape(Pscale * (freqs .^ -normb), [], 1);
    morletdata.P = morletdata.P ./ repmat( ...
        pinkspec, 1, size(morletdata.P, 2) );
end
%morletdata.ntrigs = numwaves;
%morletdata.trials = trials;
%morletdata.trialslabel = lfp_getTrialsLabel(trials, 'rule');
%morletdata.timepts = timepts(marginpts+1 : end-marginpts);
morletdata.freqs = freqs;
morletdata.win = window;
%morletdata.filenames = lfp_FileNames(filenums);
%morletdata.sessionnames = lfp_SessionNames;
%morletdata.align = lfp_AlignmentRef;
morletdata.evtavg2flag = [];
morletdata.evtbounds = [];
morletdata.multitrigflag = [];
morletdata.norefOKflag = [];
morletdata.trigfuncflag = false;

