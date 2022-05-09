function A = envwave(data, varargin)
%NOW TO REPLACE "smoothwin.m" functions used throughout lfp analysis
%1/2019, taken from DJG lfp_waveenv:
%A = lfp_waveenv(filenum)  Function intended for use with lfp_createWave.
% Constructs an envelope fitted to the peaks of the waveform obtained by
% rectifying the wave in <filenum>.
%OPTIONS
% 'mask', mask - <mask> is a logical vector with the same number of
%   elements as lfp_SamplesUnits{filenum}, or a logical array the same size
%   as lfp_SamplesUnits{filenum}.  A zero is included in the time series to
%   be splined at the last false value in mask before mask goes true, and
%   at the first false value where it returns to false.
% method - any of the methods acceptable to interp1 except for 'v5cubic'
%   can be invoked here, namely <method> = 'nearest', 'linear', 'spline',
%   'pchip', 'cubic'. Default is 'spline'.  'nearest' produces a
%   rectangular wave.  'spline' is the epitome of "good continuation", but
%   that means it overshoots.  'linear' fits the extrema exactly, but is
%   jagged.  'cubic' is a compromise between 'linear' and 'spline':  being
%   piecewise, it fits the extrema exactly without overshoot, and being
%   cubic gives it better continuation than 'linear', but you can still see
%   the discontinuities in curvature at the joints between segments.
% 'sign' - restricts fitting to those peaks where the underlying signal is
%   of the same sign as the peak, i.e. maxima must be above zero and minima
%   must be below zero.

%$Rev: 193 $
%$Date: 2011-01-12 19:03:58 -0500 (Wed, 12 Jan 2011) $
%$Author: dgibson $

mask = false(0,1);      %WHAT IS THIS FOR?
method = 'spline';
signflag = true;       %USE TRUE FOR SQ DAT
%USE SIGNFLAG = FALSE FOR JUST FILTERED DATA NOT YET SQUARED
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case {'nearest', 'linear', 'spline', 'pchip', 'cubic'}
            method = varargin{argnum};
        case 'mask'
            argnum = argnum + 1;
            mask = varargin{argnum};
        case 'nosign'
            signflag = false;
        otherwise
            error('lfp_waveenv:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if isempty(mask)
    wave = reshape(data(:), 1, []);
else
    zeros2addidx = find( ~mask(1:end-1) & mask(2:end) ...
        | ~mask(2:end) & mask(1:end-1) );
    wave = zeros(1, numel(data));
    wave(mask) = data(mask);
end
pospksamp = dg_findpks(wave);
negpksamp = dg_findpks(-wave);
if signflag
    pospksamp(wave(pospksamp) <=0 ) = [];
    negpksamp(wave(negpksamp) >=0 ) = [];
end
peaks = union(pospksamp, negpksamp);
if isempty(peaks)
    warning('lfp_waveenv:nopeaks', ...
        'There are no extrema in the waveform');
    A=nan;
    return;
end
if abs(wave(1)) > abs(wave(peaks(1)))
    startval = abs(wave(1));
else
    startval = abs(wave(peaks(1)));
end
if abs(wave(end)) > abs(wave(peaks(end)))
    endval = abs(wave(end));
else
    endval = abs(wave(peaks(end)));
end
xvals = [0 peaks length(wave)+1];
yvals = [startval abs(wave(peaks)) endval];
if isempty(mask)
    A = interp1(xvals, yvals, 1:length(wave), method);
else
    if ~isempty(intersect(xvals, zeros2addidx))
        badvals = ismember(xvals, zeros2addidx);
        xvals(badvals) = [];
        yvals(badvals) = [];
    end
    vals = [ xvals' yvals'
        zeros2addidx' zeros(length(zeros2addidx),1) ];
    vals = sortrows(vals);
    A = interp1(vals(:,1), vals(:,2), 1:length(wave), method);
end


end
