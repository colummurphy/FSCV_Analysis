function xingidx = dg_fitFSCVart(ts, samples, thresh)
%INPUTS
% ts: timestamps of all samples, same size as <samples>.
% samples: array of sampled values, one per timestamp.
% thresh: positive-going threshold to trigger interpolation.
%OUTPUTS
% xingidx: the index into <samples> of the first element above thresh in
%   each identified FSCV artifact.
%NOTES
%   The fitting of expected artifact times to actual artifact times is done
% cycle by cycle to allow for plus or minus 2% drift over the course of
% each cycle.  The first cycle is fitted by matching the initial 1/10th of
% threshold crossings against a perfectly periodic template.
%   The period of the FSCV artifact is assumed to be at least 50 ms.
%   A simple linear interpolator is provided in dg_rmFSCVart.

%$Rev: 269 $
%$Date: 2019-10-07 17:40:04 -0400 (Mon, 07 Oct 2019) $
%$Author: dgibson $

tol = 0.02; % fractional tolerance allowed for drift per cycle
sampleperiod = median(diff(ts));
minsamps = floor(0.05 / sampleperiod); % minimum conceivable FSCV period

% Find threshold crossings ("xings"), period, and start of first period.
% <xingidx> points to the first sample that is above <thresh> in each
% putative FSCV artifact.
isxing = [ false
    reshape( samples(2:end) >= thresh & samples(1:end-1) < thresh, ...
    [], 1 ) ];
xingidx = reshape(find(isxing), [], 1 );
xingdiffs = diff(xingidx);
period = round(median(xingdiffs(xingdiffs >= minsamps)));
tolpts = round(tol * period);
% Find the best alignment of a perfectly periodic template to the actual
% xings. Template should be 10 cycles long:
templen = 10 * period;
template = false(templen, 1);
template(1:period:end) = true;
fracmatch = zeros(period, 1);
% The first offset should -<tolpts>, or as close as we can get without
% going off the beginning of <isxing>:
firstoffset = max(-tolpts, 2 - xingidx(1));
offsets = firstoffset + (0 : period - 1);
for k = 1:length(offsets)
    startmatch = xingidx(1) + offsets(k);
    if startmatch + length(template) - 1 > length(isxing)
        error('fitFSCVart:length', ...
            'The data series is too short to accommodate the timing template.');
    end
    fracmatch(k) = sum( template ...
        & isxing(startmatch : startmatch + length(template) - 1) ) ...
        / sum(template);
end
if ~any(fracmatch > 0.1)
    error('fitFSCVart:fracmatch', ...
        'Beginning tenth of waveform is too noisy to fit.');
end
[~, ix] = max(fracmatch);
startT1idx = xingidx(1) + offsets(ix); % index into <samples>.
[~, bestidx2] = min(abs(xingidx - startT1idx));
% <refidx> always points at the last confirmed periodic threshold crossing
% in <samples>; it is an index into <samples>:
refidx = xingidx(bestidx2);
% "idx2" is meant to denote an index into an index, i.e. <xingidx2> is an
% index into <xingidx>.  <xingidx2> always points at the next xing in
% <xingidx> to be either confirmed or rejected:
xingidx2 = bestidx2 + 1;

%
% Eliminate threshold crossings that are off by more than <tolpts> from the
% expected index. 
%
isbadxing = false(size(xingidx)); % will get set <true> for bad xings
isbadxing(1 : bestidx2 - 1) = true;
while refidx <= numel(samples)
    % <targidx> is where we expect to find the next valid xing:
    targidx = refidx + period;
    %
    % Sort xings into good and bad until the first one is bad because it's
    % too late.
    %
    if xingidx2 > numel(xingidx)
        break
    end
    numgood = 0;
    % Mark the xings bad if they're too early:
    while xingidx2 <= numel(xingidx) && xingidx(xingidx2) < targidx - tolpts
        isbadxing(xingidx2) = true;
        xingidx2 = xingidx2 + 1;
    end
    if xingidx2 > numel(xingidx)
        break
    end
    % Count as good xings if they're not too late:
    while xingidx2 <= numel(xingidx) && xingidx(xingidx2) <= targidx + tolpts
        numgood = numgood + 1;
        xingidx2 = xingidx2 + 1;
    end
    if xingidx2 > numel(xingidx)
        break
    end
    % <xingidx2> now points at the first xing after the ones that match
    % <targidx>. We need to get rid of the extra "good" xings if there were
    % more than one, and update <refidx>, which is now stale:
    if numgood > 1
        goodxings = xingidx(xingidx2 - (1:numgood));
        [~, ix] = min(abs(goodxings - targidx));
        bogusxingsidx2 = xingidx2 - setdiff(1:numgood, ix);
        isbadxing(bogusxingsidx2) = true;
        refidx = xingidx(xingidx2 - ix);
    elseif numgood == 1
        refidx = xingidx(xingidx2 - 1);
    else
        % <numgood>  must be zero, indicating that a xing is missing.  We
        % will simply pretend that we found it, but note that as the number
        % of missing artifacts increases, it becomes increasingly likely
        % that the next artifact will be off by more than 2%, in which case
        % the entire remainder of the file will fail to match.
        refidx = targidx;
    end
    % <refidx> now points at the sample identified as the last confirmed
    % xing in <xingidx>. 
end
xingidx(isbadxing) = [];

% Insert estimates for any threshold crossings that are missing.  Missing
% xings are identified on the basis of successive xings that are
% approximately two or more periods apart.
missingxingidx = [];
ptsbetween = diff(xingidx);
% <longidx2> points into <xingidx> at the xing before a long inter-xing
% interval.
longidx2 = find(ptsbetween > period * (1 + 2 * tol)); 
for longidx3 = 1:length(longidx2)
    numperiods = round(ptsbetween(longidx2(longidx3)) / period);
    ptsperxing = ptsbetween(longidx2(longidx3)) / numperiods; % float!
    for k = 1 : numperiods - 1
        missingxingidx(end+1, 1) = xingidx(longidx2(longidx3)) + ...
            round(k * ptsperxing); %#ok<AGROW>
    end
end
xingidx = sort([xingidx; missingxingidx]);