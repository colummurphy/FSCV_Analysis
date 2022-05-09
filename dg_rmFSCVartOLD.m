function samples = dg_rmFSCVart(samples, artidx, nbefore, nafter)
%INPUTS
% samples: array of samples, one per timestamp.
% artidx: index into samples of all samples that have been identified with
%   an instance of an FSCV artifact, e.g. validated threshold crossing
%   samples.
% nbefore: number of samples to interpolate before the first sample that is
%   at or above threshold.
% nafter: number of samples to interpolate after the first sample that is
%   at or above threshold.
%OUTPUTS
% samples: a copy of the input <samples> with FSCV artifacts interpolate
%   away.
%NOTES
% This is a trivial addendum to dg_fitFSCVart.

%$Rev: 268 $
%$Date: 2019-10-04 18:45:46 -0400 (Fri, 04 Oct 2019) $
%$Author: dgibson $

for artnum = 1:length(artidx)
    idx = (-nbefore:nafter) + artidx(artnum);
    samples(idx) = linspace( ...
        samples(idx(1)), samples(idx(end)), length(idx) );
end

    