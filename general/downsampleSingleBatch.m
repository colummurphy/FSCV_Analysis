function subsamples=downsampleSingleBatch(rawsamples,N)
%taken from dg_downsample that Dan created
        % Entire file is already in memory; downsample it.
        % compute sample period <sampleT>:
        fc = 0.9/N;
        [z, p, k] = butter(4, fc, 'low');
        [sos,g]=zp2sos(z,p,k);
        h2=dfilt.df2sos(sos,g);
        % Lowpass filter <rawsamples> in place.
        % Keeps shape of <rawsamples>.
        samplesize = size(rawsamples);
        rawsamples = filter(h2, rawsamples(end:-1:1));
        rawsamples = reshape(filter(h2, rawsamples(end:-1:1)), samplesize);
        subsamples = reshape(rawsamples(1:N:end), 1, []);
       % allsubTS = reshape(allrawTS(1:N:end), 1, []);
end