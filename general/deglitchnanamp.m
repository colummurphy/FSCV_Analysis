function samples = deglitchnanamp(samples, thresh,padpts)
%faster than Dan's superdeglitch
%deglitch based on absolute signal not differential 06/24/2018
%nan signals in this period
%take differential of signal, and remove signals when the differential is
%high
%samples = dg_superdeglitch(samples, thresh, maxpts)
%HNS 03/13/2018
%Two thresholds needed, one for differential and one for actual signal
%to cut off high amplitude FSCV artifacts & diff signal for spike
%transitions
%Thresh input is an array of 2 numbers [thres_ref thres_samples]
%changed from using differential of data, to absolute threshold
%INPUTS
%ref_samples is the reference samples to take artifact time points from
%(ie. high pass filtered data)
% samples: sample data of any size or shape to interpolate down
% thresh: minimum departure from previous value that qualifies as a glitch.
% maxpts: the maximum number of points that can qualify as a glitch and
%   thus get interpolated away.
%OUTPUT
%NOTES
%  Returns the same sample data that was passed in, modified by linearly
%  interpolating glitches.  A k-point glitch is defined as a series of k
%  points that differ from a preceding point "a" by more than <thresh>, but
%  which is followed by a point "b" that differs from the predecessor point
%  "a" by less than or equal to <thresh>. It is explicitly NOT required
%  that the k members of the k-point glitch be within <thresh> of each
%  other, only that they be more than <thresh> from "a".  This makes it
%  possible to remove random noise bursts.  The interpolation extends from
%  the last point before the glitch through (including) the first point
%  after the glitch, so all glitch points get replaced by the
%  interpolation.
%   NaN values do not qualify as glitches, but they also don't
%  qualify as endpoints of a sequence to interpolate, and will prevent any
%  glitches that contain them from being removed.
%   This version has been completely rewritten to use less memory.
%  However, it does still have some memory issues, to wit:  the amount of
%  additional memory it requires is about twice the size of <samples>,
%  which means that the amount of memory allocated to the Matlab running it
%  will be at least three times the size of <samples> plus the size of a
%  freshly started Matlab (which is usually around 2 GB).

%$Rev: 214 $ $Date: 2015-03-26 00:09:24 -0400 (Thu, 26 Mar 2015) $ $Author:
%dgibson $

origsize = size(samples);
samples = reshape(samples,[],1);
difsamp = abs(diff(samples));   %ie differential of data to look for change in slope
abssamp=abs(samples);
thresh2=mean(difsamp)*7;     %arbitrarily set threshold for given signal based on mean differential

detectedtrigs=find(abssamp>thresh); 
%detectedtrigs2=find(difsamp>thresh2); 

%allTrigs=sort([detectedtrigs; detectedtrigs2]);     %merge timestamps found from both thresholds
allTrigs=detectedtrigs;     %merge timestamps found from both thresholds

allTrigs2=diff(allTrigs);       %detect when break in continuous timestamps
uniqueArtifacts=find(allTrigs2~=1); %unique event id's of allTrigsUnique
winNeg=1;   winPos=1;

for idx=1:length(uniqueArtifacts)
    if idx==1
        winNeg=allTrigs(1)-padpts;
        winPos=allTrigs(uniqueArtifacts(1))+padpts;
        %indicate window limits for interpolation
    else
        %initialize window from time point of unique artifact trigger
        idPos=allTrigs(uniqueArtifacts(idx));
        idNeg=allTrigs(uniqueArtifacts(idx));
        idxLocal=uniqueArtifacts(idx)+1;
        %find discontinuity in positive direction for window end
        while idxLocal<=length(abssamp)
            idPos=idPos+1;
            if allTrigs(idxLocal)~=idPos
                break
            else
            idxLocal=idxLocal+1;
            end
        end
        idxLocal=uniqueArtifacts(idx)-1;
        while idxLocal>=1
            idNeg=idNeg-1;
            if allTrigs(idxLocal)~=idNeg
                break
            else
            idxLocal=idxLocal-1;
            end
        end
        winNeg=idNeg-padpts;
        winPos=idPos+padpts;
        if winPos<=(length(abssamp)) && winNeg >=1
           % samples(winNeg:winPos) = linspace(samples(winNeg), ...
            %                    samples(winPos), winPos - winNeg + 1);
             samples(winNeg:winPos) =nan;                          
        end
    end
end

samples = reshape(samples, origsize);

