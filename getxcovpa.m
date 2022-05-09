function xvals=getxcovpa(xdata,tsx,targval,varargin)
%get xcov parameters
%xdata = xcov data with trials x xcov
%tsx = timepoints for xcov
%targval options = minlag, maxlag, mincoef,maxcoef...
%mincoefatlag (min coef at mean lag for min to constrain lag time points)
%xvals output targval values

argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case ''
    end
    argnum=argnum+1;
end

absxpeak=nan(1,size(xdata,1));
for ii=1:size(xdata,1)
    xpeak=nan;
    xpeakid=nan;
    if contains(targval,'min')
        [xpeak,xpeakid]=min(xdata(ii,:));
    else
        [xpeak,xpeakid]=max(xdata(ii,:));
    end
    if isnan(xpeak) || xpeak==0 || isempty(xpeakid) || isnan(xpeakid)
        absxpeak(ii)=nan;
    else
        if contains(targval,'lag')
            %get lag (ts)
            absxpeak(ii)=tsx(xpeakid);
        else
            %get coeff
            absxpeak(ii)=xdata(ii,xpeakid);
        end
    end
end

if contains(targval,'at')
    absxpeakorig=absxpeak;
    absxpeak=nan(1,size(xdata,1));
    %get values at mean min or mean max lag
    %absxpeak corresponds to lag values since 'lag' in targval
    meanlag=nanmean(absxpeakorig);
    difftsx=tsx-meanlag;
    [matchtsx,matchtsxid]=min(abs(difftsx));        %get closest  ts to rounded meanlag
    lagid=matchtsxid;
    for ii=1:size(xdata,1)
        if ~isnan(absxpeakorig)
            absxpeak(ii)=xdata(ii,lagid);
        else
            absxpeak(ii)=nan;
        end
    end
end
if contains(targval,'ard')
    absxpeakorig=absxpeak;
    absxpeak=nan(1,size(xdata,1));
    %get values around mean min or mean max lag
    %absxpeak corresponds to lag values since 'lag' in targval
    meanlag=nanmean(absxpeakorig);
    difftsx=tsx-meanlag;
    [matchtsx,matchtsxid]=min(abs(difftsx));        %get closest  ts to rounded meanlag
    lagid=matchtsxid;
    lagids=lagid-5:lagid+5;       %  +/- 3 samples (ie 600 ms window)
    lagids=lagids(lagids>0 & lagids<=length(tsx));
    for ii=1:size(xdata,1)
        if ~isnan(absxpeakorig)
            absxpeak(ii)=min(xdata(ii,lagids));
        else
            absxpeak(ii)=nan;
        end
    end
end

xvals=absxpeak;

end