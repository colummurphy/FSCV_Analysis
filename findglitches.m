function glitchids=findglitches(rawdata,glitchThres,varargin)
%pass raw data to accurately find glitches (HF spikes)
glitchids=[];
argnum=1;
smoothwidth=4;
bgmean=mean(rawdata,2);
refmatrix=repmat(bgmean,1,size(rawdata,2));   %tile BGavg vector for entire meas matrix span
rawdata=rawdata-refmatrix;      %raw data is subtracted
saturationthres=8;     %add 07/04/2018 absolute current change threshold
satwidth=[-5 30];    %atleast 2 seconds needed for saturated current to restablilize
saturationthres=[];
if isempty(saturationthres)
    saturationthres=mean(mean(abs(refmatrix),2)).*.015;
end
rodentflag=0;       %long recording
while argnum <= length(varargin)
    switch varargin{argnum}
       case 'glitchwidth'
            %width of nanned signal padding
            argnum=argnum+1;
            smoothwidth=varargin{argnum};
            case 'satwidth'
            %width of large saturation current -/+ samples
            argnum=argnum+1;
            satwidth=varargin{argnum};
            case 'satthres'
            %thres large saturation current 
            argnum=argnum+1;
            saturationthres=varargin{argnum};
        case 'rodent'
            rodentflag=1;
    end
       argnum = argnum + 1;
end
%get differential along cvs (I vs V that is) 
diffrawcvs=diff(rawdata,1,1);   
%maxdiff=max(abs(diffrawcvs),[],1);
    maxdiff=max(abs(diffrawcvs),[],1)-min(max(abs(diffrawcvs),[],1));   %fixed 10/2/2018
tbb=find(maxdiff>glitchThres);
if rodentflag
tbb=find(maxdiff>glitchThres+2);
end
idxglitch=tbb;
idxglitch=unique(idxglitch);  %remove duplicates
%get raw signal along time for v-selected
diffIrawt=diff(rawdata,2,1); %take out hf glitches along time by looking at unfiltered data
%maxdifft=max(abs(diffIrawt),[],1);
maxdifft=max(abs(diffIrawt),[],1)-mean(max(abs(diffIrawt),[],1));  %fixed 10/2/2018
%get differential of this signal to detect HF glitches
rawglitch=find(maxdifft>glitchThres);
if ~isempty(idxglitch)
    idxglitch=unique([idxglitch rawglitch]);
else
    idxglitch=rawglitch;
end
%for each detected sample add +/- smoothwidth samples around
for ii=1:length(idxglitch)
    idsaroundglitch=idxglitch(ii)-smoothwidth:idxglitch(ii)+smoothwidth;
    %07/04/2018 very large current changes (ie large movements)
    %want to take out signal more widely, more of forward current 
    if abs(mean(rawdata(:,idxglitch(ii))))>saturationthres
            idsaroundglitch=idxglitch(ii)+satwidth(1):idxglitch(ii)+satwidth(2);           
    end
    if ii==1
        glitchids=idsaroundglitch;
    else
        glitchids=[glitchids idsaroundglitch];
    end
end
glitchids=sort(unique(glitchids));    
glitchids=glitchids(glitchids>0 & ...
    glitchids<size(rawdata,2)); %only keep ids in original window
end

