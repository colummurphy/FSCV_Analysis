function tr_raster()
persistent trlists
%Plot raster plots for sel chs and options

pathsave=[trlists.path.home 'trout' filesep];   %Save path
if ~isfolder(pathsave)
    mkdir(pathsave);
end
datavar='fscv'; %Processed variable (struct) to open in trlists.processed
dach=[];
betach=[];
eventaln='display_target';
outcomeTS=30;   %30s default alignment (sample 301 for fscv, assuming first TS=0s) for all trials as done in extractionsession for FSCV/ephys
sorttrs=0;      %Default no sort, by trial order
sortrts=0;      %Sort by rts
sortwin=[0 4];
sortwin={'targ','outcome'};
win=[-2 6];     %Plot window
plotz=0;
ttypes={'big','left'};        %condition types grouped eg. {{'big','left'},{'small','left','aftersm'}}, first arg must be "big" small, targetbreak or fixbreak type
plotevents={'display_fix','display_target','reward_big','reward_small','break_target','break_fix'};%events to plot with markers if present in trial
errortr=0;          %Flag don't plot error trials = 0, plot error = 1, plot only error =2
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}   
        case 'processed'
            %Get "processed" data in trlists (specific field) (e.g. pulse)
            argnum=argnum+1;
            datavar=varargin{argnum};   %Field to open in trlists.processed
        case 'da'
            argnum=argnum+1;
            dach=varargin{argnum};     %Da chs to plot
        case 'beta'
            %Beta to plot (must be processed var in trlists)
            argnum=argnum+1;
            betach=varargin{argnum};   
        case 'event'
            %Alignment event
            argnum=argnum+1;    
            eventaln=varargin{argnum};         %event period, eg. fix, targ, outcome       
        case 'sort'
            %sort by signal quantity, define ('avg' over win, 'peak' in)
            sorttrs=1;
            argnum=argnum+1;
            sorttype=varargin{argnum};      %'avg' or 'peak'
            if contains(sorttype,'rts')
                sortrts=1;
            end
        case 'sortwin'
            %Only if sort selected, window to define avg or peak
            argnum=argnum+1;
            sortwin=varargin{argnum};       %Can be numeric or 2 events
        case 'win'
            %Plot window +/- s
            argnum=argnum+1;
            win=varargin{argnum};
        case 'plotz'
            plotz=1;           
        case 'plotmarks'
            plotmarks=1;            %Enable plot of marks
            argnum=argnum+1;
            marks=varargin{argnum};     %Marks to plot for each trial (TS)        
        case 'ttypes'
            argnum=argnum+1;
            ttypes=varargin{argnum};        %condition types grouped eg. {{'big','left'},{'small','left','aftersm'}}
        case 'ploterror'
            %Don't plot error trials = 0, plot error= 1, plot only error=2
            argnum=argnum+1;
            errortr=varargin{argnum};
        case 'plotmax'
            plotmax=1;
        case 'plotmin'
            plotmin=1;
        case 'plotmean'
            plotmean=1;
        case 'plotimmean'
            plotimmean=1;
        case 'plotabs'
            plotabs=1;
        case 'pad'
            argnum=argnum+1;
            pad=varargin{argnum};
        case 'smoothlfp'
            smoothlfp=1;
            argnum=argnum+1;
            smoothpoints=varargin{argnum};
        case 'winb'
            %win for beh
            argnum=argnum+1;
            winb=varargin{argnum};
        case 'rrstd'
            plotrrstd=1;
    end
    argnum=argnum+1;
end

trlist=trlists.trlist;

%Get Trial #'s for plotting from trlists.trlist
trids=find(contains({trlist.type},ttypes{1}));
%Open good/bad trials from FSCV, if FSCV data, then sel the act channel
goodtrs=[];
if strcmp(datavar,'fscv')
    goodtrs=find([trlists.fscv{dach}.good]);
    trids=intersect(trids,goodtrs);
else
    %Any other data variable take good trials from fscv (or of all chs)
    %Get union of all chs
    goodtrslist=[];
    for ii=1:length(trlists.fscv)
        goodtrslist(ii,:)=[trlists.fscv{ii}.good];
    end
    goodtrs2=sum(goodtrslist);
    goodtrs=find(goodtrs2~=0);
    trids=intersect(trids,goodtrs);
end


%Get plotting events and alignment of signals (TS's for each trial)
%Since multiple trials in each trial data, need to get the events closest
%and before trlist.ts (i.e. outcome ts)
eventcodes=trlists.param.eventcodes;
evtcode=getEvtID(eventaln,eventcodes);     %Get numeric code for event name

%Create array matrix TS's w/ columns refer to plotevents and rows to trids
evmat=getEvtMat(plotevents,trlist,eventcodes,'trialids',trids);     %Get event matrix for all trials


%Define plotting window in samples (based on sample rate), define uniform
%time stamps for all trials



 %more than one time point
        imagetrials=image(axa{ip}, pdata,'cdatamapping','scaled');
        set(axa{ip},'YDir','reverse')        %flip y-axis values so first trial on top
        artTime=isnan(pdata);   %find artifact points (nan periods)
        artTime=abs(artTime-1);         %make alpha data mask by inverting 1/0's
        artTime2=artTime;
        maskGray=artTime2==0;             %find Zero indices representing artifact mask
        maskGray=maskGray*.15;            %make gray rather than white default by making non-zero
        artTime=artTime+maskGray;
        set(imagetrials, 'AlphaData', artTime);    