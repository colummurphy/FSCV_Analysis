function [axa,trids,tbt_winData,tbt_selectevmatTS]=tr_raster(trlists,varargin)
%User provides channel and data varaible to plot 
% call for plotting FSCV DA ITI example: [axa,tridsP,dataP]=tr_raster(trlists,'da',2,'ttypes',{{'big','small'},{'post','break'}},'win',[-12,4],'event','display_fix');
% for plotting Nlx lick: [~,tridsL,dataL]=tr_raster(trlists,'nlx','lick','ttypes',{{'big','small'},{'post','big'}},'win',[-12,4],'event','display_fix');
%persistent trlists
%Plot raster plots for sel chs and options
%incllude output for timestamps for events relative to alignemnt event HNS
%06/12/2022

fontsize=10;
figpos=[50,50,700,900];
numplots=1;
[figsess,axa]=setupFig(figpos,numplots);
pathsave=[trlists.path.home 'trout' filesep];   %Save path
sessid=trlists.sess;    %Session id
if ~isfolder(pathsave)
    mkdir(pathsave);
end
%datavar='fscv'; %Processed variable (struct) to open in trlists.processed
signaltype='FSCV';   %Signal source type: FSCV or Nlx
sitename='';
bgwinSampSize=3;    %Size of window for averaging signal for background signal for FSCV
dach=[];
betach=[];
nlxname='pulse';
eventaln='display_target';
outcomeTS=30;   %30s default alignment (sample 301 for fscv, assuming first TS=0s) for all trials as done in extractionsession for FSCV/ephys
sorttrs=0;      %Default no sort, by trial order
sortrts=0;      %Sort by rts
sortwin=[0 4];
sortwin={'targ','outcome'};
win=[-2 6];     %Plot window
plotz=0;
ttypes={'big','left'};        %Condition types first arg must be main condition (big, small, targetbreak,fixbreak), the rest can be target side (left or right) {'big','left'} 
plotevents={'display_fix','display_target','start_target','reward_big','reward_small','break_target','break_fix'};%events to plot with markers if present in trial
pEvtSize=15;    %Marker size for events
pEvtColors=cool(length(plotevents)-1);  %Plot markers colors
cscale=[];  %Color scale
normscale=1;    %Nomralize scale
errortr=0;          %Flag don't plot error trials = 0, plot error = 1, plot only error =2
sorttype='';        %Sort type
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}   
        case 'da'
            argnum=argnum+1;
            dach=varargin{argnum};     %Da chs to plot
            signaltype='fscv';
            datavar='fscv';
        case 'nlx'
            argnum=argnum+1;
            nlxname=varargin{argnum};  %Nlx name (e.g. pulse, 'pl1-p5', etc.)
            signaltype='nlx';
            datavar='nlx';
        case 'cscale'
            argnum=argnum+1;
            cscale=varargin{argnum};    %Cscale for clim 
        case 'event'
            %Alignment event
            argnum=argnum+1;    
            eventaln=varargin{argnum};         %event period, eg. fix, targ, outcome       
        case 'sort'
            %sort by signal quantity, define ('avg' over win, 'peak' in)
            argnum=argnum+1;
            sorttype=varargin{argnum};      %'avg' or 'peak' or 'rt_fix' or rt_target
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

%Get Trial #'s for plotting from trlists.trlist based on supplied ttypes
%and signaltype
trids=[];
if strcmp(datavar,'fscv')
   trids=getTrialIDs(trlists,ttypes,signaltype,'dach',dach);
    %trids=getTrialIDs(trlists,ttypes,signaltype,'dach',dach,'ignoregood');%TESTING ONLY
   
else
    signaltype='nlx';
    trids=getTrialIDs(trlists,ttypes,signaltype);
 %   trids=getTrialIDs(trlists,ttypes,signaltype,'ignoregood');;%TESTING ONLY
end
   % trids=trids(~ismember(trids,[964 965 966 967]));;%TESTING ONLY

%Get plotting events and alignment of signals (TS's for each trial)
%Since multiple trials in each trial data, need to get the events closest
%and before trlist.ts (i.e. outcome ts)
eventcodes=trlists.param.eventcodes;
%evtcode=getEvtID(eventaln,eventcodes);     %Get numeric code for specified event name

%Get sort type variables, if beahvioral varaible call function to get/calc
%target behavior
origtrids=trids;    %Store orig unsorted trids;
if ~isempty(sorttype)
    behsig=getBehVar(trlists,sorttype,trids); %Get behavioral signal given trids and type of beh (sorttype)    
    %Sort trids based on behsig
    [sortedbehsig,sortid]=sort(behsig,'ascend');
    trids=trids(sortid);%Sorted trids;
end

%Create array matrix TS's w/ columns refer to plotevents and rows to trids
%Use FSCV timestamps (or FSCV sample ids, evmatfids) rather than NLx timestamps for FSCV signals as to not
%accumulate errors in rounding since already rounded during FSCV extract
%(probably not accumulating errors, but since there are 40-50 ms periods of
%time precision that can occur between target onset and reward, better to
%follow rounding done during original extraction than possibility of
%rounding off differently here and making another error
[evmatA,evmatfidsA]=getEvtMatAllNum(plotevents,trlist,eventcodes,'trialids',trids);     %Get event matrix for all trials, all events between cur outcome and prev outcome, 3d numeric
[evmat,evmatfids]=getEvtMat(plotevents,trlist,eventcodes,'trialids',trids);     %Get event matrix for all trials

%Get data trial by trial into a matrix array tbt_Data for the chosen trids
rate=10;    %Default fscv sample rate = 10 Hz
chnum=[];   %nlx ch
tbt_Data=[];    %trial by trial data for fscv or nlx
if strcmp(signaltype,'fscv')
    tbt_Data=vertcat(trlists.fscv{dach}(trids).da);
        site=trlists.fscvsites(find(ismember([trlists.fscvsites.ch],dach)));
        sitename=site.site;
else
    %nlx data
    chnum=find(contains(trlists.nlxnames,nlxname));%Find ch to get nlx data based on specified name
    tbt_Data=vertcat(trlists.nlx{chnum}{trids});%Convert cell to vertical double array, each row a trial
    rate=round(getSampleRate(trlist(trids(1)).NlxTS));
    if rem(rate,10)~=0
        %If not divisible by 10, something wrong with rate calculation...
        STOP;%Call debug function
        error('rate calculation incorrect; not divisible by 10');
    end
    site=trlists.nlxsites(chnum);
    sitename=site.site;
end

%Define plotting window tbt_winSamps in samples (based on signal rate), define uniform
%TSs for all selected trid trials. 
% Also BG subtract FSCV data to aln evt
%   Also get TS's (mapped to TS_plot and relative to eventaln) for all plotevents
TS_plot=win(1):1/rate:win(2);%TS to plot in seconds, relative to alignment event
winIDsRel=win(1)*rate:win(2)*rate;%Relative samples for plotting window relative to aln evt
alnID=find(contains(plotevents,eventaln));  %ID for alignment event
tbt_winSamps=[];%Sample points for each trial of data to plot for user defined window
tbt_evmatTSplot=[]; %Event time stamps relative to eventaln and TS_plot for plotting markers
%[tbt_winSamps,tbt_evmatTSplot]=getEvtAlnTBT()
if strcmp(signaltype,'fscv')
    alnFsamps=evmatfids(:,alnID);    %FSCV ids for alignment event for each trial for trids
    if length(alnID)>1
        %2 or more event IDs specified, choose the numeric or latest (max) one
        alnFsamps2=max(alnFsamps,[],2);
        alnFsamps=alnFsamps2;
    end
    tbt_winFsamps=repmat(winIDsRel,size(evmatfids,1),1)+alnFsamps;%Sample ids for each trial and window
    tbt_winSamps=tbt_winFsamps;
    %Get re-aligned TS's for marker plotting of events
    evmatS_subaln=evmatfidsA-alnFsamps;  %Event TS samples - alignment TS samples, 3d array num
    evmatTS_subaln=evmatS_subaln./rate;  %Event TS for all plotevents in TS_plot domain
    tbt_evmatTSplot=evmatTS_subaln;
    %Background subtract FSCV data
    for it=1:size(tbt_Data,1)
        bgwinids=alnFsamps(it)-bgwinSampSize:alnFsamps(it)-1;
        bgwinids(bgwinids<1)=[];
        bgwinids(bgwinids>size(tbt_Data,2))=[];
        if isempty(bgwinids)
            STOP;%Call debug function
            error('bg ids beyond boundaries of data');
        end
        bgmean=nanmean(tbt_Data(it,bgwinids));        %Mean background signal around aln evt
        tbt_Data(it,:)=tbt_Data(it,:)-bgmean;   %Subtract mean background signal
    end
else
    %Nlx data
    alnTSNlx=evmat(:,alnID);   %NLX ts's for aln event for each trial for trids
    if length(alnID)>1
        %2 or more event IDs specified, choose the numeric or latest (max) one
        alnTSNlx2=max(alnTSNlx,[],2);
        alnTSNlx=alnTSNlx2;
    end
    nlx_TSs=vertcat(trlist(trids).NlxTS);%All Nlx TS's for each selected trid trial
    tbt_alnNlxSamps=findNearestTS(alnTSNlx,nlx_TSs,1/rate/2);%Find the nlx sample for each row of data that most closely matches provided alignment event  TS (Note difference comes from downsampling), 50% margin of error
    tbt_winNlxSamps=repmat(winIDsRel,size(evmat,1),1)+tbt_alnNlxSamps;%Sample ids for each trial and window for Nlx
    tbt_winSamps=tbt_winNlxSamps;
    %Get re-aligned TS's for marker plotting of events
    evmatTS_subaln=evmatA-alnTSNlx;  %Event TS - alignment TS, 3d array numeric
    tbt_evmatTSplot=evmatTS_subaln;     %Event TS for all plotevents in TS_plot domain
end


%Process data if phys (e.g. lick, pulse, pupil), before windowing
if strcmp(signaltype,'nlx') && contains(nlxname,'lick')
    %Lick low pass filter (as used at MIT) trial by trial, really does very
    %little...may be remove.
    ldata=tbt_Data; %Save original data;
    for it=1:size(tbt_Data,1)
        tbt_Data(it,:)=filterLFP(tbt_Data(it,:),rate,[0 100]);     %Low pass filter fc = 100 Hz for lick  --> Does very little
        %datemp=smoothwin(datemp,0.1);   %Smoothing with Hanning function--> NOt sure how this ever worked for setTrialAx
    end
end

if strcmp(signaltype,'nlx') && contains(nlxname,'eyed')
    %Pupil diameter, in original setTrialAx, appears that eyed was
    %inverted..not done here yet.
    %NaN any blinks in data otherwise artificially can increase diameter
    ed_data=tbt_Data;
    blink_thres=2.5e-3;%Absolute peak threshold to remove eye blinks (usually hits rail)
    blink_pad=30;%Samples for padding around blink removal
    for it=1:size(tbt_Data,1)
       tbt_Data(it,:)=deglitchnanamp(tbt_Data(it,:),blink_thres,30);
    end
end

if strcmp(signaltype,'nlx') && contains(nlxname,'eyex')
       %Get pupil velocity from eye x traces (eye y not valid since task just
    %calibrated for left/right mvmts in Cleo/Patra)
    %eyedist=[];
   % eyeusacs=[];
    eyev=[];
    blink_pad=30;
    smoothlength=20;
    eyex=tbt_Data;%Store original data;
    meaneye=nanmean(tbt_Data,2);
    meane=nanmean(meaneye);
    stdeye=nanstd(tbt_Data,[],2);
    meanstd=nanmean(stdeye);
    threseye=abs(meane)+3*meanstd;%THreshold to remove blinks (not sure why different from eyed)
    for itrial=1:size(tbt_Data,1)
        %REMOVE BLINKS
        tbt_Data(itrial,:)=deglitchnanamp(tbt_Data(itrial,:),threseye,30); %Produce nan at blinks and large/fast transients 
        smootheye=smoothwin(tbt_Data(itrial,:),smoothlength);
        eyevel=diff(smootheye);   %Get eye velocity from differentiating x so = x/samps
       % eyevel(isnan(eyevel))=0;
        %abseyevel=abs(eyevel);
        %eyedisttemp=cumtrapz(abseyevel);%Get distance from integrating eye velocity, WHY?? JUst use x
       % thres=nanmean(abseyevel);
     %   maxlim=nanstd(abseyevel);
      %  [pks,locs,w,p] = findpeaks(abseyevel,'minpeakwidth',...
      %      round(samplespersec*.01),'minpeakdistance',round(samplespersec*.05)...
      %      ,'MinPeakprominence',thres);
      %  sacslogic=zeros(1,length(datatemp));
      %  locsbelowlim=locs(pks<=maxlim);
       % sacslogic(locsbelowlim)=1;
       % eyeusacs(itrial,:)=sacslogic;
       % eyedist(itrial,:)=eyedisttemp;
        eyev(itrial,:)=eyevel;
    end
    tbt_Data=eyex;
    %HAVE TO ADD NAN TO LAST COLUMN TO MAINTAIN SAME AMOUNT OF DATA POINTS!
    
end


%Windowed plotting data just in defined window from original tbt_Data
tbt_winData=zeros(size(tbt_winSamps));
for it=1:size(tbt_Data,1)
    tbt_winData(it,:)=tbt_Data(it,tbt_winSamps(it,:));
end




%Plot raster trial by trial data tbt_winData and format
%Display everything except aln evt
marka=[];%Event names excluding aln evt
tbt_selectevmatTS=[];%tbt_evmatTSplot excluding aln evt
if ~isempty(plotevents)
    marka=plotevents(~contains(plotevents,eventaln));       %Marker names except current aln evt
    tbt_selectevmatTS=tbt_evmatTSplot(:,find(~contains(plotevents,eventaln)),:);
end
marka=plotevents; tbt_selectevmatTS=tbt_evmatTSplot; %DISPLAY EVERYTHING INCLUDING ALN EVT
alltypes=[ttypes{:}];
plotname=[signaltype ' | Site ' sitename ' | Session ' sessid char(10)...
    'Conditions: ' [alltypes{:}] ' | 0 Aln: ' eventaln ' ' char(10)...
    'Sort: ' sorttype];
%Raster plotting function call
if ~isempty(tbt_winData)
axa{1}=plotTBTRaster(TS_plot,tbt_winData,'axes',axa{1},'title',plotname,...
    'cscale',cscale,'markers',{marka,tbt_selectevmatTS});   
end

