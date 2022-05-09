function [lfpdata, mdl]=multiregress(sessnum,trlists,datms,lfptms,binfos,plotparam,varargin)
%call from corsmultiple
%4/2020, output relevant parameters for tabulating
%supply xinfo for manual lfp win's and trial #'s comparison
%noticed differences in R, between lfptm targimwin & 0-1 win xinfo maybe
%because more trials in lfptm than xinfo
outdata={};
sessiontypes={'bigreward','smallreward'};
condtype='all';
sessid=num2str(sessnum);
event='targ';
trlistid=find(strcmp({trlists.sessid},sessid));
trlist=trlists(trlistid).list;
numtrials=length(trlist);
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,'multiregress',filesep,sessid);
if ~isdir(savepath)
mkdir(savepath);
end
sessnums=plotparam.sessnums;
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites);
uniquesites=unique({sites(1:end).site});
[dapair,lfppair]=getsitepairs(targdasites);
targses=find(strcmp({trialgrps.sessid},sessid));
trialinfo=trialgrps(targses).trialinfo;

%initialize trlist
for ix=1:length(trlist)
    trlist(ix).da=nan;
    trlist(ix).lfp=nan;   
     trlist(ix).hr=nan;   
          trlist(ix).eye=nan; 
          trlist(ix).eyev=nan;   
           trlist(ix).eyedist=nan;  
     trlist(ix).lick=nan;   
     trlist(ix).rr=nan;   
     trlist(ix).rrstd=nan;   
     trlist(ix).rmssd=nan;   
end

figpos=[50,50,1400,900];
xoff=100;
yoff=100;
mar=25;
axsize=[400,250];
fontsize=15;
axpos={};
axa={};
triallim=[];
xdata=1:numtrials;
bigtype=[];
smalltype=[];
targtype=[];
fixtype=[];
marksize=[50,50,20,30];
mark={'.','.','o','x'};
plotvar={};
argnum=1;
siteda='';
sitelfp='';
bflag=0;
plotrt=0;
plothr=0;
ploteye=0;
linfos={};
xinfo={};
xbinfos={};
hfsametype=0;
metric='targpeak';
metriclfp='targwin';
win = [0 4];
offset=[0 0];
winb= [.1 4];
rate=10;            %sample rate for xinfo data
%metriclfp='targpeakabs';
plotnum=0;
plothist=[];
plottypes=[];
wine=[.2 .8];
winhr=[0 4];
winl=[0 2];
dataout={};
relfp=0;
noplot=0;
nodata=0;
xvarflag=0;
daonly=0;
lfponly=0;
negcortarg=0;
otherplotvars={'type','side','trialnum'};
multicor=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dax'
            %values from xinfo
            argnum=argnum+1;
            siteda=varargin{argnum};
            argnum=argnum+1;
            xinfo=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='da';
        case 'datm'
            %values from datm, supplied instead of xinfo
            argnum=argnum+1;
            siteda=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='da';
        case 'lfpx'
            argnum=argnum+1;
            sitelfp=varargin{argnum};            
            argnum=argnum+1;
            xinfo=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='lfp';
            relfp=1;
        case 'lfptm'
            %values from betatm
            argnum=argnum+1;
            sitelfp=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='lfp';
        case 'behx'
            %values from xbinfos
            argnum=argnum+1;
            xbinfos=varargin{argnum};
        case 'metric'
            %for datm values
            argnum=argnum+1;
            metric=varargin{argnum};
        case 'metriclfp'
            %for betatm values
            argnum=argnum+1;
            metriclfp=varargin{argnum};
        case 'trialseq'
            plotnum=plotnum+1;
            plotvars{plotnum}='trialnum';
        case 'rt'
            plotrt=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='trt';
        case 'frt'
            plotrt=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='frt';
        case 'hr'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='hr';
        case 'rr'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rr';
        case 'rrstd'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rrstd';
        case 'rmssd'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rmssd';
        case 'lick'
            plotnum=plotnum+1;
            plotvars{plotnum}='lick';
        case 'eye'
            plotnum=plotnum+1;
            plotvars{plotnum}='eye';
        case 'regressvars'
            %struct of variables to include as addional regressors, already
            %in trlist, e.g. 'type', 'side', 'trialnum'
            argnum=argnum+1;
            otherplotvars=varargin{argnum};
        case 'types'
        case 'eyev'
            %max eye velocity within ewin (1 s)            
            plotnum=plotnum+1;
            plotvars{plotnum}='eyev';
        case 'eyedist'
            %max of integral of abs eye velocity within ewin (1 s)            
            plotnum=plotnum+1;
            plotvars{plotnum}='eyedist';
        case 'event'
            argnum=argnum+1;
            event=varargin{argnum};    
        case 'hf'
            argnum=argnum+1;
            plothist=varargin{argnum};
        case 'hfsametype'
            hfsametype=1;       %constrict types of futre/past trials
        case 'types'
            argnum=argnum+1;
            plottypes=varargin{argnum};
        case 'condition'
            argnum=argnum+1;
            condtype=varargin{argnum};
        case 'bwin'
            argnum=argnum+1;
            winb=varargin{argnum};
        case 'win'
            argnum=argnum+1;
            win=varargin{argnum};
            metriclfp=['win' num2str(win(1)) '_' num2str(win(2))];
        case 'offset'
            argnum=argnum+1;
            offset=varargin{argnum};
            metriclfp=['win' num2str(win(1)+offset(1)) '_' num2str(win(2)+offset(2))];
        case 'lwin'
            argnum=argnum+1;
            winl=varargin{argnum};
        case 'hrwin'
            argnum=argnum+1;
            winhr=varargin{argnum};
        case 'ewin'
            argnum=argnum+1;
            wine=varargin{argnum};
        case 'noplot'
            noplot=1;
        case 'xinfo'
            %supply xinfo for trial #'s comparison
            argnum=argnum+1;
            xinfo=varargin{argnum};
        case 'nodata'
            %don't save data in struct just # trials
            nodata=1;
        case 'xvar'
            %use metriclfp variable to load variable in xinfos, ignore win
            %use with lfpx argument
            xvarflag=1;
            win=[0 0];
            offset=[0 0];
        case 'daonly'
            %only da and beh
            daonly=1;
        case 'lfponly'
            %only lfp and beh
            lfponly=1;
        case 'negcortarg'
            %select data if negative sig corr between DA vs LFP & perform
            %fitlm to get p values for each variable including behvs
            negcortarg=1;
        case 'multicor'
            %independent serial correlations of da vs beta, da vs behs & beta vs behs
            %then look at corr between residuals of da vs beta VS da vs
            %behs, same for beta, to explain whether correlation da vs beta
            %mediated by corr to beh
            %and do multiregress as well to see what variable explains
            %variance in beta best
            multicor=1;
    end
    argnum=argnum+1;
end   
if multicor
savepath=fullfile(plotparam.savepath,'multicortarg',filesep,sessid);
if ~isdir(savepath)
mkdir(savepath);
end
end
if any(win~=0) && relfp && ~xvarflag
    metriclfp=['win' num2str(win(1)+offset(1)) '_' num2str(win(2)+offset(2))];
end
%set up plots
if isempty(plothist)
plothist=zeros(1,length(plotvars));       %1 means 1 trial in future value
end
if isempty(plottypes)
plottypes=repmat({'all'},1,length(plotvars));
end

numplots=sum([length(plotvars)-1:-1:1]);        % # possible unique combinations correlations, for 3 variables, this is da vs beta, da vs trial, beta vs trials, num vars being correlated to single var + remaining vars after first iteration
labels=[];
titletext=['sess ' sessid ' | ' event ' | '  plottypes{1} ' | ' condtype ' | ' siteda ' ' metric ' | ' sitelfp ' ' metriclfp ];  
for ix=1:length(plotvars)
    if ix<length(plotvars)
    labels=[labels plotvars{ix} '_'];
    else
    labels=[labels plotvars{ix}];
    end        
end
labels=[labels '_' event '_' plottypes{1} '_' condtype ];
if ~isempty(siteda)
    labels=[labels '_' metric '_' siteda];
end
if ~isempty(sitelfp)
    if ~relfp
    labels=[labels '_' metriclfp '_' sitelfp];
    elseif relfp && any(win~=0) && ~xvarflag
        labels=[labels '_win_' num2str((win(1)+offset(1))*rate) '_' num2str((win(2)+offset(2))*rate) '_' sitelfp];
    elseif relfp && xvarflag
        labels=[labels '_' metriclfp '_' sitelfp];
    end  
end
if ~isempty(winhr)
    labels=[labels '_winhr_' num2str(winhr(1)*rate) '_' num2str(winhr(2)*rate)];
        titletext=[titletext ' | ' 'winhr' num2str(winhr(1)) '-' num2str(winhr(2)) 's' ];  
end
if ~isempty(wine)
    labels=[labels '_wine_' num2str(wine(1)*rate) '_' num2str(wine(2)*rate)];
    titletext=[titletext ' | ' 'wine' num2str(wine(1)) '-' num2str(wine(2)) 's' ];  
end
if ~isempty(winl)
    labels=[labels '_winl_' num2str(winl(1)*rate) '_' num2str(winl(2)*rate)];
    titletext=[titletext ' | ' 'winl' num2str(winl(1)) '-' num2str(winl(2)) 's' ];  
end
if isempty(wine)
    wine=winb;
end
if isempty(winhr)
    winhr=winb;
end
if isempty(winl)
    winl=winb;
end
%compile data
%first get good trials list for trial type big/small/targ
maplist={};     %listing of rows mapped by sel trials; ie. seltrials{1}(1) --> maplist{1}(1) = targ row for specific type & trial
seltrials={};
dtm={};
lfptm=[];
xtargda=[];
xtarglfp=[];
goodtrials=[];
labx='';
laby='';
for itype=1:length(sessiontypes)
    if strcmp(condtype,'all')
    goodtrialgrp=find(contains(trialinfo(itype).names,'all')==1);
    incswgrp=find(strcmp(trialinfo(itype).names,'switch')==1);
    goodtrials=sort(unique([trialinfo(itype).nums{goodtrialgrp} trialinfo(itype).nums{incswgrp}]));  
    else
    goodtrialgrp=find(contains(trialinfo(itype).names,condtype)==1);
    goodtrials=sort(unique([trialinfo(itype).nums{goodtrialgrp}]));          
    end
    seltrials{itype}=goodtrials;
    tlabelsgood=goodtrials+99;              %convert to trial labels ie 100,101,etc.
    typelistids=find(strcmp({trlist.type},sessiontypes{itype}));
    maplist(itype).trlist=find(ismember([trlist.id],tlabelsgood) & ...
        strcmp({trlist.type},sessiontypes{itype}));
    maplist(itype).trialnums=goodtrials;
    maplist(itype).type=sessiontypes{itype};    
    
%get hr  from binfos
hr=[];
rr=[];
rrstd=[];
rmssd=[];
hrtrials=[];
if isempty(xbinfos)
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},event));            %event matters since where calculated
    hr=binfos(btarg).pulse;
    hrtrials=binfos(btarg).seltrials;
else 
    %get hr from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'pulse') &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));            %event matters since where calculated
    if ~isempty(targrows)
    temp=xbinfos(targrows(1)).daall.lfptracesaln;
    hrtrials=xbinfos(targrows(1)).daall.trials;
        alnidx=xbinfos(targrows(1)).daall.mididx;
            hhwin=alnidx+winhr(1)*rate:alnidx+winhr(2)*rate;
    hrwin=temp(:,hhwin); 
    hr=nanmean(hrwin,2);
   rrtemp=1./temp.*60; 
   rrwin=rrtemp(:,hhwin);
   rr=nanmean(rrwin,2);
    rrstd=nanstd(rrwin,[],2).*1000; %in milliseconds
    %outofrange=find(rrstd<.01);
   % rrstd(outofrange)=nan;
    rmssd=rms(diff(rrtemp(:,hhwin),1,2),2);    
    end
    
end
for ida=1:length(hrtrials)
    tid=find(ismember(maplist(itype).trialnums,hrtrials(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).hr=hr(ida);
        if ~isempty(xbinfos)
        trlist(maplist(itype).trlist(tid)).rr=rr(ida);
        trlist(maplist(itype).trlist(tid)).rrstd=rrstd(ida);
        trlist(maplist(itype).trlist(tid)).rmssd=rmssd(ida);
        end
    end
end

%get lick  from binfos
lick=[];
licktrials=[];
if isempty(xbinfos)
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},event));            %event matters since where calculated
    lick=binfos(btarg).lickpost;
    licktrials=binfos(btarg).seltrials;  
else
    %get  from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'lick') &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));            %event matters since where calculated
        if ~isempty(targrows)

    temp=xbinfos(targrows(1)).daall.lfptracesaln;
    licktrials=xbinfos(targrows(1)).daall.trials;
            alnidx=xbinfos(targrows(1)).daall.mididx;

    licks=temp(:,alnidx+winl(1)*rate:alnidx+winl(2)*rate); 
    lick=nanmean(licks,2);
        end
end
for ida=1:length(hrtrials)
    tid=find(ismember(maplist(itype).trialnums,licktrials(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).lick=lick(ida);
    end
end 
    
%eye diameter, only right side valid
eye=[];
etrials=[];
eyevmax=[];
eyedist=[];
etrialsall=[];
btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
    strcmp({binfos.sessionid},sessid) &...
    strcmp({binfos.evt},event)); 
if isempty(xbinfos)    
           %event matters since where calculated
    eye=-binfos(btarg).reyed;
    etrials=binfos(btarg).seltrialsr;    
else
    %get  from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'eye') &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));   
    %pupil diameter, only right side, since opened door on left
            %many times   
        if ~isempty(targrows)
            etrialsall=xbinfos(targrows(1)).daall.trials;                     
            etrialsr=binfos(btarg).seltrialsr;      %right side only
            etrialsel=find(ismember(etrialsall,etrialsr));
            etrials=intersect(etrialsall,etrialsr);
            temp=xbinfos(targrows(1)).daall.lfptracesaln(etrialsel,:);
            alnidx=xbinfos(targrows(1)).daall.mididx;
            eewin=[alnidx+wine(1)*rate:alnidx+wine(2)*rate];
            eyewin=-temp(:,eewin); 
            eye=nanmean(eyewin,2);            
            %Z SCORE EYE 072019 from plottracesel
            glitchwidth=2;
            basewin=[alnidx-7*rate:alnidx+4*rate];          %fix baseline window for eye below
            %zscore eye wrt baseline as done in olivier et al 2014 frontier
            %remove blinks
            extwin=[alnidx-3*rate:alnidx+4*rate];
            relwin=intersect(eewin,extwin)-extwin(1)+1;
            datawin=-temp(:,extwin); %flip so original sign
            datawin=deglitchinterpsimple(datawin,3.5,glitchwidth);
            basedata=-temp(:,basewin);   %normalizing to pre-fix and targ window allows normalization trial to trial to look at faster changes within 1 sec 06/04 without phase trend eg. 67/68/69
            basedata=-basedata;
            for ix=1:size(basedata,1)
                [xx, xb]=rmoutliers(basedata(ix,:));
                basedata(ix,xb)=nan;
            end
            baseline=nanstd(basedata,[],2);          %normalized to each trial so independent of phase changes
            dataz=(datawin-nanmean(basedata,2))./baseline;
            datawin=dataz;
            eyewin=dataz(:,relwin);     %user selected windo
            eye=nanmean(eyewin,2);
        end
        %max velocity within 1 s of targ
        %eye distance integral of velocity
        targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
            strcmp({xbinfos.sitelfp},'eyev') &...
            strcmp({xbinfos.sessionid},sessid) &...
            strcmp({xbinfos.event},event)); 
       if ~isempty(targrows)
            etrialsall=xbinfos(targrows(1)).daall.trials;   
            temp=xbinfos(targrows(1)).daall.lfptracesaln;
            alnidx=xbinfos(targrows(1)).daall.mididx;
            eewin=[alnidx+wine(1)*rate:alnidx+wine(2)*rate];
            eyewin=-temp(:,eewin); 
            eyevmax=max(abs(eyewin),[],2);    
            eyevcum=cumtrapz(abs(eyewin),2);
            eyedist=max(eyevcum,[],2);
        end
end
for ida=1:length(etrials)
    tid=find(ismember(maplist(itype).trialnums,etrials(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).eye=eye(ida);
    end
end
for ida=1:length(etrialsall)
    tid=find(ismember(maplist(itype).trialnums,etrialsall(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).eyev=eyevmax(ida);
        trlist(maplist(itype).trlist(tid)).eyedist=eyedist(ida);
    end
end

%get signal values 
    dtm{itype}=[];    
    %datm/betatm values
    %get ch # for corresponding session, based on site label
    ich=0;
    lch=0;
    if ~isempty(siteda)
        %targeted da values from datm
    targsi=find([sites.sessnum]==sessnum & ...
    strcmp({sites.probeid},siteda));
    ich=sites(targsi(1)).ch;        %ch # defined by fscv recording ch
    dtm{itype}=datms{targses}{itype}{ich};           %get all da targeted values for current trial type big/small
    targtrials=dtm{itype}.trialnums;
            if ~isempty(xinfo)
                %compare trial #'s, get trial #'s from xinfo
                 targrows=find((strcmp({xinfo.sitelfp},sitelfp) & ...
                        strcmp({xinfo.siteda},siteda) & ...
                        strcmp({xinfo.event},event)) & ...
                        strcmp({xinfo.sessionid},sessid) &...
                        contains({xinfo.sessiontype},sessiontypes{itype})==1);   
                 xtrials=xinfo(targrows(1)).daall.trials;
                 targtrials=intersect(dtm{itype}.trialnums,xtrials);
            end
   % targtrials=find(ismember(dtm{itype}.trialnums,maplist(itype).trialnums));
    datatemp=getfield(dtm{itype},metric);
    targtrialids=find(ismember(dtm{itype}.trialnums,targtrials));
    dadata=datatemp(targtrialids);  
            for ida=1:length(targtrials)
                tid=find(ismember(maplist(itype).trialnums,targtrials(ida)));
                if ~isempty(tid)
                    trlist(maplist(itype).trlist(tid)).da=dadata(ida);
                end
            end
    end

    if ~isempty(sitelfp) 
        %targeted beta values from betatm
        if ~relfp
            %get precalc values in lfptms
            lfptm=lfptms{targses}{itype};
            for ii=1:length(lfptm)
                if strcmp(lfptm{ii}.site,sitelfp)
                    lch=ii;
                end
            end
            ltm{itype}=lfptm{lch};           %get all da targeted values for current trial type big/small
            targtrials=ltm{itype}.trialnums;
            if ~isempty(xinfo)
                %compare trial #'s, get trial #'s from xinfo
                 targrows=find((strcmp({xinfo.sitelfp},sitelfp) & ...
                        strcmp({xinfo.siteda},siteda) & ...
                        strcmp({xinfo.event},event)) & ...
                        strcmp({xinfo.sessionid},sessid) &...
                        contains({xinfo.sessiontype},sessiontypes{itype})==1);   
                 xtrials=xinfo(targrows(1)).daall.trials;
                 targtrials=intersect(ltm{itype}.trialnums,xtrials);
            end
          %  targtrials=find(ismember(ltm{itype}.trialnums,maplist(itype).trialnums));
            ltemp=getfield(ltm{itype},metriclfp);
                targtrialids=find(ismember(ltm{itype}.trialnums,targtrials));
            ldata=ltemp(targtrialids);
            for ida=1:length(targtrials)
                tid=find(ismember(maplist(itype).trialnums,targtrials(ida)));
                if ~isempty(tid)
                    trlist(maplist(itype).trlist(tid)).lfp=ldata(ida);
                end
            end
        elseif relfp && ~isempty(xinfo)
            %else recalculate from xinfos and defined win
            emptyflag=0;
            targrows=find((strcmp({xinfo.sitelfp},sitelfp) & ...
                contains({xinfo.siteda},siteda) & ...
            strcmp({xinfo.event},event)) & ...
            strcmp({xinfo.sessionid},sessid) &...
            contains({xinfo.sessiontype},sessiontypes{itype})==1);   
            if isempty(targrows)
                emptyflag=1;
                disp(['ch ' sitelfp ' not found in xinfo']);
            end
            if length(targrows)>1
                disp(['ch ' sitelfp ' & event ' event ' more than one instance in xinfo']);
                if isempty(siteda)
                    %choose lfp channel with more good trials, don't care
                    %about da pairing, only lfp signal used
                    trialsintargrows={xinfo(targrows).goodtrials};
                    maxid=targrows(1);
                    for itt=2:length(trialsintargrows)
                        if length(trialsintargrows{itt})>length(trialsintargrows{itt-1})
                            maxid=targrows(itt);
                        end
                    end
                    targrows=maxid;
                else                    
                    dbstop in plotcortrials2 at 563 if length(targrows)>1
                end
            end
            xtrials=xinfo(targrows(1)).daall.trials;
            alnidx=xinfo(targrows(1)).daall.mididx;
            lfpvals=[];
         %   targtrials=find(ismember(xtrials,maplist(itype).trialnums));
            if any(win~=0)
                %recalculate average based on user supplied window
            xdata=xinfo(targrows(1)).daall.lfptracesaln;
           % lfpwinids=[alnidx+win(1)*rate+offset(1)*rate:alnidx+win(2)*rate+offset(2)*rate];
           % %windowing used previously < 4/2020 for manual [0 1] wins
                     %   lfpwinids=[alnidx+win(1)*rate+1:alnidx+win(2)*rate+1];      %windowing used in targimwin generated in setTrialAx in plotsession
            lfpwinids=[alnidx+win(1)*rate+offset(1)*rate:alnidx+win(2)*rate+offset(2)*rate];
            lfpvals=nanmean(xdata(:,lfpwinids),2)';
            end
            if xvarflag
                %open xinfos variable as defined by metriclfp
                lfpvals=getfield(xinfo(targrows(1)).daall,metriclfp);
            end
            for ida=1:length(xtrials)
                tid=find(ismember(maplist(itype).trialnums,xtrials(ida)));
                if ~isempty(tid)
                    trlist(maplist(itype).trlist(tid)).lfp=lfpvals(ida);
                end
            end
        end
    end   
    
end
bigtype=double(strcmp({trlist.type},'bigreward'));
bigtype(bigtype==0)=nan;
smalltype=double(strcmp({trlist.type},'smallreward'));
smalltype(smalltype==0)=nan;
targtype=double(strcmp({trlist.type},'targetbreak'));
targtype(targtype==0)=nan;
fixtype=double(strcmp({trlist.type},'fixbreak'));
fixtype(fixtype==0)=nan;
sidetype=[trlist.side];
iax=1;
trialids=find(bigtype==1 | smalltype==1);
bigtrials=contains({trlist(trialids).type},'big');
list2=[];
list2.side=[trlist(trialids).side];              % 0 for left and 1 for right
list2.type=bigtrials;                            %0 for small and 1 for big
list2.trialnum=trialids;                         %trial index # for all successful/error trials

prevrew=contains({trlist(trialids(trialids>1)-1).type},'big');
if trialids(1)==1 && length(prevrew)~=length(trialids)
    %need to include first trial as 0, no reward
    prevrew=[0 prevrew];
end
prevsuc=contains({trlist(trialids(trialids>1)-1).type},'break');
if trialids(1)==1 && length(prevsuc)~=length(trialids)
    %need to include first trial as 0, no reward
    prevsuc=[0 prevsuc];
end
prevdni=double(prevsuc);        %do not include terms where previous trial was error for big vs small comparison
prevdni(prevdni==1)=nan;
list2.prevrew=double(prevrew)+prevdni;   %reward history, 0 if small previous reward, 1 if previously big reward, nan if error (ie remove trial)
list2.prevsuc=~prevsuc;               %performance history, 0 if previous error, 1 if previously big/small reward
%Only use one or the other history variable, if both, take out all data,
%because nan observations (rows) completely removed

Xvar=[];
Xvarnorm=[];
Xvarlabels={};
count=1;
origplotvars=plotvars;
for ix=1:length(origplotvars)
    %data matrix of regressors, except lfp
    if ~strcmp(origplotvars{ix},'lfp')
        datax=[trlist(trialids).(origplotvars{ix})]';
        badtids=find(isoutlier(datax)); %based on median outlier removal rather than mean 04/25/2020
        datax(badtids)=nan; %set outlier values to nan
        Xvar(:,count)=datax;
        datax=(datax-nanmean(datax))./ nanstd(datax); %standardize values (zscore) to get standardized beta coefficients from regress to evaluate how much influence each independent variable on y
        Xvarnorm(:,count)=datax;
        Xvarlabels=[Xvarlabels origplotvars{ix}];
        count=count+1;
    end
end
if ~isempty(otherplotvars)
    %include other regressors (e.g. side, type, trialnum, prevrew, etc.)
    plotvars=[plotvars otherplotvars];
    iid=size(Xvar,2)+1;
    for ix=1:length(otherplotvars)
        datax=[list2.(otherplotvars{ix})]';
        Xvar(:,iid)=datax;
        datax=(datax-nanmean(datax))./ nanstd(datax); %standardize values (zscore) to get standardized beta coefficients from regress to evaluate how much influence each independent variable on y
        Xvarnorm(:,iid)=datax;
        Xvarlabels=[Xvarlabels otherplotvars{ix}];
        iid=iid+1;
    end
end
%Xvar=[ones(size(Xvar,1),1) Xvar];
%Xvarnorm=[ones(size(Xvar,1),1) Xvarnorm];
%Xvarlabels=['intersect' Xvarlabels];

%don't need ones vector for fitlm, automatic
dadata=Xvar(:,1);
Yvar=[trlist(trialids).lfp]';
badtids=find(isoutlier(Yvar));
Yvar(badtids)=nan;
Yvarnorm=(Yvar-nanmean(Yvar))./ nanstd(Yvar);
Yvarlabel='lfp';
nonnanids=find(~isnan(dadata) & ~isnan(Yvar));
[r,psig]=corr(dadata(nonnanids),Yvar(nonnanids));       %find if da vs beta sig corr neg, assume DA first vecotr in Xvar
%mdl=fitlm(Xvar,Yvar,'VarNames',[Xvarlabels(1:end) Yvarlabel]);
mdl={};
%if r<0 && psig<0.05
mdl=fitlm(Xvarnorm,Yvarnorm,'VarNames',[Xvarlabels(1:end) Yvarlabel]);          %multiregression analysis with normalized values
p_mdl_betas=[mdl.Coefficients.pValue(2:end)];       %store beta coefficients p values what we care about if significant or not
coef_mdl_betas=[mdl.Coefficients.Estimate(2:end)];
dabeta=fitlm(Yvarnorm,Xvarnorm(:,1)); 
p_dabeta=dabeta.Coefficients.pValue(2);
b_dabeta=dabeta.Coefficients.Estimate(2);
res_dabeta=dabeta.Residuals.Raw;
betada=fitlm(Yvarnorm,Xvarnorm(:,1)); 
p_betada=betada.Coefficients.pValue(2);
b_betada=betada.Coefficients.Estimate(2);
res_betada=betada.Residuals.Raw;
dabehlm={};
p_dabeh=[];
b_dabeh=[];
res_dabeh=[];
beh=[];
corsmulti={};
cor_res_dabeh_betabeh={};
corsmulti=setfield(corsmulti,{1},['sessnum'],sessnum);
corsmulti=setfield(corsmulti,{1},['siteda'],siteda);
corsmulti=setfield(corsmulti,{1},['sitelfp'],sitelfp);
corsmulti=setfield(corsmulti,{1},['win'],win);
corsmulti=setfield(corsmulti,{1},['Xvar'],Xvar);
corsmulti=setfield(corsmulti,{1},['Yvar'],Yvar);
corsmulti=setfield(corsmulti,{1},['Xlabel'],Xvarlabels);
corsmulti=setfield(corsmulti,{1},['Ylabel'],Yvarlabel);
corsmulti=setfield(corsmulti,{1},['mdl'],mdl);
corsmulti=setfield(corsmulti,{1},['pvals'],p_mdl_betas);
corsmulti=setfield(corsmulti,{1},['betas'],coef_mdl_betas);

for ix=1:size(Xvar,2)-1
    %regress DA vs beh
    beh(ix).beh=Xvarlabels{ix+1};
    beh(ix).dabehlm=fitlm(Xvarnorm(:,ix+1), Xvarnorm(:,1));         %DA response / dependent variable, beh independent
     %   beh(ix).dabehlm=fitlm(Xvarnorm(:,1),Xvarnorm(:,ix+1) );         %DA response / dependent variable, beh independent

    beh(ix).p_dabeh = beh(ix).dabehlm.Coefficients.pValue(2);
    beh(ix).b_dabeh= beh(ix).dabehlm.Coefficients.Estimate(2);
    beh(ix).res_dabeh = beh(ix).dabehlm.Residuals.Raw;
    %Beta vs behs 
    beh(ix).betabehlm=fitlm(Yvarnorm,Xvarnorm(:,ix+1));         %beta response / dependent variable, beh independent
      %  beh(ix).betabehlm=fitlm(Xvarnorm(:,ix+1),Yvarnorm);         %beta response / dependent variable, beh independent

    beh(ix).p_betabeh= beh(ix).betabehlm.Coefficients.pValue(2);
    beh(ix).b_betabeh = beh(ix).betabehlm.Coefficients.Estimate(2);
    beh(ix).res_betabeh = beh(ix).betabehlm.Residuals.Raw;
  %  da_subres=res_dabeta-beh(ix).res_dabeh;     %subtract residuals from da vs beh from da vs beta  WHY DO WE NEED TO SUBTRACT FROM DA VS BETA?
   % beta_subres=res_betada-beh(ix).res_betabeh ;     %subtract residuals from beta vs beh from beta vs da
    nonnanids=find(~isnan(beh(ix).res_dabeh) & ~isnan(beh(ix).res_betabeh));
    [r,psig]=corr(beh(ix).res_dabeh(nonnanids),beh(ix).res_betabeh(nonnanids)); 
    beh(ix).r_res=r;          %correlation of residuals
    beh(ix).p_res=psig;
   
end
corsmulti=setfield(corsmulti,{1},'beh',beh);
%end
%[b,bint,r,rint,stats]  = regress(Yvarnorm,Xvarnorm);

outdata.sessnum=sessnum;
outdata.siteda=siteda;
outdata.sitelfp=sitelfp;
outdata.trlist=trlist;
outdata.Xvar=Xvar;
outdata.Xvarnorm=Xvarnorm;
outdata.Xvarlabels=Xvarlabels;
outdata.Yvar=Yvar;
outdata.Yvarnorm=Yvarnorm;
outdata.Yvarlabel=Yvarlabel;
outdata.r=r;
outdata.p=psig;
%outdata.regcoefs=b;
%for ix=1:length(Xvarlabels)
%outdata.(['b' num2str(ix)])=b(ix);
%outdata.(['x' num2str(ix)])=Xvarlabels{ix};
%end
outdata.metricda=metric;
outdata.metriclfp=metriclfp;
if relfp && ~xvarflag
    outdata.metriclfp=win;
end
outdata.trialids=trialids;
%outdata.bint=bint;
%outdata.r=r;
%outdata.rint=rint;
%outdata.stats=stats;
%{
%Make DA dependent Y variable, rather than LFP
Xvar=[];
Xvarnorm=[];
Xvarlabels={};
count=1;
plotvars=origplotvars;
for ix=1:length(origplotvars)
    %data matrix of regressors, except DA
    if ~strcmp(origplotvars{ix},'da')
        datax=[trlist(trialids).(origplotvars{ix})]';
        badtids=find(isoutlier(datax)); %based on median outlier removal rather than mean 04/25/2020
        datax(badtids)=nan; %set outlier values to nan
        Xvar(:,count)=datax;
        datax=(datax-nanmean(datax))./ nanstd(datax); %standardize values (zscore) to get standardized beta coefficients from regress to evaluate how much influence each independent variable on y
        Xvarnorm(:,count)=datax;
        Xvarlabels=[Xvarlabels origplotvars{ix}];
        count=count+1;
    end
end
if ~isempty(otherplotvars)
    %include other regressors (e.g. side, type, trialnum, prevrew, etc.)
    plotvars=[plotvars otherplotvars];
    iid=size(Xvar,2)+1;
    for ix=1:length(otherplotvars)
        datax=[list2.(otherplotvars{ix})]';
        Xvar(:,iid)=datax;
        datax=(datax-nanmean(datax))./ nanstd(datax); %standardize values (zscore) to get standardized beta coefficients from regress to evaluate how much influence each independent variable on y
        Xvarnorm(:,iid)=datax;
        Xvarlabels=[Xvarlabels otherplotvars{ix}];
        iid=iid+1;
    end
end
Xvar=[ones(size(Xvar,1),1) Xvar];
Xvarnorm=[ones(size(Xvar,1),1) Xvarnorm];
Xvarlabels=['intersect' Xvarlabels];
Yvar=[trlist(trialids).da]';
badtids=find(isoutlier(Yvar));
Yvar(badtids)=nan;
Yvarnorm=(Yvar-nanmean(Yvar))./ nanstd(Yvar);
Yvarlabel='da';
[b,bint,r,rint,stats]  = regress(Yvarnorm,Xvarnorm);

dadata={};
dadata.sessnum=sessnum;
dadata.siteda=siteda;
dadata.sitelfp=sitelfp;
dadata.trlist=trlist;
dadata.Xvar=Xvar;
dadata.Xvarnorm=Xvarnorm;
dadata.Xvarlabels=Xvarlabels;
dadata.Yvar=Yvar;
dadata.Yvarnorm=Yvarnorm;
dadata.Yvarlabel=Yvarlabel;
dadata.regcoefs=b;
for ix=1:length(Xvarlabels)
dadata.(['b' num2str(ix)])=b(ix);
dadata.(['x' num2str(ix)])=Xvarlabels{ix};
end
dadata.metricda=metric;
dadata.metriclfp=metriclfp;
if relfp && ~xvarflag
    dadata.metriclfp=win;
end
dadata.trialids=trialids;
dadata.bint=bint;
dadata.r=r;
dadata.rint=rint;
dadata.stats=stats;
%}

lfpdata=[];
lfpdata=outdata;
if multicor
lfpdata=corsmulti;
end


savename=[savepath '_' labels];
save(savename,'lfpdata','mdl')


end


